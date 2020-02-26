/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: PM Larsen (Technical University of Denmark)
------------------------------------------------------------------------- */

#include "compute_icna_atom.h"
#include <mpi.h>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXNEAR 16
#define MAXCOMMON 8

enum{UNKNOWN,FCC,HCP,BCC,ICOS,OTHER};
enum{NCOMMON,NBOND,MAXBOND,MINBOND};

/* ---------------------------------------------------------------------- */

ComputeICNAAtom::ComputeICNAAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  list(NULL), nearest(NULL), nnearest(NULL), pattern(NULL)
{
  if (narg != 4) error->all(FLERR,"Illegal compute icna/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  double cutoff = force->numeric(FLERR,arg[3]);
  if (cutoff < 0.0) error->all(FLERR,"Illegal compute icna/atom command");
  cutsq = cutoff*cutoff;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeICNAAtom::~ComputeICNAAtom()
{
  memory->destroy(nearest);
  memory->destroy(nnearest);
  memory->destroy(pattern);
}

/* ---------------------------------------------------------------------- */

void ComputeICNAAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute icna/atom requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute icna/atom cutoff is longer than pairwise cutoff");

  // cannot use neighbor->cutneighmax b/c neighbor has not yet been init

  if (2.0*sqrt(cutsq) > force->pair->cutforce + neighbor->skin &&
      comm->me == 0)
    error->warning(FLERR,"Compute icna/atom cutoff may be too large to find "
                   "ghost atom neighbors");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"icna/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute icna/atom defined");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeICNAAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

static int get_neighbors(void* list, int atom_index, double (*nbrs)[3])
{
  int i = list->ilist[atom_index];
  int jlist = list->firstneigh[i];
  int num_available = list->numneigh[i];

  double sizes[

  int n = 0;
  for (int jj=0;jj<num_available;jj++) {
    int j = jlist[jj];
    if (mask[j] & groupbit) {
        j &= NEIGHMASK;

        double dx = x[i][0] - x[j][0];
        double dy = x[i][1] - x[j][1];
        double dz = x[i][2] - x[j][2];
        double rsq = dx * dx + dy * dy + dz * dz;
        n++;
    }
  }

  return num_found;
}

static double squared_distance(double* a, double* b)
{
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;    
}

static int get_neighbours(int *numneigh, int **firstneigh, double **x,
                          size_t atom_index, int num, double (*nbr_pos)[3])
{
  double *pos = x[atom_index];

  int *jlist = NULL;
  int jnum = 0;
  jlist = data->firstneigh[atom_index];
  jnum = data->numneigh[atom_index];

  std::vector<nbr_t> nbr_order;
  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];
    j &= NEIGHMASK;
    if (j == atom_index)
      continue;

    nbr_t nbr = {j, square_distanced(pos, x[j])};
    nbr_order.push_back(nbr);
  }

  std::sort(nbr_order.begin(), nbr_order.end(), &sorthelper_compare);
  int num_nbrs = std::min(num, (int)nbr_order.size());

  for (int jj = 0; jj < num_nbrs; jj++) {
    int j = nbr_order[jj].index;
    nbr_pos[jj][0] = x[j][0] - pos[0];
    nbr_pos[jj][1] = x[j][1] - pos[1];
    nbr_pos[jj][2] = x[j][2] - pos[2];
  }

  return num_nbrs;
}

void ComputeICNAAtom::compute_peratom()
{
  int i,j,k,ii,jj,kk,m,n,inum,jnum,inear,jnear;
  int firstflag,ncommon,nbonds,maxbonds,minbonds;
  int nfcc,nhcp,nbcc4,nbcc6,nico,cj,ck,cl,cm;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int icna[MAXNEAR][4],onenearest[MAXNEAR];
  int common[MAXCOMMON],bonds[MAXCOMMON];
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;

  invoked_peratom = update->ntimestep;

  // grow arrays if necessary

  if (atom->nmax > nmax) {
    //memory->destroy(nearest);
    //memory->destroy(nnearest);
    memory->destroy(pattern);
    nmax = atom->nmax;

    //memory->create(nearest,nmax,MAXNEAR,"icna:nearest");
    //memory->create(nnearest,nmax,"icna:nnearest");
    memory->create(pattern,nmax,"icna:icna_pattern");
    vector_atom = pattern;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // find the neigbours of each atom within cutoff using full neighbor list
  // nearest[] = atom indices of nearest neighbors, up to MAXNEAR
  // do this for all atoms, not just compute group
  // since ICNA calculation requires neighbors of neighbors

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // compute ICNA for each atom in group
  // only performed if # of nearest neighbors = 12 or 14 (fcc,hcp)

  int nerror = 0;
  for (ii = 0; ii < inum; ii++) {
    int atom_index = ilist[ii];

    if (!(mask[atom_index] & groupbit)) {
      pattern[atom_index] = UNKNOWN;
      continue;
    }

    #define ICNA_MAX_NBRS 17
    pattern[atom_index] = OTHER;
    double nbrs[ICNA_MAX_NBRS][3];
    int num_found = get_neighbours(numneigh, firstneigh, x, atom_index, ICNA_MAX_NBRS, nbrs);
    if (num_found < 12) {
        continue;
    }

    double distances[ICNA_MAX_NBRS][ICNA_MAX_NBRS];
    for (int j=0;j<num_found;j++) {
        distances[j][j] = INFINITY;
        for (int i=j+1;i<num_found;i++) {
            distances[i, j] = squared_distance(nbrs[i], nbrs[j]);
        }
    }

    // loop over near neighbors of I to build icna data structure
    // icna[k][NCOMMON] = # of common neighbors of I with each of its neighs
    // icna[k][NBONDS] = # of bonds between those common neighbors
    // icna[k][MAXBOND] = max # of bonds of any common neighbor
    // icna[k][MINBOND] = min # of bonds of any common neighbor

    for (m = 0; m < nnearest[i]; m++) {
      j = nearest[i][m];

      // common = list of neighbors common to atom I and atom J
      // if J is an owned atom, use its near neighbor list to find them
      // if J is a ghost atom, use full neighbor list of I to find them
      // in latter case, must exclude J from I's neighbor list

      if (j < nlocal) {
        firstflag = 1;
        ncommon = 0;
        for (inear = 0; inear < nnearest[i]; inear++)
          for (jnear = 0; jnear < nnearest[j]; jnear++)
            if (nearest[i][inear] == nearest[j][jnear]) {
              if (ncommon < MAXCOMMON) common[ncommon++] = nearest[i][inear];
              else if (firstflag) {
                nerror++;
                firstflag = 0;
              }
            }

      } else {
        xtmp = x[j][0];
        ytmp = x[j][1];
        ztmp = x[j][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        n = 0;
        for (kk = 0; kk < jnum; kk++) {
          k = jlist[kk];
          k &= NEIGHMASK;
          if (k == j) continue;

          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq) {
            if (n < MAXNEAR) onenearest[n++] = k;
            else break;
          }
        }

        firstflag = 1;
        ncommon = 0;
        for (inear = 0; inear < nnearest[i]; inear++)
          for (jnear = 0; (jnear < n) && (n < MAXNEAR); jnear++)
            if (nearest[i][inear] == onenearest[jnear]) {
              if (ncommon < MAXCOMMON) common[ncommon++] = nearest[i][inear];
              else if (firstflag) {
                nerror++;
                firstflag = 0;
              }
            }
      }

      icna[m][NCOMMON] = ncommon;

      // calculate total # of bonds between common neighbor atoms
      // also max and min # of common atoms any common atom is bonded to
      // bond = pair of atoms within cutoff

      for (n = 0; n < ncommon; n++) bonds[n] = 0;

      nbonds = 0;
      for (jj = 0; jj < ncommon-1; jj++) {
        j = common[jj];
        xtmp = x[j][0];
        ytmp = x[j][1];
        ztmp = x[j][2];
        for (kk = jj+1; kk < ncommon; kk++) {
          k = common[kk];
          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq) {
            nbonds++;
            bonds[jj]++;
            bonds[kk]++;
          }
        }
      }

      icna[m][NBOND] = nbonds;

      maxbonds = 0;
      minbonds = MAXCOMMON;
      for (n = 0; n < ncommon; n++) {
        maxbonds = MAX(bonds[n],maxbonds);
        minbonds = MIN(bonds[n],minbonds);
      }
      icna[m][MAXBOND] = maxbonds;
      icna[m][MINBOND] = minbonds;
    }

    // detect ICNA pattern of the atom
    pattern[i] = OTHER;
    nfcc = nhcp = nbcc4 = nbcc6 = nico = 0;

    if (nnearest[i] == 12) {
      for (inear = 0; inear < 12; inear++) {
        cj = icna[inear][NCOMMON];
        ck = icna[inear][NBOND];
        cl = icna[inear][MAXBOND];
        cm = icna[inear][MINBOND];
        if (cj == 4 && ck == 2 && cl == 1 && cm == 1) nfcc++;
        else if (cj == 4 && ck == 2 && cl == 2 && cm == 0) nhcp++;
        else if (cj == 5 && ck == 5 && cl == 2 && cm == 2) nico++;
      }
      if (nfcc == 12) pattern[i] = FCC;
      else if (nfcc == 6 && nhcp == 6) pattern[i] = HCP;
      else if (nico == 12) pattern[i] = ICOS;

    } else if (nnearest[i] == 14) {
      for (inear = 0; inear < 14; inear++) {
        cj = icna[inear][NCOMMON];
        ck = icna[inear][NBOND];
        cl = icna[inear][MAXBOND];
        cm = icna[inear][MINBOND];
        if (cj == 4 && ck == 4 && cl == 2 && cm == 2) nbcc4++;
        else if (cj == 6 && ck == 6 && cl == 2 && cm == 2) nbcc6++;
      }
      if (nbcc4 == 6 && nbcc6 == 8) pattern[i] = BCC;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeICNAAtom::memory_usage()
{
  double bytes = nmax * sizeof(int);
  bytes += nmax * MAXNEAR * sizeof(int);
  bytes += nmax * sizeof(double);
  return bytes;
}
