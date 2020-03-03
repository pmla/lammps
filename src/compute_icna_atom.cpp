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
#include "citeme.h"
#include <mpi.h>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
#include <cassert>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

static const char cite_user_ptm_package[] =
    "i-CNA package:\n\n"
    "@Article{larsen2020revisiting,\n"
    " author={Larsen, Peter Mahler},\n"
    " title={Revisiting the Common Neighbour Analysis and the Centrosymmetry Parameter},\n"
    " journal={Modelling~Simul.~Mater.~Sci.~Eng.},\n"
    " year={2020},\n"
    " number={XXX},\n"
    " volume={XXX},\n"
    " pages={XXX},\n"
    " DOI = {XXX}"
    "}\n\n";


enum{UNKNOWN,FCC,HCP,BCC,ICOS,OTHER};

#define MAX_NEIGHBORS 14

/* ---------------------------------------------------------------------- */

ComputeICNAAtom::ComputeICNAAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  list(NULL), pattern(NULL)
{
  if (narg != 3)
    error->all(FLERR,"Illegal compute icna/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeICNAAtom::~ComputeICNAAtom()
{
  memory->destroy(pattern);
}

/* ---------------------------------------------------------------------- */

void ComputeICNAAtom::init()
{
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

/// Pair of neighbor atoms that form a bond (bit-wise storage).
typedef unsigned int CNAPairBond;

/**
 * A bit-flag array indicating which pairs of neighbors are bonded
 * and which are not.
 */
struct NeighborBondArray
{
	/// Two-dimensional bit array that stores the bonds between neighbors.
	unsigned int neighborArray[32];

	/// Resets all bits.
	NeighborBondArray() {
		memset(neighborArray, 0, sizeof(neighborArray));
	}

	/// Returns whether two nearest neighbors have a bond between them.
	inline bool neighborBond(int neighborIndex1, int neighborIndex2) const {
		assert(neighborIndex1 < 32);
		assert(neighborIndex2 < 32);
		return (neighborArray[neighborIndex1] & (1<<neighborIndex2));
	}

	/// Sets whether two nearest neighbors have a bond between them.
	inline void setNeighborBond(int neighborIndex1, int neighborIndex2, bool bonded) {
		assert(neighborIndex1 < 32);
		assert(neighborIndex2 < 32);
		if(bonded) {
			neighborArray[neighborIndex1] |= (1<<neighborIndex2);
			neighborArray[neighborIndex2] |= (1<<neighborIndex1);
		}
		else {
			neighborArray[neighborIndex1] &= ~(1<<neighborIndex2);
			neighborArray[neighborIndex2] &= ~(1<<neighborIndex1);
		}
	}
};


/******************************************************************************
* Find all atoms that are nearest neighbors of the given pair of atoms.
******************************************************************************/
int findCommonNeighbors(const NeighborBondArray& neighborArray, int neighborIndex, unsigned int& commonNeighbors, int numNeighbors)
{
	commonNeighbors = neighborArray.neighborArray[neighborIndex];
#ifndef Q_CC_MSVC
	// Count the number of bits set in neighbor bit-field.
	return __builtin_popcount(commonNeighbors);
#else
	// Count the number of bits set in neighbor bit-field.
	unsigned int v = commonNeighbors - ((commonNeighbors >> 1) & 0x55555555);
	v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
	return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
#endif
}

/******************************************************************************
* Finds all bonds between common nearest neighbors.
******************************************************************************/
int findNeighborBonds(const NeighborBondArray& neighborArray, unsigned int commonNeighbors, int numNeighbors, CNAPairBond* neighborBonds)
{
	int numBonds = 0;

	unsigned int nib[32];
	int nibn = 0;
	unsigned int ni1b = 1;
	for(int ni1 = 0; ni1 < numNeighbors; ni1++, ni1b <<= 1) {
		if(commonNeighbors & ni1b) {
			unsigned int b = commonNeighbors & neighborArray.neighborArray[ni1];
			for(int n = 0; n < nibn; n++) {
				if(b & nib[n]) {
					neighborBonds[numBonds++] = ni1b | nib[n];
				}
			}
			nib[nibn++] = ni1b;
		}
	}
	return numBonds;
}

/******************************************************************************
* Find all chains of bonds.
******************************************************************************/
static int getAdjacentBonds(unsigned int atom, CNAPairBond* bondsToProcess, int& numBonds, unsigned int& atomsToProcess, unsigned int& atomsProcessed)
{
    int adjacentBonds = 0;
	for(int b = numBonds - 1; b >= 0; b--) {
		if(atom & *bondsToProcess) {
            ++adjacentBonds;
   			atomsToProcess |= *bondsToProcess & (~atomsProcessed);
   			memmove(bondsToProcess, bondsToProcess + 1, sizeof(CNAPairBond) * b);
   			numBonds--;
		}
		else ++bondsToProcess;
	}
	return adjacentBonds;
}

/******************************************************************************
* Find all chains of bonds between common neighbors and determine the length
* of the longest continuous chain.
******************************************************************************/
int calcMaxChainLength(CNAPairBond* neighborBonds, int numBonds)
{
    // Group the common bonds into clusters.
	int maxChainLength = 0;
	while(numBonds) {
        // Make a new cluster starting with the first remaining bond to be processed.
		numBonds--;
        unsigned int atomsToProcess = neighborBonds[numBonds];
        unsigned int atomsProcessed = 0;
		int clusterSize = 1;
        do {
#ifndef Q_CC_MSVC
        	// Determine the number of trailing 0-bits in atomsToProcess, starting at the least significant bit position.
			int nextAtomIndex = __builtin_ctz(atomsToProcess);
#else
			unsigned long nextAtomIndex;
			_BitScanForward(&nextAtomIndex, atomsToProcess);
			assert(nextAtomIndex >= 0 && nextAtomIndex < 32);
#endif
			unsigned int nextAtom = 1 << nextAtomIndex;
        	atomsProcessed |= nextAtom;
			atomsToProcess &= ~nextAtom;
			clusterSize += getAdjacentBonds(nextAtom, neighborBonds, numBonds, atomsToProcess, atomsProcessed);
		}
        while(atomsToProcess);
        if(clusterSize > maxChainLength)
        	maxChainLength = clusterSize;
	}
	return maxChainLength;
}

int analyzeSmallSignature(NeighborBondArray& neighborArray)
{
	int nn = 12;
	int n421 = 0;
	int n422 = 0;
	int n555 = 0;
	for(int ni = 0; ni < nn; ni++) {

		// Determine number of neighbors the two atoms have in common.
		unsigned int commonNeighbors;
		int numCommonNeighbors = findCommonNeighbors(neighborArray, ni, commonNeighbors, nn);
		if(numCommonNeighbors != 4 && numCommonNeighbors != 5)
			break;

		// Determine the number of bonds among the common neighbors.
		CNAPairBond neighborBonds[MAX_NEIGHBORS*MAX_NEIGHBORS];
		int numNeighborBonds = findNeighborBonds(neighborArray, commonNeighbors, nn, neighborBonds);
		if(numNeighborBonds != 2 && numNeighborBonds != 5)
			break;

		// Determine the number of bonds in the longest continuous chain.
		int maxChainLength = calcMaxChainLength(neighborBonds, numNeighborBonds);
		if(numCommonNeighbors == 4 && numNeighborBonds == 2) {
			if(maxChainLength == 1) n421++;
			else if(maxChainLength == 2) n422++;
			else break;
		}
		else if(numCommonNeighbors == 5 && numNeighborBonds == 5 && maxChainLength == 5) n555++;
		else break;
	}
	if(n421 == 12) return FCC;
	else if(n421 == 6 && n422 == 6) return HCP;
	else if(n555 == 12) return ICOS;
	return OTHER;
}

int analyzeLargeSignature(NeighborBondArray& neighborArray)
{
	int nn = 14;
	int n444 = 0;
	int n666 = 0;
	for(int ni = 0; ni < nn; ni++) {

		// Determine number of neighbors the two atoms have in common.
		unsigned int commonNeighbors;
		int numCommonNeighbors = findCommonNeighbors(neighborArray, ni, commonNeighbors, nn);
		if(numCommonNeighbors != 4 && numCommonNeighbors != 6)
			break;

		// Determine the number of bonds among the common neighbors.
		CNAPairBond neighborBonds[MAX_NEIGHBORS*MAX_NEIGHBORS];
		int numNeighborBonds = findNeighborBonds(neighborArray, commonNeighbors, nn, neighborBonds);
		if(numNeighborBonds != 4 && numNeighborBonds != 6)
			break;

		// Determine the number of bonds in the longest continuous chain.
		int maxChainLength = calcMaxChainLength(neighborBonds, numNeighborBonds);
		if(numCommonNeighbors == 4 && numNeighborBonds == 4 && maxChainLength == 4) n444++;
		else if(numCommonNeighbors == 6 && numNeighborBonds == 6 && maxChainLength == 6) n666++;
		else break;
	}
	if(n666 == 8 && n444 == 6) return BCC;
	return OTHER;
}


typedef struct {
  int index;
  double d;
} nbr_t;

static bool sorthelper_compare(nbr_t const &a, nbr_t const &b) {
  return a.d < b.d;
}

static double squared_distance(double* a, double* b)
{
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;    
}

static int get_neighbours(int *numneigh, int **firstneigh, double **x,
                          size_t atom_index, int num, double (*nbr_pos)[3], double* lengths)
{
  double *pos = x[atom_index];

  int *jlist = NULL;
  int jnum = 0;
  jlist = firstneigh[atom_index];
  jnum = numneigh[atom_index];

  std::vector<nbr_t> nbr_order;
  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];
    j &= NEIGHMASK;
    if (j == atom_index)
      continue;

    nbr_t nbr = {j, squared_distance(pos, x[j])};
    nbr_order.push_back(nbr);
  }

  std::sort(nbr_order.begin(), nbr_order.end(), &sorthelper_compare);
  int num_nbrs = std::min(num, (int)nbr_order.size());

  for (int jj = 0; jj < num_nbrs; jj++) {
    int j = nbr_order[jj].index;
    nbr_pos[jj][0] = x[j][0] - pos[0];
    nbr_pos[jj][1] = x[j][1] - pos[1];
    nbr_pos[jj][2] = x[j][2] - pos[2];

    lengths[jj] = nbr_order[jj].d;
  }

  return num_nbrs;
}

enum GraphEdgeType {
  NONE,
  SHORT,
  LONG
};

struct GraphEdge {
  GraphEdge(int _i, int _j, double _length, int _edgeType)
    : i(_i), j(_j), length(_length), edgeType(_edgeType) {}

  int i = 0;
  int j = 0;
  double length = 0;
  int edgeType = 0;
  GraphEdge* nextShort = NULL;
  GraphEdge* nextLong = NULL;
};

/******************************************************************************
* Builds an edge list sorted by length
******************************************************************************/
class EdgeIterator {
public:
  EdgeIterator(int nn, double (*neighborVectors)[3], double shortThreshold, double longThreshold) {

    if (nn < 12) shortThreshold = 0;
    if (nn < 14) longThreshold = 0;

    // End points are the shortest edges lengths which exceed their respective thresholds.
    GraphEdge shortEnd(-1, -1, INFINITY, SHORT);
    GraphEdge longEnd(-1, -1, INFINITY, LONG);

    // Find edges which will make up intervals.
    for (int i=0;i<nn;i++) {
      for (int j=i+1;j<nn;j++) {
        double length = sqrt(squared_distance(neighborVectors[i], neighborVectors[j]));

        int edgeType = NONE;
        if (i < 12 && j < 12 && length < shortThreshold) {
          edgeType |= SHORT;
        }
        if (length < longThreshold) {
          edgeType |= LONG;
        }

        if (edgeType == NONE) {
          if (length < longEnd.length) {
            longEnd = GraphEdge(i, j, length, LONG);
          }
          else if (length < shortEnd.length) {
            shortEnd = GraphEdge(i, j, length, SHORT);
          }
        }
        else {
          edges.push_back(GraphEdge(i, j, length, edgeType));
        }
      }
    }

    // Sort edges by length to create intervals.
    std::sort(edges.begin(), edges.end(), &edge_compare);

    if (shortEnd.i != -1) {
      edges.push_back(shortEnd);
    }
    if (longEnd.i != -1) {
      edges.push_back(longEnd);
    }

    // Create two paths through intervals: short and long.
    for (int i=edges.size() - 1;i>=0;i--) {
      if (edges[i].edgeType & SHORT) {
        edges[i].nextShort = nextShort;
        nextShort = &edges[i];
      }
      if (edges[i].edgeType & LONG) {
        edges[i].nextLong = nextLong;
        nextLong = &edges[i];
      }
    }
  }

  static bool edge_compare(GraphEdge const &a, GraphEdge const &b) {
    return a.length < b.length;
  }

  std::vector< GraphEdge > edges;
  GraphEdge* nextLong = NULL;
  GraphEdge* nextShort = NULL;
};

static int analyze_atom(int *numneigh, int **firstneigh, double **positions,
                          size_t atom_index)
{
	// Find nearest neighbors of current atom.
    #define ICNA_MAX_NBRS 17
    double nbrs[ICNA_MAX_NBRS][3];
    double neighborLengths[ICNA_MAX_NBRS];
    int num_found = get_neighbours(numneigh, firstneigh, positions, atom_index,
                                   ICNA_MAX_NBRS, nbrs, neighborLengths);

	// Determine which structure types to search for.
	bool analyzeShort = num_found >= 12;
	bool analyzeLong = num_found >= 14;
	if (analyzeLong) num_found = 14;
	else if (analyzeShort) num_found = 12;
	else return OTHER;

	// We will set the threshold for interval start points two thirds of the way between
	// the first and second neighbor shells.
	const double x = 2.0f / 3.0f;
	const double fraction = ((1 - x) * 1 + x * sqrt(2));

	// Calculate length thresholds and local scaling factors.
	double shortLengthThreshold = 0, longLengthThreshold = 0;

	if (analyzeShort) {
		int num = 12;
		double shortLocalScaling = 0;
		for(int i = 0; i < num; i++)
			shortLocalScaling += neighborLengths[i];
		shortLocalScaling /= num;
		shortLengthThreshold = fraction * shortLocalScaling;
	}
	if (analyzeLong) {
		int num = 14;
		double longLocalScaling = 0;
		for(int i = 0; i < 8; i++)
			longLocalScaling += neighborLengths[i] / sqrt(3.0f / 4.0f);
		for(int i = 8; i < num; i++)
			longLocalScaling += neighborLengths[i];
		longLocalScaling /= num;
		longLengthThreshold = fraction * longLocalScaling;
	}

	// Use interval width to resolve ambiguities in traditional CNA classification
	double bestIntervalWidth = 0;
	int bestType = OTHER;

	EdgeIterator it = EdgeIterator(num_found, nbrs, shortLengthThreshold, longLengthThreshold);

	/////////// 12 neighbors ///////////
	if(analyzeShort) {
		const int num = 12; //Number of neighbors to analyze for FCC, HCP and Icosahedral atoms
		int n4 = 0, n5 = 0;
		int coordinations[num] = {0};
		NeighborBondArray neighborArray;

		GraphEdge* edge = it.nextShort;
		GraphEdge* next = edge != NULL ? edge->nextShort : NULL;
		while (next != NULL) {
			coordinations[edge->i]++;
			coordinations[edge->j]++;
			neighborArray.setNeighborBond(edge->i, edge->j, true);

			if (coordinations[edge->i] == 4) n4++;
			if (coordinations[edge->i] == 5) {n4--; n5++;}
			if (coordinations[edge->i] > 5) break;

			if (coordinations[edge->j] == 4) n4++;
			if (coordinations[edge->j] == 5) {n4--; n5++;}
			if (coordinations[edge->j] > 5) break;

			if (n4 == num || n5 == num) {
				// Coordination numbers are correct - perform traditional CNA
				int type = analyzeSmallSignature(neighborArray);
				if (type != OTHER) {
					double intervalWidth = next->length - edge->length;
					if (intervalWidth > bestIntervalWidth) {
						bestIntervalWidth = intervalWidth;
						bestType = type;
					}
				}
			}

			edge = next;
			next = next->nextShort;
		}
	}

	/////////// 14 neighbors ///////////
	if(analyzeLong) {
		const int num = 14; //Number of neighbors to analyze for BCC atoms
		int n4 = 0, n6 = 0;
		int coordinations[num] = {0};
		NeighborBondArray neighborArray;

		GraphEdge* edge = it.nextLong;
		GraphEdge* next = edge != NULL ? edge->nextLong : NULL;
		while (next != NULL) {
			coordinations[edge->i]++;
			coordinations[edge->j]++;
			neighborArray.setNeighborBond(edge->i, edge->j, true);

			if (coordinations[edge->i] == 4) n4++;
			if (coordinations[edge->i] == 5) n4--;
			if (coordinations[edge->i] == 6) n6++;
			if (coordinations[edge->i] > 6) break;

			if (coordinations[edge->j] == 4) n4++;
			if (coordinations[edge->j] == 5) n4--;
			if (coordinations[edge->j] == 6) n6++;
			if (coordinations[edge->j] > 6) break;

			if (n4 == 6 && n6 == 8) {
				// Coordination numbers are correct - perform traditional CNA
				int type = analyzeLargeSignature(neighborArray);
				if (type != OTHER) {
					double intervalWidth = next->length - edge->length;
					if (intervalWidth > bestIntervalWidth) {
						bestIntervalWidth = intervalWidth;
						bestType = type;
					}
				}
			}

			edge = next;
			next = next->nextLong;
		}
	}

	return bestType;
}

void ComputeICNAAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow arrays if necessary
  if (atom->nmax > nmax) {
    memory->destroy(pattern);
    nmax = atom->nmax;
    memory->create(pattern,nmax,"icna:icna_pattern");
    vector_atom = pattern;
  }

  // invoke full neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // compute ICNA for each atom in group
  int nerror = 0;
  for (int ii = 0; ii < inum; ii++) {
    int atom_index = ilist[ii];

    if (!(mask[atom_index] & groupbit)) {
      pattern[atom_index] = UNKNOWN;
    }
    else {
      pattern[atom_index] = analyze_atom(numneigh, firstneigh, x, atom_index);
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeICNAAtom::memory_usage()
{
  double bytes = nmax * sizeof(int);
  bytes += nmax * sizeof(double);
  return bytes;
}

