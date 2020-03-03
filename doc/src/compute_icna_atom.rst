.. index:: compute icna/atom

compute icna/atom command
=========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID icna/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* icna/atom = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all icna/atom

Description
"""""""""""

Define a computation that calculates the i-CNA (Interval Common Neighbor
Analysis) pattern for each atom in the group.  In solid-state systems
the i-CNA pattern is a useful measure of the local crystal structure
around an atom.  The i-CNA method is a parameter-free adaptation of the
:doc:`conventional CNA method <compute_cna_atom>`; it determines a per-atom
cutoff by exhaustive analysis of the coordination intervals. The method is
described in :ref:`(Larsen) <Larsen>`.

Currently, there are five kinds of i-CNA patterns LAMMPS recognizes:

* fcc = 1
* hcp = 2
* bcc = 3
* icosahedral = 4
* unknown = 5

The value of the i-CNA pattern will be 0 for atoms not in the specified
compute group.  Note that normally a i-CNA calculation should only be
performed on mono-component systems.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (e.g. each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently or to have multiple compute/dump commands, each with a
*icna/atom* style.

**Output info:**

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-atom vector values will be a number from 0 to 5, as explained
above.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute cna/atom <compute_cna_atom>`

:doc:`compute centro/atom <compute_centro_atom>`

**Default:** none


----------


.. _LarsenICNA:



**(Larsen)** PM Larsen, XXX, XXX, XXX (2020).
