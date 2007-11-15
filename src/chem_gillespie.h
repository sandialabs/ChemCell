/* ----------------------------------------------------------------------
   ChemCell - Cell simulator for particle diffusion and reactions
   Steve Plimpton (sjplimp@sandia.gov), Alex Slepoy (aslepoy@sandia.gov)
   Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/chemcell.html

   Copyright (2004) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level ChemCell directory.
------------------------------------------------------------------------- */

#ifndef CHEM_GILLESPIE_H
#define CHEM_GILLESPIE_H

#include "chem.h"

class ChemGillespie : public Chem {
 public:
  ChemGillespie();
  ~ChemGillespie();
  void init();
  void reactions();

 private:

  // tree stores summed propensities as binary tree
  // tree is stored contiguously, nodes followed by leaves
  // leaves are linear list of propensity values, padded with 0's to 2^n
  // each node is sum of 2 children
  // tree[0] is the sum of all propensities
  // left child is odd index, right child is even index
  // parent of tree node i = (i-1)/2
  // children of tree node i = 2*i + 1, 2*(i+1)
  // offset = always 2^n - 1 = loc where 1st leaf starts
  // unused tree values (both leaves and nodes) must be init to 0

  int allocated;

  int nreactions;            // # of reactions
  int *nreactant;            // # of reactants for each reaction
  int **reactants;           // i,j = jth reactant of ith reaction
  int *nproduct;             // # of products for each reaction
  int **products;            // i,j = jth product of ith reaction
  double *rate;              // input rate for each reaction
  int *pcount;               // count of particles of each species

  int *ndepends;             // # of reactions that depend on each reaction
  int **depends;             // i,j = jth reaction that depends on ith reaction

  double *tree;              // tree structure of summed propensities
  double *propensity;        // propensity of each reaction
  int offset;                // loc in tree where 1st leaf starts

  double factor_zero;        // conversion factor for different reactions
  double factor_dual;

  void set(int, double);
  int find(double);
  void sum();
  void build_dependency_graph();
  double compute_propensity(int);
  void free_arrays();
};

#endif
