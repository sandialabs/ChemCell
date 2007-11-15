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

#include "math.h"
#include "string.h"
#include "chem_gillespie.h"
#include "simulator.h"
#include "particle.h"
#include "react.h"
#include "random.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#define AVOGADRO 6.023e23

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

ChemGillespie::ChemGillespie() : Chem()
{
  char *str = "gillespie";
  int n = strlen(str) + 1;
  style = new char[n];
  strcpy(style,str);

  allocated = 0;
  volume = 0.0;
  nbinbin = nbinpair = ndist = noverlap = 0;
}

/* ---------------------------------------------------------------------- */

ChemGillespie::~ChemGillespie()
{
  if (allocated) free_arrays();
}

/* ----------------------------------------------------------------------
   init a simulation
------------------------------------------------------------------------- */

void ChemGillespie::init()
{
  // memory allocation

  if (allocated) free_arrays();
  allocated = 1;

  if (volume == 0.0) error->all("Must set volume for run style gillespie");

  factor_zero = AVOGADRO * volume;
  factor_dual = 1.0 / (AVOGADRO * volume);

  // grab needed values from other classes

  nreactions = react->nreactions;
  nreactant = react->nreactant;
  reactants = react->reactants;
  nproduct = react->nproduct;
  products = react->products;
  rate = react->rate;
  pcount = particle->pcount;

  ndepends = new int[nreactions];
  build_dependency_graph();

  // n = value such that 2^n >= nreactions

  int n = 0;
  while ((1 << n) < nreactions) n++;

  // create tree of length ntotal
  // init all leaves to 0.0

  int ntotal = 2*(1 << n) - 1;
  tree = new double[ntotal];
  offset = (1 << n) - 1;
  for (int i = offset; i < ntotal; i++) tree[i] = 0.0;

  // compute initial propensity for each reaction
  // sum the tree initially

  propensity = &tree[offset];
  for (int m = 0; m < nreactions; m++) propensity[m] = compute_propensity(m);
  sum();

  // zero reaction counts

  rcount = new int[nreactions];
  for (int m = 0; m < nreactions; m++) rcount[m] = 0;

  doneflag = 0;
}

/* ----------------------------------------------------------------------
   perform a single reaction and update data structures
------------------------------------------------------------------------- */

void ChemGillespie::reactions()
{
  int i,m,n;
  double r1,r2,delta,value;

  if (tree[0] > 0.0) {
    r1 = random->gillespie();
    delta = -1.0/tree[0] * log(r1);
  } else {
    doneflag = 1;
    return;
  }
  simulator->ctime += delta;

  r2 = random->gillespie();
  m = find(r2*tree[0]);

  rcount[m]++;
  for (i = 0; i < nreactant[m]; i++) pcount[reactants[m][i]]--;
  for (i = 0; i < nproduct[m]; i++) pcount[products[m][i]]++;
  
  for (i = 0; i < ndepends[m]; i++) {
    n = depends[m][i];
    value = compute_propensity(n);
    set(n,value);
  }
}

/* ----------------------------------------------------------------------
   set propensity[i] to value
   recompute sum tree for all its ancestors
------------------------------------------------------------------------- */

void ChemGillespie::set(int i, double value)
{
  int parent,sibling;

  propensity[i] = value;

  // i walks tree from leaf to root, summing children at each step
  // left child is odd index, right child is even index

  i += offset;
  while (i > 0) {
    if (i % 2) sibling = i + 1;
    else sibling = i - 1;
    parent = (i-1)/2;
    tree[parent] = tree[i] + tree[sibling];
    i = parent;
  }
}

/* ----------------------------------------------------------------------
   value = uniform RN from 0 to tree[0]
   return index (0 to M-1) of propensity bin it falls in
------------------------------------------------------------------------- */

int ChemGillespie::find(double value)
{
  int i,leftchild;

  // i walks tree from root to appropriate leaf
  // value is modified when right branch of tree is traversed

  i = 0;
  while (i < offset) {
    leftchild = 2*i + 1;
    if (value <= tree[leftchild]) i = leftchild;
    else {
      value -= tree[leftchild];
      i = leftchild + 1;
    }
  }
  return i - offset;
}

/* ----------------------------------------------------------------------
   sum entire tree, all nodes are computed
------------------------------------------------------------------------- */

void ChemGillespie::sum()
{
  int child1,child2;

  for (int parent = offset-1; parent >= 0; parent--) {
    child1 = 2*parent + 1;
    child2 = 2*parent + 2;
    tree[parent] = tree[child1] + tree[child2];
  }
}

/* ----------------------------------------------------------------------
   build dependency graph for entire set of reactions
   reaction N depends on M if a reactant of N is a reactant or product of M
------------------------------------------------------------------------- */

void ChemGillespie::build_dependency_graph()
{
  int i,j,k,m,n,mspecies,nspecies;

  // count the dependencies in flag array:
  // loop over reactants & products of each reaction
  // for each reaction m, mspecies = its reactants and products
  // for each species, loop over reactants of all reactions n
  // if a match, then set flag[n] since n is in dependency list of m

  int *flag = new int[nreactions];

  for (m = 0; m < nreactions; m++) {
    for (n = 0; n < nreactions; n++) flag[n] = 0;

    for (i = 0; i < nreactant[m]; i++) {
      mspecies = reactants[m][i];
      for (n = 0; n < nreactions; n++) {
	for (j = 0; j < nreactant[n]; j++) {
	  nspecies = reactants[n][j];
	  if (mspecies == nspecies) flag[n] = 1;
	}
      }
    }

    for (i = 0; i < nproduct[m]; i++) {
      mspecies = products[m][i];
      for (n = 0; n < nreactions; n++) {
	for (j = 0; j < nreactant[n]; j++) {
	  nspecies = reactants[n][j];
	  if (mspecies == nspecies) flag[n] = 1;
	}
      }
    }

    ndepends[m] = 0;
    for (n = 0; n < nreactions; n++) if (flag[n]) ndepends[m]++;
  }

  delete [] flag;

  // allocate depends array, 2nd dim is max of ndepends[]

  int nmax = 0;
  for (m = 0; m < nreactions; m++)
    nmax = MAX(nmax,ndepends[m]);

  depends = memory->create_2d_int_array(nreactions,nmax,"gillespie:depends");

  // zero the dependencies

  for (m = 0; m < nreactions; m++) ndepends[m] = 0;

  // store the dependencies via same loops as before
  // k loop insures dependency was not already stored

  for (m = 0; m < nreactions; m++) ndepends[m] = 0;

  for (m = 0; m < nreactions; m++) {
    for (i = 0; i < nreactant[m]; i++) {
      mspecies = reactants[m][i];
      for (n = 0; n < nreactions; n++) {
	for (j = 0; j < nreactant[n]; j++) {
	  nspecies = reactants[n][j];
	  if (mspecies == nspecies) {
	    for (k = 0; k < ndepends[m]; k++)
	      if (n == depends[m][k]) break;
	    if (k == ndepends[m]) depends[m][ndepends[m]++] = n;
	  }
	}
      }
    }

    for (i = 0; i < nproduct[m]; i++) {
      mspecies = products[m][i];
      for (n = 0; n < nreactions; n++) {
	for (j = 0; j < nreactant[n]; j++) {
	  nspecies = reactants[n][j];
	  if (mspecies == nspecies) {
	    for (k = 0; k < ndepends[m]; k++)
	      if (n == depends[m][k]) break;
	    if (k == ndepends[m]) depends[m][ndepends[m]++] = n;
	  }
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute propensity of a single reaction:
   for zero reaction: propensity = (Avogadro*Volume) * rate
   for mono reaction: propensity = count * rate
   for dual reaction: propensity = count1 * count2 * rate / (Avogadro*Volume)
   for dual reaction: cut in half if reactants are same species
------------------------------------------------------------------------- */

double ChemGillespie::compute_propensity(int m)
{
  double p;
  if (nreactant[m] == 0) p = factor_zero * rate[m];
  else if (nreactant[m] == 1) p = pcount[reactants[m][0]] * rate[m];
  else {
    if (reactants[m][0] == reactants[m][1]) 
      p = 0.5 * factor_dual * pcount[reactants[m][0]] * 
	(pcount[reactants[m][1]] - 1) * rate[m];
    else
      p = factor_dual * pcount[reactants[m][0]] * 
	pcount[reactants[m][1]] * rate[m];
  }
  return p;
}

/* ----------------------------------------------------------------------
   free arrays used by Gillspie solver
------------------------------------------------------------------------- */

void ChemGillespie::free_arrays()
{
  delete [] ndepends;
  memory->destroy_2d_int_array(depends);
  delete [] tree;
  delete [] rcount;
}
