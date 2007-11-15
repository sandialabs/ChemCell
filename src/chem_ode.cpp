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
#include "chem_ode.h"
#include "simulator.h"
#include "particle.h"
#include "react.h"
#include "random.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#define AVOGADRO 6.023e23

/* ---------------------------------------------------------------------- */

ChemODE::ChemODE() : Chem()
{
  char *str = "ode";
  int n = strlen(str) + 1;
  style = new char[n];
  strcpy(style,str);

  allocated = 0;
  volume = 0.0;
  nbinbin = nbinpair = ndist = noverlap = 0;
}

/* ---------------------------------------------------------------------- */

ChemODE::~ChemODE()
{
  if (allocated) free_arrays();
}

/* ----------------------------------------------------------------------
   init a simulation
------------------------------------------------------------------------- */

void ChemODE::init()
{
  // memory allocation

  if (allocated) free_arrays();
  allocated = 1;

  // volume only needed for output as particle count instead of concentration

  if (volume == 0.0) error->all("Must set volume for run style ode");

  // grab needed values from other classes

  nreactions = react->nreactions;
  nreactant = react->nreactant;
  reactants = react->reactants;
  nproduct = react->nproduct;
  products = react->products;
  rate = react->rate;
  wreactant = react->wreactant;
  wproduct = react->wproduct;

  nspecies = particle->nspecies;
  ccopy = new double[nspecies];
  ccount = particle->ccount;

  dt = simulator->dt;

  // initialize reaction counters

  rcount = new int[nreactions];
  for (int ireact = 0; ireact < nreactions; ireact++) rcount[ireact] = 0;
}

/* ----------------------------------------------------------------------
   perform each reaction for a timestep and update data structures
   delta change in each species for a single reaction:
   for zero reaction: dt * rate
   for mono reaction: dt * concentration * rate
   for dual reaction: dt * conc1 * conc2
   for dual reaction: cut in half if reactants are same species
------------------------------------------------------------------------- */

void ChemODE::reactions()
{
  int i,j;
  double delta;

  // make copy of concentration to update from previous timestep consistently

  for (i = 0; i < nspecies; i++) ccopy[i] = ccount[i];

  // loop over reactions
  // delta = amount to change reactants and products by in dt
  // update all reactants and products by delta

  for (i = 0; i < nreactions; i++) {
    if (nreactant[i] == 0) delta = rate[i] * dt;
    else if (nreactant[i] == 1) delta = ccopy[reactants[i][0]] * rate[i] * dt;
    else {
      if (reactants[i][0] == reactants[i][1]) 
	delta = 0.5 * ccopy[reactants[i][0]] * 
	  ccopy[reactants[i][1]] * rate[i] * dt;
      else
	delta = ccopy[reactants[i][0]] * 
	  ccopy[reactants[i][1]] * rate[i] * dt;
    }

    for (j = 0; j < nreactant[i]; j++)
      ccount[reactants[i][j]] -= delta * wreactant[i][j];
    for (j = 0; j < nproduct[i]; j++)
      ccount[products[i][j]] += delta * wproduct[i][j];
  }
}

/* ----------------------------------------------------------------------
   free arrays used by ODE solver
------------------------------------------------------------------------- */

void ChemODE::free_arrays()
{
  delete [] ccopy;
  delete [] rcount;
}
