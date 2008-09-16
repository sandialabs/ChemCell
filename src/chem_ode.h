/* ----------------------------------------------------------------------
   ChemCell - Cell simulator for particle diffusion and reactions
   Steve Plimpton (sjplimp@sandia.gov)
   Alex Slepoy (alexander.slepoy@nnsa.doe.gov)
   Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/chemcell.html

   Copyright (2004) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level ChemCell directory.
------------------------------------------------------------------------- */

#ifndef CHEM_ODE_H
#define CHEM_ODE_H

#include "chem.h"

class ChemODE : public Chem {
 public:
  ChemODE();
  ~ChemODE();
  void init();
  void reactions();

 private:
  int nreactions;            // # of reactions
  int *nreactant;            // # of reactants for each reaction
  int **reactants;           // i,j = jth reactant of ith reaction
  int *nproduct;             // # of products for each reaction
  int **products;            // i,j = jth product of ith reaction
  double *rate;              // input rate for each reaction
  double **wreactant;        // volume weights for j products of ith reaction
  double **wproduct;         // volume weights for j reactants of ith reaction

  int nspecies;              // # of species

  double *ccount;            // species concentration from particle class
  double *ccopy;             // work copy
  
  int allocated;
  double dt;
  double factor_zero,factor_dual;

  void free_arrays();
};

#endif
