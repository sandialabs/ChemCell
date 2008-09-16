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

#ifndef CHEM_ODE_RK_H
#define CHEM_ODE_RK_H

#include "chem.h"

class ChemODERK : public Chem {
 public:
  ChemODERK();
  ~ChemODERK();
  void init();
  void reactions();
  void f(double, double *, double *);

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

  double abserr;
  int flag;
  int i_step;
  int n_step;
  double relerr;
  double t;
  double t_out;
  
  int allocated;
  double dt;
  double factor_zero,factor_dual;
  
  double d_epsilon(void);
  double d_max(double, double);
  double d_min(double, double);
  double d_sign(double);
  
  void fehl_d(int , 
	      double *, double, double, double *, double *, 
	      double *, double *, double *, double *, double *);
  int rkf45_d(int, double *, double *, double *, 
	      double, double *, double, int);

  void timestamp (void);
  void free_arrays();
};

#endif
