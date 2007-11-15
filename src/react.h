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

#ifndef REACT_H
#define REACT_H

#include "system.h"

class React : public System {
 public:
  int nreactions;                // # of user defined reactions
  char **name;                   // name of each reaction (react-ID)

  int *nreactant;                // nreactant[I] = # of reactants of reaction I
                                 //   always 1 or 2
  int **reactants;               // reactants[I][J] = particle species of Jth
				 //   reactant of reaction I

  int *nproduct;                 // nproduct[I] = # of products of reaction I
  int **products;                // products[I][J] = particle species of Jth
				 //   product of reaction I

  double *rate;                  // rate[I] = input rate for reaction I

  int **setwhich;                // setwhich[I][J] = requested location
                                 //   (def,1,2,1/2) of Jth of I
  int **setstyle;                // setstyle[I][J] = requested location style
                                 //   (def,at,near) of Jth of I
  int **setdir;                  // setdir[I][J] = requested location direction
                                 //   (def,in,out,in/out) of Jth of I
  double *setdist;               // setdist[I] = requested cutoff for react I
  double *setprob;               // setprob[I] = requested prob for react I

  double **wreactant,**wproduct; // volume weightings for non-spatial sims
                                 // w[I][J] = w of Jth product of reaction I

  int **npair;                   // npair[I][J] = how many pair reactions
                                 //   particle species I,J take part in
  int ***pairlist;               // pairlist[I][J][K] = list (1:K-1) of
                                 //   reactions between species I,J

  int *nmono;                    // nmono[I] = how many mono reactions
                                 //   particle species I takes part in
  int **monolist;                // monolist[I][K] = list (1:K-1) of
                                 //   mono reactions for species I

  int **locwhich;                // locwhich[I][J] = location (ONE,TWO,ONETWO)
                                 //   of Jth product of reaction I
  int **locstyle;                // locstyle[I][J] = location style (AT,NEAR)
                                 //   of Jth product of reaction I
  int **locdir;                  // locdir[I][J] = location direction from surf
                                 //   (INSIDE,OUTSIDE,INOUT) of Jth of I

  React();
  ~React();
  void add(int, char **);
  void delete_reaction(char *);
  void modify(int, char **);
  void init();
  int find(char *);
  int create_product(int, int, int, int);

 private:
  int allocated;

  void invalid(int);
  void free_arrays();
};

#endif
