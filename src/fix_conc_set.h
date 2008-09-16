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

#ifndef FIX_CONC_SET_H
#define FIX_CONC_SET_H

#include "fix.h"

class FixConcSet : public Fix {
 public:
  FixConcSet(int, char **);
  ~FixConcSet();
  int setmask();
  void init();
  void initial();

 private:
  int nevery;           // how often to enforce specification
  int nsp;              // # of species
  int *list;            // indices of species

  int itime;            // time index to use for setting species probs
  int ntimes;           // # of time course entries from file
  double *times;        // time stamp for each time course entry
  double **values;      // concentration values for each time and species
};

#endif
