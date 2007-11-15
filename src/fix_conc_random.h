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

#ifndef FIX_CONC_RANDOM_H
#define FIX_CONC_RANDOM_H

#include "fix.h"

class FixConcRandom : public Fix {
 public:
  FixConcRandom(int, char **);
  ~FixConcRandom();
  int setmask();
  void init();
  void initial();

 private:
  int nevery;           // how often to enforce specification
  int ninput;           // # of input species
  int *inlist;          // indices of input species
  int noutput;          // # of output species
  int *outlist;         // indices of output species
  int *spflag;          // 0 if species is not an input, 1 if it is

  int itime;            // time index to use for setting species probs
  int ntimes;           // # of time course entries from file
  double *times;        // time stamp for each time course entry
  double **probs;       // prob values for each time course entry
                        // stored as cummulative prob for each noutput species
                        // hence plist[*][noutput-1] = 1.0
};

#endif
