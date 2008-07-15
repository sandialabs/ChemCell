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

#ifndef FIX_DNA_TOGGLE_H
#define FIX_DNA_TOGGLE_H

#include "fix.h"

class FixDNAToggle : public Fix {
 public:
  FixDNAToggle(int, char **);
  ~FixDNAToggle();
  int setmask();
  void init();
  void initial();

 private:
  int nevery;            // how often to adjust rates

  int nlist;             // number of rates to adjust
  int *list;             // reaction indices of rates to adjust
  double *rate_initial;  // initial unadjusted rates
  double half;           // concentration at which rate is cut in half
  double volscale;       // scale factor on species concentration for stoch
  int ispecies;          // species index which affects rate dynamically

  double *rate;
  char *species;
  char **reactions;
};

#endif