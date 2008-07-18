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
  void cleanup();

 private:
  int nevery;            // how often to adjust rates
  double kon,koff;
  double ktranscription,kconstitutive;

  int ieqrna,ieqdna;
  int idna,ibind;
  char *dnaspecies,*bindspecies;
  char *eqrna,*eqdna;
  double rate_rna,rate_dna;

  double *rate;
  int *pcount;
  double *ccount;
};

#endif
