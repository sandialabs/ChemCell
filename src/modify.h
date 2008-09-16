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

#ifndef MODIFY_H
#define MODIFY_H

#include "stdio.h"
#include "system.h"

class Fix;

class Modify : public System {
 public:
  int nfix;
  int maxfix;
  int n_initial,n_final,n_cleanup;
  Fix **fix;                 // list of fixes

  Modify();
  ~Modify();
  void init();
  void initial();
  void final();
  void cleanup();

  void add_fix(int, char **);
  void delete_fix(char *);

 private:
  int *fmask;                // bit mask of when fix is applied
  int *list_initial;         // lists of fixes to apply at different times
  int *list_final;
  int *list_cleanup;

  void list_init(int, int &, int *&);
};

#endif
