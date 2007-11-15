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

#ifndef FIX_H
#define FIX_H

#include "system.h"

class Fix : public System {
 public:
  char *id,*style;

  int INITIAL;                   // mask settings

  Fix(int, char **);
  virtual ~Fix();
  virtual int setmask() = 0;
  virtual void init() {}
  virtual void initial() {}
};

#endif
