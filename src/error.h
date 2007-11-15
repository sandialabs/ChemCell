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

#ifndef ERROR_H
#define ERROR_H

#include "stdio.h"
#include "system.h"

class Error : public System {
 public:
  Error() {}
  ~Error() {}

  void universe_all(char *);
  void universe_one(char *);

  void all(char *);
  void one(char *);
  void warning(char *);
};

#endif
