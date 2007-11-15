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

#ifndef FINISH_H
#define FINISH_H

#include "system.h"

class Finish : public System {
 public:
  Finish() {}
  ~Finish() {}
  void end();

 private:
  int me,nprocs;

  void memory_usage();
  void stats(int, double *, double *, double *, double *, int, int *);
};

#endif
