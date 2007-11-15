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

#ifndef DOMAIN_H
#define DOMAIN_H

#include "system.h"

class Domain : public System {
 public:
  double xlo,xhi,ylo,yhi,zlo,zhi;    // bounding box around entire domain
  double xsize,ysize,zsize;          // size of bounding box
  int xperiodic,yperiodic,zperiodic; // 1 if periodic, 0 if nonperiodic
  int setflag;                       // 0 if global domain not set, 1 if set

  Domain();
  ~Domain() {}
  void init();
  void global(int, char **);
  void boundary(int, char **);
  int inside(double *);
};

#endif
