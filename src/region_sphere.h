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

#ifndef REGION_SPHERE_H
#define REGION_SPHERE_H

#include "region.h"

class RegionSphere : public Region {
 public:
  RegionSphere(int, char **);
  ~RegionSphere() {}
  
  void bbox(double *, double *);
  int inside(double *);
  int hex_intersect(double *, double *);
  bool line_intersect(double *, double *, double *, double *, double &, int &);
  void compute_normal(double *, double *);
  void move2d(double *, double *, double *, double, double *);
  double distance(double *);

 private:
  double xc,yc,zc;
  double r;
};

#endif
