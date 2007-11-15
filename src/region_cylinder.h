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

#ifndef REGION_CYLINDER_H
#define REGION_CYLINDER_H

#include "region.h"

class RegionCylinder : public Region {
 public:
  friend class Surf;

  RegionCylinder(int, char **);
  ~RegionCylinder() {}
  
  void bbox(double *, double *);
  int inside(double *);
  int hex_intersect(double *, double *);
  bool line_intersect(double *, double *, double *, double *, double &, int &);
  void compute_normal(double *, double *);
  void move2d(double *, double *, double *, double, double *);
  double distance(double *);

 private:
  int axis;
  double c1,c2;
  int c1dim,c2dim;
  double r;
};

#endif
