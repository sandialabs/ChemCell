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

#ifndef REGION_H
#define REGION_H

#include "system.h"

class Region : public System {
 public:
  char *style;

  Region(int, char **);
  virtual ~Region();

  virtual void bbox(double *, double *) = 0;
  virtual int inside(double *) = 0;
  virtual int hex_intersect(double *, double *) = 0;
  virtual bool line_intersect(double *, double *, double *,
			      double *, double &, int &) = 0;
  virtual void compute_normal(double *, double *) = 0;
  virtual void move2d(double *, double *, double *, double, double *) = 0;
  virtual double distance(double *) = 0;
};

#endif
