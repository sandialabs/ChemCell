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

#include "math.h"
#include "string.h"
#include "region_cylinder.h"
#include "domain.h"
#include "geometry.h"
#include "error.h"

enum {NONE,OUTSIDE,INSIDE,BOTH};        // matches geometry.cpp

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

RegionCylinder::RegionCylinder(int narg, char **arg) : Region(narg, arg)
{
  if (narg != 5) error->all("Invalid region arguments");

  if (strcmp(arg[1],"x") == 0) {
    axis = 0;
    c1dim = 1;
    c2dim = 2;
  } else if (strcmp(arg[1],"y") == 0) {
    axis = 1;
    c1dim = 0;
    c2dim = 2;
  } else if (strcmp(arg[1],"z") == 0) {
    axis = 2;
    c1dim = 0;
    c2dim = 1;
  } else error->all("Invalid region arguments");

  c1 = atof(arg[2]);
  c2 = atof(arg[3]);
  r = atof(arg[4]);
}

/* ----------------------------------------------------------------------
   bounding box of cylinder
------------------------------------------------------------------------- */

void RegionCylinder::bbox(double *lo, double *hi)
{
  if (axis == 0) {
    lo[0] = domain->xlo;
    hi[0] = domain->xhi;
    lo[1] = c1 - r;
    hi[1] = c1 + r;
    lo[2] = c2 - r;
    hi[2] = c2 + r;
  } else if (axis == 1) {
    lo[0] = c1 - r;
    hi[0] = c1 + r;
    lo[1] = domain->ylo;
    hi[1] = domain->yhi;
    lo[2] = c2 - r;
    hi[2] = c2 + r;
  } else if (axis == 2) {
    lo[0] = c1 - r;
    hi[0] = c1 + r;
    lo[1] = c2 - r;
    hi[1] = c2 + r;
    lo[2] = domain->zlo;
    hi[2] = domain->zhi;
  }
}

/* ----------------------------------------------------------------------
   check if point x is inside of cylinder
   pt on cylinder surface is outside
------------------------------------------------------------------------- */

int RegionCylinder::inside(double *x)
{
  double d1 = x[c1dim] - c1;
  double d2 = x[c2dim] - c2;
  if (d1*d1 + d2*d2 < r*r) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if hex defined by lo/hi pts is intersected by cylinder surface
------------------------------------------------------------------------- */

int RegionCylinder::hex_intersect(double *lo, double *hi)
{
  // close = distance from cylinder center to nearest pt in hex volume
  // far = distance from cylinder center to farthest pt in hex volume

  double d1lo,d2lo,d1hi,d2hi;

  if (c1 < lo[c1dim]) d1lo = lo[c1dim] - c1;
  else if (c1 > hi[c1dim]) d1lo = c1 - hi[0];
  else d1lo = 0.0;
  if (c1 < lo[c2dim]) d2lo = lo[c2dim] - c2;
  else if (c2 > hi[c2dim]) d2lo = c2 - hi[c2dim];
  else d2lo = 0.0;

  d1hi = MAX(hi[c1dim]-c1,c1-lo[c1dim]);
  d2hi = MAX(hi[c2dim]-c2,c2-lo[c2dim]);

  double closesq = d1lo*d1lo + d2lo*d2lo;
  double farsq = d1hi*d1hi + d2hi*d2hi;

  // if close is inside cyl and far is outside cyl, then yes intersect
  //   use <= and >= in case close or far is on cylinder surf
  // if close/far are both inside or both outside cyl, then no intersect

  if (closesq <= r*r && farsq >= r*r) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if line segment intersects cylinder surface
   if so, return intersection point and normal at that point,
     parametric distance along line segment and side (in/out) it hits
------------------------------------------------------------------------- */

bool RegionCylinder::line_intersect(double *start, double *stop, double *point,
				    double *normal, double &param, int &side)
{
  // intersection of line with cylinder surface
  // ctr = cylinder center, ss = stop - start, cs = start - ctr
  // pt on line = x = start + alpha * ss
  // cyl surface is (x - ctr)^2 = r^2
  // combine 2 eqs and solve resulting quadratic for param
  // A param^2 + B param + C = 0
  // A = ss^2, B = 2 sc dot ss, C = sc^2 - r^2

  double ctr[3],ss[3],cs[3];
  ctr[0] = c1; ctr[1] = c2; ctr[2] = 0.0;

  subtract(ctr,start,cs);
  subtract(start,stop,ss);
  cs[2] = ss[2] = 0.0;               // set 3rd dimension = 0.0

  double A = dot(ss,ss);
  double B = 2.0 * dot(cs,ss);
  double C = dot(cs,cs) - r*r;

  // bsq4ac < 0 means line does not intersect cylinder
  // param = smallest positive of 2 solutions
  // if both params <= 0 or smallest positive > 1, then no intersection
  // demand param > 0 so that start pt on cylinder surface is not a collision
  // this allows a particle to hit cylinder multiple times in 1 move

  double bsq4ac = B*B - 4.0*A*C;
  if (bsq4ac < 0.0) return false;
  double param1 = (-B - sqrt(bsq4ac)) / (2.0*A);
  double param2 = (-B + sqrt(bsq4ac)) / (2.0*A);
  if (param1 <= 0.0 && param2 <= 0.0) return false;
  if (param1 > 0.0 && param2 > 0.0) param = MIN(param1,param2);
  else if (param1 > 0.0) param = param1;
  else param = param2;
  if (param > 1.0) return false;

  point[0] = start[0] + param * ss[0];
  point[1] = start[1] + param * ss[1];
  point[2] = start[2] + param * ss[2];
  compute_normal(point,normal);

  // set side to INSIDE or OUTSIDE based on starting point

  double rstart = sqrt(dot(cs,cs));
  if (rstart < r) side = INSIDE;
  else side = OUTSIDE;

  return true;
}

/* ----------------------------------------------------------------------
   compute unit normal to cylinder at pt x
   assume x is on cylinder surface
------------------------------------------------------------------------- */

void RegionCylinder::compute_normal(double *x, double *normal)
{
  normal[axis] = 0.0;
  normal[c1dim] = (x[c1dim] - c1) / r;
  normal[c2dim] = (x[c2dim] - c2) / r;
}

/* ----------------------------------------------------------------------
   compute xnew on cylinder due to move in dir of distance from xold
------------------------------------------------------------------------- */

void RegionCylinder::move2d(double *xold, double *normal, double *dir,
			    double distance, double *xnew)
{
  xnew[0] = xold[0] + distance*dir[0];
  xnew[1] = xold[1] + distance*dir[1];
  xnew[2] = xold[2] + distance*dir[2];

  // force xnew to be on cylinder surface

  double d1 = xnew[c1dim] - c1;
  double d2 = xnew[c2dim] - c2;
  double length = sqrt(d1*d1 + d2*d2);

  xnew[c1dim] = c1 + d1/length;
  xnew[c2dim] = c2 + d2/length;
}

/* ----------------------------------------------------------------------
   compute and return distance from infinite cylinder to point x
------------------------------------------------------------------------- */

double RegionCylinder::distance(double *x)
{
  double d1 = x[c1dim] - c1;
  double d2 = x[c2dim] - c2;
  double rpoint = sqrt(d1*d1 + d2*d2);
  double dist = fabs(rpoint - r);
  return dist;
}
