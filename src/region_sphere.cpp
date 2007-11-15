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
#include "region_sphere.h"
#include "geometry.h"
#include "error.h"

enum {NONE,OUTSIDE,INSIDE,BOTH};        // matches geometry.cpp

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

RegionSphere::RegionSphere(int narg, char **arg) : Region(narg, arg)
{
  if (narg != 5) error->all("Invalid region arguments");

  xc = atof(arg[1]);
  yc = atof(arg[2]);
  zc = atof(arg[3]);
  r = atof(arg[4]);
}

/* ----------------------------------------------------------------------
   bounding box of sphere
------------------------------------------------------------------------- */

void RegionSphere::bbox(double *lo, double *hi)
{
  lo[0] = xc - r;
  lo[1] = yc - r;
  lo[2] = zc - r;
  hi[0] = xc + r;
  hi[1] = yc + r;
  hi[2] = zc + r;
}

/* ----------------------------------------------------------------------
   check if point x is inside of sphere
   pt on sphere surface is outside
------------------------------------------------------------------------- */

int RegionSphere::inside(double *x)
{
  double dx = x[0] - xc;
  double dy = x[1] - yc;
  double dz = x[2] - zc;

  if (dx*dx + dy*dy + dz*dz < r*r) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if hex defined by lo/hi pts is intersected by sphere surface
------------------------------------------------------------------------- */

int RegionSphere::hex_intersect(double *lo, double *hi)
{
  // close = distance from sphere center to nearest pt in hex volume
  // far = distance from sphere center to farthest pt in hex volume

  double dxlo,dylo,dzlo,dxhi,dyhi,dzhi;

  if (xc < lo[0]) dxlo = lo[0] - xc;
  else if (xc > hi[0]) dxlo = xc - hi[0];
  else dxlo = 0.0;
  if (yc < lo[1]) dylo = lo[1] - yc;
  else if (yc > hi[1]) dylo = yc - hi[1];
  else dylo = 0.0;
  if (zc < lo[2]) dzlo = lo[2] - zc;
  else if (zc > hi[2]) dzlo = zc - hi[2];
  else dzlo = 0.0;

  dxhi = MAX(hi[0]-xc,xc-lo[0]);
  dyhi = MAX(hi[1]-yc,yc-lo[1]);
  dzhi = MAX(hi[2]-zc,zc-lo[2]);

  double closesq = dxlo*dxlo + dylo*dylo + dzlo*dzlo;
  double farsq = dxhi*dxhi + dyhi*dyhi + dzhi*dzhi;

  // if close is inside sph and far is outside sph, then yes intersect
  //   use <= and >= in case close or far is on sphere surf
  // if close/far are both inside or both outside sph, then no intersect

  if (closesq <= r*r && farsq >= r*r) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if line segment intersects sphere surface
   if so, return intersection point and normal at that point,
     parametric distance along line segment and side (in/out) it hits
------------------------------------------------------------------------- */

bool RegionSphere::line_intersect(double *start, double *stop, double *point,
				  double *normal, double &param, int &side)
{
  // intersection of line with sphere surface
  // ctr = sphere center, ss = stop - start, cs = start - ctr
  // pt on line = x = start + alpha * ss
  // sph surface is (x - ctr)^2 = r^2
  // combine 2 eqs and solve resulting quadratic for param
  // A param^2 + B param + C = 0
  // A = ss^2, B = 2 sc dot ss, C = sc^2 - r^2

  double ctr[3],ss[3],cs[3];
  ctr[0] = xc; ctr[1] = yc; ctr[2] = zc;

  subtract(ctr,start,cs);
  subtract(start,stop,ss);
  double A = dot(ss,ss);
  double B = 2.0 * dot(cs,ss);
  double C = dot(cs,cs) - r*r;

  // bsq4ac < 0 means line does not intersect sphere
  // param = smallest positive of 2 solutions
  // if both params <= 0 or smallest positive > 1, then no intersection
  // demand param > 0 so that start pt on sphere surface is not a collision
  // this allows a particle to hit sphere multiple times in 1 move

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
   compute unit normal to sphere at pt x
   assume x is on sphere surface
------------------------------------------------------------------------- */

void RegionSphere::compute_normal(double *x, double *normal)
{
  normal[0] = (x[0] - xc) / r;
  normal[1] = (x[1] - yc) / r;
  normal[2] = (x[2] - zc) / r;
}

/* ----------------------------------------------------------------------
   compute xnew on sphere due to move in dir of distance from xold
------------------------------------------------------------------------- */

void RegionSphere::move2d(double *xold, double *normal, double *dir,
			  double distance, double *xnew)
{
  xnew[0] = xold[0] + distance*dir[0];
  xnew[1] = xold[1] + distance*dir[1];
  xnew[2] = xold[2] + distance*dir[2];

  // force xnew to be on sphere surface

  double dx = xnew[0] - xc;
  double dy = xnew[1] - yc;
  double dz = xnew[2] - zc;
  double length = sqrt(dx*dx + dy*dy + dz*dz);

  xnew[0] = xc + dx/length;
  xnew[1] = yc + dy/length;
  xnew[2] = zc + dz/length;
}

/* ----------------------------------------------------------------------
   compute and return distance from sphere to point x
------------------------------------------------------------------------- */

double RegionSphere::distance(double *x)
{
  double dx = x[0] - xc;
  double dy = x[1] - yc;
  double dz = x[2] - zc;
  double rpoint = sqrt(dx*dx + dy*dy + dz*dz);
  double dist = fabs(rpoint - r);
  return dist;
}
