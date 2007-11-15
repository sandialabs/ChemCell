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
#include "region_box.h"
#include "geometry.h"
#include "error.h"

enum {NONE,OUTSIDE,INSIDE,BOTH};        // matches geometry.cpp

/* ---------------------------------------------------------------------- */

RegionBox::RegionBox(int narg, char **arg) : Region(narg, arg)
{
  if (narg != 7) error->all("Invalid region arguments");

  xlo = atof(arg[1]);
  ylo = atof(arg[2]);
  zlo = atof(arg[3]);
  xhi = atof(arg[4]);
  yhi = atof(arg[5]);
  zhi = atof(arg[6]);

  if (xlo >= xhi || ylo >= yhi || zlo >= zhi)
    error->all("Illegal region arguments");
}

/* ----------------------------------------------------------------------
   bounding box of box
------------------------------------------------------------------------- */

void RegionBox::bbox(double *lo, double *hi)
{
  lo[0] = xlo;
  lo[1] = ylo;
  lo[2] = zlo;
  hi[0] = xhi;
  hi[1] = yhi;
  hi[2] = zhi;
}

/* ----------------------------------------------------------------------
   check if point x is inside of box
   pt on box surface is outside
------------------------------------------------------------------------- */

int RegionBox::inside(double *x)
{
  if (x[0] > xlo && x[0] < xhi && x[1] > ylo && x[1] < yhi &&
      x[2] > zlo && x[2] < zhi) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check if hex defined by lo/hi pts is intersected by box surface
------------------------------------------------------------------------- */

int RegionBox::hex_intersect(double *lo, double *hi)
{
  // if hex is outside box, no intersection

  if (xhi < lo[0]) return 0;
  if (xlo > hi[0]) return 0;
  if (yhi < lo[1]) return 0;
  if (ylo > hi[1]) return 0;
  if (zhi < lo[2]) return 0;
  if (zlo > hi[2]) return 0;

  // if hex is fully inside box, no intersection
  // else surface must intersect

  if (xlo < lo[0] && xhi > hi[0] && ylo < lo[1] && yhi > hi[1] && 
      zlo < lo[2] && zhi > hi[2]) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check if line segment intersects box surface
   if so, return intersection point and normal at that point,
     parametric distance along line segment and side (in/out) it hits
------------------------------------------------------------------------- */

bool RegionBox::line_intersect(double *start, double *stop, double *point,
			       double *normal, double &param, int &side)
{
  // ss = vector from start to stop
  
  double ss[3];
  subtract(start,stop,ss);

  // compute intersection point with each of 6 planes
  // check if point is inside face of box corresponding to that plane
  // store smallest parametric distance
  // demand alpha > 0 so that start pt on box surface is not a collision
  // this allows a particle to hit box multiple times in 1 move

  double alpha,ptx,pty,ptz;
  param = 2.0;

  if (start[0] != stop[0]) {
    alpha = (xlo - start[0]) / (stop[0] - start[0]);
    if (alpha > 0.0 and alpha < param) {
      pty = start[1] + alpha * ss[1];
      ptz = start[2] + alpha * ss[2];
      if (pty >= ylo && pty <= yhi && ptz >= zlo && ptz <= zhi)
	param = alpha;
    }
    alpha = (xhi - start[0]) / (stop[0] - start[0]);
    if (alpha > 0.0 and alpha < param) {
      pty = start[1] + alpha * ss[1];
      ptz = start[2] + alpha * ss[2];
      if (pty >= ylo && pty <= yhi && ptz >= zlo && ptz <= zhi)
	param = alpha;
    }
  }

  if (start[1] != stop[1]) {
    alpha = (ylo - start[1]) / (stop[1] - start[1]);
    if (alpha > 0.0 and alpha < param) {
      ptx = start[0] + alpha * ss[0];
      ptz = start[2] + alpha * ss[2];
      if (ptx >= xlo && ptx <= xhi && ptz >= zlo && ptz <= zhi)
	param = alpha;
    }
    alpha = (yhi - start[1]) / (stop[1] - start[1]);
    if (alpha > 0.0 and alpha < param) {
      ptx = start[0] + alpha * ss[0];
      ptz = start[2] + alpha * ss[2];
      if (ptx >= xlo && ptx <= xhi && ptz >= zlo && ptz <= zhi)
	param = alpha;
    }
  }

  if (start[2] != stop[2]) {
    alpha = (zlo - start[2]) / (stop[2] - start[2]);
    if (alpha > 0.0 and alpha < param) {
      ptx = start[0] + alpha * ss[0];
      pty = start[1] + alpha * ss[1];
      if (ptx >= xlo && ptx <= xhi && pty >= ylo && pty <= yhi)
	param = alpha;
    }
    alpha = (zhi - start[2]) / (stop[2] - start[2]);
    if (alpha > 0.0 and alpha < param) {
      ptx = start[0] + alpha * ss[0];
      pty = start[1] + alpha * ss[1];
      if (ptx >= xlo && ptx <= xhi && pty >= ylo && pty <= yhi)
	param = alpha;
    }
  }

  if (param <= 1.0) {
    point[0] = start[0] + param * ss[0];
    point[1] = start[1] + param * ss[1];
    point[2] = start[2] + param * ss[2];
    compute_normal(point,normal);

    // set side to INSIDE or OUTSIDE based on starting point

    if (start[0] > xlo && start[0] < xhi && start[1] > ylo && start[1] < yhi &&
	start[2] > zlo && start[2] < zhi) side = INSIDE;
    else side = OUTSIDE;

    return true;
  }
  return false;
}

/* ----------------------------------------------------------------------
   compute unit normal to box at pt x
   assume x is on box surface
------------------------------------------------------------------------- */

void RegionBox::compute_normal(double *x, double *normal)
{
  normal[0] = normal[1] = normal[2] = 0.0;

  if (x[0] == xlo) normal[0] = -1.0;
  else if (x[0] == xhi) normal[0] = 1.0;
  else if (x[1] == ylo) normal[1] = -1.0;
  else if (x[1] == yhi) normal[1] = 1.0;
  else if (x[2] == zlo) normal[2] = -1.0;
  else if (x[2] == zhi) normal[2] = 1.0;
}

/* ----------------------------------------------------------------------
   compute xnew on box surface due to move in dir of distance from xold
   have to wrap around box edges
------------------------------------------------------------------------- */

void RegionBox::move2d(double *xold, double *normal, double *dir,
		       double distance, double *xnew)
{
  xnew[0] = xold[0] + distance * dir[0];
  xnew[1] = xold[1] + distance * dir[1];
  xnew[2] = xold[2] + distance * dir[2];
}

/* ----------------------------------------------------------------------
   compute and return distance from box to point x
   compute distances to 6 planes and return smallest
------------------------------------------------------------------------- */

double RegionBox::distance(double *x)
{
  double dist = fabs(x[0] - xlo);
  if (fabs(x[0] - xhi) < dist) dist = fabs(x[0] - xhi);
  if (fabs(x[1] - ylo) < dist) dist = fabs(x[1] - ylo);
  if (fabs(x[1] - yhi) < dist) dist = fabs(x[1] - yhi);
  if (fabs(x[2] - zlo) < dist) dist = fabs(x[2] - zlo);
  if (fabs(x[2] - zhi) < dist) dist = fabs(x[2] - zhi);
  return dist;
}
