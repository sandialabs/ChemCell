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

#include "math.h"
#include "region_plane.h"
#include "domain.h"
#include "geometry.h"
#include "error.h"

enum {NONE,OUTSIDE,INSIDE,BOTH};        // matches geometry.cpp

#define BIG 1.0e20
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

RegionPlane::RegionPlane(int narg, char **arg) : Region(narg, arg)
{
  if (narg != 7) error->all("Invalid region arguments");

  ctr[0] = atof(arg[1]);
  ctr[1] = atof(arg[2]);
  ctr[2] = atof(arg[3]);
  norm[0] = atof(arg[4]);
  norm[1] = atof(arg[5]);
  norm[2] = atof(arg[6]);

  double length = sqrt(dot(norm,norm));
  if (length == 0.0) error->all("Invalid region arguments");
  norm[0] /= length;
  norm[1] /= length;
  norm[2] /= length;
}

/* ----------------------------------------------------------------------
   bounding box of plane
   plane is infinite, but bbox is limited to interior of global domain
------------------------------------------------------------------------- */

void RegionPlane::bbox(double *lo, double *hi)
{
  lo[0] = lo[1] = lo[2] = BIG;
  hi[0] = hi[1] = hi[2] = -BIG;
  
  double xlo = domain->xlo;
  double ylo = domain->ylo;
  double zlo = domain->zlo;
  double xhi = domain->xhi;
  double yhi = domain->yhi;
  double zhi = domain->zhi;

  // intersect 12 segments of global domain box with plane
  // shrink bounding box each time there is an intersection

  double start[3],stop[3],point[3];
  double param;
  int side;

  // segments in x direction

  start[0] = xlo; start[1] = ylo; start[2] = zlo;
  stop[0]  = xhi; stop[1]  = ylo; stop[2]  = zlo;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  start[0] = xlo; start[1] = yhi; start[2] = zlo;
  stop[0]  = xhi; stop[1]  = yhi; stop[2]  = zlo;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  start[0] = xlo; start[1] = yhi; start[2] = zhi;
  stop[0]  = xhi; stop[1]  = yhi; stop[2]  = zhi;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  start[0] = xlo; start[1] = ylo; start[2] = zhi;
  stop[0]  = xhi; stop[1]  = ylo; stop[2]  = zhi;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  // segments in y direction

  start[0] = xlo; start[1] = ylo; start[2] = zlo;
  stop[0]  = xlo; stop[1]  = yhi; stop[2]  = zlo;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  start[0] = xhi; start[1] = ylo; start[2] = zlo;
  stop[0]  = xhi; stop[1]  = yhi; stop[2]  = zlo;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  start[0] = xhi; start[1] = ylo; start[2] = zhi;
  stop[0]  = xhi; stop[1]  = yhi; stop[2]  = zhi;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  start[0] = xlo; start[1] = ylo; start[2] = zhi;
  stop[0]  = xlo; stop[1]  = yhi; stop[2]  = zhi;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  // segments in z direction

  start[0] = xlo; start[1] = ylo; start[2] = zlo;
  stop[0]  = xlo; stop[1]  = ylo; stop[2]  = zhi;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  start[0] = xhi; start[1] = ylo; start[2] = zlo;
  stop[0]  = xhi; stop[1]  = ylo; stop[2]  = zhi;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  start[0] = xhi; start[1] = yhi; start[2] = zhi;
  stop[0]  = xhi; stop[1]  = yhi; stop[2]  = zhi;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);

  start[0] = xlo; start[1] = yhi; start[2] = zhi;
  stop[0]  = xlo; stop[1]  = yhi; stop[2]  = zhi;
  if (line_intersect(start,stop,point,norm,param,side))
    shrink_bbox(lo,hi,point);
}

/* ----------------------------------------------------------------------
   check if point x is inside of plane
   plane normal points to outside
   pt on plane surface is outside
------------------------------------------------------------------------- */

int RegionPlane::inside(double *x)
{
  // if dot product of (x - ctr) and normal < 0 then x is inside
  // else outside

  double vec[3];
  subtract(ctr,x,vec);
  if (dot(vec,norm) < 0.0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if hex defined by lo/hi pts is intersected by plane
------------------------------------------------------------------------- */

int RegionPlane::hex_intersect(double *lo, double *hi)
{
  // if all 8 hex pts are on same side of plane, no intersection
  // side is determined by dot product = -1,0,1

  int sum = dotside(ctr,norm,lo[0],lo[1],lo[2]);
  sum += dotside(ctr,norm,hi[0],lo[1],lo[2]);
  sum += dotside(ctr,norm,lo[0],hi[1],lo[2]);
  sum += dotside(ctr,norm,hi[0],hi[1],lo[2]);
  sum += dotside(ctr,norm,lo[0],lo[1],hi[2]);
  sum += dotside(ctr,norm,hi[0],lo[1],hi[2]);
  sum += dotside(ctr,norm,lo[0],hi[1],hi[2]);
  sum += dotside(ctr,norm,hi[0],hi[1],hi[2]);
  
  if (sum == 8 || sum == -8) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check if line segment intersects plane
   if so, return intersection point and normal at that point,
     parametric distance along line segment and side (in/out) it hits
------------------------------------------------------------------------- */

bool RegionPlane::line_intersect(double *start, double *stop, double *point,
				 double *normal, double &param, int &side)
{
  double vec[3],start2stop[3];

  // if start,stop are on same side of plane, no intersection
  // if start,stop are both in plane, no intersection

  subtract(ctr,start,vec);
  double dotstart = dot(norm,vec);
  subtract(ctr,stop,vec);
  double dotstop = dot(norm,vec);

  if (dotstart < 0.0 && dotstop < 0.0) return false;
  if (dotstart > 0.0 && dotstop > 0.0) return false;
  if (dotstart == 0.0 && dotstop == 0.0) return false;

  // param = parametric distance from start to stop
  //   at which plane is intersected
  // param = must be > 0 and <= 1, else no intersection
  // demand param > 0 so that start pt on plane surface is not a collision

  subtract(start,ctr,vec);
  subtract(start,stop,start2stop);
  param = dot(norm,vec) / dot(norm,start2stop);
  if (param <= 0.0 || param > 1.0) return false;

  // point = intersection pt with plane of triangle

  point[0] = start[0] + param * start2stop[0];
  point[1] = start[1] + param * start2stop[1];
  point[2] = start[2] + param * start2stop[2];

  // set side to INSIDE or OUTSIDE based on starting point

  if (dotstart < 0.0) side = INSIDE;
  else side = OUTSIDE;

  normal[0] = norm[0];
  normal[1] = norm[1];
  normal[2] = norm[2];

  return true;
}

/* ----------------------------------------------------------------------
   return unit normal of plane
------------------------------------------------------------------------- */

void RegionPlane::compute_normal(double *x, double *normal)
{
  normal[0] = norm[0];
  normal[1] = norm[1];
  normal[2] = norm[2];
}

/* ----------------------------------------------------------------------
   compute xnew on plane due to move in dir of distance from xold
------------------------------------------------------------------------- */

void RegionPlane::move2d(double *xold, double *normal, double *dir,
			 double distance, double *xnew)
{
  xnew[0] = xold[0] + distance*dir[0];
  xnew[1] = xold[1] + distance*dir[1];
  xnew[2] = xold[2] + distance*dir[2];

  // force xnew to be in plane

  push_pt_to_plane(xnew,ctr,norm);
}

/* ----------------------------------------------------------------------
   union of 2 bounding boxes to smaller one
------------------------------------------------------------------------- */

void RegionPlane::shrink_bbox(double *lo, double *hi, double *x)
{
  lo[0] = MIN(lo[0],x[0]);
  lo[1] = MIN(lo[1],x[1]);
  lo[2] = MIN(lo[2],x[2]);
  hi[0] = MAX(hi[0],x[0]);
  hi[1] = MAX(hi[1],x[1]);
  hi[2] = MAX(hi[2],x[2]);
}

/* ----------------------------------------------------------------------
   compute and return distance from plane to point x
   xc = x - c
   xc dot norm = len(xc) * cos(theta)
   sin(theta) = distance / len(xc)
------------------------------------------------------------------------- */

double RegionPlane::distance(double *x)
{
  double xc[3];
  xc[0] = x[0] - ctr[0];
  xc[0] = x[1] - ctr[1];
  xc[0] = x[2] - ctr[2];
  double xclength = sqrt(dot(xc,xc));
  double xcdotnorm = dot(xc,norm);
  double costheta = xcdotnorm / xclength;
  double sintheta = sqrt(1.0 - costheta*costheta);
  double dist = sintheta * xclength;
  //printf("CCC %g %g: %g %g: %g\n",xclength,xcdotnorm,costheta,sintheta,dist);
  return dist;
}
