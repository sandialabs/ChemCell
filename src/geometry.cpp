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

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "geometry.h"

enum {NONE,OUTSIDE,INSIDE,BOTH};

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ----------------------------------------------------------------------
   detect intersection between a line segment and a triangle
   intersection is defined as any line segment pt (including end pts)
     in common with any triangle pt (interior, edge, vertex)
   one exception is if both line end pts are in plane of triangle
     then is not an intersection
   v0,v1,v2 = 3 vertices of triangle
   norm = unit vector normal to triangle plane
     pointing OUTSIDE via right-hand rule
   start,stop = end points of directed line segment
   return TRUE if there is an intersection, else FALSE
   if TRUE also return:
     point = pt of intersection
     param = intersection pt is this fraction along segment (0-1 inclusive)
     side = OUTSIDE or INSIDE (enum value)
------------------------------------------------------------------------- */

bool tri_line_intersect(double *v0, double *v1, double *v2, double *norm,
			double *start, double *stop,
			double *point, double &param, int &side)
{
  double vec[3],start2stop[3],edge[3],pvec[3],xproduct[3];

  // if start,stop are on same side of triangle, no intersection
  // if start,stop are both in plane of triangle, no intersection

  subtract(v0,start,vec);
  double dotstart = dot(norm,vec);
  subtract(v0,stop,vec);
  double dotstop = dot(norm,vec);

  if (dotstart < 0.0 && dotstop < 0.0) return false;
  if (dotstart > 0.0 && dotstop > 0.0) return false;
  if (dotstart == 0.0 && dotstop == 0.0) return false;

  // param = parametric distance from start to stop
  //   at which tri plane is intersected
  // param = must be 0.0 to 1.0 inclusive, else no intersection

  subtract(start,v0,vec);
  subtract(start,stop,start2stop);
  param = dot(norm,vec) / dot(norm,start2stop);
  if (param < 0.0 || param > 1.0) return false;

  // point = intersection pt with plane of triangle

  point[0] = start[0] + param * start2stop[0];
  point[1] = start[1] + param * start2stop[1];
  point[2] = start[2] + param * start2stop[2];

  // test if intersection pt is inside triangle
  // edge = edge vector of triangle
  // pvec = vector from vertex to intersection point
  // xproduct = cross prodcut of edge with pvec
  // if dot product of xproduct with norm < 0.0 for any of 3 edges,
  //   intersection point is outside tri

  subtract(v0,v1,edge);
  subtract(v0,point,pvec);
  cross(edge,pvec,xproduct);
  if (dot(xproduct,norm) < 0.0) return false;

  subtract(v1,v2,edge);
  subtract(v1,point,pvec);
  cross(edge,pvec,xproduct);
  if (dot(xproduct,norm) < 0.0) return false;

  subtract(v2,v0,edge);
  subtract(v2,point,pvec);
  cross(edge,pvec,xproduct);
  if (dot(xproduct,norm) < 0.0) return false;

  // there is a valid intersection inside the triangle
  // set side to INSIDE or OUTSIDE
  // if start point is inside or outside then is INSIDE or OUTSIDE
  // if particle started on surface, then is opposite of stop point

  if (dotstart < 0.0) side = INSIDE;
  else if (dotstart > 0.0) side = OUTSIDE;
  else if (dotstop > 0.0) side = INSIDE;
  else if (dotstop < 0.0) side = OUTSIDE;
  return true;
}

// debug version

bool tri_line_intersect2(double *v0, double *v1, double *v2, double *norm,
			double *start, double *stop,
			double *point, double &param, int &side)
{
  double vec[3],start2stop[3],edge[3],pvec[3],xproduct[3];

  // if start,stop are on same side of triangle, no intersection
  // if start,stop are both in plane of triangle, no intersection

  subtract(v0,start,vec);
  double dotstart = dot(norm,vec);
  subtract(v0,stop,vec);
  double dotstop = dot(norm,vec);

  printf("AAA %g %g\n",dotstart,dotstop);
  printf("NORM %g %g %g\n",norm[0],norm[1],norm[2]);

  if (dotstart < 0.0 && dotstop < 0.0) return false;
  if (dotstart > 0.0 && dotstop > 0.0) return false;
  if (dotstart == 0.0 && dotstop == 0.0) return false;

  // param = parametric distance from start to stop
  //   at which tri plane is intersected
  // param = must be 0.0 to 1.0 inclusive, else no intersection

  subtract(start,v0,vec);
  subtract(start,stop,start2stop);
  param = dot(norm,vec) / dot(norm,start2stop);
  if (param < 0.0 || param > 1.0) return false;

  printf("PARAM %g\n",param);

  // point = intersection pt with plane of triangle

  point[0] = start[0] + param * start2stop[0];
  point[1] = start[1] + param * start2stop[1];
  point[2] = start[2] + param * start2stop[2];

  printf("INT POINT %g %g %g\n",point[0],point[1],point[2]);

  // test if intersection pt is inside triangle
  // edge = edge vector of triangle
  // pvec = vector from vertex to intersection point
  // xproduct = cross prodcut of edge with pvec
  // if dot product of xproduct with norm < 0.0 for any of 3 edges,
  //   intersection point is outside tri

  subtract(v0,v1,edge);
  subtract(v0,point,pvec);
  cross(edge,pvec,xproduct);
  if (dot(xproduct,norm) < 0.0) return false;

  printf("Pass 1\n");

  subtract(v1,v2,edge);
  subtract(v1,point,pvec);
  cross(edge,pvec,xproduct);
  if (dot(xproduct,norm) < 0.0) return false;

  printf("Pass 2\n");

  subtract(v2,v0,edge);
  subtract(v2,point,pvec);
  cross(edge,pvec,xproduct);
  if (dot(xproduct,norm) < 0.0) return false;

  printf("Pass 3\n");

  // there is a valid intersection inside the triangle
  // set side to INSIDE or OUTSIDE
  // if start point is inside or outside then is INSIDE or OUTSIDE
  // if particle started on surface, then is opposite of stop point

  if (dotstart < 0.0) side = INSIDE;
  else if (dotstart > 0.0) side = OUTSIDE;
  else if (dotstop > 0.0) side = INSIDE;
  else if (dotstop < 0.0) side = OUTSIDE;
  return true;
}

/* ----------------------------------------------------------------------
   compute intersection of a triangle with a hex cell
   intersection is defined as
     any triangle pt (interior, edge, vertex) in common with
     any hex pt (interior, face, edge, vertex)
   v0,v1,v2 and norm = 3 vertices of triangle and unit normal vec
   lo,hi = opposite corner pts of hex
   return 1 if intersection, else 0
------------------------------------------------------------------------- */

int tri_hex_intersect(double *v0, double *v1, double *v2, double *norm,
		      double *lo, double *hi)
{
  double xlo,ylo,zlo,xhi,yhi,zhi,sum;
  double b[3],e[3],h0[3],h1[3],h2[3],h3[3],n[3],point[3];
  double param;
  int side;

  xlo = lo[0];
  ylo = lo[1];
  zlo = lo[2];
  xhi = hi[0];
  yhi = hi[1];
  zhi = hi[2];

  // if all 8 hex pts are on same side of tri plane, no intersection
  // side is determined by dot product = -1,0,1

  sum = dotside(v0,norm,xlo,ylo,zlo);
  sum += dotside(v0,norm,xhi,ylo,zlo);
  sum += dotside(v0,norm,xlo,yhi,zlo);
  sum += dotside(v0,norm,xhi,yhi,zlo);
  sum += dotside(v0,norm,xlo,ylo,zhi);
  sum += dotside(v0,norm,xhi,ylo,zhi);
  sum += dotside(v0,norm,xlo,yhi,zhi);
  sum += dotside(v0,norm,xhi,yhi,zhi);
  
  if (sum == 8 || sum == -8) return 0;
	
  // if any of 3 tri vertices are inside hex, intersection
  // use <= and >= so touching hex surface is same as inside it

  if (v0[0] >= xlo && v0[0] <= xhi && v0[1] >= ylo && v0[1] <= yhi &&
      v0[2] >= zlo && v0[2] <= zhi) return 1;

  if (v1[0] >= xlo && v1[0] <= xhi && v1[1] >= ylo && v1[1] <= yhi &&
      v1[2] >= zlo && v1[2] <= zhi) return 1;

  if (v2[0] >= xlo && v2[0] <= xhi && v2[1] >= ylo && v2[1] <= yhi &&
      v2[2] >= zlo && v2[2] <= zhi) return 1;

  // test 12 hex edges for intersection with tri
  // b,e = begin/end of hex edge line segment

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = ylo;   e[2] = zlo;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zlo;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zhi;
  e[0] = xhi;   e[1] = ylo;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zhi;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xlo;   e[1] = yhi;   e[2] = zlo;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zlo;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zhi;
  e[0] = xlo;   e[1] = yhi;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zhi;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xlo;   e[1] = ylo;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = ylo;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zlo;
  e[0] = xlo;   e[1] = yhi;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xhi;   b[1] = yhi;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  // test 3 tri edges for intersection with 6 faces of hex
  // h0,h1,h2,h3 = 4 corner pts of hex face
  // n = normal to xyz faces, depends on vertex ordering
  // each face is treated as 2 triangles -> 6 tests per face
  
  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xlo;  h1[1] = yhi;  h1[2] = zlo;
  h2[0] = xlo;  h2[1] = yhi;  h2[2] = zhi;
  h3[0] = xlo;  h3[1] = ylo;  h3[2] = zhi;
  n[0]  = 1.0;  n[1]  = 0.0;  n[2]  = 0.0;
  
  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  h0[0] = h1[0] = h2[0] = h3[0] = xhi;

  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xhi;  h1[1] = ylo;  h1[2] = zlo;
  h2[0] = xhi;  h2[1] = ylo;  h2[2] = zhi;
  h3[0] = xlo;  h3[1] = ylo;  h3[2] = zhi;
  n[0]  = 0.0;  n[1]  = -1.0;  n[2]  = 0.0;
  
  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  h0[1] = h1[1] = h2[1] = h3[1] = yhi;

  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xhi;  h1[1] = ylo;  h1[2] = zlo;
  h2[0] = xhi;  h2[1] = yhi;  h2[2] = zlo;
  h3[0] = xlo;  h3[1] = yhi;  h3[2] = zlo;
  n[0]  = 0.0;  n[1]  = 0.0;  n[2]  = 1.0;
  
  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  h0[2] = h1[2] = h2[2] = h3[2] = zhi;
  
  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  return 0;
}

/* ----------------------------------------------------------------------
   compute distance from x to tri defined by vertices v0,v1,v2
   project x to plane of triangle
   if projected point is not in triangle, return false
   if it is, return true and projected distance
------------------------------------------------------------------------- */

bool point_tri_distance(double *x, double *v0, double *v1, double *v2,
			double *norm, double &distance)
{
  double vec[3],point[3],edge[3],pvec[3],xproduct[3];

  // point = intersection pt with plane of triangle

  subtract(x,v0,vec);
  double param = dot(norm,vec);
  point[0] = x[0] + param * norm[0];
  point[1] = x[1] + param * norm[1];
  point[2] = x[2] + param * norm[2];
    
  // test if intersection pt is inside triangle
  // edge = edge vector of triangle
  // pvec = vector from vertex to intersection point
  // xproduct = cross prodcut of edge with pvec
  // if dot product of xproduct with norm < 0.0 for any of 3 edges,
  //   intersection point is outside tri

  subtract(v0,v1,edge);
  subtract(v0,point,pvec);
  cross(edge,pvec,xproduct);
  if (dot(xproduct,norm) < 0.0) return false;

  subtract(v1,v2,edge);
  subtract(v1,point,pvec);
  cross(edge,pvec,xproduct);
  if (dot(xproduct,norm) < 0.0) return false;

  subtract(v2,v0,edge);
  subtract(v2,point,pvec);
  cross(edge,pvec,xproduct);
  if (dot(xproduct,norm) < 0.0) return false;

  // there is a valid intersection inside the triangle
  // compute distance from x to point

  double dx = point[0] - x[0];
  double dy = point[1] - x[1];
  double dz = point[2] - x[2];
  distance = sqrt(dx*dx + dy*dy + dz*dz);
  return true;
}

/* ----------------------------------------------------------------------
   move distance from xold in direction dir on a triangle surface
   triangle is defined by vertices v0,v1,v2 and unit-normal norm
   xold is assumed to be inside triangle and near its plane
   dir is assumed to be in plane of triangle
   if final pt is inside triangle:
     return true
     xnew = final pt
   if final pt is outside triangle:
     return false
     xnew = collision pt with boundary of triangle
     moved = distance moved
     flag = 1,2,3 for 3 edges if collision pt is on an edge
     flag = 4,5,6 for 3 vertices if collision pt is a vertex
   return false and flag = 0 if error occurs
------------------------------------------------------------------------- */

bool point_tri_move(double *v0, double *v1, double *v2, double *norm,
		    double *xold, double *dir, double distance,
		    double *xnew, double &moved, int &edge)
{
  // xnew = new pt via full move

  xnew[0] = xold[0] + distance * dir[0];
  xnew[1] = xold[1] + distance * dir[1];
  xnew[2] = xold[2] + distance * dir[2];

  // force xnew to be in plane of triangle

  push_pt_to_plane(xnew,v0,norm);

  // enorm123 = 3 vecs normal to each edge of tri
  // are in plane of tri, pointing towards center of tri
  // enorms are NOT unit vectors

  double delta[3];
  double enorm1[3],enorm2[3],enorm3[3];

  subtract(v0,v1,delta);
  cross(norm,delta,enorm1);
  subtract(v1,v2,delta);
  cross(norm,delta,enorm2);
  subtract(v2,v0,delta);
  cross(norm,delta,enorm3);

  // done if xnew is in triangle

  if (pt_in_triangle(xnew,v0,v1,v2,enorm1,enorm2,enorm3)) return true;

  // compute parametric distance from xold to xnew when collision with
  //   each of 3 edges occurs
  // no collision possible if dot product with edge normal is >= 0.0
  // error if no parametric distance is <= 1.0

  double param1,param2,param3;
  param1 = param2 = param3 = 2.0;
  subtract(xold,xnew,delta);

  if (dot(delta,enorm1) < 0.0)
    param1 = vec_plane_intersect(xold,xnew,v0,enorm1);
  if (dot(delta,enorm2) < 0.0)
    param2 = vec_plane_intersect(xold,xnew,v1,enorm2);
  if (dot(delta,enorm3) < 0.0)
    param3 = vec_plane_intersect(xold,xnew,v2,enorm3);

  double param = MIN(param1,MIN(param2,param3));
  if (param == 2.0) {
    edge = 0;
    return false;
  }

  // xnew = collision pt with nearest edge
  // moved = distance to xnew

  xnew[0] = xold[0] + param*delta[0];
  xnew[1] = xold[1] + param*delta[1];
  xnew[2] = xold[2] + param*delta[2];

  moved = param * distance;

  // if param is unique, set edge variable to matching edge (1,2,3)
  // if 2 params are the same, set edge variable to matching vertex (4,5,6)

  if (param == param1) edge = 1;
  else if (param == param2) edge = 2;
  else edge = 3;

  if (edge == 1 && param1 == param3) {
    edge = 4;
    xnew[0] = v0[0];
    xnew[1] = v0[1];
    xnew[2] = v0[2];
  } else if (edge == 2 && param1 == param2) {
    edge = 5;
    xnew[0] = v1[0];
    xnew[1] = v1[1];
    xnew[2] = v1[2];
  } else if (edge == 3 && param2 == param3) {
    edge = 6;
    xnew[0] = v2[0];
    xnew[1] = v2[1];
    xnew[2] = v2[2];
  }

  return false;
}

/* ----------------------------------------------------------------------
   v3 = v2 - v1
------------------------------------------------------------------------- */

void subtract(double *v1, double *v2, double *v3)
{
  v3[0] = v2[0] - v1[0];
  v3[1] = v2[1] - v1[1];
  v3[2] = v2[2] - v1[2];
}

/* ----------------------------------------------------------------------
   return dotproduct = v1 dot v2
------------------------------------------------------------------------- */

double dot(double *v1, double *v2)
{
  return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

/* ----------------------------------------------------------------------
   v3 = v1 x v2
------------------------------------------------------------------------- */

void cross(double *v1, double *v2, double *v3)
{
  v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
  v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
  v3[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

/* ----------------------------------------------------------------------
   determine which side of plane the point x,y,z is on
   plane is defined by vertex pt v and unit normal vec
   return -1,0,1 for below,on,above plane
------------------------------------------------------------------------- */

int dotside(double *v, double *norm, double x, double y, double z)
{
  double vec[3];
  vec[0] = x - v[0];
  vec[1] = y - v[1];
  vec[2] = z - v[2];

  double dotproduct = dot(norm,vec);
  if (dotproduct < 0.0) return -1;
  else if (dotproduct > 0.0) return 1;
  else return 0;
}

/* ----------------------------------------------------------------------
   force pt to be on the plane defined by x and unit-length norm
   return pt projected onto the plane
------------------------------------------------------------------------- */

void push_pt_to_plane(double *pt, double *x, double *norm)
{
  double diff[3],offset[3];
  subtract(x,pt,diff);
  double dotproduct = dot(diff,norm);

  offset[0] = dotproduct * norm[0];
  offset[1] = dotproduct * norm[1];
  offset[2] = dotproduct * norm[2];
  subtract(offset,pt,pt);
}

/* ----------------------------------------------------------------------
   determine if pt is in triangle
   triangle has 3 vertices v1,v2,v3 and
     3 non-unit edge normals norm1,norm2,norm3 pointing into triangle center
   return 0 if pt is outside triangle, 1 if inside
------------------------------------------------------------------------- */

int pt_in_triangle(double *pt, double *v1, double *v2, double *v3,
		   double *norm1, double *norm2, double *norm3)
{
  // if (pt - vertex) dotted into edge normal is negative, then outside tri
  
  double diff[3];
  subtract(v1,pt,diff);
  if (dot(diff,norm1) < 0.0) return 0;
  subtract(v2,pt,diff);
  if (dot(diff,norm2) < 0.0) return 0;
  subtract(v3,pt,diff);
  if (dot(diff,norm3) < 0.0) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   compute parametric distance from xold to xnew at which the line
     segment intersects the plane containing v with normal norm
   assume denominator dot product is non-zero
------------------------------------------------------------------------- */

double vec_plane_intersect(double *xold, double *xnew, double *v, double *norm)
{
  // solve for parametric distance t
  // variable x = xold + t (xnew - xold)
  // solve for x that satisfies norm dot (x-v) = 0

  double diff[3];
  subtract(v,xold,diff);
  double numerator = dot(norm,diff);
  subtract(xold,xnew,diff);
  double denominator = dot(norm,diff);
  
  double t = - numerator / denominator;
  return t;
}

/* ----------------------------------------------------------------------
   reflect the vector from x1 to x2 off of unit normal
   return new x2
   new vector from x1 to x2 is of same length but points in reflected dir
   A = vec from x1 to x2, N = unit normal vec
   vec from x1 to new x2 = A - 2(A dot N)N
   new x2 = x2 - 2(A dot N)N
------------------------------------------------------------------------- */

void reflect_normal(double *x1, double *x2, double *normal)
{
  double dot = (x2[0]-x1[0]) * normal[0] + (x2[1]-x1[1]) * normal[1] + 
    (x2[2]-x1[2]) * normal[2];
  x2[0] -= 2.0*dot * normal[0];
  x2[1] -= 2.0*dot * normal[1];
  x2[2] -= 2.0*dot * normal[2];
}

/* ----------------------------------------------------------------------
   reflect the unit vector dir off the edge from x1 to x2
   return new unit dir
   E = vec from x1 to x2 (not unit length)
   Esq = (length E)^2
   A = old unit dir
   new dir = 2(A dot E)E/Esq - A
   new dir should be unit-length, but normalize to make sure
------------------------------------------------------------------------- */

void reflect_edge(double *x1, double *x2, double *dir)
{
  double ex = x2[0] - x1[0];
  double ey = x2[1] - x1[1];
  double ez = x2[2] - x1[2];

  double aedot = dir[0]*ex + dir[1]*ey + dir[2]*ez;
  double esqinv = 1.0 / (ex*ex + ey*ey + ez*ez);

  dir[0] = 2.0*aedot*ex*esqinv - dir[0];
  dir[1] = 2.0*aedot*ey*esqinv - dir[1];
  dir[2] = 2.0*aedot*ez*esqinv - dir[2];

  double lensq = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
  double scale = 1.0/sqrt(lensq);
  dir[0] *= scale;
  dir[1] *= scale;
  dir[2] *= scale;
}

/* ----------------------------------------------------------------------
   bend the vector dir as it goes from a triangle with norm1 to one with norm2
   return new unit dir
   dir, norm1, norm2 are all unit vectors
   formula:  p' = p + sin(b) (a x p) + [1 - cos(b)] (a x (a x p)) 
     where p' = new dir, p = old dir, b = angle of rotation
     a = axis of rotation = unit vec in edge direction norm1 x norm2
     if norm1 and norm2 are in same dir, just return p' = p
   new dir should be unit-length, but normalize to make sure
------------------------------------------------------------------------- */

void bend2d(double *norm1, double *norm2, double *dir)
{
  double a[3],ap[3],aap[3];

  // A = unit edge vector around which rotation of angle b takes place
  // A = N1 x N2
  // sin(b) = length of (N1 x N2)
  
  a[0] = norm1[1]*norm2[2] - norm1[2]*norm2[1];
  a[1] = norm1[2]*norm2[0] - norm1[0]*norm2[2];
  a[2] = norm1[0]*norm2[1] - norm1[1]*norm2[0];

  double lensq = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
  if (lensq == 0.0) return;

  double sinb = sqrt(lensq);
  double scale = 1.0/sinb;
  a[0] *= scale;
  a[1] *= scale;
  a[2] *= scale;

  // cos(b) = N1 dot N2

  double cosb = norm1[0]*norm2[0] + norm1[1]*norm2[1] + norm1[2]*norm2[2];
  double onecosb = 1.0 - cosb;
  
  // AP = A x initial direction P
  // AAP = A x AP

  ap[0] = a[1]*dir[2] - a[2]*dir[1];
  ap[1] = a[2]*dir[0] - a[0]*dir[2];
  ap[2] = a[0]*dir[1] - a[1]*dir[0];

  aap[0] = a[1]*ap[2] - a[2]*ap[1];
  aap[1] = a[2]*ap[0] - a[0]*ap[2];
  aap[2] = a[0]*ap[1] - a[1]*ap[0];

  // new direction = initial direction + sin(b) AP + (1-cos(b)) AAP

  dir[0] += sinb*ap[0] + onecosb*aap[0];
  dir[1] += sinb*ap[1] + onecosb*aap[1];
  dir[2] += sinb*ap[2] + onecosb*aap[2];

  lensq = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
  scale = 1.0/sqrt(lensq);
  dir[0] *= scale;
  dir[1] *= scale;
  dir[2] *= scale;
}

/* ----------------------------------------------------------------------
   update pt x by moving distance in direction dir (unit vector)
   dir must be of finite length
   return new position as pt y
------------------------------------------------------------------------- */

void push_off(double *x, double *dir, double distance, double *y)
{
  y[0] = x[0] + distance*dir[0];
  y[1] = x[1] + distance*dir[1];
  y[2] = x[2] + distance*dir[2];
}
