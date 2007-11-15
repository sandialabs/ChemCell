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

#ifndef GEOMETRY_H
#define GEOMETRY_H

bool tri_line_intersect(double *v0, double *v1, double *v2, double *norm,
			double *start, double *stop,
			double *intersect, double &param, int &flagface);

// debug routine
bool tri_line_intersect2(double *v0, double *v1, double *v2, double *norm,
			double *start, double *stop,
			double *intersect, double &param, int &flagface);

int tri_hex_intersect(double *v0, double *v1, double *v2, double *norm,
		      double *lo, double *hi);

bool point_tri_distance(double *x, double *v0, double *v1, double *v2,
			double *norm, double &distance);

bool point_tri_move(double *v1, double *v2, double *v3, double *norm,
		    double *start, double *dir, double distance,
		    double *end, double &moved, int &index);

void subtract(double *, double *, double *);
double dot(double *, double *);
void cross(double *, double *, double *);
int dotside(double *, double *, double, double, double);
void push_pt_to_plane(double *,  double *, double *);
int pt_in_triangle(double *, double *, double *, double *,
		   double *, double *, double *);
double vec_plane_intersect(double *, double *, double *, double *);
void reflect_normal(double *, double *, double *);
void reflect_edge(double *, double *, double *);
void bend2d(double *, double *, double *);
void push_off(double *, double *, double, double *);

#endif
