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

#ifndef SURF_H
#define SURF_H

#include "system.h"

class Hash;
class Region;

class Surf : public System {
 public:
  struct OneTri {
    double normal[3];   // unit normal to tri plane, points outward
    int isurf;          // index of surf this tri belongs to
    int vert[3];        // index of 3 vertices
    int connect[3];     // index of tri connected to on each side (-1 if none)
  };

  int nsurf;            // # of surfaces (triangulated or region)
  char **name;          // surf-ID of each surf

  int nregion;          // # of regions
  Region **rlist;       // list of region ptrs
  int *region2surf;     // r2s[i] = which surf region I is
  int *surf2region;     // s2r[i] = which region surf I is (-1 = tri surf)

  int ntri;             // # of global triangles
  OneTri *tlist;        // list of triangles

  int nvert;            // # of global vertices for triangles
  double **vlist;       // list of vertices

  Surf();
  ~Surf();
  void init();
  int find(char *);                     // which surf matches ID string
  void read_triangles(int, char **);    // read triangulated surf and vertices
  void add_region(int, char **);        // add a region to rlist
  int memory_usage();

 private:
  int me,nprocs;

  int add_surf(char *);                 // add a surface to lists
  void tri_normal(int);                 // compute triangle unit normal
};

#endif
