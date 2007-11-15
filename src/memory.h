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

#ifndef MEMORY_H
#define MEMORY_H

#include "system.h"

class Memory : public System {
 public:
  Memory() {}
  ~Memory() {}

  void *smalloc(int n, char *);
  void sfree(void *);
  void *srealloc(void *, int n, char *name);

  double *create_1d_double_array(int, int, char *);
  void destroy_1d_double_array(double *, int);
  
  double **create_2d_double_array(int, int, char *);
  void destroy_2d_double_array(double **);
  double **grow_2d_double_array(double **, int, int, char *);

  int **create_2d_int_array(int, int, char *);
  void destroy_2d_int_array(int **);
  int **grow_2d_int_array(int **, int, int, char *);

  double **create_2d_double_array(int, int, int, char *);
  void destroy_2d_double_array(double **, int);

  double ***create_3d_double_array(int, int, int, char *);
  void destroy_3d_double_array(double ***);
  double ***grow_3d_double_array(double ***, int, int, int, char *);

  double ***create_3d_double_array(int, int, int, int, char *);
  void destroy_3d_double_array(double ***, int);

  double ***create_3d_double_array(int, int, int, int, int, int, char *);
  void destroy_3d_double_array(double ***, int, int, int);

  int ***create_3d_int_array(int, int, int, char *);
  void destroy_3d_int_array(int ***);
};

#endif
