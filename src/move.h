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

#ifndef MOVE_H
#define MOVE_H

#include "system.h"

class Move : public System {
  friend class MoveTest;

 public:
  int nmove;                  // total # of moves performed on this proc
  int ncheck,nreflect,nnear,nstick,nfar,nthru;     // collision stats
  int checkflag;
  int debug_proc,debug_step,debug_index;

  struct Values {
    double in_prob[5];
    double out_prob[5];
    int in_species[5];
    int out_species[5];
    int check,checklo,checkhi;
    int plo,phi;
  };
  struct Part {
    Values *surf;
  };
  Part *part;                 // particle/surface interaction settings

  Move();
  ~Move();
  void set_style(int, char **);
  void set_permeable(int, char **);
  void set_check(int, char **);
  void grow_partsurf(int, int);
  void default_partsurf(int, int);
  void init();
  void motion();
  double maxbin(double);
  void check();

 private:
  int me;
  int allocated;
  int npartspec,npartsurf;
  int geomflag,sampleflag;
  int *exlist;

  double *scale;
  double *displace3d,*displace2d;
  double *maxmove;           // max move for a species in any dimension
  double maxradsq;           // max radius^2 of sampling sphere
  double maxradsq_one;       // max radius^2 of sampling sphere plus one

  typedef void (Move::*FnPtr3d)(int, int *, double *);
  FnPtr3d delta3d;
  void delta_square_uniform(int, int *, double *);
  void delta_cube_uniform(int, int *, double *);
  void delta_square_brownian(int, int *, double *);
  void delta_cube_brownian(int, int *, double *);
  void delta_circle_uniform(int, int *, double *);
  void delta_sphere_uniform(int, int *, double *);
  void delta_circle_brownian(int, int *, double *);
  void delta_sphere_brownian(int, int *, double *);

  typedef void (Move::*FnPtr2d)(int, int *, double *, double *, double *);
  FnPtr2d delta2d;
  void delta_2d_uniform(int, int *, double *, double *, double *);
  void delta_2d_brownian(int, int *, double *, double *, double *);

  double erfcc(double);
  void free_arrays();
};

#endif
