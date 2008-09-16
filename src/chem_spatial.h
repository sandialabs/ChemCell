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

#ifndef CHEM_SPATIAL_H
#define CHEM_SPATIAL_H

#include "chem.h"
#include "particle.h"
#include "grid.h"

class ChemSpatial : public Chem {
 public:
  ChemSpatial();
  ~ChemSpatial();
  void init();
  void reactions();
  void setup_stencil();

 private:
  int me;
  int allocated;

  struct Migrate {          // particle data to copy or migrate a particle
    double x[3];
    int species;
    int seed;
    int ibin;
    int itri;
    int flag;
  };

  Migrate *buf1,*buf2,*buf3;
  int *proclist;
  int size1,size2,size3,sizeproc;

  int *nreactant;
  int **reactants;
  int *nproduct;
  int **products;
  int *nmono;
  int **monolist;
  int **npair;
  int ***pairlist;

  int ncolor;                    // # of bin colors
  int *firstcolor;               // first local bin of each color
  int nstencil;                  // # of bin pairings in one color
  int *stencil1,*stencil2;       // list of bin pairings

  double **distsq;               // distsq[I][J] = reaction cutoff distance
                                 //   (squared) for all reacts between sp I,J

  double ***pairprob;            // pairprob[I][J][K] = cummulative probability
                                 //   (1:K-1) of reactions between species I,J
  double **pairprobsum;          // pairprobsum[I][J] = summed probability
                                 //   for all reactions between species I,J
  double **monoprob;             // monoprob[I][K] = cummulative probability
                                 //   (1:K-1) of mono reactions for species I
  double *monoprobsum;           // monoprobsum[I] = summed probability
                                 //   for all mono reactions for species I

  void unpack(int, Migrate *);
  void fill_rm(Particle::OnePart *, Grid::Migrate *, int, int *);
  void fill_rc(Particle::OnePart *, Grid::Migrate *, int, int *);
  int match(Migrate *, int);

  void setup_colors();
  int whichcolor(int);
  void free_arrays();
};

#endif
