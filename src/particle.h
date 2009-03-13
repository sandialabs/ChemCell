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

#ifndef PARTICLE_H
#define PARTICLE_H

#include "system.h"
#include "grid.h"

class Particle : public System {
 public:
  struct OnePart {
    double x[3];              // coords of particle
    int species;              // species of particle
    int ibin;                 // which local grid bin particle is in
    int itri;                 // which triangle part is on, -1 if 3d
    int seed;                 // RNG seed for this particle
    int next;                 // next particle in bin, -1 if no more
    int flag;                 // flag for this particle
  };

  int nlocal;               // # of particles owned by this proc
  int nghost;               // # of ghost particles on this proc
  int ntotal;               // sum of owned and ghost particles on this proc
  int nghost_last;          // # of ghost particles before compacting to 0
  OnePart *plist;           // list of particles

  int nspecies;             // # of unique species
  int maxspecies;           // max # of species arrays can hold
  char **name;              // root name of each species (species-ID)
  double *diffusivity;      // diffusivity of each species
  int *dimension;           // 2/3 if species is a 2d/3d diffuser

  int nalias;               // total # of alias strings
  int maxalias;             // max # of alias strings array can hold
  char **alias;             // list of alias strings
  int *alias2name;          // which species name each alias maps to

  int *pcount;              // integer count of each species
  double *ccount;           // concentration of each species (molarity)
  double *scount;           // concentration of each species for stats output
  int size_restart;         // # of quantities per particle to save for restart

  Particle();
  ~Particle();
  void init();
  void read(int, char **);
  int add(int, double, double, double);
  int find(char *);
  int add_species(char *);
  void set_species(int, char **);
  void set_dimension(int, char **);
  void set_diffusion(int, char **);
  void set_count(int, char **);
  void compute_count(int);
  void migrate();
  void migrate2();
  void ghost_acquire();
  void link();
  void unlink(int);
  void compact();
  int match(char *, char *);
  int allmatch(char *, int *);
  int pack_restart(double *);
  void unpack_restart(int, double *);
  int memory_usage();

 private:
  int me;
  int maxpart;              // max # of particles plist can hold

  struct Migrate {          // particle data to copy or migrate a particle
    double x[3];
    int species;
    int seed;
    int ibin;
    int itri;
  };

  Migrate *buf1,*buf2,*buf3;
  int *proclist;
  int size1,size2,size3,sizeproc;

  void add_alias(char *, int);
  void unpack(int, Migrate *, int);
  void fill_pm(int, Grid::Migrate *, int *);
  void fill_pc(int, Grid::Migrate *, int *);
  void check_dimension(int);
};

#endif
