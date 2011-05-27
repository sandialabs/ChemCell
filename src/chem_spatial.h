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
  friend class Particle;

 public:
  ChemSpatial();
  ~ChemSpatial();
  void init();
  void reactions();
  void create();
  void dynamic();
  double maxbin();
  int memory_usage();

 private:
  int me;
  int allocated;
  int nlink;                // # of particles in linked list of react bins

  struct Migrate {          // particle data to copy or migrate a particle
    double x[3];            // particle coord, remapped for PBC
    int species;            // species
    int seed;               // seed
    int ibin;               // local grid bin of particle in receiver domain
    int itri;               // which triangle particle is on
    int flag;               // REACTANT or PRODUCT particle
  };

  Migrate *buf1,*buf2,*buf3;
  int *proclist;
  int size1,size2,size3,sizeproc;

  int *nreactant;               // local copies of react ptrs
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

  struct ReactBin {              // reaction bins for finding reaction partners
    int id;                      // global ID of bin for owned, -1 for ghost
    int ghost;                   // 0 = owned, 1 = upwind ghost, 2 = far ghost
    int igridbin;                // index of local grid bin this rbin is inside
    int nparts;                  // # of particles in reaction bin
    int first;                   // index of 1st particle in bin (-1 if none)
    int next;                    // next bin of same color (-1 if last)
    int flag;                    // flag for this bin
  };

  ReactBin *rblist;              // list of reaction bins (3d array underneath)
                                 // these bins tile all local grid bins
                                 // including full tiling of ghost grid bins

  // global reaction bin counts not including ghosts
  // global IDs numbered from 0 to gbins-1 (x varies fast, y middle, z slow)
  // numbered in each dim from 0 to gbinxyz-1 inclusive

  int grbins;                    // global reaction bin count
  int grbinx,grbiny,grbinz;      // global reaction bin count in each dir
  int nperx,npery,nperz;         // reaction bins per grid bin

  // local reaction bin counts including ghosts
  // local IDs numbered from 0 to nbins-1 (x varies fast, y middle, z slow)
  // numbered in each dim from 0 to nbinxyz-1 inclusive
  // first and last nper in each dim are ghost since rbins tile grid bins

  int nrbins;                    // local rbin count (with ghosts)
  int nrbinx,nrbiny,nrbinz;      // local rbin count in each dir (with ghosts)
  int maxrbin;                   // max size of rblist

  // offset between local and global reaction bins
  // global = local + offset

  int xoffset,yoffset,zoffset;

  int **bbounds;    // bbounds[m][] = first/last global reaction bin in
                    //                Mth local grid bin
                    // 0/1 = xlo/xhi, 2/3 = ylo/yhi, 4/5 = zlo/zhi
  int maxnbin;      // max size of 1st dim of bbounds

  double xbinsize,ybinsize,zbinsize;   // reaction bin size in each dir
  double xbininv,ybininv,zbininv;      // inverse bin size in each dir
  double xorigin,yorigin,zorigin;      // xyz pt at which global rbins start

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

  void fill_rm(Particle::OnePart *, Grid::Migrate *, int, int *);
  void fill_rc(Particle::OnePart *, Grid::Migrate *, int, int *);
  void unpack(int, Migrate *);
  int match(Migrate *, int);

  void topology();
  void setup_stencil();
  void setup_colors();
  int whichcolor(int);
  void link();
  void unlink();
  int whichlocal(int, double *);
  void sort();
  static int compare(const void *, const void *);

  void free_arrays();
};

#endif
