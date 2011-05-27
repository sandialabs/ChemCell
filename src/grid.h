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

#ifndef GRID_H
#define GRID_H

#include "stdio.h"
#include "system.h"

class Grid : public System {
 public:
  struct Index {        // linked list of index objects
    int index;          // object index
    Index *ptr;         // ptr to next object in list (NULL if last)
  };

  struct Migrate {      // linked list of migrate actions for parts in this bin
    int proc;           // proc to send to (if me, just copy)
    int ibin;           // local bin ID of particle in receiving location
    int rcheck;         // 1 if should check reaction bins, 0 if no
    int ilo,jlo,klo;    // 1 if recv bin is downwind from send bin in that dir
    int pbc;            // 0 if no PBC remapping, 1 if yes
    int xpbc,ypbc,zpbc; // -1,0,1 box lengths to add to particle coords
    Migrate *ptr;       // ptr to next object in list (NULL if last)
  };

  struct OneBin {       // grid bins for geometry and particle tracking
    int id;             // global ID of bin (0 to gbins-1)
                        // if global ghost and periodic, ID of PBC image
                        // if global ghost and non-periodic, id = -1
    int proc;           // owner of bin (0 to nprocs-1)
                        // if global ghost and periodic, owner of PBC image
                        // if global ghost and non-periodic, proc = -1
    int ghost;          // 0 if local owned, 1 if local ghost
    int nmigrate;       // # of migrate actions for this bin's particles
    Migrate *migrate;   // 1st migrate action
    int ntri;           // # of triangles in bin
    Index *tri;         // 1st triangle
    int nparts;         // # of particles in bin
    int first;          // index of 1st particle in bin (-1 if none)
  };

  OneBin *blist;        // list of bins (3d array underneath)

  // global grid bin counts not including ghosts
  // global IDs numbered from 0 to gbins-1 (x varies fast, y middle, z slow)
  // numbered in each dim from 0 to gbinxyz-1 inclusive

  int gbins;                    // global bin count
  int gbinx,gbiny,gbinz;        // global bin count in each dir
  int nperx,npery,nperz;        // reaction bin count per grid bin

  // local grid bin counts including ghosts
  // local IDs numbered from 0 to nbins-1 (x varies fast, y middle, z slow)
  // numbered in each dim from 0 to nbinxyz-1 inclusive (ghost = 0,nbinxyz-1)

  int nbins;                    // local bin count (with ghosts)
  int nbinx,nbiny,nbinz;        // local bin count in each dir (with ghosts)

  // offset between local and global bins
  // global = local + offset
  // local is 0 for leftmost ghost, global is 0 for leftmost actual bin
  // offset for leftmost proc is -1 since its ghost is global ghost

  int xoffset,yoffset,zoffset;

  double xbinsize,ybinsize,zbinsize;   // bin size in each dir
  double xbininv,ybininv,zbininv;      // inverse bin size in each dir
  double xbinhalf,ybinhalf,zbinhalf;   // half bin size
  double xorigin,yorigin,zorigin;      // xyz pt at which global bins start

  int setflag;              // 0 if bins are not setup, 1 if setup
  int decomp_style;         // current decomp is (0) 3d or (1) RCB
  int balance_style;        // what kind of balancing to perform
  int nbalance;             // balance every this many steps
  double threshhold;        // imbalance factor threshhold
  int next_balance;         // next timestep at which to perform balancing
  int ncount;               // # of times balancer has been invoked
  int weightflag;           // 0 if no particle weights in RCB, 1 if yes
  int balance_file;         // 0 if no file write, 1 if yes
  FILE *fpbal;              // file pointer to balance file

  Grid();
  ~Grid();
  void init();
  int check();
  void dynamic();
  void create(int, char **);
  int whichlocal(double *);
  int whichglobal(double *, int *, int *, int *);
  int global2local(int);
  int global2local(int, int, int);
  int local2global(int);
  void local2global(int, int *, int *, int *);
  void local2local(int, int *, int *, int *);
  int local2one(int, int, int);
  void coord_center(int, double *, double *, double *);
  void coord_center(int, int, int, double *, double *, double *);
  double coord_cut(int, int);
  void coord_corners(int, double *, double *, double *,
		     double *, double *, double *);
  void coord_corners(int, int, int, double *, double *);
  void set_balance(int, char **);
  void rebalance();
  int box_overlap(double *, double *,
		  int &, int &, int &, int &, int &, int &);
  void add_tri(int, int);
  void delete_tri(Index *);
  void delete_migrate(Migrate *);
  void link();
  int memory_usage();

 private:
  int me,nprocs;
  int *list;                   // big enough to hold all procs

  int nprocx,nprocy,nprocz;    // procs in each dim for 3d decomp
  int *x2proc,*y2proc,*z2proc; // map from global grid to proc as 3 1d decomps

  struct BinAsk {
    int requester;             // proc that made request
    int id;                    // global bin ID being requested
    int ibin;                  // local bin ID on the requester
    int flag;                  // style of request: surf, particle, both
  };
  struct SurfReply {
    int itri;                  // index of surf object (tri or region)
    int ibin;                  // local bin ID on the requester
  };
  struct PartReply {
    double x[3];               // particle coord
    int species;               // particle species
    int seed;                  // particle seed
    int itri;                  // triangle index (-1 if none)
    int ibin;                  // local bin ID on the requester
  };

  char *buf1,*buf2,*buf3,*buf4;
  int *proclist1,*proclist2;
  int size1,size2,size3,size4;
  int sizeproc1,sizeproc2;
  int nichunks,nmchunks;
  Index **ichunks;
  Migrate **mchunks;
  Index *ifreelist;
  Migrate *mfreelist;
  Migrate possible[8];
  int **bininfo;

  void fill_ba(int, int, int, int *);
  void fill_sr(Index *, int, int, int *);
  void fill_pr(int, int, int, int *);

  void global(int, char **);
  void topology();
  void decomp3d();
  void decomp1d(int, int, int, int *, int *);
  void proc3d(int, int, int, int, int, 
	      int *, int *, int *, int *, int *, int *);
  int owner(int, int, int);
  void write_balance(double, double, double, double, double, double);
  Grid::Index *index_request();
  void index_return(Index *);
  Grid::Migrate *migrate_request();
  void migrate_return(Migrate *);
  void migrate_possible(int, int, int, int, int, int, int);
};

#endif
