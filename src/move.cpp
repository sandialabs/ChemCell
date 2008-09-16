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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "move.h"
#include "grid.h"
#include "particle.h"
#include "surf.h"
#include "region.h"
#include "domain.h"
#include "time.h"
#include "random.h"
#include "simulator.h"
#include "geometry.h"
#include "memory.h"
#include "error.h"

enum {CUBE,SPHERE,SQUARE,CIRCLE};
enum {UNIFORM,BROWNIAN};
enum {NONE,OUTSIDE,INSIDE,BOTH};        // matches geometry.cpp
enum {REFLECT,NEAR,STICK,FAR,THRU};

#define MAXITER 25
#define PI 3.1415926
#define NSLICE 8192
#define TOLERANCE 1.0e-6
#define EPSILON 1.0e-6
#define DELTA 1.0e-4
#define BIGINT 1 << 30

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// comment out if don't want DEBUG enabled

//#define DEBUG_MOVE

/* ---------------------------------------------------------------------- */

Move::Move()
{
  MPI_Comm_rank(world,&me);

  allocated = 0;

  geomflag = CUBE;
  sampleflag = BROWNIAN;

  scale = NULL;
  displace3d = NULL;
  displace2d = NULL;
  maxmove = NULL;

  exlist = NULL;

  npartspec = npartsurf = 0;
  part = NULL;

  debug_proc = debug_step = debug_index = -1;
}

/* ---------------------------------------------------------------------- */

Move::~Move()
{
  if (allocated) free_arrays();

  for (int i = 0; i < npartspec; i++)
    memory->sfree(part[i].surf);
  memory->sfree(part);
}

/* ----------------------------------------------------------------------
   set move style
------------------------------------------------------------------------- */

void Move::set_style(int narg, char **arg)
{
  for (int iarg = 0; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"cube") == 0) geomflag = CUBE;
    else if (strcmp(arg[iarg],"sphere") == 0) geomflag = SPHERE;
    else if (strcmp(arg[iarg],"square") == 0) geomflag = SQUARE;
    else if (strcmp(arg[iarg],"circle") == 0) geomflag = CIRCLE;
    else if (strcmp(arg[iarg],"uniform") == 0) sampleflag = UNIFORM;
    else if (strcmp(arg[iarg],"brownian") == 0) sampleflag = BROWNIAN;
    else error->all("Illegal move_style command");
  }
}

/* ----------------------------------------------------------------------
   set in/out checking for a species/surface pair
------------------------------------------------------------------------- */

void Move::set_check(int narg, char **arg)
{
  if (narg < 3) error->all("Illegal check command");

  // species check

  int ispecies = particle->find(arg[0]);
  if (ispecies == -1) {
    char str[128];
    sprintf(str,"Unknown species %s in check command",arg[0]);
    error->all(str);
  }
  if (particle->dimension[ispecies] == 2)
    error->all("Cannot check for a 2d species");

  // surface check

  int isurf = surf->find(arg[1]);
  if (isurf == -1) {
    char str[128];
    sprintf(str,"Unknown surface %s in check command",arg[1]);
    error->all(str);
  }

  // side check

  int side;
  if (strcmp(arg[2],"out") == 0) side = OUTSIDE;
  else if (strcmp(arg[2],"in") == 0) side = INSIDE;
  else if (strcmp(arg[2],"none") == 0) side = NONE;
  else error->all("Illegal check command");

  // check = which side or NONE
  // lo/hi = timestep values if args are set, else 0 to BIGINT

  part[ispecies].surf[isurf].check = side;
  if (side != NONE) {
    part[ispecies].surf[isurf].checklo = 0;
    part[ispecies].surf[isurf].checkhi = BIGINT;
    part[ispecies].surf[isurf].plo = -1;
    part[ispecies].surf[isurf].phi = BIGINT;
    if (narg >= 5) {
      part[ispecies].surf[isurf].checklo = atoi(arg[3]);
      part[ispecies].surf[isurf].checkhi = atoi(arg[4]);
    }
    if (narg == 7) {
      part[ispecies].surf[isurf].plo = atoi(arg[5]);
      part[ispecies].surf[isurf].phi = atoi(arg[6]);
    }
  }
}

/* ----------------------------------------------------------------------
   set permeability for a species/surface pair
------------------------------------------------------------------------- */

void Move::set_permeable(int narg, char **arg)
{
  if (narg < 3) error->all("Illegal permeable command");

  // species check

  int ispecies = particle->find(arg[0]);
  if (ispecies == -1) {
    char str[128];
    sprintf(str,"Unknown species %s in permeable command",arg[0]);
    error->all(str);
  }
  if (particle->dimension[ispecies] == 2)
    error->all("Cannot set permeability for a 2d species");

  // surface check

  int isurf = surf->find(arg[1]);
  if (isurf == -1) {
    char str[128];
    sprintf(str,"Unknown surface %s in permeable command",arg[1]);
    error->all(str);
  }

  // side check

  int side;
  if (strcmp(arg[2],"out") == 0) side = OUTSIDE;
  else if (strcmp(arg[2],"in") == 0) side = INSIDE;
  else if (strcmp(arg[2],"both") == 0) side = BOTH;
  else error->all("Illegal permeable command");

  // default values for probabilities and species changes

  double rprob,nprob,sprob,fprob,tprob;
  rprob = 1.0;
  nprob = sprob = fprob = tprob = 0.0;
  int rsp,nsp,ssp,fsp,tsp;
  rsp = nsp = ssp = fsp = tsp = ispecies;

  // loop over all keyword/value pairs

  int iarg = 3;
  while (iarg < narg) {
    if (narg < iarg+1) error->all("Illegal permeable command");
    if (strcmp(arg[iarg],"reflect") == 0) rprob = atof(arg[iarg+1]);
    else if (strcmp(arg[iarg],"near") == 0) nprob = atof(arg[iarg+1]);
    else if (strcmp(arg[iarg],"stick") == 0) sprob = atof(arg[iarg+1]);
    else if (strcmp(arg[iarg],"far") == 0) fprob = atof(arg[iarg+1]);
    else if (strcmp(arg[iarg],"thru") == 0) tprob = atof(arg[iarg+1]);
    else if (strcmp(arg[iarg],"rsp") == 0) {
      rsp = particle->find(arg[iarg+1]);
      if (rsp == -1) {
	char str[128];
	sprintf(str,"Unknown species %s in permeable command",arg[iarg+1]);
	error->all(str);
      }
    } else if (strcmp(arg[iarg],"nsp") == 0) {
      nsp = particle->find(arg[iarg+1]);
      if (nsp == -1) {
	char str[128];
	sprintf(str,"Unknown species %s in permeable command",arg[iarg+1]);
	error->all(str);
      }
    } else if (strcmp(arg[iarg],"ssp") == 0) {
      ssp = particle->find(arg[iarg+1]);
      if (ssp == -1) {
	char str[128];
	sprintf(str,"Unknown species %s in permeable command",arg[iarg+1]);
	error->all(str);
      }
    } else if (strcmp(arg[iarg],"fsp") == 0) {
      fsp = particle->find(arg[iarg+1]);
      if (fsp == -1) {
	char str[128];
	sprintf(str,"Unknown species %s in permeable command",arg[iarg+1]);
	error->all(str);
      }
    } else if (strcmp(arg[iarg],"tsp") == 0) {
      tsp = particle->find(arg[iarg+1]);
      if (tsp == -1) {
	char str[128];
	sprintf(str,"Unknown species %s in permeable command",arg[iarg+1]);
	error->all(str);
      }
    } else error->all("Illegal permeable command");
    iarg += 2;
  }

  // check that probabilities sum to 1.0
  // check that stick species is 2d if stick prob is non-zero

  double sum = rprob + nprob + sprob + fprob + tprob;
  if (sum < 1.0-EPSILON || sum > 1.0+EPSILON) 
    error->all("Permeable probabilities do not sum to 1.0");

  if (sprob > 0.0 && particle->dimension[ssp] == 3)
    error->all("Permeable stick species is not 2d");

  // store parameters in permeability data struct

  if (side == INSIDE || side == BOTH) {
    part[ispecies].surf[isurf].in_prob[REFLECT] = rprob;
    part[ispecies].surf[isurf].in_prob[NEAR] = rprob + nprob;
    part[ispecies].surf[isurf].in_prob[STICK] = rprob + nprob + sprob;
    part[ispecies].surf[isurf].in_prob[FAR] = rprob + nprob + sprob + fprob;
    part[ispecies].surf[isurf].in_prob[THRU] = 1.0;
    part[ispecies].surf[isurf].in_species[REFLECT] = rsp;
    part[ispecies].surf[isurf].in_species[NEAR] = nsp;
    part[ispecies].surf[isurf].in_species[STICK] = ssp;
    part[ispecies].surf[isurf].in_species[FAR] = fsp;
    part[ispecies].surf[isurf].in_species[THRU] = tsp;
  } else if (side == OUTSIDE || side == BOTH) {
    part[ispecies].surf[isurf].out_prob[REFLECT] = rprob;
    part[ispecies].surf[isurf].out_prob[NEAR] = rprob + nprob;
    part[ispecies].surf[isurf].out_prob[STICK] = rprob + nprob + sprob;
    part[ispecies].surf[isurf].out_prob[FAR] = rprob + nprob + sprob + fprob;
    part[ispecies].surf[isurf].out_prob[THRU] = 1.0;
    part[ispecies].surf[isurf].out_species[REFLECT] = rsp;
    part[ispecies].surf[isurf].out_species[NEAR] = nsp;
    part[ispecies].surf[isurf].out_species[STICK] = ssp;
    part[ispecies].surf[isurf].out_species[FAR] = fsp;
    part[ispecies].surf[isurf].out_species[THRU] = tsp;
  }
}

/* ----------------------------------------------------------------------
   grow the particle/surface data struct
   one of the 2 args will be larger than the current 2d array of Part
------------------------------------------------------------------------- */

void Move::grow_partsurf(int nspecies, int nsurf)
{
  if (nspecies > npartspec) {
    part = (Part *) memory->srealloc(part,nspecies*sizeof(Part),"move:part");
    for (int i = npartspec; i < nspecies; i++)
      part[i].surf = 
	(Values *) memory->smalloc(nsurf*sizeof(Values),"move:part");
    npartspec = nspecies;
  }
  if (nsurf > npartsurf) {
    for (int i = 0; i < nspecies; i++)
      part[i].surf = 
	(Values *) memory->srealloc(part[i].surf,nsurf*sizeof(Values),
				       "move:part");
    npartsurf = nsurf;
  }
}

/* ----------------------------------------------------------------------
   set defaults in the particle/surface data struct
   one of the 2 args will be -1, the other will be >= 0
------------------------------------------------------------------------- */

void Move::default_partsurf(int ispecies, int isurf)
{
  if (ispecies >= 0)
    for (int i = 0; i < npartsurf; i++) {
      part[ispecies].surf[i].in_prob[REFLECT] = 1.0;
      part[ispecies].surf[i].in_prob[NEAR] = 1.0;
      part[ispecies].surf[i].in_prob[STICK] = 1.0;
      part[ispecies].surf[i].in_prob[FAR] = 1.0;
      part[ispecies].surf[i].in_prob[THRU] = 1.0;
      part[ispecies].surf[i].in_species[REFLECT] = ispecies;
      part[ispecies].surf[i].in_species[NEAR] = ispecies;
      part[ispecies].surf[i].in_species[STICK] = ispecies;
      part[ispecies].surf[i].in_species[FAR] = ispecies;
      part[ispecies].surf[i].in_species[THRU] = ispecies;
      part[ispecies].surf[i].out_prob[REFLECT] = 1.0;
      part[ispecies].surf[i].out_prob[NEAR] = 1.0;
      part[ispecies].surf[i].out_prob[STICK] = 1.0;
      part[ispecies].surf[i].out_prob[FAR] = 1.0;
      part[ispecies].surf[i].out_prob[THRU] = 1.0;
      part[ispecies].surf[i].out_species[REFLECT] = ispecies;
      part[ispecies].surf[i].out_species[NEAR] = ispecies;
      part[ispecies].surf[i].out_species[STICK] = ispecies;
      part[ispecies].surf[i].out_species[FAR] = ispecies;
      part[ispecies].surf[i].out_species[THRU] = ispecies;

      part[ispecies].surf[i].check = NONE;
    }
  if (isurf >= 0)
    for (int i = 0; i < npartspec; i++) {
      part[i].surf[isurf].in_prob[REFLECT] = 1.0;
      part[i].surf[isurf].in_prob[NEAR] = 1.0;
      part[i].surf[isurf].in_prob[STICK] = 1.0;
      part[i].surf[isurf].in_prob[FAR] = 1.0;
      part[i].surf[isurf].in_prob[THRU] = 1.0;
      part[i].surf[isurf].in_species[REFLECT] = i;
      part[i].surf[isurf].in_species[NEAR] = i;
      part[i].surf[isurf].in_species[STICK] = i;
      part[i].surf[isurf].in_species[FAR] = i;
      part[i].surf[isurf].in_species[THRU] = i;
      part[i].surf[isurf].out_prob[REFLECT] = 1.0;
      part[i].surf[isurf].out_prob[NEAR] = 1.0;
      part[i].surf[isurf].out_prob[STICK] = 1.0;
      part[i].surf[isurf].out_prob[FAR] = 1.0;
      part[i].surf[isurf].out_prob[THRU] = 1.0;
      part[i].surf[isurf].out_species[REFLECT] = i;
      part[i].surf[isurf].out_species[NEAR] = i;
      part[i].surf[isurf].out_species[STICK] = i;
      part[i].surf[isurf].out_species[FAR] = i;
      part[i].surf[isurf].out_species[THRU] = i;

      part[i].surf[isurf].check = NONE;
    }
}

/* ----------------------------------------------------------------------
   setup for motion
------------------------------------------------------------------------- */

void Move::init()
{
  int ispecies,isurf,jspecies;
  double prob;

  if (allocated) free_arrays();
  allocated = 1;

  int nspecies = particle->nspecies;
  int nsurf = surf->nsurf;

  // checkflag = 1 if any checks are to be made for particles vs surfaces

  checkflag = 0;
  for (ispecies = 0; ispecies < nspecies; ispecies++)
    for (isurf = 0; isurf < nsurf; isurf++)
      if (part[ispecies].surf[isurf].check != NONE) checkflag = 1;
  
  // if debug requested, error if not enabled

#ifndef DEBUG_MOVE
  if (debug_proc >= 0)
    error->all("Code is not compiled with move debug option");
#endif

  // check that all permeability STICK species are 2d

  for (ispecies = 0; ispecies < nspecies; ispecies++)
    for (isurf = 0; isurf < nsurf; isurf++) {
      prob = part[ispecies].surf[isurf].in_prob[STICK] -
	part[ispecies].surf[isurf].in_prob[NEAR];
      jspecies = part[ispecies].surf[isurf].in_species[STICK];
      if (prob > 0.0 && particle->dimension[jspecies] != 2)
	error->all("Sticking species is not 2d");
      prob = part[ispecies].surf[isurf].out_prob[STICK] -
	part[ispecies].surf[isurf].out_prob[NEAR];
      jspecies = part[ispecies].surf[isurf].out_species[STICK];
      if (prob > 0.0 && particle->dimension[jspecies] != 2)
	error->all("Sticking species is not 2d");
    }

  // set function pointers for correct move styles

  if (geomflag == CUBE && sampleflag == UNIFORM)
    delta3d = &Move::delta_cube_uniform;
  else if (geomflag == CUBE && sampleflag == BROWNIAN)
    delta3d = &Move::delta_cube_brownian;
  else if (geomflag == SPHERE && sampleflag == UNIFORM)
    delta3d = &Move::delta_sphere_uniform;
  else if (geomflag == SPHERE && sampleflag == BROWNIAN)
    delta3d = &Move::delta_sphere_brownian;
  else if (geomflag == SQUARE && sampleflag == UNIFORM)
    delta3d = &Move::delta_square_uniform;
  else if (geomflag == SQUARE && sampleflag == BROWNIAN)
    delta3d = &Move::delta_square_brownian;
  else if (geomflag == CIRCLE && sampleflag == UNIFORM)
    delta3d = &Move::delta_square_uniform;
  else if (geomflag == CIRCLE && sampleflag == BROWNIAN)
    delta3d = &Move::delta_square_brownian;
  else error->all("No matching move style parameters");

  if (sampleflag == UNIFORM)
    delta2d = &Move::delta_2d_uniform;
  else if (sampleflag == BROWNIAN)
    delta2d = &Move::delta_2d_brownian;

  // setup move parameters for each species
  // scale converts to diffusion factor for a particular species
  // displace is used to store pre-computed samplings for Brownian
  // maxmove is max distance moved which must account for all dimensions
  //   since a surface relection can change particle direction

  double dt = simulator->dt;
  scale = new double[nspecies];
  displace3d = new double[NSLICE];
  displace2d = new double[NSLICE];
  maxmove = new double[nspecies];

  // sample uniformly
  // 3d samples box extent from -0.5 to 0.5
  // 2d samples radius from 0 to 1

  if (sampleflag == UNIFORM) {
    for (ispecies = 0; ispecies < nspecies; ispecies++) {
      if (particle->dimension[ispecies] == 3) {
	scale[ispecies] = 32.0 / (3.0 * sqrt(PI)) * 
	  sqrt(1.0e8*particle->diffusivity[ispecies] * dt);
	maxmove[ispecies] = scale[ispecies] * 0.5;
	if (geomflag == CUBE) maxmove[ispecies] *= sqrt(3.0);
	else if (geomflag == SQUARE) maxmove[ispecies] *= sqrt(2.0);
      } else {
	scale[ispecies] = 1.5 * sqrt(PI) * 
	  sqrt(1.0e8*particle->diffusivity[ispecies] * dt);
	maxmove[ispecies] = scale[ispecies];
      }
    }
  }

  // sample in Brownian sense
  // 3d samples 1d Gaussians from 1/NSLICE to (NSLICE-1)/NSLICE
  // 2d samples 2d cummulative distribution from 
  //   1 / 2*NSLICE to (2*NSLICE-1) / 2*NSLICE

  if (sampleflag == BROWNIAN) {
    for (ispecies = 0; ispecies < nspecies; ispecies++) {

      if (particle->dimension[ispecies] == 3) {

	scale[ispecies] = 
	  sqrt(4.0 * 1.0e8*particle->diffusivity[ispecies] * dt);

	double lo,hi,guess,fraction,target,difference;
	for (int i = 1; i < NSLICE; i += 2) {
	  fraction = 1.0*i / NSLICE;
	  target = 1.0 - fraction;
	  lo = 0.0;
	  hi = 1.0e6;
	  guess = 0.5 * (lo+hi);
	  difference = erfcc(guess) - target;
	  while (fabs(difference) > TOLERANCE) {
	    if (difference < 0.0) hi = guess;
	    else lo = guess;
	    guess = 0.5 * (lo+hi);
	    difference = erfcc(guess) - target;
	  }
	  displace3d[NSLICE/2 + (i-1)/2] = guess;
	  displace3d[NSLICE/2 - (i+1)/2] = -guess;
	}
	
	maxradsq = displace3d[NSLICE-1] * displace3d[NSLICE-1];
	maxradsq_one = maxradsq + 1.0;

	maxmove[ispecies] = scale[ispecies] * displace3d[NSLICE-1];
	if (geomflag == CUBE) maxmove[ispecies] *= sqrt(3.0);
	else if (geomflag == SQUARE) maxmove[ispecies] *= sqrt(2.0);

      } else {

	scale[ispecies] = 
	  sqrt(4.0 * 1.0e8*particle->diffusivity[ispecies] * dt);

	double fraction;
	for (int i = 0; i < NSLICE; i++) {
	  fraction = (2.0*i+1) / (2.0*NSLICE);
	  displace2d[i] = sqrt(-log(1.0-fraction));
	}

	maxmove[ispecies] = scale[ispecies] * displace2d[NSLICE-1];
      }
    }
  }

  // print move statistics

  int max3dspec = -1;
  int max2dspec = -1;
  double max3ddiff = 0.0;
  double max2ddiff = 0.0;

  for (ispecies = 0; ispecies < nspecies; ispecies++) {
    if (particle->dimension[ispecies] == 3) {
      if (particle->diffusivity[ispecies] > max3ddiff) {
	max3dspec = ispecies;
	max3ddiff = particle->diffusivity[ispecies];
      }
    } else {
      if (particle->diffusivity[ispecies] > max2ddiff) {
	max2dspec = ispecies;
	max2ddiff = particle->diffusivity[ispecies];
      }
    }
  }

  double ave3d,ave2d;

  if (max3dspec >= 0)
    ave3d = 4.0/sqrt(PI) * sqrt(1.0e8*particle->diffusivity[max3dspec] * dt);
  if (max2dspec >= 0)
    ave2d = 4.0/sqrt(PI) * sqrt(1.0e8*particle->diffusivity[max2dspec] * dt);
  
  if (me == 0) {
    fprintf(screen,"Move:\n");
    if (max3dspec >= 0)
      fprintf(screen,"  ave/max move of 3d particle with largest D = %g %g\n",
	      ave3d,maxmove[max3dspec]);
    if (max2dspec >= 0)
      fprintf(screen,"  ave/max move of 2d particle with largest D = %g %g\n",
	      ave2d,maxmove[max2dspec]);
    fprintf(logfile,"Move:\n");
    if (max3dspec >= 0)
      fprintf(logfile,"  ave/max move of 3d particle with largest D = %g %g\n",
	      ave3d,maxmove[max3dspec]);
    if (max2dspec >= 0)
      fprintf(logfile,"  ave/max move of 2d particle with largest D = %g %g\n",
	      ave2d,maxmove[max2dspec]);
  }
 
  // error if move distance exceeds bin size

  if (max3dspec >= 0) {
    if (maxmove[max3dspec] > grid->xbinsize || 
	maxmove[max3dspec] > grid->ybinsize || 
	maxmove[max3dspec] > grid->zbinsize)
      error->all("Max 3d move distance > bin size");
  }
  if (max2dspec >= 0) {
    if (maxmove[max2dspec] > grid->xbinsize || 
	maxmove[max2dspec] > grid->ybinsize || 
	maxmove[max2dspec] > grid->zbinsize)
      error->all("Max 2d move distance > bin size");
  }

  // set size of collision-with-triangle lists to 2 * maxperbin
  // max = scenario where starts at vertex and transmits thru vertex

  int nbins = grid->nbins;
  Grid::OneBin *blist = grid->blist;

  int maxperbin = 0;
  for (int i = 0; i < nbins; i++) maxperbin = MAX(maxperbin,blist[i].ntri);
  exlist = new int[2*maxperbin];

  // init counters

  nmove = 0;
  ncheck = 0;
  nreflect = nnear = nstick = nfar = nthru = 0;
}

/* ----------------------------------------------------------------------
   move each particle I own
------------------------------------------------------------------------- */

void Move::motion()
{
  double *xold,*v0,*v1,*v2;
  double delta[3],xnew[3],xborder[3],xcollide[3];
  double xfirst[3],dir[3],regnormal[3],whichnormal[3];
  double param,minparam,rn,distance,moved,cut,mincut;
  double *problist,*normal;
  int ibin,jbin,species,style,ncollide,iconnect,niterate,flag;
  int itri,jtri,whichside,whichtri,side,iregion;
  int isurf,edge,nexclude,iexclude,ix,iy,iz,jx,jy,jz,whichdim,nextbin;
  int *speclist;

  Grid::OneBin *blist = grid->blist;
  Surf::OneTri *tlist = surf->tlist;
  double **vlist = surf->vlist;
  Region **rlist = surf->rlist;
  Particle::OnePart *plist = particle->plist;
  int nlocal = particle->nlocal;
  double *diffusivity = particle->diffusivity;
  int *dimension = particle->dimension;
  int *region2surf = surf->region2surf;

  // increment counter

  nmove += nlocal;

  // loop over owned particles

  for (int i = 0; i < nlocal; i++) {

    // if diffusion constant = 0.0, no movement
    // just update RN seed so reactions will be randomized

    species = plist[i].species;
    if (diffusivity[species] == 0.0) {
      random->move(&plist[i].seed);
      continue;
    }

    xold = plist[i].x;
    ibin = plist[i].ibin;
    itri = plist[i].itri;

    // 3d particles

    if (dimension[species] == 3) {

      // delta = 3d displacement
      // jbin = bin for xnew
      
      (this->*delta3d)(species,&plist[i].seed,delta);
      xnew[0] = xold[0] + delta[0];
      xnew[1] = xold[1] + delta[1];
      xnew[2] = xold[2] + delta[2];
      jbin = grid->whichlocal(xnew);

      // loop on 3d move among triangles/regions:
      //   enter with xold, ibin, xnew, jbin, species, exclude list
      //   exit with xnew, jbin, species, itri (if now 2d)
      //   exclude list is list of tri IDs to skip for collision
      //   regions are not excluded b/c can be hit multiple times in 1 move
      //     e.g. grazing along inwardly curved surface
      // while (1)
      //   if ibin = jbin:
      //     xborder = xnew
      //   else:
      //     find face of ibin the move hits first
      //       if hits edge or vertex of ibin, any face is OK
      //     xborder = intersection point with ibin face
      //     nextbin = bin on opposite side of that face
      //   if any triangles or region indices in ibin:
      //     loop over them and test for collisions with xold->xborder segment
      //     use exclude list to skip a tri in exclude list
      //     store collision with smallest param
      //     if collision:
      //       random # with permeability settings -> style of collsion
      //       change species according to style
      //       if REFLECT:
      //         xold = collision pt
      //         xnew = reflect_normal() off triangle/region normal
      //         if smallest param > 0, clear exclude list (can be 0 on vertex)
      //         if triangle collided with, add to exclude list
      //         if region collided with:
      //           xold = push_off() along normal by EPS in correct dir
      //           this prevents instant collision with region on next iter
      //           re-compute ibin of xold
      //         compute jbin of xnew
      //         back to beginning of loop
      //       if NEAR:
      //         xnew = push_off() along tri/region normal in correct dir
      //         compute jbin of xnew
      //         move is done, break from loop
      //       if STICK:
      //         xnew = collision pt
      //         set itri for particle
      //         compute jbin of xnew
      //         move is done, break from loop
      //       if FAR:
      //         xnew = push_off() along tri/region normal in correct dir
      //         compute jbin of xnew
      //         move is done, break from loop
      //       if THRU: 
      //         xold = collision pt
      //         xnew, jbin = unchanged
      //         if smallest param > 0, clear exclude list (can be 0 on vertex)
      //         if triangle collided with, add to exclude list
      //         if region collided with:
      //           xold = push_off() along normal by EPS in correct dir
      //           this prevents instant collision with region on next iter
      //           re-compute ibin of xold
      //         back to beginning of loop
      //   if ibin != jbin:
      //     xold = xborder
      //     clear exclude list
      //     ibin = nextbin
      //     back to beginning of loop
      //   else ibin = jbin, move is done, break from loop

      niterate = nexclude = 0;

      while (1) {

	niterate++;
	if (niterate > MAXITER) {
	  char str[128];
	  sprintf(str,"Part,spec %d,%d on step %d exceeds MAXITER",
		  i,plist[i].species,simulator->ntimestep);
	  error->one(str);
	}

#ifdef DEBUG_MOVE
	if (me == debug_proc && simulator->ntimestep == debug_step && 
	    i == debug_index) {
	  printf("DEBUG iteration %d\n",niterate);
	  printf("  binning: %d\n",grid->whichlocal(xold));
	  printf("  xold: %g %g %g\n",xold[0],xold[1],xold[2]);
	  printf("  xnew: %g %g %g\n",xnew[0],xnew[1],xnew[2]);
	  printf("  ibin: %d, jbin: %d, species %d\n",ibin,jbin,species);
	  printf("  nexclude %d\n",nexclude);
	  for (iexclude = 0; iexclude < nexclude; iexclude++)
	    printf("    %d",exlist[iexclude]);
	  if (nexclude) printf("\n");
	}
#endif

	if (ibin == jbin) {
	  xborder[0] = xnew[0];
	  xborder[1] = xnew[1];
	  xborder[2] = xnew[2];
	  nextbin = ibin;

	} else {
	  grid->local2global(ibin,&ix,&iy,&iz);
	  grid->local2global(jbin,&jx,&jy,&jz);
	  minparam = 2.0;

	  if (ix != jx) {
	    if (ix < jx) cut = grid->coord_cut(0,jx);
	    else cut = grid->coord_cut(0,ix);
	    param = (cut-xold[0]) / (xnew[0]-xold[0]);
	    if (param < minparam) {
	      minparam = param;
	      whichdim = 0;
	      mincut = cut;
	    }
	  }
	  if (iy != jy) {
	    if (iy < jy) cut = grid->coord_cut(1,jy);
	    else cut = grid->coord_cut(1,iy);
	    param = (cut-xold[1]) / (xnew[1]-xold[1]);
	    if (param < minparam) {
	      minparam = param;
	      whichdim = 1;
	      mincut = cut;
	    }
	  }
	  if (iz != jz) {
	    if (iz < jz) cut = grid->coord_cut(2,jz);
	    else cut = grid->coord_cut(2,iz);
	    param = (cut-xold[2]) / (xnew[2]-xold[2]);
	    if (param < minparam) {
	      minparam = param;
	      whichdim = 2;
	      mincut = cut;
	    }
	  }
	  
	  if (whichdim == 0) {
	    xborder[0] = mincut;
	    xborder[1] = xold[1] + minparam * (xnew[1]-xold[1]);
	    xborder[2] = xold[2] + minparam * (xnew[2]-xold[2]);
	    nextbin = grid->global2local(jx,iy,iz);
	  } else if (whichdim == 1) {
	    xborder[0] = xold[0] + minparam * (xnew[0]-xold[0]);
	    xborder[1] = mincut;
	    xborder[2] = xold[2] + minparam * (xnew[2]-xold[2]);
	    nextbin = grid->global2local(ix,jy,iz);
	  } else {
	    xborder[0] = xold[0] + minparam * (xnew[0]-xold[0]);
	    xborder[1] = xold[1] + minparam * (xnew[1]-xold[1]);
	    xborder[2] = mincut;
	    nextbin = grid->global2local(ix,iy,jz);
	  }
	  
#ifdef DEBUG_MOVE
	  if (me == debug_proc && simulator->ntimestep == debug_step && 
	      i == debug_index) {
	    printf("  xborder: %g %g %g\n",xborder[0],xborder[1],xborder[2]);
	    printf("  bin status\n");
	    printf("    ibin,glob ix,iy,iz: %d, %d %d %d\n",ibin,ix,iy,iz);
	    printf("    jbin,glob jx,jy,jz: %d, %d %d %d\n",jbin,jx,jy,jz);
	    printf("    nextbin,loc nbinx,nbiny,nbinz: %d, %d %d %d\n",
		   nextbin,grid->nbinx,grid->nbiny,grid->nbinz);
	    printf("    whichdim,minparam,mincut: %d %g %g\n",
		   whichdim,minparam,mincut);
	  }
#endif

	  // error check which should never be needed

	  if (nextbin < 0 || nextbin >= grid->nbins)
	    error->one("Move to next bin out of range");
	}

	if (blist[ibin].ntri) {
	  ncollide = 0;
	  minparam = 2.0;

	  Grid::Index *ptr = blist[ibin].tri;
	  while (ptr) {
	    jtri = ptr->index;
	    ptr = ptr->ptr;

#ifdef DEBUG_MOVE
	    if (me == debug_proc && simulator->ntimestep == debug_step && 
		i == debug_index) {
	      printf("  triangle check:\n");
	      if (jtri >= 0) {
		printf("    tri, verts: %d: %d %d %d\n",
		       jtri+1,tlist[jtri].vert[0]+1,
		       tlist[jtri].vert[1]+1,tlist[jtri].vert[2]+1);
	      } else
		printf("    local region: %d %s\n",
		       -jtri-1,surf->name[region2surf[-jtri-1]]);
	    }
#endif

	    if (nexclude > 0) {
	      for (iexclude = 0; iexclude < nexclude; iexclude++)
		if (jtri == exlist[iexclude]) break;
	      if (iexclude < nexclude) continue;
	    }

#ifdef DEBUG_MOVE
	    if (me == debug_proc && simulator->ntimestep == debug_step && 
		i == debug_index) {
	      printf("      not excluded\n");
	    }
#endif

	    ncheck++;
	    flag = 0;
	    if (jtri >= 0) {
	      v0 = vlist[tlist[jtri].vert[0]];
	      v1 = vlist[tlist[jtri].vert[1]];
	      v2 = vlist[tlist[jtri].vert[2]];
	      normal = tlist[jtri].normal;
	      flag = tri_line_intersect(v0,v1,v2,normal,
					xold,xborder,xcollide,param,side);
	    } else {
	      normal = regnormal;
	      flag = rlist[-jtri-1]->line_intersect(xold,xborder,xcollide,
						    normal,param,side);
	    }

#ifdef DEBUG_MOVE
	    if (me == debug_proc && simulator->ntimestep == debug_step && 
		i == debug_index) {
	      if (!flag) printf("      no collision\n");
	      else {
		printf("      collision: param,whichside: %g %d\n",
		       param,side);
		printf("                 pt: %g %g %g\n",
		       xcollide[0],xcollide[1],xcollide[2]);
		printf("                 norm: %g %g %g\n",
		       normal[0],normal[1],normal[2]);
	      }
	    }
#endif

	    if (flag && param < minparam) {
	      ncollide = 1;
	      minparam = param;
	      whichtri = jtri;
	      whichside = side;
	      xfirst[0] = xcollide[0];
	      xfirst[1] = xcollide[1];
	      xfirst[2] = xcollide[2];
	      whichnormal[0] = normal[0];
	      whichnormal[1] = normal[1];
	      whichnormal[2] = normal[2];
	    }
	  }

	  if (ncollide) {
	    if (whichtri >= 0) isurf = tlist[whichtri].isurf;
	    else {
	      iregion = -whichtri - 1;
	      isurf = region2surf[iregion];
	    }

	    if (whichside == OUTSIDE) {
	      speclist = part[species].surf[isurf].out_species;
	      problist = part[species].surf[isurf].out_prob;
	    } else {
	      speclist = part[species].surf[isurf].in_species;
	      problist = part[species].surf[isurf].in_prob;
	    }

	    rn = random->move(&plist[i].seed);

	    if (rn < problist[0]) style = REFLECT;
	    else if (rn < problist[1]) style = NEAR;
	    else if (rn < problist[2]) style = STICK;
	    else if (rn < problist[3]) style = FAR;
	    else style = THRU;

	    species = speclist[style];

	    if (style == REFLECT) {
	      nreflect++;
	      xold[0] = xfirst[0];
	      xold[1] = xfirst[1];
	      xold[2] = xfirst[2];
	      reflect_normal(xold,xnew,whichnormal);
	      if (minparam > 0.0) nexclude = 0;
	      if (whichtri >= 0) exlist[nexclude++] = whichtri;
	      else {
		if (whichside == OUTSIDE) 
		  push_off(xfirst,whichnormal,EPSILON,xold);
		else
		  push_off(xfirst,whichnormal,-EPSILON,xold);
		ibin = grid->whichlocal(xold);
	      }
	      jbin = grid->whichlocal(xnew);
	      continue;

	    } else if (style == NEAR) {
	      nnear++;
	      if (whichside == OUTSIDE) 
		push_off(xfirst,whichnormal,EPSILON,xnew);
	      else
		push_off(xfirst,whichnormal,-EPSILON,xnew);
	      jbin = grid->whichlocal(xnew);
	      break;

	    } else if (style == STICK) {
	      nstick++;
	      xnew[0] = xfirst[0];
	      xnew[1] = xfirst[1];
	      xnew[2] = xfirst[2];
	      itri = whichtri;
	      jbin = grid->whichlocal(xnew);
	      break;

	    } else if (style == FAR) {
	      nfar++;
	      if (whichside == OUTSIDE)
		push_off(xfirst,whichnormal,-EPSILON,xnew);
	      else
		push_off(xfirst,whichnormal,EPSILON,xnew);
	      jbin = grid->whichlocal(xnew);
	      break;

	    } else {
	      nthru++;
	      xold[0] = xfirst[0];
	      xold[1] = xfirst[1];
	      xold[2] = xfirst[2];
	      if (minparam > 0.0) nexclude = 0;
	      if (whichtri >= 0) exlist[nexclude++] = whichtri;
	      else {
		if (whichside == OUTSIDE) 
		  push_off(xfirst,whichnormal,-EPSILON,xold);
		else
		  push_off(xfirst,whichnormal,EPSILON,xold);
		ibin = grid->whichlocal(xold);
	      }
	      continue;
	    }
	  }
	}

	if (ibin != jbin) {
	  xold[0] = xborder[0];
	  xold[1] = xborder[1];
	  xold[2] = xborder[2];
	  nexclude = 0;
	  ibin = nextbin;
	} else break;
      }

    } else {

      // 2d particles

      // normal = normal to plane of motion
      // distance/dir = distance/dir to move in plane of triangle/region itri
      
      if (itri >= 0) normal = tlist[itri].normal;
      else {
	normal = regnormal;
	rlist[-itri-1]->compute_normal(xold,normal);
      }
      (this->*delta2d)(species,&plist[i].seed,normal,&distance,dir);

      // loop on 2d move on triangles:
      //   enter with xold, itri, distance, dir to move
      //   exit with xnew, jbin, updated itri
      // point_tri_move for each iteration:
      //   computes xnew as interior point or collision with edge/vertex
      //   if interior:
      //     move is done, compute jbin and exit
      //   if partial move hitting an edge:
      //     update distance
      //     xold = collision pt (xnew)
      //     if hit connected edge:
      //       itri = adjacent triangle
      //       dir = bend2d() on normals of old and new triangles
      //       normal = normal of new adjacent triangle
      //     else if hit unconnected edge:
      //       dir = reflect_edge() of 2 vertices of edge
      //       itri and normal are unchanged
      //   if partial move hitting a vertex:
      //     error for now
      //     should try again and exclude current edge via exclude list?

      if (itri >= 0) {
	niterate = 0;

	while (1) {
	  niterate++;
	  if (niterate > MAXITER) {
	    char str[128];
	    sprintf(str,"Part,spec %d,%d on step %d exceeds MAXITER",
		    i,plist[i].species,simulator->ntimestep);
	    error->one(str);
	  }

#ifdef DEBUG_MOVE
	  if (me == debug_proc && simulator->ntimestep == debug_step && 
	      i == debug_index) {
	    printf("DEBUG iteration %d\n",niterate);
	    printf("  xold: %g %g %g\n",xold[0],xold[1],xold[2]);
	    printf("  ibin: %d, itri: %d, species %d\n",ibin,itri,species);
	    printf("  distance, dir: %g %g %g %g\n",
		   distance,dir[0],dir[1],dir[2]);
	  }
#endif

	  v0 = vlist[tlist[itri].vert[0]];
	  v1 = vlist[tlist[itri].vert[1]];
	  v2 = vlist[tlist[itri].vert[2]];

	  if (point_tri_move(v0,v1,v2,normal,xold,dir,distance,
			     xnew,moved,edge)) {
	    jbin = grid->whichlocal(xnew);
	    break;
	  }

	  if (edge == 0) error->one("Invalid 2d particle move");
	  distance -= moved;
	  xold[0] = xnew[0];
	  xold[1] = xnew[1];
	  xold[2] = xnew[2];

#ifdef DEBUG_MOVE
	  if (me == debug_proc && simulator->ntimestep == debug_step && 
	      i == debug_index) {
	    printf("  distance, moved, edge: %g %g %d\n",distance,moved,edge);
	    printf("  xnew: %g %g %g\n",xnew[0],xnew[1],xnew[2]);
	  }
#endif

	  if (edge <= 3) {
	    iconnect = tlist[itri].connect[edge-1];
	    if (iconnect >= 0) {
	      itri = iconnect;
	      bend2d(normal,tlist[itri].normal,dir);
	      normal = tlist[itri].normal;
	    } else {
	      if (edge == 1) reflect_edge(v0,v1,dir);
	      else if (edge == 2) reflect_edge(v1,v2,dir);
	      else reflect_edge(v2,v0,dir);
	    }
	  } else error->one("2d particle move hit a triangle vertex");

#ifdef DEBUG_MOVE
	  if (me == debug_proc && simulator->ntimestep == debug_step && 
	      i == debug_index) {
	    printf("  new tri to move on: %d\n",itri);
	    printf("  new distance to move: %g\n",distance);
	    printf("  new dir to move: %g %g %g\n",dir[0],dir[1],dir[2]);
	    printf("  new xnew to move to: %g %g %g\n",xnew[0],xnew[1],xnew[2]);
	  }
#endif
	}

      // 2d move on region surface:
      //   enter with xold, distance, dir to move
      //   exit with xnew, jbin
      // move2d returns xnew

      } else {
	rlist[itri-1]->move2d(xold,normal,dir,distance,xnew);
	jbin = grid->whichlocal(xnew);
      }
    }

#ifdef DEBUG_MOVE
      if (me == debug_proc && simulator->ntimestep == debug_step && 
	  i == debug_index) {
	printf("  final xnew: %g %g %g\n",xnew[0],xnew[1],xnew[2]);
	printf("  final bin: %d\n",jbin);
      }
#endif

    // set particle's new position, bin, triangle, and species
	
    plist[i].x[0] = xnew[0];
    plist[i].x[1] = xnew[1];
    plist[i].x[2] = xnew[2];
    plist[i].ibin = jbin;
    plist[i].itri = itri;
    plist[i].species = species;
  }

  if (checkflag) check();
}

/* ----------------------------------------------------------------------
   check if my particles are inside/outside surfaces as requested
------------------------------------------------------------------------- */

void Move::check()
{
  int isurf,ispecies,nhit,i,m,flag,iregion;
  double *v0,*v1,*v2;
  double xoutside[3],xcollide[3];
  double param;
  int flagface;

  Surf::OneTri *tlist = surf->tlist;
  double **vlist = surf->vlist;
  Region **rlist = surf->rlist;
  Particle::OnePart *plist = particle->plist;
  int *surf2region = surf->surf2region;
  int *dimension = particle->dimension;
  int nlocal = particle->nlocal;

  // double loop over all surfaces and particles

  for (isurf = 0; isurf < surf->nsurf; isurf++) {
    for (i = 0; i < nlocal; i++) {

      // only check if a 3d species and this surf and species are active

      ispecies = plist[i].species;
      if (dimension[ispecies] == 2) continue;
      if (part[ispecies].surf[isurf].check == NONE) continue;
      if (part[ispecies].surf[isurf].checklo > simulator->ntimestep) continue;
      if (part[ispecies].surf[isurf].checkhi < simulator->ntimestep) continue;
      if (part[ispecies].surf[isurf].plo > i) continue;
      if (part[ispecies].surf[isurf].phi < i) continue;

      // if surface is a region, test inside/outside directly
      // if surface is triangluted, expensive loop over all triangles
      //   if # of triangle hits is odd, particle is inside, else outside
      //   skip triangles not in this surface

      flag = nhit = 0;

      if (surf2region[isurf] >= 0) {
	iregion = surf2region[isurf];
	nhit = rlist[iregion]->inside(plist[i].x);
	if (nhit == -1 && part[ispecies].surf[isurf].check == INSIDE)
	  flag = 1;
	else if (nhit == 1 && part[ispecies].surf[isurf].check == OUTSIDE)
	  flag = 1;

      } else {
	xoutside[0] = domain->xlo - 1.0;
	xoutside[1] = plist[i].x[1];
	xoutside[2] = plist[i].x[2];

	for (m = 0; m < surf->ntri; m++) {
	  if (tlist[m].isurf != isurf) continue;
	  v0 = vlist[tlist[m].vert[0]];
	  v1 = vlist[tlist[m].vert[1]];
	  v2 = vlist[tlist[m].vert[2]];
	  if (tri_line_intersect(v0,v1,v2,tlist[m].normal,plist[i].x,
				 xoutside,xcollide,param,flagface)) nhit++;
	}
    
	if (nhit % 2 == 0 && part[ispecies].surf[isurf].check == INSIDE)
	  flag = 1;
	else if (nhit % 2 == 1 && part[ispecies].surf[isurf].check == OUTSIDE)
	  flag = 1;
      }
      
      if (flag) {
	if (part[ispecies].surf[isurf].check == INSIDE)
	  error->warning("Particle is outside surface");
	else
	  error->warning("Particle is inside surface");
	printf("  Timestep %d\n",simulator->ntimestep);
	printf("  Particle species, index, hit count = %s %d %d\n",
	       particle->name[ispecies],i,nhit);
	printf("  Particle coord = %g %g %g\n",
	       plist[i].x[0],plist[i].x[1],plist[i].x[2]);
	printf("  Surface name, index = %s %d\n",surf->name[isurf],isurf);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   2d uniform move in a xy square
------------------------------------------------------------------------- */

void Move::delta_square_uniform(int ispecies, int *seed, double *delta)
{
  double prefactor = scale[ispecies];
  delta[0] = prefactor * (random->move(seed) - 0.5);
  delta[1] = prefactor * (random->move(seed) - 0.5);
  delta[2] = 0.0;
}

/* ----------------------------------------------------------------------
   3d uniform move in a cube
------------------------------------------------------------------------- */

void Move::delta_cube_uniform(int ispecies, int *seed, double *delta)
{
  double prefactor = scale[ispecies];
  delta[0] = prefactor * (random->move(seed) - 0.5);
  delta[1] = prefactor * (random->move(seed) - 0.5);
  delta[2] = prefactor * (random->move(seed) - 0.5);
}

/* ----------------------------------------------------------------------
   2d Brownian move in a square
------------------------------------------------------------------------- */

void Move::delta_square_brownian(int ispecies, int *seed, double *delta)
{
  int m;
  double prefactor = scale[ispecies];
  m = static_cast<int> (NSLICE * random->move(seed));
  delta[0] = prefactor * displace3d[m];
  m = static_cast<int> (NSLICE * random->move(seed));
  delta[1] = prefactor * displace3d[m];
  delta[2] = 0.0;
}

/* ----------------------------------------------------------------------
   3d Brownian move in a cube
------------------------------------------------------------------------- */

void Move::delta_cube_brownian(int ispecies, int *seed, double *delta)
{
  int m;
  double prefactor = scale[ispecies];
  m = static_cast<int> (NSLICE * random->move(seed));
  delta[0] = prefactor * displace3d[m];
  m = static_cast<int> (NSLICE * random->move(seed));
  delta[1] = prefactor * displace3d[m];
  m = static_cast<int> (NSLICE * random->move(seed));
  delta[2] = prefactor * displace3d[m];
}

/* ----------------------------------------------------------------------
   2d uniform move in a circle
------------------------------------------------------------------------- */

void Move::delta_circle_uniform(int ispecies, int *seed, double *delta)
{
  double x,y;
  double rsq = 1.0;
  while (rsq > 0.25) {
    x = random->move(seed) - 0.5;
    y = random->move(seed) - 0.5;
    rsq = x*x + y*y;
  }
  double prefactor = scale[ispecies];
  delta[0] = prefactor * x;
  delta[1] = prefactor * y;
  delta[2] = 0.0;
}

/* ----------------------------------------------------------------------
   3d uniform move in a sphere
------------------------------------------------------------------------- */

void Move::delta_sphere_uniform(int ispecies, int *seed, double *delta)
{
  double x,y,z;
  double rsq = 1.0;
  while (rsq > 0.25) {
    x = random->move(seed) - 0.5;
    y = random->move(seed) - 0.5;
    z = random->move(seed) - 0.5;
    rsq = x*x + y*y + z*z;
  }
  double prefactor = scale[ispecies];
  delta[0] = prefactor * x;
  delta[1] = prefactor * y;
  delta[2] = prefactor * z;
}

/* ----------------------------------------------------------------------
   2d Brownian move in a circle
------------------------------------------------------------------------- */

void Move::delta_circle_brownian(int ispecies, int *seed, double *delta)
{
  int m;
  double x,y;
  double rsq = maxradsq_one;
  while (rsq > maxradsq) {
    m = static_cast<int> (NSLICE * random->move(seed));
    x = displace3d[m];
    m = static_cast<int> (NSLICE * random->move(seed));
    y = displace3d[m];
    rsq = x*x + y*y;
  }
  double prefactor = scale[ispecies];
  delta[0] = prefactor * x;
  delta[1] = prefactor * y;
  delta[2] = 0.0;
}

/* ----------------------------------------------------------------------
   3d Brownian move in a sphere
------------------------------------------------------------------------- */

void Move::delta_sphere_brownian(int ispecies, int *seed, double *delta)
{
  int m;
  double x,y,z;
  double rsq = maxradsq_one;
  while (rsq > maxradsq) {
    m = static_cast<int> (NSLICE * random->move(seed));
    x = displace3d[m];
    m = static_cast<int> (NSLICE * random->move(seed));
    y = displace3d[m];
    m = static_cast<int> (NSLICE * random->move(seed));
    z = displace3d[m];
    rsq = x*x + y*y + z*z;
  }
  double prefactor = scale[ispecies];
  delta[0] = prefactor * x;
  delta[1] = prefactor * y;
  delta[2] = prefactor * z;
}

/* ----------------------------------------------------------------------
   2d uniform move in a circle
   generate distance to move and random unit direction
------------------------------------------------------------------------- */

void Move::delta_2d_uniform(int ispecies, int *seed,
			    double *normal, double *distance, double *dir)
{
  // distance distribution varies as square of radius

  *distance = scale[ispecies] * sqrt(random->move(seed));

  // random dir = (normal of itri) x (random 3-vec)
  // insure dir is not 0-length

  double lensq = 0.0;

  while (lensq == 0.0) {
    double rx = random->move(seed) - 0.5;
    double ry = random->move(seed) - 0.5;
    double rz = random->move(seed) - 0.5;

    dir[0] = normal[1]*rz - normal[2]*ry;
    dir[1] = normal[2]*rx - normal[0]*rz;
    dir[2] = normal[0]*ry - normal[1]*rx;

    lensq = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
  }

  // convert dir to unit vector

  double scale = 1.0/sqrt(lensq);
  dir[0] *= scale;
  dir[1] *= scale;
  dir[2] *= scale;
}

/* ----------------------------------------------------------------------
   2d Brownian move in a circle
   generate distance to move and random unit direction
------------------------------------------------------------------------- */

void Move::delta_2d_brownian(int ispecies, int *seed,
			     double *normal, double *distance, double *dir)
{
  // distance distribution from displace2d[]

  int m = static_cast<int> (NSLICE * random->move(seed));
  *distance = scale[ispecies] * displace2d[m];

  // random dir = (normal of itri) x (random 3-vec)
  // insure dir is not 0-length

  double lensq = 0.0;

  while (lensq == 0.0) {
    double rx = random->move(seed) - 0.5;
    double ry = random->move(seed) - 0.5;
    double rz = random->move(seed) - 0.5;

    dir[0] = normal[1]*rz - normal[2]*ry;
    dir[1] = normal[2]*rx - normal[0]*rz;
    dir[2] = normal[0]*ry - normal[1]*rx;

    lensq = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
  }

  // convert dir to unit vector

  double scale = 1.0/sqrt(lensq);
  dir[0] *= scale;
  dir[1] *= scale;
  dir[2] *= scale;
}

/* ----------------------------------------------------------------------
   Numerical Recipes erfc(x) function
------------------------------------------------------------------------- */

double Move::erfcc(double x)
{
  double t,z,ans;

  z = fabs(x);
  t = 1.0/(1.0+0.5*z);
  ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
	      t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
	      t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}

/* ----------------------------------------------------------------------
   compute bin size needed for a given diffusivity
   called by bin::global() which sets up bins
   assume 3d diffusion and current settings for move style and dt
   does same computation as init() for maxmove[ispecies]
------------------------------------------------------------------------- */

double Move::maxbin(double diffusivity)
{
  double max;

  if (sampleflag == UNIFORM)
    max = 0.5 * 32.0 /
      (3.0 * sqrt(PI)) * sqrt(1.0e8*diffusivity*simulator->dt);

  else {
    double lo,hi,guess,fraction,target,difference;
    fraction = 1.0*(NSLICE-1) / NSLICE;
    target = 1.0 - fraction;
    lo = 0.0;
    hi = 1.0e6;
    guess = 0.5 * (lo+hi);
    difference = erfcc(guess) - target;
    while (fabs(difference) > TOLERANCE) {
      if (difference < 0.0) hi = guess;
      else lo = guess;
      guess = 0.5 * (lo+hi);
      difference = erfcc(guess) - target;
    }
    max = sqrt(4.0*1.0e8*diffusivity*simulator->dt) * guess;
  }

  if (geomflag == CUBE) max *= sqrt(3.0);
  else if (geomflag == SQUARE) max *= sqrt(2.0);

  return max;
}

/* ----------------------------------------------------------------------
   free arrays
------------------------------------------------------------------------- */

void Move::free_arrays()
{
  if (scale) delete [] scale;
  if (displace3d) delete [] displace3d;
  if (displace2d) delete [] displace2d;
  if (maxmove) delete [] maxmove;

  if (exlist) delete [] exlist;
}
