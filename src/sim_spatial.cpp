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

#include "math.h"
#include "sim_spatial.h"
#include "move.h"
#include "particle.h"
#include "chem_spatial.h"
#include "grid.h"
#include "domain.h"
#include "modify.h"
#include "react.h"
#include "output.h"
#include "surf.h"
#include "balance.h"
#include "timer.h"

#define REACTANT 1
#define PRODUCT  2

#define AVOGADRO 6.023e23
#define PI 3.1415926
#define BIG 1.0e20

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

SimSpatial::SimSpatial(int narg, char **arg) : Simulator(narg, arg)
{
  spatial_flag = 1;
  stochastic_flag = 1;

  particle = new Particle;
  react = new React;
  chem = new ChemSpatial;
  domain = new Domain;
  grid = new Grid;
  move = new Move;
  surf = new Surf;
  modify = new Modify;
  balance = new Balance;
}

/* ---------------------------------------------------------------------- */

SimSpatial::~SimSpatial()
{
  delete particle;
  delete react;
  delete domain;
  delete grid;
  delete move;
  delete surf;
  delete chem;
  delete modify;
  delete balance;
}

/* ----------------------------------------------------------------------
   init a simulation
------------------------------------------------------------------------- */

void SimSpatial::init()
{
  // duration of simulation
  // convert elapse to nsteps if needed

  firststep = ntimestep;
  firsttime = ctime;

  if (nsteps >= 0) {
    laststep = firststep + nsteps;
    lasttime = firsttime + nsteps*dt;
  } else {
    nsteps = static_cast<int> ((elapse+0.5*dt)/dt);
    laststep = firststep + nsteps;
    lasttime = firsttime + nsteps*dt;
  }

  // init other classes

  particle->init();
  react->init();
  domain->init();
  grid->init();
  move->init();
  surf->init();
  chem->init();
  modify->init();
  timer->init();
  output->init();
}

/* ----------------------------------------------------------------------
   run a simulation
------------------------------------------------------------------------- */

void SimSpatial::run()
{
  timer->barrier_start(TIME_LOOP);

  for (int istep = 0; istep < nsteps; istep++) {
    ntimestep++;

    if (modify->n_initial) modify->initial();

    timer->stamp();
    move->motion();
    timer->stamp(TIME_MOVE);
    particle->migrate();
    timer->stamp(TIME_MIGRATE);
    chem->reactions();
    timer->stamp(TIME_REACT);

    if (modify->n_final) modify->final();

    if (grid->next_balance == ntimestep) {
      grid->dynamic();
      chem->setup_stencil();
      timer->stamp(TIME_BALANCE);
    }

    if (output->next == ntimestep) {
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }

  if (modify->n_cleanup) modify->cleanup();

  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   call this to check if RN are duplicated across particles (on 1 proc)
------------------------------------------------------------------------- */

void SimSpatial::sort(int flag)
{
  int n = particle->nlocal;
  int *a = new int[n];

  for (int i = 0; i < n; i++) a[i] = particle->plist[i].seed;
  qsort(a,n,sizeof(int),compare);

  for (int i = 0; i < n-1; i++)
    if (a[i] == a[i+1]) {
      printf("SEED match on step %d with flag %d: %d %d\n",
	     ntimestep,flag,a[i],a[i+1]);
      for (int j = 0; j < particle->nlocal; j++)
	if (particle->plist[j].seed == a[i])
	  printf("  part: %d %d %g %g %g\n",j,particle->plist[j].species,
		 particle->plist[j].x[0],particle->plist[j].x[1],
		 particle->plist[j].x[2]);
    }

  delete [] a;
}

/* ---------------------------------------------------------------------- */

int SimSpatial::compare(const void *pi, const void *pj)
{
  int i = *((int *) pi);
  int j = *((int *) pj);

  if (i < j) return -1;
  else if (i > j) return 1;
  else return 0;
}
