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
#include "sim_gillespie.h"
#include "particle.h"
#include "react.h"
#include "chem_gillespie.h"
#include "modify.h"
#include "random.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#define AVOGADRO 6.023e23

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

SimGillespie::SimGillespie(int narg, char **arg) : Simulator(narg, arg)
{
  int nprocs;
  MPI_Comm_size(world,&nprocs);
  if (nprocs > 1) error->all("Cannot use this run style in parallel");

  spatial_flag = 0;
  stochastic_flag = 1;

  particle = new Particle;
  react = new React;
  chem = new ChemGillespie;
  modify = new Modify;
}

/* ---------------------------------------------------------------------- */

SimGillespie::~SimGillespie()
{
  delete particle;
  delete react;
  delete chem;
  delete modify;
}

/* ----------------------------------------------------------------------
   init a simulation
------------------------------------------------------------------------- */

void SimGillespie::init()
{
  // duration of simulation

  firststep = ntimestep;
  firsttime = ctime;
  if (nsteps >= 0) {
    laststep = firststep + nsteps;
    lasttime = firsttime;
  } else {
    laststep = firststep;
    lasttime = firsttime + elapse;
  }

  // if doing output every delta, set stats_time to next time > current time

  if (nsteps >= 0) runflag = 1;
  else runflag = 0;
  if (output->stats_every >= 0) outflag = 1;
  else {
    outflag = 0;
    int n = static_cast<int> (simulator->ctime / output->stats_delta);
    stats_time = (n+1) * output->stats_delta;
  }

  // init other classes

  particle->init();
  react->init();
  chem->init();
  modify->init();
  timer->init();
  output->init();
}

/* ----------------------------------------------------------------------
   run a simulation
------------------------------------------------------------------------- */

void SimGillespie::run()
{
  int flag;
  int doneflag = 0;

  timer->barrier_start(TIME_LOOP);

  while (1) {
    ntimestep++;

    if (modify->n_initial) modify->initial();

    timer->stamp();
    chem->reactions();
    timer->stamp(TIME_REACT);

    if (runflag && ntimestep == laststep) doneflag = 1;
    else if (!runflag && ctime >= lasttime) doneflag = 1;
    else if (chem->doneflag) doneflag = 1;

    flag = 0;
    if (outflag && ntimestep == output->next_stats) flag = 1;
    else if (!outflag && simulator->ctime > stats_time) flag = 1;
    else if (doneflag) flag = 1;

    if (modify->n_final) modify->final();

    if (flag) {
      output->stats();
      if (outflag) output->next_stats += output->stats_every;
      else stats_time += output->stats_delta;
      timer->stamp(TIME_OUTPUT);
    }

    if (doneflag) break;
  }

  if (modify->n_cleanup) modify->cleanup();

  timer->barrier_stop(TIME_LOOP);
  nsteps = ntimestep - firststep;
}
