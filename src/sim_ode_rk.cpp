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
#include "sim_ode_rk.h"
#include "particle.h"
#include "react.h"
#include "chem_ode_rk.h"
#include "random.h"
#include "modify.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#define AVOGADRO 6.023e23

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

SimODERK::SimODERK(int narg, char **arg) : Simulator(narg, arg)
{
  int nprocs;
  MPI_Comm_size(world,&nprocs);
  if (nprocs > 1) error->all("Cannot use this run style in parallel");

  spatial_flag = 0;
  stochastic_flag = 0;

  particle = new Particle;
  react = new React;
  chem = new ChemODERK;
  modify = new Modify;
}

/* ---------------------------------------------------------------------- */

SimODERK::~SimODERK()
{
  delete particle;
  delete react;
  delete chem;
  delete modify;
}

/* ----------------------------------------------------------------------
   init a simulation
------------------------------------------------------------------------- */

void SimODERK::init()
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
  chem->init();
  modify->init();
  timer->init();
  output->init();
}

/* ----------------------------------------------------------------------
   run a simulation
------------------------------------------------------------------------- */

void SimODERK::run()
{
  timer->barrier_start(TIME_LOOP);

  for (int istep = 0; istep < nsteps; istep++) {
    ntimestep++;

    if (modify->n_initial) modify->initial();

    timer->stamp();
    chem->reactions();
    timer->stamp(TIME_REACT);

    if (output->next == ntimestep) {
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }

  if (modify->n_cleanup) modify->cleanup();

  timer->barrier_stop(TIME_LOOP);
}
