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
#include "sim_sssa.h"
#include "move.h"
#include "particle.h"
#include "chem_sssa.h"
#include "grid.h"
#include "domain.h"
#include "modify.h"
#include "random.h"
#include "react.h"
#include "output.h"
#include "surf.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

SimSpatialSSA::SimSpatialSSA(int narg, char **arg) : Simulator (narg, arg)
{
  int nprocs;
  MPI_Comm_size(world,&nprocs);
  if (nprocs > 1) error->all("Cannot use this run style in parallel");

  spatial_flag = 1;

  particle = new Particle;
  react = new React;
  chem = new ChemSpatialSSA;
  domain = new Domain;
  grid = new Grid;
  move = new Move;
  surf = new Surf;
  modify = new Modify;
}

/* ---------------------------------------------------------------------- */

SimSpatialSSA::~SimSpatialSSA()
{
  delete particle;
  delete react;
  delete domain;
  delete grid;
  delete move;
  delete surf;
  delete chem;
  delete modify;
}

/* ----------------------------------------------------------------------
   init a simulation
------------------------------------------------------------------------- */

void SimSpatialSSA::init()
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

void SimSpatialSSA::run()
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

    if (output->next == ntimestep) {
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }

  timer->barrier_stop(TIME_LOOP);
}
