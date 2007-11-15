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

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "system.h"

class Simulator : public System {
 public:
  char *style;                    // style of simulation

  int nsteps;                     // # of steps to run
  double elapse;                  // amount of time to run
  double dt;                      // timestep size

  int ntimestep;                  // current timestep
  double ctime;                   // current time
  int firststep,laststep;         // 1st & last timestep of run
  double firsttime,lasttime;      // initial and final time of run

  int spatial_flag;               // 0 if non-spatial, 1 if spatial
  int stochastic_flag;            // 0 if deterministic (ODE), 1 if stochastic

  Simulator(int, char **);
  virtual ~Simulator();
  virtual void init() = 0;
  virtual void run() = 0;
};

#endif
