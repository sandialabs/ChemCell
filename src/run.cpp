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

#include "string.h"
#include "run.h"
#include "simulator.h"
#include "finish.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

void Run::command(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal run command");
  if (!simulator) error->all("Must set run_style first");

  if (strchr(arg[0],'.') == NULL) {
    simulator->nsteps = atoi(arg[0]);
    simulator->elapse = -1.0;
  } else {
    simulator->nsteps = -1;
    simulator->elapse = atof(arg[0]);
  }

  simulator->init();
  simulator->run();
  Finish finish;
  finish.end();
}
