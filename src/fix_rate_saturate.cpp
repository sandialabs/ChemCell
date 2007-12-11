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

#include "string.h"
#include "stdlib.h"
#include "fix_rate_saturate.h"
#include "simulator.h"
#include "particle.h"
#include "random.h"
#include "memory.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

FixRateSaturate::FixRateSaturate(int narg, char **arg) : Fix(narg, arg)
{
  if (narg < 4) error->all("Illegal fix rate/saturate command");
  if (simulator->spatial_flag == 1)
    error->all("Cannot use fix rate/saturate with spatial simulations");

  nevery = atoi(arg[2]);
}

/* ---------------------------------------------------------------------- */

FixRateSaturate::~FixRateSaturate()
{
}

/* ---------------------------------------------------------------------- */

int FixRateSaturate::setmask()
{
  int mask = 0;
  mask |= INITIAL;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRateSaturate::init()
{
}

/* ---------------------------------------------------------------------- */

void FixRateSaturate::initial()
{
  if (simulator->ntimestep % nevery) return;
}
