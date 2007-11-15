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
#include "simulator.h"

/* ---------------------------------------------------------------------- */

Simulator::Simulator(int narg, char **arg)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  ntimestep = 0;
  ctime = 0.0;
  dt = 0.0;
}

/* ---------------------------------------------------------------------- */

Simulator::~Simulator()
{
  delete [] style;
}
