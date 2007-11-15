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
#include "chem.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

Chem::Chem()
{
  prob_style = 0;
  max_prob = 0.5;
  rcount = NULL;
}

Chem::~Chem()
{
  delete [] style;
}

/* ----------------------------------------------------------------------
   set probablility
------------------------------------------------------------------------- */

void Chem::set_prob(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal probability command");

  if (strcmp(arg[0],"max") == 0) {
    prob_style = 0;
    max_prob = atof(arg[1]);
  } else if (strcmp(arg[0],"diff") == 0) {
    prob_style = 1;
    diff_factor = atof(arg[1]);
  } else error->all("Illegal probability command");
}

