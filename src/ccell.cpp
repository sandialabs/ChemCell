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

// ChemCell main driver

#include "string.h"
#include "mpi.h"
#include "system.h"
#include "input.h"

#define CommandInclude
#include "style.h"
#undef CommandInclude

/* ---------------------------------------------------------------------- */

int main (int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  System sys;
  sys.open(argc,argv,MPI_COMM_WORLD);
  sys.create();
  Input *input = sys.input;

  // loop over input script commands
  // next() returns when there is a top-level command to execute

  char *command;
  while (command = input->next()) {

    if (0) break;      // dummy line to enable else-if macro expansion

#define CommandClass
#define CommandStyle(key,Class) \
    else if (strcmp(command,#key) == 0) { \
      Class key; \
      key.command(input->narg,input->arg); \
    }
#include "style.h"
#undef CommandClass

  }
  
  sys.destroy();
  sys.close();
  MPI_Finalize();
}
