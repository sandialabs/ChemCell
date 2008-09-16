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
#include "mpi.h"
#include "system.h"
#include "input.h"

/* ----------------------------------------------------------------------
   main program to drive LAMMPS
------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  System sys;
  sys.open(argc,argv,MPI_COMM_WORLD);
  sys.create();

  Input *input = sys.input;
  input->file();

  sys.destroy();
  sys.close();
  MPI_Finalize();
}
