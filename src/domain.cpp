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
#include "string.h"
#include "stdlib.h"
#include "mpi.h"
#include "domain.h"
#include "grid.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

Domain::Domain()
{
  setflag = 0;
  xperiodic = yperiodic = zperiodic = 1;
}

/* ---------------------------------------------------------------------- */

void Domain::init()
{
  if (setflag == 0)
    error->all("Must set simulation domain via global command");

  int me;
  MPI_Comm_rank(world,&me);

  if (me == 0) {
    fprintf(screen,"Domain:\n");
    fprintf(screen,"  global size in xyz = %g %g %g\n",xsize,ysize,zsize);
    fprintf(logfile,"Domain:\n");
    fprintf(logfile,"  global size in xyz = %g %g %g\n",xsize,ysize,zsize);
  }
}

/* ----------------------------------------------------------------------
   set bounds of global simulation domain
------------------------------------------------------------------------- */

void Domain::global(int narg, char **arg)
{
  if (setflag == 1) error->all("Simulation domain is already set");
  if (grid->setflag == 1) error->all("Bins are already setup");
  if (narg != 6) error->all("Illegal global command");
  setflag = 1;

  xlo = atof(arg[0]);
  ylo = atof(arg[1]);
  zlo = atof(arg[2]);
  xhi = atof(arg[3]);
  yhi = atof(arg[4]);
  zhi = atof(arg[5]);
  if (xlo >= xhi || ylo >= yhi || zlo >= zhi)
    error->all("Illegal global command");

  xsize = xhi - xlo;
  ysize = yhi - ylo;
  zsize = zhi - zlo;
}

/* ----------------------------------------------------------------------
   set periodicity of global simulation domain
------------------------------------------------------------------------- */

void Domain::boundary(int narg, char **arg)
{
  if (narg != 3) error->all("Illegal boundary command");
  if (grid->setflag == 1) error->all("Bins are already setup");

  if (strcmp(arg[0],"p") == 0) xperiodic = 1;
  else if (strcmp(arg[0],"n") == 0) xperiodic = 0;
  else error->all("Illegal boundary command");

  if (strcmp(arg[1],"p") == 0) yperiodic = 1;
  else if (strcmp(arg[1],"n") == 0) yperiodic = 0;
  else error->all("Illegal boundary command");

  if (strcmp(arg[2],"p") == 0) zperiodic = 1;
  else if (strcmp(arg[2],"n") == 0) zperiodic = 0;
  else error->all("Illegal boundary command");
}

/* ----------------------------------------------------------------------
   return 1 if pt is inside global domain, 0 if not
   pt on domain boundary is considered inside
------------------------------------------------------------------------- */

int Domain::inside(double *x)
{
  if (x[0] < xlo || x[0] > xhi) return 0;
  if (x[1] < ylo || x[1] > yhi) return 0;
  if (x[2] < zlo || x[2] > zhi) return 0;
  return 1;
}
