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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "read_restart.h"
#include "universe.h"
#include "simulator.h"
#include "particle.h"
#include "grid.h"
#include "memory.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

void ReadRestart::command(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal read_restart command");
  if (!simulator) error->all("Must set run_style first");
  if (grid->setflag == 0) error->all("Must set bins before read restart");

  // open restart file

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  if (me == 0) {
    if (screen) fprintf(screen,"Reading restart file ...\n");
    fp = fopen(arg[0],"rb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",arg[0]);
      error->one(str);
    }
  }

  // read header info

  header();

  // proc 0 reads particles from file
  // bcast to other procs
  // atoms are in chunks, one chunk per proc when restart file was written
  // each proc unpacks the particles, saving ones in it's sub-domain

  int maxbuf = 0;
  double *buf = NULL;
  int n;

  for (int iproc = 0; iproc < nprocs_file; iproc++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n > maxbuf) {
      maxbuf = n;
      delete [] buf;
      buf = new double[maxbuf];
    }

    if (me == 0) fread(buf,sizeof(double),n,fp);
    MPI_Bcast(buf,n,MPI_DOUBLE,0,world);

    particle->unpack_restart(n,buf);
  }

  delete [] buf;

  // close file

  if (me == 0) fclose(fp);

  // check that all particles were assigned to procs

  int ntotal;
  int nlocal = particle->nlocal;
  MPI_Allreduce(&nlocal,&ntotal,1,MPI_INT,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  %d particles\n",ntotal);
    if (logfile) fprintf(logfile,"  %d particles\n",ntotal);
  }

  if (ntotal != nparts) error->all("Did not assign all particles correctly");
}

/* ----------------------------------------------------------------------
   read header of restart file 
------------------------------------------------------------------------- */

void ReadRestart::header()
{
  // read flags and values until flag = -1

  int flag = read_int();
  while (flag >= 0) {

    // check restart file version, compare to ChemCell version

    if (flag == 0) {
      char *version = read_char();
      if (strcmp(version,universe->version) != 0) {
	error->warning("Restart file version does not match ChemCell version");
	if (screen) fprintf(screen,"  restart file = %s, ChemCell = %s\n",
			    version,universe->version);
      }
      delete [] version;

    // check run_style, warn if different

    } else if (flag == 1) {
      char *style = read_char();
      if (strcmp(style,simulator->style) != 0 && me == 0)
	error->warning("Run styles do not match");
      delete [] style;

    } else if (flag == 2) {
      simulator->ntimestep = read_int();

    // read nprocs_file from restart file, warn if different

    } else if (flag == 3) {
      nprocs_file = read_int();
      if (nprocs_file != nprocs && me == 0)
	error->warning("Restart file used different # of processors");

    } else if (flag == 4) {
      nparts = read_int();

    } else error->all("Invalid flag in header of restart file");

    flag = read_int();
  }
}

/* ----------------------------------------------------------------------
   read an int from restart file 
------------------------------------------------------------------------- */

int ReadRestart::read_int()
{
  int value;
  if (me == 0) fread(&value,sizeof(int),1,fp);
  MPI_Bcast(&value,1,MPI_INT,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a double from restart file 
------------------------------------------------------------------------- */

double ReadRestart::read_double()
{
  double value;
  if (me == 0) fread(&value,sizeof(double),1,fp);
  MPI_Bcast(&value,1,MPI_DOUBLE,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a char str from restart file 
------------------------------------------------------------------------- */

char *ReadRestart::read_char()
{
  int n;
  if (me == 0) fread(&n,sizeof(int),1,fp);
  MPI_Bcast(&n,1,MPI_INT,0,world);
  char *value = new char[n];
  if (me == 0) fread(value,sizeof(char),n,fp);
  MPI_Bcast(value,n,MPI_CHAR,0,world);
  return value;
}
