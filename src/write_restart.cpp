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
#include "write_restart.h"
#include "universe.h"
#include "simulator.h"
#include "particle.h"
#include "memory.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

WriteRestart::WriteRestart()
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
   called as write_restart command in input script
------------------------------------------------------------------------- */

void WriteRestart::command(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal write_restart command");
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot write restart files for non-spatial simulation");

  // insure simulation is fully setup

  simulator->init();
  write(arg[0]);
}

/* ----------------------------------------------------------------------
   called from output within timestep loop
------------------------------------------------------------------------- */

void WriteRestart::write(char *file)
{
  // nparts = sum of nlocal = value to write into restart file

  int nlocal = particle->nlocal;
  MPI_Allreduce(&nlocal,&nparts,1,MPI_INT,MPI_SUM,world);

  // open restart file

  if (me == 0) {
    fp = fopen(file,"wb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",file);
      error->one(str);
    }
  }

  // write header

  if (me == 0) header();

  // communication buffer for all my particle info
  // max_size = largest buffer needed by any proc

  int max_size;
  int send_size = nlocal * particle->size_restart;
  MPI_Allreduce(&send_size,&max_size,1,MPI_INT,MPI_MAX,world);

  double *buf;
  if (me == 0) 
    buf = (double *) 
      memory->smalloc(max_size*sizeof(double),"write_restart:buf");
  else
    buf = (double *) 
      memory->smalloc(send_size*sizeof(double),"write_restart:buf");

  // pack my particles into buf

  int n = particle->pack_restart(buf);

  // write one chunk of particles per proc to file
  // proc 0 pings each proc, receives its data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp,recv_size;
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buf,max_size,MPI_DOUBLE,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&recv_size);
      } else recv_size = send_size;

      fwrite(&recv_size,sizeof(int),1,fp);
      fwrite(buf,sizeof(double),recv_size,fp);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,send_size,MPI_DOUBLE,0,0,world);
  }

  memory->sfree(buf);
  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out problem description 
------------------------------------------------------------------------- */

void WriteRestart::header()
{
  write_char(0,universe->version);
  write_char(1,simulator->style);
  write_int(2,simulator->ntimestep);
  write_int(3,nprocs);
  write_int(4,nparts);

  // -1 flag signals end of header

  int flag = -1;
  fwrite(&flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   write_int - write a flag and an int into restart file 
------------------------------------------------------------------------- */

void WriteRestart::write_int(int flag, int value)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&value,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   write_double - write a flag and a double into restart file 
------------------------------------------------------------------------- */

void WriteRestart::write_double(int flag, double value)
{
  fwrite(&flag,sizeof(int),1,fp);
  fwrite(&value,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
   write_char - write a flag and a char str into restart file 
------------------------------------------------------------------------- */

void WriteRestart::write_char(int flag, char *value)
{
  fwrite(&flag,sizeof(int),1,fp);
  int n = strlen(value) + 1;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(value,sizeof(char),n,fp);
}
