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
#include "dump.h"
#include "domain.h"
#include "simulator.h"
#include "particle.h"
#include "surf.h"
#include "region.h"
#include "error.h"
#include "memory.h"

/* ---------------------------------------------------------------------- */

Dump::Dump(int narg, char **arg)
{
  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  if (me == 0) {
    fp = fopen(arg[2],"w");
    if (fp == NULL) error->one("Could not open dump file");
  }

  format = NULL;
  orient_flag = 0;
  maxbuf = 0;
  buf = NULL;

  // set spmask to 1 if that species is dumped
  // no species args: all species dumped
  // for each species arg (with wildcard): spwork = list of matches

  nsp = particle->nspecies;
  spmask = new int[nsp];
  int *spwork = new int[nsp];

  for (int i = 0; i < nsp; i++) spmask[i] = 0;

  if (narg == 3) {
    for (int i = 0; i < nsp; i++) spmask[i] = 1;
  } else {
    for (int iarg = 3; iarg < narg; iarg++) {
      int nmatch = particle->allmatch(arg[iarg],spwork);
      if (nmatch == 0) error->all("No particle species match dump species");
      for (int i = 0; i < nsp; i++) if (spwork[i]) spmask[i] = 1;
    }
  }

  delete [] spwork;
}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
  delete [] id;
  if (me == 0) fclose(fp);
  if (format) delete [] format;
  memory->sfree(buf);
  delete [] spmask;
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
  if (orient_flag == 0) size_one = 5;
  else size_one = 8;

  // default format depends on orient flag

  if (format) delete [] format;
  char *str;
  if (orient_flag == 0) str = "%d %d %g %g %g";
  else str = "%d %d %g %g %g %g %g %g";
  int n = strlen(str) + 3;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  // test that spmask is correct length

  if (nsp != particle->nspecies)
    error->all("Dump species count does not match particle species count");
}

/* ----------------------------------------------------------------------
   write a single snapshot to file
------------------------------------------------------------------------- */

void Dump::write()
{
  int iproc,send_size,recv_size;

  // grow communication buffer if necessary

  int nmax;
  MPI_Allreduce(&particle->nlocal,&nmax,1,MPI_INT,MPI_MAX,world);

  if (nmax*size_one > maxbuf) {
    memory->sfree(buf);
    maxbuf = nmax*size_one;
    buf = (double *) memory->smalloc(maxbuf*sizeof(double),"dump:buf");
  }

  // pack my data into buf

  send_size = pack();

  // ndump = # of quantities in this dump

  int n = send_size/size_one;
  int ndump;
  MPI_Allreduce(&n,&ndump,1,MPI_INT,MPI_SUM,world);

  // proc 0 writes header to file

  if (me == 0) write_header(ndump);
  
  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp;
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    for (iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buf,maxbuf,MPI_DOUBLE,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&recv_size);
	recv_size /= size_one;
      } else recv_size = send_size/size_one;

      write_data(recv_size,buf);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,send_size,MPI_DOUBLE,0,0,world);
  }
}

/* ----------------------------------------------------------------------
   write dump header to file, only called by proc 0
------------------------------------------------------------------------- */

void Dump::write_header(int ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",simulator->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS\n");
  fprintf(fp,"%g %g\n",domain->xlo,domain->xhi);
  fprintf(fp,"%g %g\n",domain->ylo,domain->yhi);
  fprintf(fp,"%g %g\n",domain->zlo,domain->zhi);
  fprintf(fp,"ITEM: ATOMS\n");
}

/* ----------------------------------------------------------------------
   pack a proc's particles into buf
   convert species to 1-Nsp
------------------------------------------------------------------------- */

int Dump::pack()
{
  Particle::OnePart *plist = particle->plist;
  int *dimension = particle->dimension;
  Surf::OneTri *tlist = surf->tlist;
  Region **rlist = surf->rlist;
  double normal[3];
  int nlocal = particle->nlocal;

  int m = 0;

  // pack one of 2 ways depending on orient_flag

  if (orient_flag == 0) {

    for (int i = 0; i < nlocal; i++) {
      if (spmask[plist[i].species]) {
	buf[m++] = plist[i].seed;
	buf[m++] = plist[i].species+1;
	buf[m++] = plist[i].x[0];
	buf[m++] = plist[i].x[1];
	buf[m++] = plist[i].x[2];
      }
    }

  } else {

    int itri;
    for (int i = 0; i < nlocal; i++) {
      if (spmask[plist[i].species]) {
	buf[m++] = plist[i].seed;
	buf[m++] = plist[i].species+1;
	buf[m++] = plist[i].x[0];
	buf[m++] = plist[i].x[1];
	buf[m++] = plist[i].x[2];
	if (dimension[plist[i].species] == 2) {
	  itri = plist[i].itri;
	  if (itri >= 0) {
	    buf[m++] = tlist[itri].normal[0];
	    buf[m++] = tlist[itri].normal[1];
	    buf[m++] = tlist[itri].normal[2];
	  } else {
	    rlist[-itri-1]->compute_normal(plist[i].x,normal);
	    buf[m++] = normal[0];
	    buf[m++] = normal[1];
	    buf[m++] = normal[2];
	  }
	} else {
	  buf[m++] = 0.0;
	  buf[m++] = 0.0;
	  buf[m++] = 0.0;
	}
      }
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   write a chunk of particles from buf into file
------------------------------------------------------------------------- */

void Dump::write_data(int n, double *buf)
{
  int m = 0;

  if (orient_flag == 0) {
    for (int i = 0; i < n; i++) {
      fprintf(fp,format,
	      static_cast<int> (buf[m]), static_cast<int> (buf[m+1]),
	      buf[m+2],buf[m+3],buf[m+4]);
      m += size_one;
    }
  } else {
    for (int i = 0; i < n; i++) {
      fprintf(fp,format,
	      static_cast<int> (buf[m]), static_cast<int> (buf[m+1]),
	      buf[m+2],buf[m+3],buf[m+4],buf[m+5],buf[m+6],buf[m+7]);
      m += size_one;
    }
  }
}

/* ----------------------------------------------------------------------
   modify dump parameters
------------------------------------------------------------------------- */

void Dump::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"orient") == 0) {
      if (iarg+1 >= narg) error->all("Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) orient_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) orient_flag = 0;
      else error->all("Illegal dump_modify command");
      iarg += 2;
    } else error->all("Illegal dump_modify command");
  }
}
