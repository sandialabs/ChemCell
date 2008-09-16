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
#include "stdlib.h"
#include "fix_conc_set.h"
#include "simulator.h"
#include "particle.h"
#include "random.h"
#include "memory.h"
#include "error.h"

#define CHUNK 16
#define MAXLINE 128

/* ---------------------------------------------------------------------- */

FixConcSet::FixConcSet(int narg, char **arg) : Fix(narg, arg)
{
  if (narg < 4) error->all("Illegal fix conc/set command");
  if (simulator->spatial_flag == 1)
    error->all("Cannot use fix conc/set with spatial simulations");

  nevery = atoi(arg[2]);

  // parse species list

  int isp;
  int istart = 4;
  if (narg < istart+1) error->all("Illegal fix conc/set command");
  nsp = atoi(arg[istart]);
  if (narg < istart+1+nsp) error->all("Illegal fix conc/set command");

  list = new int[nsp];
  for (int iarg = istart+1; iarg <= istart+nsp; iarg++) {
    isp = particle->find(arg[iarg]);
    if (isp < 0) error->all("Invalid species ID in fix conc/set command");
    list[iarg-istart-1] = isp;
  }

  // read and bcast time course file
  // store time/values in times and values

  int me;
  MPI_Comm_rank(world,&me);
  int maxtimes;

  if (me == 0) {
    ntimes = 0;
    maxtimes = CHUNK;
    times = (double *) memory->smalloc(maxtimes*sizeof(double),
				       "fix_conc:times");
    values = memory->create_2d_double_array(maxtimes,nsp,"fix_conc:values");
    char *line = new char[MAXLINE];
    char *word;

    FILE *fp = fopen(arg[3],"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix conc/set file %s",arg[3]);
      error->one(str);
    }

    while (fgets(line,MAXLINE,fp)) {
      if (ntimes == maxtimes) {
	maxtimes += CHUNK;
	times = (double *) memory->srealloc(times,maxtimes*sizeof(double),
					    "fix_conc:times");
	values = memory->grow_2d_double_array(values,maxtimes,nsp,
					      "fix_conc:values");
      }
      word = strtok(line," \t\n");
      if (word == NULL) continue;
      if (word[0] == '#') continue;
      times[ntimes] = atof(word);
      for (int i = 0; i < nsp; i++) {
	word = strtok(NULL," \t\n");
	values[ntimes][i] = atof(word);
      }
      ntimes++;
    }
    fclose(fp);
    delete [] line;
  }

  MPI_Bcast(&maxtimes,1,MPI_INT,0,world);
  if (me) {
    times = (double *) memory->smalloc(maxtimes*sizeof(double),
				       "fix_conc:times");
    values = memory->create_2d_double_array(maxtimes,nsp,"fix_conc:values");
  }
  MPI_Bcast(&ntimes,1,MPI_INT,0,world);
  MPI_Bcast(times,ntimes,MPI_DOUBLE,0,world);
  MPI_Bcast(values[0],ntimes*nsp,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

FixConcSet::~FixConcSet()
{
  delete [] list;

  memory->sfree(times);
  memory->destroy_2d_double_array(values);
}

/* ---------------------------------------------------------------------- */

int FixConcSet::setmask()
{
  int mask = 0;
  mask |= INITIAL;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixConcSet::init()
{
  itime = -1;
}

/* ---------------------------------------------------------------------- */

void FixConcSet::initial()
{
  if (simulator->ntimestep % nevery) return;

  // itime = time index to use for setting species concentrations
  // don't allow itime to increment to ntimes

  double ctime = simulator->firsttime +
    (simulator->ntimestep-simulator->firststep) * simulator->dt;
  while (itime < ntimes-1 && ctime >= times[itime+1]) itime++;
  if (itime < 0) return;

  // reset pcount (Gillespie) or ccount (ODE) for species in list

  int *pcount = particle->pcount;
  double *ccount = particle->ccount;
  int nspecies = particle->nspecies;

  int ispecies;
  for (int isp = 0; isp < nsp; isp++) {
    ispecies = list[isp];
    if (simulator->stochastic_flag)
      pcount[ispecies] = static_cast<int> (values[itime][isp]);
    else
      ccount[ispecies] = values[itime][isp];
  }
}
