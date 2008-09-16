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
#include "fix_conc_random.h"
#include "simulator.h"
#include "particle.h"
#include "random.h"
#include "memory.h"
#include "error.h"

#define CHUNK 16
#define MAXLINE 128

/* ---------------------------------------------------------------------- */

FixConcRandom::FixConcRandom(int narg, char **arg) : Fix(narg, arg)
{
  if (narg < 4) error->all("Illegal fix conc/random command");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use fix conc/random with non-spatial simulations");

  nevery = atoi(arg[2]);

  // parse input species list

  int isp;
  int istart = 4;
  if (narg < istart+1) error->all("Illegal fix conc/random command");
  ninput = atoi(arg[istart]);
  if (narg < istart+1+ninput) error->all("Illegal fix conc command");

  inlist = new int[ninput];
  for (int iarg = istart+1; iarg <= istart+ninput; iarg++) {
    isp = particle->find(arg[iarg]);
    if (isp < 0) error->all("Invalid species ID in fix conc/random command");
    inlist[iarg-istart-1] = isp;
  }

  // parse output species list

  istart += ninput + 1;
  if (narg < istart+1) error->all("Illegal fix conc/random command");
  noutput = atoi(arg[istart]);
  if (narg < istart+1+noutput) error->all("Illegal fix conc/random command");

  noutput = atoi(arg[istart]);
  outlist = new int[noutput];
  for (int iarg = istart+1; iarg <= istart+noutput; iarg++) {
    isp = particle->find(arg[iarg]);
    if (isp < 0) error->all("Invalid species ID in fix conc/random command");
    outlist[iarg-istart-1] = isp;
  }

  // read and bcast time course file
  // store time/prob values in times and probs

  int me;
  MPI_Comm_rank(world,&me);
  int maxtimes;

  if (me == 0) {
    ntimes = 0;
    maxtimes = CHUNK;
    times = (double *) memory->smalloc(maxtimes*sizeof(double),
				       "fix_conc:times");
    probs = memory->create_2d_double_array(maxtimes,noutput,
					   "fix_conc:probs");
    char *line = new char[MAXLINE];
    char *word;

    FILE *fp = fopen(arg[3],"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix conc/random file %s",arg[3]);
      error->one(str);
    }

    while (fgets(line,MAXLINE,fp)) {
      if (ntimes == maxtimes) {
	maxtimes += CHUNK;
	times = (double *) memory->srealloc(times,maxtimes*sizeof(double),
					    "fix_conc:times");
	probs = memory->grow_2d_double_array(probs,maxtimes,noutput,
					     "fix_conc:probs");
      }
      word = strtok(line," \t\n");
      if (word == NULL) continue;
      if (word[0] == '#') continue;
      times[ntimes] = atof(word);
      for (int i = 0; i < noutput-1; i++) {
	word = strtok(NULL," \t\n");
	probs[ntimes][i] = atof(word);
      }
      probs[ntimes][noutput-1] = 1.0;
      ntimes++;
    }
    fclose(fp);
    delete [] line;
  }

  MPI_Bcast(&maxtimes,1,MPI_INT,0,world);
  if (me) {
    times = (double *) memory->smalloc(maxtimes*sizeof(double),
				       "fix_conc:times");
    probs = memory->create_2d_double_array(maxtimes,noutput,
					   "fix_conc:probs");
  }
  MPI_Bcast(&ntimes,1,MPI_INT,0,world);
  MPI_Bcast(times,ntimes,MPI_DOUBLE,0,world);
  MPI_Bcast(probs[0],ntimes*noutput,MPI_DOUBLE,0,world);

  spflag = NULL;
}

/* ---------------------------------------------------------------------- */

FixConcRandom::~FixConcRandom()
{
  delete [] inlist;
  delete [] outlist;
  if (spflag) delete [] spflag;

  memory->sfree(times);
  memory->destroy_2d_double_array(probs);
}

/* ---------------------------------------------------------------------- */

int FixConcRandom::setmask()
{
  int mask = 0;
  mask |= INITIAL;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixConcRandom::init()
{
  // setup spflag of input species based on current species list

  if (spflag) delete [] spflag;

  int n = particle->nspecies;
  spflag = new int[n];
  for (int i = 0; i < n; i++) spflag[i] = 0;
  for (int i = 0; i < ninput; i++) spflag[inlist[i]] = 1;

  itime = -1;
}

/* ---------------------------------------------------------------------- */

void FixConcRandom::initial()
{
  if (simulator->ntimestep % nevery) return;

  // itime = time index to use for setting species proportions
  // don't allow itime to increment to ntimes

  double ctime = simulator->firsttime +
    (simulator->ntimestep-simulator->firststep) * simulator->dt;
  while (itime < ntimes-1 && ctime >= times[itime+1]) itime++;
  if (itime < 0) return;

  // set particle species if it is flagged in splist
  // particle seed generate RNG for picking new species from prob dist
  // don't change seed since are not checking particles in global order

  Particle::OnePart *plist = particle->plist;
  int nlocal = particle->nlocal;
  int j,seed;
  double rnd;
  
  for (int i = 0; i < nlocal; i++) {
    if (spflag[plist[i].species]) {
      seed = plist[i].seed;
      rnd = random->move(&seed);
      for (j = 0; j < noutput; j++) if (rnd < probs[itime][j]) break;
      plist[i].species = outlist[j];
      }
  }
}
