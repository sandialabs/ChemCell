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
#include "mpi.h"
#include "move_test.h"
#include "simulator.h"
#include "particle.h"
#include "move.h"
#include "error.h"

#define EPSILON 1.0e-6
#define PI 3.1415926

/* ---------------------------------------------------------------------- */

void MoveTest::command(int narg, char **arg)
{
  if (narg != 5) error->all("Illegal move_test command");
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use move_test command with non-spatial simulation");
  if (simulator->dt == 0.0)
    error->all("Cannot do move_test with timestep = 0");

  // parse command args

  int ispecies = particle->find(arg[0]);
  if (ispecies == -1) error->all("Unknown species in move_test command");

  int ntest = atoi(arg[1]);
  int nhisto = atoi(arg[2]);
  if (ntest <= 0 || nhisto <= 0) error->all("Illegal move_test command");

  int seed = atoi(arg[3]);
  char *file = arg[4];

  // initialize the simulation

  simulator->init();

  // only proc 0 performs a move test

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) return;

  // create histogram vec and bin size
  // maxdist = maximum distance a particle can move

  int *histo = new int[nhisto];
  for (int i = 0; i < nhisto; i++) histo[i] = 0;
  double maxdist = move->maxmove[ispecies];
  double binsize = (maxdist + EPSILON) / nhisto;
  double invbin = 1.0/binsize;

  // ideal = Brownian average distance for 3d or 2d for this species

  double ideal;
  if (particle->dimension[ispecies] == 3)
    ideal = 4.0/sqrt(PI) * 
      sqrt(1.0e8*particle->diffusivity[ispecies] * simulator->dt);
  else 
    ideal = sqrt(PI) * 
      sqrt(1.0e8*particle->diffusivity[ispecies] * simulator->dt);

  // loop over N moves for 3d or 2d particle

  double delta[3],dir[3],normal[3],distance,dave;
  int ihisto;
  normal[0] = 1.0;
  normal[1] = normal[2] = 0.0;

  double time1 = MPI_Wtime();

  if (particle->dimension[ispecies] == 3) {
    for (int i = 0; i < ntest; i++) {
      (move->*(move->delta3d))(ispecies,&seed,delta);
      distance = 
	sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
      dave += distance;
      ihisto = static_cast<int> (distance*invbin);
      histo[ihisto]++;
    }

  } else {
    for (int i = 0; i < ntest; i++) {
      (move->*(move->delta2d))(ispecies,&seed,normal,&distance,dir);
      dave += distance;
      ihisto = static_cast<int> (distance*invbin);
      histo[ihisto]++;
    }
  }

  double time2 = MPI_Wtime();

  if (screen) {
    fprintf(screen,"Move test of species %s:\n",particle->name[ispecies]);
    fprintf(screen,"Dimension = %d\n",particle->dimension[ispecies]);
    fprintf(screen,"Diff coeff = %g\n",particle->diffusivity[ispecies]);
    fprintf(screen,"Timestep = %g\n",simulator->dt);
    fprintf(screen,"Ideal Ave distance = %g\n",ideal);
    fprintf(screen,"Ave distance moved = %g\n",dave/ntest);
    fprintf(screen,"Max distance moved = %g\n",maxdist);
    fprintf(screen,"Number of moves = %d\n",ntest);
    fprintf(screen,"CPU Time = %g\n",time2-time1);
  }

  if (logfile) {
    fprintf(logfile,"Move test of species %s:\n",particle->name[ispecies]);
    fprintf(logfile,"Dimension = %d\n",particle->dimension[ispecies]);
    fprintf(logfile,"Diff coeff = %g\n",particle->diffusivity[ispecies]);
    fprintf(logfile,"Timestep = %g\n",simulator->dt);
    fprintf(logfile,"Ideal Ave distance = %g\n",ideal);
    fprintf(logfile,"Ave distance moved = %g\n",dave/ntest);
    fprintf(logfile,"Max distance moved = %g\n",maxdist);
    fprintf(logfile,"Number of moves = %d\n",ntest);
    fprintf(logfile,"CPU Time = %g\n",time2-time1);
  }

  // write histogram distribution into file

  FILE *fp = fopen(file,"w");
  if (fp == NULL) error->one("Could not open move_test file");

  int count;
  for (int i = 0; i < nhisto; i++) {
    count += histo[i];
    fprintf(fp,"%g %d\n",(i+0.5)*binsize,histo[i]);
  }

  fclose(fp);

  // free local memory

  delete [] histo;
}
