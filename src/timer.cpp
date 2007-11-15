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

#include "mpi.h"
#include "timer.h"
#include "memory.h"

/* ---------------------------------------------------------------------- */

Timer::Timer()
{
  array = (double *) memory->smalloc(TIME_N*sizeof(double),"array");
}

/* ---------------------------------------------------------------------- */

Timer::~Timer()
{
  memory->sfree(array);
}

/* ---------------------------------------------------------------------- */

void Timer::init()
{
  for (int i = 0; i < TIME_N; i++) array[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void Timer::stamp()
{
  previous_time = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::stamp(int which)
{
  double current_time = MPI_Wtime();
  array[which] += current_time - previous_time;
  previous_time = current_time;
}

/* ---------------------------------------------------------------------- */

void Timer::sub_stamp()
{
  previous_time_sub = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::sub_stamp(int which)
{
  double current_time = MPI_Wtime();
  array[which] += current_time - previous_time_sub;
  previous_time_sub = current_time;
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_start(int which)
{
  MPI_Barrier(world);
  array[which] = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_stop(int which)
{
  MPI_Barrier(world);
  double current_time = MPI_Wtime();
  array[which] = current_time - array[which];
}

/* ---------------------------------------------------------------------- */

double Timer::elapsed(int which)
{
  double current_time = MPI_Wtime();
  return (current_time - array[which]);
}
