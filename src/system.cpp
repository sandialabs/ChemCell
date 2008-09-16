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

// parent class for all of ChemCell

#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "system.h"
#include "simulator.h"
#include "universe.h"
#include "input.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "random.h"
#include "timer.h"

// set all static ptrs in System to NULL

Memory *System::memory = NULL;
Error *System::error = NULL;
Universe *System::universe = NULL;
Input *System::input = NULL;

Balance *System::balance = NULL;
Chem *System::chem = NULL;
Domain *System::domain = NULL;
Grid *System::grid = NULL;
Modify *System::modify = NULL;
Move *System::move = NULL;
Output *System::output = NULL;
Particle *System::particle = NULL;
Random *System::random = NULL;
React *System::react = NULL;
Simulator *System::simulator = NULL;
Surf *System::surf = NULL;
Timer *System::timer = NULL;

MPI_Comm System::world = 0;
FILE *System::infile = NULL;
FILE *System::screen = NULL;
FILE *System::logfile = NULL;

/* ----------------------------------------------------------------------
   allocate fundamental classes (memory, error, universe, input)
   parse input switches
   initialize communicators, screen & logfile output
   input is allocated at end after MPI info is setup
------------------------------------------------------------------------- */

void System::open(int narg, char **arg, MPI_Comm communicator)
{
  memory = new Memory;
  error = new Error;
  universe = new Universe(communicator);

  // parse input switches

  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;
  int iarg = 1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"-partition") == 0) {
      if (iarg+1 > narg) 
	error->universe_all("Invalid command-line argument");
      iarg++;
      while (iarg < narg && arg[iarg][0] != '-') {
	universe->add_world(arg[iarg]);
	iarg++;
      }
    } else if (strcmp(arg[iarg],"-in") == 0) {
      if (iarg+2 > narg) 
	error->universe_all("Invalid command-line argument");
      inflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-screen") == 0) {
      if (iarg+2 > narg) 
	error->universe_all("Invalid command-line argument");
      screenflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-log") == 0) {
      if (iarg+2 > narg) 
	error->universe_all("Invalid command-line argument");
      logflag = iarg + 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"-var") == 0) {
      if (iarg+3 > narg) 
	error->universe_all("Invalid command-line argument");
      iarg += 3;
    } else if (strcmp(arg[iarg],"-echo") == 0) {
      if (iarg+2 > narg) 
	error->universe_all("Invalid command-line argument");
      iarg += 2;
    } else error->universe_all("Invalid command-line argument");
  }

  // if procs was not a command-line switch, universe is one world w/ all procs

  if (universe->nworlds == 0) universe->add_world(NULL);

  // sum of procs in all worlds must equal total # of procs

  if (!universe->consistent())
    error->universe_all("Processor partitions are inconsistent");

  // multiple-world universe must define input file

  if (universe->nworlds > 1 && inflag == 0)
    error->universe_all("Must use -in switch with multiple partitions");

  // set universe screen and logfile

  if (universe->me == 0) {
    if (screenflag == 0)
      universe->uscreen = stdout;
    else if (strcmp(arg[screenflag],"none") == 0)
      universe->uscreen = NULL;
    else {
      universe->uscreen = fopen(arg[screenflag],"w");
      if (universe->uscreen == NULL) 
	error->universe_one("Could not open universe screen file");
    }
    if (logflag == 0) {
      universe->ulogfile = fopen("log.ccell","w");
      if (universe->ulogfile == NULL) 
	error->universe_one("Could not open log.ccell");
    } else if (strcmp(arg[logflag],"none") == 0)
      universe->ulogfile = NULL;
    else {
      universe->ulogfile = fopen(arg[logflag],"w");
      if (universe->ulogfile == NULL) 
	error->universe_one("Could not open universe log file");
    }
  }

  if (universe->me > 0) {
    if (screenflag == 0) universe->uscreen = stdout;
    else universe->uscreen = NULL;
    universe->ulogfile = NULL;
  }

  // universe is single world
  // inherit settings from universe
  // set world screen, logfile, communicator, infile
  // open input script if from file

  if (universe->nworlds == 1) {
    screen = universe->uscreen;
    logfile = universe->ulogfile;
    world = universe->uworld;
    infile = NULL;

    if (universe->me == 0) {
      if (inflag == 0) infile = stdin;
      else infile = fopen(arg[inflag],"r");
      if (infile == NULL) error->one("Could not open input script");
    }

    if (universe->me == 0) {
      if (screen) fprintf(screen,"ChemCell (%s)\n",universe->version);
      if (logfile) fprintf(logfile,"ChemCell (%s)\n",universe->version);
    }

  // universe is multiple worlds
  // split into separate communicators
  // set world screen, logfile, communicator, infile
  // open input script

  } else {
    int me;
    MPI_Comm_split(universe->uworld,universe->iworld,0,&world);
    MPI_Comm_rank(world,&me);

    if (me == 0) {
      if (screenflag == 0) {
	char str[32];
	sprintf(str,"screen.%d",universe->iworld);
	screen = fopen(str,"w");
	if (screen == NULL) error->one("Could not open screen file");
      } else if (strcmp(arg[screenflag],"none") == 0)
	screen = NULL;
      else {
	char str[32];
	sprintf(str,"%s.%d",arg[screenflag],universe->iworld);
	screen = fopen(str,"w");
	if (screen == NULL) error->one("Could not open screen file");
      }
    } else screen = NULL;
    
    if (me == 0) {
      if (logflag == 0) {
	char str[32];
	sprintf(str,"log.ccell.%d",universe->iworld);
	logfile = fopen(str,"w");
	if (logfile == NULL) error->one("Could not open logfile");
      } else if (strcmp(arg[logflag],"none") == 0)
	logfile = NULL;
      else {
	char str[32];
	sprintf(str,"%s.%d",arg[logflag],universe->iworld);
	logfile = fopen(str,"w");
	if (logfile == NULL) error->one("Could not open logfile");
      }
    } else logfile = NULL;
    
    if (me == 0) {
      infile = fopen(arg[inflag],"r");
      if (infile == NULL) error->one("Could not open input script");
    } else infile = NULL;
    
    // screen and logfile messages for universe and world
    
    if (universe->me == 0) {
      if (universe->uscreen) {
	fprintf(universe->uscreen,"ChemCell (%s)\n",universe->version);
	fprintf(universe->uscreen,"Running on %d partitions of processors\n",
		universe->nworlds);
      }
      if (universe->ulogfile) {
	fprintf(universe->ulogfile,"ChemCell (%s)\n",universe->version);
	fprintf(universe->ulogfile,"Running on %d partitions of processors\n",
		universe->nworlds);
      }
    }

    if (me == 0) {
      if (screen) {
	fprintf(screen,"ChemCell (%s)\n",universe->version);
	fprintf(screen,"Processor partition = %d\n",universe->iworld);
      }
      if (logfile) {
	fprintf(logfile,"ChemCell (%s)\n",universe->version);
	fprintf(logfile,"Processor partition = %d\n",universe->iworld);
      }
    }
  }

  // allocate input class now that MPI is fully setup

  input = new Input(narg,arg);
}

/* ----------------------------------------------------------------------
   allocate single instance of top-level classes
   fundamental classes are allocated in open()
   simulator is allocated via input script, it allocates other classes
------------------------------------------------------------------------- */

void System::create()
{
  output = new Output;
  random = new Random;
  timer = new Timer;
}

/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in close()
   simulator deletes other classes
------------------------------------------------------------------------- */

void System::destroy()
{
  if (simulator) delete simulator;
  simulator = NULL;

  delete output;
  delete random;
  delete timer;
}

/* ----------------------------------------------------------------------
   shutdown system
   close screen and log files in world and universe
   output files were already closed in system::destroy()
   delete fundamental classes
------------------------------------------------------------------------- */

void System::close()
{
  if (universe->nworlds > 1) {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
  }
  if (universe->ulogfile) fclose(universe->ulogfile);

  if (world != universe->uworld) MPI_Comm_free(&world);

  delete input;
  delete universe;
  delete error;
  delete memory;
}
