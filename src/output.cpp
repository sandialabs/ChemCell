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

#include "stdio.h"
#include "string.h"
#include "mpi.h"
#include "output.h"
#include "simulator.h"
#include "particle.h"
#include "dump.h"
#include "write_restart.h"
#include "memory.h"
#include "error.h"

#define DELTA 1

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

enum {COUNT,MOLARITY,UMOLARITY,NMOLARITY};    // same as in particle.cpp

/* ---------------------------------------------------------------------- */

Output::Output() {
  stats_every = 0;
  stats_delta = -1.0;
  nsp_stats = 0;
  spname = (char **) memory->smalloc(nsp_stats*sizeof(char *),"output:spname");
  stat2species = NULL;
  stats_units = COUNT;
  format = NULL;
  format_user = NULL;

  ndump = 0;
  maxdump = 0;
  next_dump = NULL;
  last_dump = NULL;
  dump_every = NULL;
  dump = NULL;

  restart = NULL;
  restart1 = restart2 = NULL;
  restart_every = 0;

  MPI_Comm_rank(world,&me);
}

/* ---------------------------------------------------------------------- */

Output::~Output()
{
  for (int isp = 0; isp < nsp_stats; isp++) delete [] spname[isp];
  memory->sfree(spname);

  if (stat2species) delete [] stat2species;
  delete [] format;
  delete [] format_user;

  memory->sfree(next_dump);
  memory->sfree(last_dump);
  memory->sfree(dump_every);
  for (int i = 0; i < ndump; i++) delete dump[i];
  memory->sfree(dump);

  delete restart;
  delete [] restart1;
  delete [] restart2;
}

/* ----------------------------------------------------------------------
   init output counters
------------------------------------------------------------------------- */

void Output::init()
{
  if (simulator->spatial_flag && stats_units != COUNT)
    error->all("Must print stats in count units for spatial simulation");

  if (format) delete [] format;
  if (format_user == NULL && stats_units == COUNT) {
    char *str = " %d";
    int n = strlen(str) + 1;
    format = new char[n];
    strcpy(format,str);
  } else if (format_user == NULL) {
    char *str = " %g";
    int n = strlen(str) + 1;
    format = new char[n];
    strcpy(format,str);
  } else {
    int n = strlen(format_user) + 2;
    format = new char[n];
    strcpy(format," ");
    strcat(format,format_user);
  }

  int ntimestep = simulator->ntimestep;

  // init species list for stats

  if (stat2species) delete [] stat2species;
  nstats = particle->nspecies;
  if (nsp_stats) nstats = nsp_stats;
  stat2species = new int[nstats];

  int isp,jsp;

  if (nsp_stats == 0) {
    for (isp = 0; isp < particle->nspecies; isp++) stat2species[isp] = isp;
  } else {
    for (isp = 0; isp < nsp_stats; isp++) {
      jsp = particle->find(spname[isp]);
      if (jsp == -1) {
	char str[128];
	sprintf(str,"Invalid species %s in stats command",spname[isp]);
	error->all(str);
      }
      stat2species[isp] = jsp;
    }
  }

  // always do stats at start of run
  // if stats_every = -1 and Gillespie (dt = 0)
  //   do nothing, Gillespie does it's own calls to stats()
  // if stats_every = -1 and non-Gillespie (dt > 0)
  //   use stats_delta to set next_stats
  // insure next stats is a multiple of every
  // always do stats on last timestep, even if every = 0

  if (me == 0) stats_header();
  stats();

  if (simulator->dt == 0.0) {
    if (stats_every == -1) next_stats = 0;
    else if (stats_every == 0) next_stats = 0;
    else {
      next_stats = ntimestep + stats_every;
      if (ntimestep % stats_every != 0)
	next_stats = (ntimestep/stats_every)*stats_every + stats_every;
    }
  } else {
    if (stats_every == -1)
      stats_every = static_cast<int> 
	((stats_delta+0.5*simulator->dt) / simulator->dt);
    next_stats = ntimestep + stats_every;
    if (stats_every == 0) next_stats = simulator->laststep;
    if (stats_every && ntimestep % stats_every != 0)
      next_stats = (ntimestep/stats_every)*stats_every + stats_every;
    next_stats = MIN(next_stats,simulator->laststep);
  }

  // init all dumps

  if (ndump && simulator->spatial_flag == 0)
    error->all("Cannot write dump files for non-spatial simulation");

  int idump;
  for (idump = 0; idump < ndump; idump++) dump[idump]->init();
  
  // dump at start of run only if dump has never been done before
  // set next dump timestep
  // insure next dump is a multiple of every

  for (int idump = 0; idump < ndump; idump++)
    if (last_dump[idump] == -1) {
      dump[idump]->write();
      last_dump[idump] = ntimestep;
      next_dump[idump] = ntimestep + dump_every[idump];
      if (ntimestep % dump_every[idump] != 0)
	next_dump[idump] = (ntimestep/dump_every[idump])*dump_every[idump] + 
	  dump_every[idump];
    }

  // do not write a restart file at start of run
  // set next restart timestep
  // if every = 0, set to last + 1
  // insure next restart is a multiple of every

  if (restart_every && simulator->spatial_flag == 0)
    error->all("Cannot write restart files for non-spatial simulation");

  next_restart = ntimestep + restart_every;
  if (restart_every == 0) next_restart = simulator->laststep + 1;
  if (restart_every && ntimestep % restart_every != 0)
    next_restart = (ntimestep/restart_every)*restart_every + restart_every;
  // next = next timestep any output will be done

  next = MIN(next_stats,next_restart);
  for (idump = 0; idump < ndump; idump++) next = MIN(next,next_dump[idump]);
}

/* ----------------------------------------------------------------------
   all output for this timestep
   Gillespie style never calls this routine, it does own stats() call
------------------------------------------------------------------------- */

void Output::write(int ntimestep)
{
  // do stats if next_stats = ntimestep and set next stats timestep
  // always do stats on last timestep of run

  if (next_stats == ntimestep) {
    simulator->ctime = simulator->firsttime + 
      (ntimestep-simulator->firststep)*simulator->dt;
    stats();
    next_stats += stats_every;
    next_stats = MIN(next_stats,simulator->laststep);
  }

  // decide whether any dump is due and set next dump timestep

  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (next_dump[idump] == ntimestep) {
      dump[idump]->write();
      next_dump[idump] += dump_every[idump];
    }

  // restart output and set next restart timestep

  if (next_restart == ntimestep) {
    if (restart_toggle == 0) {
      char file[128];
      sprintf(file,"%s.%d",restart1,ntimestep);
      restart->write(file);
    } else if (restart_toggle == 1) {
      restart->write(restart1);
      restart_toggle = 2;
    } else if (restart_toggle == 2) {
      restart->write(restart2);
      restart_toggle = 1;
    }
    next_restart += restart_every;
  }

  // next = next timestep any output will be done

  next = MIN(next_stats,next_restart);
  for (idump = 0; idump < ndump; idump++) next = MIN(next,next_dump[idump]);
}

/* ----------------------------------------------------------------------
   user stats command
------------------------------------------------------------------------- */

void Output::set_stats(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal stats command");

  if (strchr(arg[0],'.') == NULL) {
    stats_every = atoi(arg[0]);
    stats_delta = -1.0;
    if (stats_every < 0) error->all("Illegal stats command");
  } else {
    stats_every = -1;
    stats_delta = atof(arg[0]);
    if (stats_delta <= 0.0) error->all("Illegal stats command");
  }

  // free old species list

  for (int isp = 0; isp < nsp_stats; isp++) delete [] spname[isp];
  memory->sfree(spname);

  // set ptrs to species names in particle class

  char **name = particle->name;
  int nalias = particle->nalias;
  char **alias = particle->alias;
  int *alias2name = particle->alias2name;

  // nsp_stats = # of species in new species list
  // each arg with a wildcard is expanded to a list of names

  nsp_stats = 0;
  for (int iarg = 1; iarg < narg; iarg++) {
    if (strchr(arg[iarg],'*') == NULL) nsp_stats++;
    else {
      if (strchr(arg[iarg],'*') != strrchr(arg[iarg],'*'))
	error->all("Two or more wildcard * in stats species ID");
      for (int ialias = 0; ialias < nalias; ialias++)
	if (particle->match(arg[iarg],alias[ialias])) nsp_stats++;
    }
  }

  // allocate new list

  spname = (char **) memory->smalloc(nsp_stats*sizeof(char *),"output:spname");

  // fill new list with species names
  // each arg with a wildcard is expanded to a list of root species names

  int nlen;
  nsp_stats = 0;
  for (int iarg = 1; iarg < narg; iarg++) {
    if (strchr(arg[iarg],'*') == NULL) {
      nlen = strlen(arg[iarg]) + 1;
      spname[nsp_stats] = new char[nlen];
      strcpy(spname[nsp_stats],arg[iarg]);
      nsp_stats++;
    } else {
      for (int ialias = 0; ialias < nalias; ialias++)
	if (particle->match(arg[iarg],alias[ialias])) {
	  nlen = strlen(name[alias2name[ialias]]) + 1;
	  spname[nsp_stats] = new char[nlen];
	  strcpy(spname[nsp_stats],name[alias2name[ialias]]);
	  nsp_stats++;
	}
    }
  }
}

/* ----------------------------------------------------------------------
   modify stats parameters
------------------------------------------------------------------------- */

void Output::stats_modify(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal stats_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+1 >= narg) error->all("Illegal stats_modify command");
      if (strcmp(arg[iarg+1],"count") == 0) stats_units = COUNT;
      else if (strcmp(arg[iarg+1],"molarity") == 0) stats_units = MOLARITY;
      else if (strcmp(arg[iarg+1],"um") == 0) stats_units = UMOLARITY;
      else if (strcmp(arg[iarg+1],"nm") == 0) stats_units = NMOLARITY;
      else error->all("Illegal stats_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+1 >= narg) error->all("Illegal stats_modify command");
      delete [] format_user;
      int n = strlen(arg[iarg+1]) + 1;
      format_user = new char[n];
      strcpy(format_user,arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal stats_modify command");
  }
}

/* ----------------------------------------------------------------------
   print header for stats
   only called by proc 0
------------------------------------------------------------------------- */

void Output::stats_header()
{
  fprintf(screen,"Step Time ");
  fprintf(logfile,"Step Time ");

  if (nsp_stats == 0)
    for (int i = 0; i < particle->nspecies; i++) {
      if (screen) fprintf(screen,"%s ",particle->name[i]);
      if (logfile) fprintf(logfile,"%s ",particle->name[i]);
    }
  else
    for (int i = 0; i < nsp_stats; i++) {
      if (screen) fprintf(screen,"%s ",particle->name[stat2species[i]]);
      if (logfile) fprintf(logfile,"%s ",particle->name[stat2species[i]]);
    }

  if (screen) fprintf(screen,"\n");
  if (logfile) fprintf(logfile,"\n");
}

/* ----------------------------------------------------------------------
   print particle counts to screen & log file
------------------------------------------------------------------------- */

void Output::stats()
{
  int i;

  particle->compute_count(stats_units);

  // print stats to screen and log file

  // kludge for Red Storm timing to avoid last print-out
  // if (simulator->ntimestep == simulator->laststep) return;

  if (stats_units == COUNT) {
    int *pcount = particle->pcount;
    if (me == 0) {
      if (screen) {
	fprintf(screen,"%d %g",simulator->ntimestep,simulator->ctime);
	for (i = 0; i < nstats; i++)
	  fprintf(screen,format,pcount[stat2species[i]]);
	fprintf(screen,"\n");
      }
      if (logfile) {
	fprintf(logfile,"%d %g",simulator->ntimestep,simulator->ctime);
	for (i = 0; i < nstats; i++)
	  fprintf(logfile,format,pcount[stat2species[i]]);
	fprintf(logfile,"\n");
      }
    }

  } else {
    double *scount = particle->scount;
    if (me == 0) {
      if (screen) {
	fprintf(screen,"%d %g",simulator->ntimestep,simulator->ctime);
	for (i = 0; i < nstats; i++)
	  fprintf(screen,format,scount[stat2species[i]]);
	fprintf(screen,"\n");
      }
      if (logfile) {
	fprintf(logfile,"%d %g",simulator->ntimestep,simulator->ctime);
	for (i = 0; i < nstats; i++)
	  fprintf(logfile,format,scount[stat2species[i]]);
	fprintf(logfile,"\n");
      }
    }
  }
}

/* ----------------------------------------------------------------------
   add a Dump to list of Dumps
------------------------------------------------------------------------- */

void Output::add_dump(int narg, char **arg)
{
  if (narg < 3) error->all("Illegal dump command");

  // error checks

  for (int idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0) error->all("Reuse of dump ID");
  if (atoi(arg[1]) <= 0) error->all("Invalid dump frequency");

  // extend Dump lists if necessary

  if (ndump == maxdump) {
    maxdump += DELTA;
    dump = (Dump **)
      memory->srealloc(dump,maxdump*sizeof(Dump *),"output:dump");
    dump_every = (int *)
      memory->srealloc(dump_every,maxdump*sizeof(int *),"output:dump_every");
    next_dump = (int *)
      memory->srealloc(next_dump,maxdump*sizeof(int *),"output:next_dump");
    last_dump = (int *)
      memory->srealloc(last_dump,maxdump*sizeof(int *),"output:last_dump");
  }

  // create the Dump

  dump[ndump] = new Dump(narg,arg);
  dump_every[ndump] = atoi(arg[1]);
  last_dump[ndump] = -1;
  ndump++;
}

/* ----------------------------------------------------------------------
   modify parameters of a Dump
------------------------------------------------------------------------- */

void Output::dump_modify(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal dump_modify command");

  // find which dump it is

  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0) break;
  if (idump == ndump) error->all("Cound not find dump_modify ID");

  dump[idump]->modify_params(narg-1,&arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Dump from list of Dumps
------------------------------------------------------------------------- */

void Output::delete_dump(char *id)
{
  // find which dump it is and delete it

  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (strcmp(id,dump[idump]->id) == 0) break;
  if (idump == ndump) error->all("Could not find undump ID");

  delete dump[idump];

  // move other dumps down in list one slot

  for (int i = idump+1; i < ndump; i++) {
    dump[i-1] = dump[i];
    dump_every[i-1] = dump_every[i];
    next_dump[i-1] = next_dump[i];
    last_dump[i-1] = last_dump[i];
  }
  ndump--;
}

/* ----------------------------------------------------------------------
   setup restart capability 
------------------------------------------------------------------------- */

void Output::create_restart(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal restart command");

  if (restart) delete restart;
  delete [] restart1;
  delete [] restart2;

  restart_every = atoi(arg[0]);
  if (restart_every == 0) {
    if (narg != 1) error->all("Illegal restart command");
    return;
  }

  restart = new WriteRestart;

  int n = strlen(arg[1]) + 1;
  restart1 = new char[n];
  strcpy(restart1,arg[1]);

  if (narg == 2) {
    restart2 = NULL;
    restart_toggle = 0;
  } else if (narg == 3) {
    n = strlen(arg[2]) + 1;
    restart2 = new char[n];
    strcpy(restart2,arg[2]);
    restart_toggle = 1;
  } else error->all("Illegal restart command");
}
