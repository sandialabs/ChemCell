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

#ifndef OUTPUT_H
#define OUTPUT_H

#include "system.h"

class Dump;
class WriteRestart;

class Output : public System {
 public:
  int next;                          // next timestep for any kind of output
  int next_stats;                    // next timestep to do stats on
  int stats_every;                   // stats output every this many steps
  double stats_delta;                // stats output every this delta time

  Output();
  ~Output();
  void init();
  void stats();
  void set_stats(int, char **);      // setup screen stats
  void stats_modify(int, char **);
  void add_dump(int, char **);       // add a Dump to Dump list
  void dump_modify(int, char **);
  void delete_dump(char *);          // delete a Dump from Dump list
  void write(int);                   // output for current timestep
  void create_restart(int, char **);

 private:
  int me;

  int nsp_stats;               // # of species to print stats on, 0 if all
  int nstats;                  // # of species to print stats on, true # if all
  char **spname;               // names of species to print stats on
  int *stat2species;           // ptr to species index for each stat
  int stats_units;             // unit setting for stats output
  char *format;                // format string for stats output
  char *format_user;           // user-specified stats format

  int next_dump_any;           // next timestep for any Dump
  int ndump;                   // # of Dumps defined
  int maxdump;                 // max size of Dump list
  Dump **dump;                 // list of defined Dumps
  int *dump_every;             // output of each Dump every this many steps
  int *next_dump;              // next timestep to do each Dump
  int *last_dump;              // last timestep each Dump was done
 
  int next_restart;            // next timestep to write a restart file
  int restart_every;           // write a restart file every this many steps
  int restart_toggle;          // 0 if use restart1 as prefix
                               // 1 if use restart1 as file, 2 for restart2
  char *restart1,*restart2;    // names for restart files
  WriteRestart *restart;       // Restart output

  void stats_header();
};

#endif
