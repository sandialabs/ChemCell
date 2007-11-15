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

/* ----------------------------------------------------------------------
   parent class for all of ChemCell
   all other classes inherit from this
   contains static ptrs to single instance of other classes
   contains MPI communicator and file handles for my world of procs
------------------------------------------------------------------------- */

#ifndef SYSTEM_H
#define SYSTEM_H

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"

class Universe;
class Error;
class Memory;
class Input;

class Balance;
class Chem;
class Domain;
class Grid;
class Modify;
class Move;
class Output;
class Particle;
class Random;
class React;
class Simulator;
class Surf;
class Timer;

class System {
 public:
  static Memory *memory;          // memory allocation functions
  static Error *error;            // error handling
  static Universe *universe;      // universe of processors
  static Input *input;            // input script processing

  static Balance *balance;        // load-balancer
  static Chem *chem;              // reaction chemistry
  static Domain *domain;          // global domain
  static Grid *grid;              // volumetric grid
  static Modify *modify;          // fixes
  static Move *move;              // particle motion
  static Output *output;          // screen and file output
  static Particle *particle;      // particles
  static Random *random;          // random # generators
  static React *react;            // chemical reactions
  static Simulator *simulator;    // top-level simulation
  static Surf *surf;              // surfaces
  static Timer *timer;            // CPU timing info

  static MPI_Comm world;          // communicator for my world of procs
  static FILE *infile;            // infile for my world
  static FILE *screen;            // screen output for my world
  static FILE *logfile;           // logfile for my world

  void open(int, char **, MPI_Comm);
  void create();
  void destroy();
  void close();
};

#endif
