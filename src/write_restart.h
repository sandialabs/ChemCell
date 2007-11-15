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

#ifndef WRITE_RESTART_H
#define WRITE_RESTART_H

#include "stdio.h"
#include "system.h"

class WriteRestart : public System {
 public:
  WriteRestart();
  ~WriteRestart() {}
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  FILE *fp;
  int nparts;

  void header();

  void write_int(int, int);
  void write_double(int, double);
  void write_char(int, char *);
};

#endif
