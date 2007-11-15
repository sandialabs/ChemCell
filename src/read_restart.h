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

#ifndef READ_RESTART_H
#define READ_RESTART_H

#include "stdio.h"
#include "system.h"

class ReadRestart : public System {
 public:
  ReadRestart() {}
  ~ReadRestart() {}
  void command(int, char **);

 private:
  int me,nprocs;
  FILE *fp;
  int nprocs_file,nparts;

  void header();

  int read_int();
  double read_double();
  char *read_char();
};

#endif
