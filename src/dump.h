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

#ifndef DUMP_H
#define DUMP_H

#include "stdio.h"
#include "system.h"

class Dump : public System {
 public:
  char *id;                  // user-defined dump ID

  Dump(int, char **);
  ~Dump();
  void init();
  void write();
  void modify_params(int, char **);

 private:
  int me,nprocs;             // proc info
  int size_one;              // # of datums in one line of output
  int maxbuf;                // size of buf
  double *buf;               // memory for atom quantities
  FILE *fp;                  // file to write dump to
  char *format;              // format string to write one line of output
  int orient_flag;           // whether to print surface particle orientation

  int nsp;                   // size of species mask array, 0 if all species
  int *spmask;               // mask flag for each species: 1 = dump, 0 = no

  void write_header(int);
  int pack();
  void write_data(int, double *);
};

#endif
