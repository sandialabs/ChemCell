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

#ifndef CHEM_H
#define CHEM_H

#include "system.h"

class Chem : public System {
 public:
  char *style;

  double volume;                  // system volume for Gillespie model
  int nbinbin;                    // # of bin-bin interactions
  int nbinpair;                   // # of pairwise checks between bins
  int ndist;                      // # of distance checks
  int noverlap;                   // # of reaction overlaps
  int *rcount;                    // count of reactions of each kind
  int doneflag;                   // Gillespie flag for no more reactions
  int sortflag;                   // spatial flag for sorting particles

  int prob_style;
  double max_prob;
  double diff_factor;

  Chem();
  virtual ~Chem();
  virtual void init() = 0;
  virtual void reactions() = 0;
  virtual void create() {}
  virtual void dynamic() {}
  virtual double maxbin() {return 0.0;}
  virtual int memory_usage() {return 0;};

  void set_prob(int, char **);
};

#endif
