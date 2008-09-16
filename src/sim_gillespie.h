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

#ifndef SIM_GILLESPIE_H
#define SIM_GILLESPIE_H

#include "simulator.h"

class SimGillespie : public Simulator {
 public:
  SimGillespie(int, char **);
  ~SimGillespie();
  void init();
  void run();

 private:
  int runflag,outflag;
  double stats_time;
};

#endif
