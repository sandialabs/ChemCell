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

#ifndef SIM_SPATIAL_H
#define SIM_SPATIAL_H

#include "simulator.h"

class SimSpatial : public Simulator {
 public:
  SimSpatial(int, char **);
  ~SimSpatial();
  void init();
  void run();

  void sort(int);
  static int compare(const void *, const void *);
};

#endif
