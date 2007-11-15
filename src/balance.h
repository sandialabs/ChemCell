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

#ifndef BALANCE_H
#define BALANCE_H

#include "system.h"
#include "zoltan.h"

class Balance : public System {

public:
  Zoltan_Struct *lb;           // current bin decomposition
  Zoltan_Struct *lb_new;       // new bin decomposition
  
  Balance();
  ~Balance();

  // load-balancing methods for bins

  void RCB_setup();
  void RCB(int, double *, double *, double *, double *, double *, double *);
  void swap();

  // Zoltan callbacks for bins

  static int geometry_dimension(void *, int *);
  static int bin_number(void *, int *);
  static void bin_coords(void *, int, int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR,
			 int, double *, int *);
  static void bin_wts_bin(void *, int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR,
			  int, float *, int *);
  static void bin_wts_particle(void *, int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR,
			       int, float *, int *);

 private:
  int me;
};

#endif
