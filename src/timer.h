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

#ifndef TIMER_H
#define TIMER_H

#include "system.h"

#define TIME_LOOP       0
#define TIME_MOVE       1
#define TIME_MIGRATE    2
#define TIME_REACT      3
#define TIME_REACT_COMM 4
#define TIME_OUTPUT     5
#define TIME_BALANCE    6

#define TIME_N          7

class Timer : public System {
 public:
  double *array;

  Timer();
  ~Timer();
  void init();
  void stamp();
  void stamp(int);
  void sub_stamp();
  void sub_stamp(int);
  void barrier_start(int);
  void barrier_stop(int);
  double elapsed(int);

 private:
  double previous_time;
  double previous_time_sub;
};

#endif
