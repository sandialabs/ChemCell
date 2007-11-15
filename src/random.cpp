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

#include "stdio.h"
#include "math.h"
#include "random.h"
#include "time.h"

#define IM 2147483647
#define AM (1.0/IM)

#define IA 16807
#define IQ 127773
#define IR 2836

#define IA2 48271
#define IQ2 44488
#define IR2 3399

#define IA3 69621
#define IQ3 30845
#define IR3 23902

#define IA4 40692
#define IQ4 52744
#define IR4 3791
#define IM4 2147483399

/* ---------------------------------------------------------------------- */

Random::Random()
{
  time_t seconds;
  time(&seconds);
  seed = (unsigned int) seconds;

  // seed = 123456;
}

/* ----------------------------------------------------------------------
   RNG called for each particle read from file
   uses and updates internal seed (which can be set by user)
   RNG is Park/Miller with 2nd set of constants
   return int seed to be stored with read-in particle
------------------------------------------------------------------------- */

int Random::read()
{
  int k = seed/IQ2;
  seed = IA2*(seed-k*IQ2) - IR2*k;
  if (seed < 0) seed += IM;
  return seed;
}

/* ----------------------------------------------------------------------
   RNG called for particle moves
   uses and updates particle's seed
   RNG is Park/Miller with 1st set of constants
   return uniform RN from 0.0 to 1.0
------------------------------------------------------------------------- */

double Random::move(int *pseed)
{
  int iseed = *pseed;
  int k = iseed/IQ;
  iseed = IA*(iseed-k*IQ) - IR*k;
  if (iseed < 0) iseed += IM;
  double ans = AM*iseed;
  *pseed = iseed;
  return ans;
}

/* ----------------------------------------------------------------------
   RNG called for particle reactions
   uses combination of 2 reactant's seeds and does not update them
     combines 2 seeds by adding them as unsigned ints to avoid overflow
   RNG is Park/Miller with 3rd set of constants
   return uniform RN from 0.0 to 1.0
------------------------------------------------------------------------- */

double Random::react(int seed1, int seed2)
{
  unsigned int useed1 = seed1;
  unsigned int useed2 = seed2;
  int iseed = (useed1+useed2) % IM;
  int k = iseed/IQ3;
  iseed = IA3*(iseed-k*IQ3) - IR3*k;
  if (iseed < 0) iseed += IM;
  double ans = AM*iseed;
  return ans;
}

/* ----------------------------------------------------------------------
   RNG called for each product created in a reaction
   uses and updates 1st reactant's seed
     which will be discarded after all reactants are created
   RNG is L'Ecuyer with 4th set of constants (Num Recip p 271)
     must be a different period than move() (Park/Miller)
     else will get multiple particles with same seed
     happens when 2 products of same parent, 1st of them later produces prod
     grandchild and its uncle will have same seed
     even though moves with different RNG happened in between!
   return int iseed to be stored with product particle
------------------------------------------------------------------------- */

int Random::product(int *pseed)
{
  int iseed = *pseed;
  int k = iseed/IQ4;
  iseed = IA4*(iseed-k*IQ4) - IR4*k;
  if (iseed < 0) iseed += IM4;
  *pseed = iseed;
  return iseed;
}

/* ----------------------------------------------------------------------
   RNG called for Gillespie simulation
   uses and updates internal seed (which can be set by user)
   RNG is standard Park/Miller with 1st set of constants
   return uniform RN from 0.0 to 1.0
------------------------------------------------------------------------- */

double Random::gillespie()
{
  int k = seed/IQ;
  seed = IA*(seed-k*IQ) - IR*k;
  if (seed < 0) seed += IM;
  double ans = AM*seed;
  return ans;
}
