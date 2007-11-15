// parent class for all of ChemCreate
// all other classes inherit from this
// contains static ptrs to single instance of other core classes

#ifndef SYSTEM_H
#define SYSTEM_H

#include <math.h>

class Random;
class Distribution;
class Particle;

class System {
 public:
  static Random *random;      // random numbers
  static Distribution *distribution; //distribution
  static Particle *particle; //particle class
  void create();
  void destroy();
};

#endif
