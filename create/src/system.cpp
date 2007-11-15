#include "system.h"

#include <string>
#include <Distribution.h>
#include <random.h>
#include <Particle.h>
#include <control.h>

/* set all static ptrs in System to NULL */

Random *System::random = NULL;
Distribution *System::distribution = NULL;
Particle *System::particle = NULL;

/* allocate single instance of System classes */

void System::create()
{
  string dummy = "dummy";
  double dumm[3];
  bool dum = true;

  dumm[0]=0;dumm[1]=0;dumm[2]=0;

  random = new Random;
  distribution = new Distribution;
  particle = new Particle(dumm,dummy,dum);

}

/* delete single instance of all System classes */

void System::destroy()
{
  delete random;
  delete distribution;
  delete particle;

}
