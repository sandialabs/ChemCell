#ifndef Particle_H
#define Particle_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <system.h>
#include <map>

using namespace std;

class Particle : public System
{
 public:
 Particle(double [3], string &, bool);
 ~Particle();
 double crd[3];
 string type;
 bool membrane;
 string surfname;
 bool output_flag;

 private:


};

#endif
