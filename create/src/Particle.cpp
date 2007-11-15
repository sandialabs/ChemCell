#include <Particle.h>
#include <system.h>
#include <random.h>
#include <Surface.h>
#include <iostream>

Particle::Particle(double crd_in[3], string & type_in, bool dim)
{
  type = type_in;
  crd[0] = crd_in[0];crd[1] = crd_in[1];crd[2] = crd_in[2];
  membrane = dim;
  output_flag = false ;
}
Particle::~Particle(){}
