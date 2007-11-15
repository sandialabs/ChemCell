#ifndef Distribution_H
#define Distribution_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

//#include <Surface.h>
//#include <random.h>

#include <system.h>

class Surface;

using namespace std;

class Distribution : public System {
 public:
  Distribution();
  ~Distribution();

  void param_set(double[3],double[2]); /*uniform*/
  void param_set(double[3],double[2],double); /*exponential*/
  void param_set(double[3],double[2],double,double); /*gaussian*/
  void param_set(Surface *); /*2d uniform*/
  void Distribution::param_set(Surface *,int,double);/*2d clumps*/

 void sample_particle(double[3]);

 private:

 //type (0 uniform, 1 exponential, 2 gaussian)
 int dist_type;
 //origin
 double origin[3];
 //outer radius
 double rmax;
 //inner radius
 double rmin;
 //exponential
 double exponent;
 //gaussian
 double mu,sigma;
 //2d surface
 Surface *surf_p;

 vector <double>centers;
 int clumps;
 
};


#endif
