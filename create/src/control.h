#ifndef CONTROL_H
#define CONTROL_H

#include <string>
#include <map>
#include <system.h>
#include <Surface.h>

using namespace std;

struct Parttype{string name; int num;bool twod;};

class control: public System
{
public:

  control (void);
  ~control (void);
  int execute (char *);


private:


  string script_file;

  string surface_file;

  string output_file;

  SurfaceMap surfaces;

  vector <Surface *> surface_vect;

  vector <Particle *> particles;

  vector <Parttype *> types;

  Parttype *add_particle_type(string &);

  int pdb_out(ofstream &, int);

  int write_particle_file(ofstream &,int);
  int write_particle_file(ofstream &,Parttype *, int);
  void write_file(string &);

};

#endif // CONTROL_H
