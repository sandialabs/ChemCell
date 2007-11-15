
#ifndef Surface_H
#define Surface_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <system.h>
#include <map>

using namespace std;

class Surface : public System
{
 public:
 Surface(string &,string &,string &);
 Surface(string &);
 ~Surface();

 string name;
 double origin[3];
 bool origin_set;

 struct vertex{double crd[3];};
 struct triangle
 {
   vertex *vert[3];
   int iv[3];
   double area;
   int edgtr[3];
   int side[3];
   double xbound[2];
   double ybound[2];
   double zbound[2];
 };
 struct edge {triangle *trg[2];int itrg[2];int side[2];};

 vector<vertex *> vertices;
 vector<triangle *> triangles;
 vector<edge *> edges;
 
 vector <double> cumulative_area;

 void read_surface(string &,string &);

 void add_surface(Surface *);

 void set_vertex(double , double , double );
 void set_vertex(vertex *);
 void set_triangle(vertex *,vertex *,vertex *,int,int,int);
 void set_triangle(triangle *, int, int);
 void set_xyzbounds();
 void set_edge(int,triangle *, int, int,triangle *, int);
 void set_edge(edge *,int,int);

 void translate(double [3]);
 void scale(double [3]);
 void rotate(double [3]);

 void set_area();
 void calculate_bounding();
 double bounding_box[6];

 void set_outside_point();
 double outside[3];

 double triangle_area(triangle *);
 int find_triangle_bin(double);
 void populate_triangle(triangle *, double[3]);
 void populate_triangle(triangle *, double[3], double);
 void detect_edges();
 bool edge_compare(int, int, int, int);
 bool edge_set(int, int, int, int);

 void find_center(double [3]);

 int pdb_out(ofstream &, int);
 int write_file(ofstream &);

 bool surface_set; 


 bool inside(double, double, double);
 bool inside(double[3]);

 double dot_product(double [3],double [3]);
 void cross_product(double [3],double [3],double [3]);
 void subtract_vectors(double [3],double [3],double [3]);
 void add_vectors(double [3],double [3],double [3]);

 private:

 void clear();
 void triangle_out(triangle *);
 double tot_area;

};

typedef map<string, Surface*> SurfaceMap;
typedef SurfaceMap::value_type SurfacePair;

#endif
