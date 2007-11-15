#include <string>
#include <Surface.h>
#include <system.h>
#include <random.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <ios>
#include <iomanip>
#include <geometry.h>
#include <vector>


Surface::Surface(string & name_in, string & filename_in,string &filesurface)
{
  name = name_in;
  clear();
  read_surface(filename_in,filesurface);
}

Surface::Surface(string & name_in)
{
  name = name_in;
  clear();
}

// Private methods

void Surface::clear()
{
    vertices.clear();
    triangles.clear();
    edges.clear();
    tot_area = 0.0;

    origin[0]=0.0;origin[1]=0.0;origin[2]=0.0;
    origin_set = false;

    surface_set = false;
    return;
}

Surface::~Surface()
{
  vector <edge *>::iterator edge_i;
  for (edge_i=edges.begin() ; edge_i != edges.end() ; edge_i++) 
    {
      (*edge_i)->~edge();
      delete (*edge_i);
    }
  vector <triangle *>::iterator triangle_i;
  for (triangle_i=triangles.begin() ; triangle_i != triangles.end() ; triangle_i++) 
    {
      (*triangle_i)->~triangle();
      delete (*triangle_i);
    }
  vector <vertex *>::iterator vertex_i;
  for (vertex_i=vertices.begin() ; vertex_i != vertices.end() ; vertex_i++) 
    {
      (*vertex_i)->~vertex();
      delete (*vertex_i);
    }
  vertices.clear();
  edges.clear();
  triangles.clear();
}

void Surface::add_surface(Surface * psurf_in)
{

  int old_nvrt = vertices.size();
  int old_ntri = triangles.size();


  vector <vertex *>::iterator vertex_i;

  for (vertex_i=psurf_in->vertices.begin() ; 
       vertex_i != psurf_in->vertices.end() ; vertex_i++) 

    set_vertex(*vertex_i);

  vector <triangle *>::iterator triangle_i;
  for (triangle_i=psurf_in->triangles.begin() ; 
       triangle_i != psurf_in->triangles.end() ; triangle_i++) 
    set_triangle(*triangle_i,old_nvrt, old_ntri);

  vector <edge *>::iterator edge_i;
  for (edge_i=psurf_in->edges.begin() ; 
       edge_i != psurf_in->edges.end() ; edge_i++) 
    set_edge(*edge_i,old_nvrt,old_ntri);

  psurf_in->vertices.clear();
  psurf_in->edges.clear();
  psurf_in->triangles.clear();

   surface_set = true;
//   detect_edges();
   calculate_bounding();
}

void Surface::read_surface(string & filename,string & filesurface)
{
  string line, input_line, token="",token1="", surface_name_in="";
  double x,y,z;
  int a,b,c,d,e,f;
  int num_vertices, num_triangles,num_edges;
  string::size_type p;
  int i;

  ifstream inp (filename.data());
  if (!inp) {
    cout << "Unable to open surface file: " << filename << "\n";
  }

  istringstream lin;
  surface_name_in = "";

  while(surface_name_in!=filesurface && !inp.eof())
    {
      while(token!="ChemCell" && token1!="surface"&& !inp.eof())
	{

	  token = "";token1="";
	  getline (inp, input_line);
	  lin.clear();
	  lin.str(input_line.data());
	  lin >> token;
	  lin >> token1;
	}
      getline (inp, input_line);
      lin.clear();
      lin.str(input_line.data());
      lin >> surface_name_in;

    }
  if(inp.eof()) 
    {
      cout << "Surface "<<name<<" not found."<<endl;
    }
  else
    {
      lin >> num_vertices;
      lin >> num_triangles;
    }

  while(token != "1")
     {
      getline (inp, input_line);
      lin.clear();
      lin.str(input_line.data());
      lin >> token;
     }

  for (i=0;i<num_vertices;i++)
    {
      lin >> x; lin >> y; lin >> z;
      getline (inp, input_line);
      lin.clear();
      lin.str(input_line.data());
      lin >> token;
      set_vertex(x,y,z);
      //      cout << "vertex: "<< x << " "<< y << " " << z << endl;
    }
  while(token != "1")
     {
      getline (inp, input_line);
      lin.clear();
      lin.str(input_line.data());
      lin >> token;
     }

  for (i=0;i<num_triangles;i++)
    {
      lin >> a; lin >> b; lin >> c;
      getline (inp, input_line);
      lin.clear();
      lin.str(input_line.data());
      lin >> token;
      set_triangle(vertices[a-1],vertices[b-1],vertices[c-1],a,b,c);
      //      cout << "triangle: "<< a << " "<< b << " " << c << endl;
    }
  while(token != "1")
     {
      getline (inp, input_line);
      lin.clear();
      lin.str(input_line.data());
      lin >> token;
     }

  for (i=0;i<num_triangles;i++)
    {
      lin >> a; lin >> b; lin >> c;lin >> d; lin >> e; lin >> f;
      set_edge(i+1,triangles[i],1,a,triangles[a-1],b);
      set_edge(i+1,triangles[i],2,c,triangles[c-1],d);
      set_edge(i+1,triangles[i],3,e,triangles[e-1],f);
      getline (inp, input_line);
      lin.clear();
      lin.str(input_line.data());
      lin >> token;
      //       cout << "edge: "<< a << " "<< b << " " << c 
      //	    << " "<< d <<" "<< e << " " << f <<endl;
    }

  
   inp.close();
   surface_set = true;
//   detect_edges();
   calculate_bounding();
   //   set_outside_point();

}
void Surface::set_vertex(double x, double y, double z)
{
  vertex* vrt = new vertex;
  vrt->crd[0]=x;
  vrt->crd[1]=y;
  vrt->crd[2]=z;
  vertices.push_back(vrt);
}

void Surface::set_vertex(vertex * vertin)
{
  vertices.push_back(vertin);
}
void Surface::set_xyzbounds()
{
  vector<triangle *> ::iterator triangle_i; 
  vertex * vert_p;

  for(triangle_i=triangles.begin();triangle_i!=triangles.end();triangle_i++)
    {
      //set xy bounds
      (*triangle_i)->xbound[0] = 10000.0;
      (*triangle_i)->xbound[1] = -10000.0;
      (*triangle_i)->ybound[0] = 10000.0;
      (*triangle_i)->ybound[1] = -10000.0;
      (*triangle_i)->zbound[0] = 10000.0;
      (*triangle_i)->zbound[1] = -10000.0;
      
      for (int i = 0; i < 3;i++)
	{
	  vert_p = (*triangle_i)->vert[i];
	  if (vert_p->crd[0]<(*triangle_i)->xbound[0]) 
	    (*triangle_i)->xbound[0] = vert_p->crd[0];
	  if (vert_p->crd[0]>(*triangle_i)->xbound[1]) 
	    (*triangle_i)->xbound[1] = vert_p->crd[0];
	  if (vert_p->crd[1]<(*triangle_i)->ybound[0]) 
	    (*triangle_i)->ybound[0] = vert_p->crd[1];
	  if (vert_p->crd[1]>(*triangle_i)->ybound[1]) 
	    (*triangle_i)->ybound[1] = vert_p->crd[1];
	  if (vert_p->crd[2]<(*triangle_i)->zbound[0]) 
	    (*triangle_i)->zbound[0] = vert_p->crd[2];
	  if (vert_p->crd[2]>(*triangle_i)->zbound[1]) 
	    (*triangle_i)->zbound[1] = vert_p->crd[2];
	}
    }
}
void Surface::set_triangle(vertex *a, vertex *b, vertex *c,int ia,int ib,int ic)
{


  triangle* trg = new triangle;
  trg->iv[0] = ia;
  trg->iv[1] = ib;
  trg->iv[2] = ic;
  trg->vert[0]=a;
  trg->vert[1]=b;
  trg->vert[2]=c;

  trg->edgtr[0]=0; trg->edgtr[1]=0; trg->edgtr[2]=0;
  trg->side[0]=0; trg->side[1]=0; trg->side[2]=0;

  triangles.push_back(trg);
  //  triangle_out(trg);
}

void Surface::set_triangle(triangle * trgin, int old_nvrt, int old_ntri)
{
  triangles.push_back(trgin);
  trgin->iv[0] += old_nvrt;
  trgin->iv[1] += old_nvrt;
  trgin->iv[2] += old_nvrt;

  trgin->edgtr[0] += old_ntri;
  trgin->edgtr[1] += old_ntri;
  trgin->edgtr[2] += old_ntri;

}

void Surface::set_edge(int ia,triangle *a, int b, int ic,triangle *c, int d)
{
  edge* edg = new edge;
  edg->itrg[0] = ia;
  edg->itrg[1] = ic;
  edg->trg[0]=a;
  edg->side[0]=b-1;
  edg->trg[1]=c;
  edg->side[1]=d-1;
  edges.push_back(edg);

  a->edgtr[b-1]=ic; c->edgtr[d-1]=ia;
  a->side[b-1]=d; c->side[d-1]=b;
}

void Surface::set_edge(edge * edgin,int old_nvrt, int old_ntri)
{
  edges.push_back(edgin);
  edgin->itrg[0] += old_ntri;
  edgin->itrg[1] += old_ntri; 
}

void Surface::detect_edges()
{
  int i, j, k,l,n;

  n = triangles.size();

  for (i=0;i<n-1;i++)
    for (j=i+1;j<n;j++)
      for (k=0;k<3;k++)
	for (l=0;l<3;l++)   
	  {
	    if(edge_compare(i,k,j,l)&&!edge_set(i,k,j,l)) 
	      {
		set_edge(i+1,triangles[i],k+1,j+1,triangles[j],l+1);
		//    cout << "edge: "<<i<<" "<<k<<" "<<j<<" "<<l<< endl;
	      }
	  }
}

bool Surface::edge_set(int tri1,int side1, int tri2,int side2)
{
  bool test = false;

  triangle *tr1 = triangles[tri1];
  triangle *tr2 = triangles[tri2];
  
  if(tr1->edgtr[side1]>0&&tr2->edgtr[side2]>0) test = true;

  return test;

}

bool Surface::edge_compare(int tri1,int side1, int tri2,int side2)
{
  int v1, v2, v11, v12, v21, v22;

  triangle *tr1 = triangles[tri1];
  triangle *tr2 = triangles[tri2];

  //side 0 (0,1); side 1 (1,2); side 2 (2,0)
  v1 = side1 %3; v2 = (side1+1) %3;
  v11 = tr1->iv[v1]; v12 = tr1->iv[v2];
  v1 = side2 %3; v2 = (side2+1) %3;
  v21 = tr2->iv[v1]; v22 = tr2->iv[v2];

  if((v11==v21&&v12==v22)||(v11==v22&&v12==v21)) return true;
  else return false;

}

void Surface::set_area()
{
  vector<triangle *> ::iterator triangle_i;
  double area;

  tot_area = 0;
  cumulative_area.clear();

  for(triangle_i=triangles.begin();triangle_i!=triangles.end();triangle_i++)
    {
      area = triangle_area(*triangle_i);
      (*triangle_i)->area = area;
      tot_area += area;
      cumulative_area.push_back(tot_area);
    }
}

double Surface::triangle_area(triangle *trg)
{
  double v1[3], v2[3], v3[3];
  double vec1[3], vec2[3], vec3[3];
  double area;

  v1[0] = trg->vert[0]->crd[0];
  v1[1] = trg->vert[0]->crd[1];
  v1[2] = trg->vert[0]->crd[2];  
  v2[0] = trg->vert[1]->crd[0];
  v2[1] = trg->vert[1]->crd[1];
  v2[2] = trg->vert[1]->crd[2];  
  v3[0] = trg->vert[2]->crd[0];
  v3[1] = trg->vert[2]->crd[1];
  v3[2] = trg->vert[2]->crd[2];  

  subtract_vectors(v2,v1,vec1);
  subtract_vectors(v3,v1,vec2);
  cross_product(vec1,vec2,vec3);
  area = sqrt(dot_product(vec3,vec3))/2.;

  return area;
}
void Surface::find_center(double cntr[3])
{
  vector<vertex *> ::iterator vertex_i;

  cntr[0] = cntr[1] = cntr[2] = 0;

  for(vertex_i=vertices.begin();vertex_i!=vertices.end();vertex_i++)
    {
      cntr[0] += (*vertex_i)->crd[0];
      cntr[1] += (*vertex_i)->crd[1];
      cntr[2] += (*vertex_i)->crd[2];
    }
  cntr[0] = cntr[0]/vertices.size();
  cntr[1] = cntr[1]/vertices.size();
  cntr[2] = cntr[2]/vertices.size();
}

void Surface::translate(double trans[3])
{
  vector<vertex *> ::iterator vertex_i;

  if (origin_set)
    {
      origin[0] += trans[0];
      origin[1] += trans[1]; 
      origin[2] += trans[2];
    }


  for(vertex_i=vertices.begin();vertex_i!=vertices.end();vertex_i++)
    {
      (*vertex_i)->crd[0] = (*vertex_i)->crd[0] + trans[0];
      (*vertex_i)->crd[1] = (*vertex_i)->crd[1] + trans[1];     
      (*vertex_i)->crd[2] = (*vertex_i)->crd[2] + trans[2];     
    }
  calculate_bounding();

}
void Surface::scale(double scale[3])
{
  vector<vertex *> ::iterator vertex_i;
  double center1[3],center2[3],differ[3];

  if (origin_set)
    {
      center1[0] = origin[0];
      center1[1] = origin[1];
      center1[2] = origin[2];
    }
  else
      find_center(center1);

  for(vertex_i=vertices.begin();vertex_i!=vertices.end();vertex_i++)
    {
      (*vertex_i)->crd[0] = (*vertex_i)->crd[0]*scale[0];
      (*vertex_i)->crd[1] = (*vertex_i)->crd[1]*scale[1];     
      (*vertex_i)->crd[2] = (*vertex_i)->crd[2]*scale[2];     
    }
  if (origin_set)
    {
      center2[0] = origin[0]*scale[0];
      center2[1] = origin[1]*scale[1];
      center2[2] = origin[2]*scale[2];
    }
  else
    find_center(center2);

  subtract_vectors(center1,center2,differ);
  translate(differ);
  calculate_bounding();

}
void Surface::rotate(double angle[3])
{
  double sth,sph,sps,cth,cph,cps;
  double rot[3][3];
  vector<vertex *> ::iterator vertex_i;
  double center[3];
  double tvec[3];

  if (origin_set)
    {
      center[0] = origin[0];
      center[1] = origin[1];
      center[2] = origin[2];
    }
  else
    find_center(center);

  center[0] = -center[0];center[1] = -center[1];center[2] = -center[2];
  translate(center);

  sth = sin(angle[0]);cth = cos(angle[0]);
  sph = sin(angle[1]);cph = cos(angle[1]);
  sps = sin(angle[2]);cps = cos(angle[2]);

  rot[0][0] = cph*cth*cps - sph*sps;
  rot[0][1] = -cph*cth*sps - sph*cps;
  rot[0][2] = cph*sth;

  rot[1][0] = sph*cth*cps + cph*sps;
  rot[1][1] = -sph*cth*sps + cph*cps;
  rot[1][2] = sph*sth;

  rot[2][0] = -sth*cps;
  rot[2][1] = sth*sps;
  rot[2][2] = cth;

  for(vertex_i=vertices.begin();vertex_i!=vertices.end();vertex_i++)
    {
      tvec[0] = (*vertex_i)->crd[0];
      tvec[1] = (*vertex_i)->crd[1];
      tvec[2] = (*vertex_i)->crd[2];
      (*vertex_i)->crd[0] = tvec[0]*rot[0][0] + tvec[1]*rot[0][1] + 
	tvec[2]*rot[0][2];
      (*vertex_i)->crd[1] = tvec[0]*rot[1][0] + tvec[1]*rot[1][1] + 
	tvec[2]*rot[1][2];
      (*vertex_i)->crd[2] = tvec[0]*rot[2][0] + tvec[1]*rot[2][1] + 
	tvec[2]*rot[2][2];
    }

  center[0] = -center[0];center[1] = -center[1];center[2] = -center[2];
  translate(center);
  calculate_bounding();

}

int Surface::find_triangle_bin(double uni_rnd)
{
  bool no;
  int lo, hi;
  double stop = uni_rnd * tot_area;
  int search;
  int i;

  for (i=0;i<triangles.size();i++)
    {
      if(cumulative_area[i]>stop) break;
    }

//   no = true;
//   lo = 0;
//   hi = triangles.size()-1;
//   search = hi+lo/2;

//   while(no)
//     {
//       if (cumulative_area[search] < stop)
// 	{
// 	  lo = search;
// 	  search = (hi+lo)/2; 
// 	}
//       if  (cumulative_area[search] > stop)
// 	{
// 	  hi = search;
// 	  search = (hi+lo)/2; 
// 	}
//       if(hi-search<1) no = false;
//     }

//  cout << "rand: "<< uni_rnd <<" bin: "<< i << endl;
  return i;

}
void Surface::populate_triangle(triangle * trg, double coord[3])
{
  double v1[3], v2[3], v3[3];
  double vec1[3], vec2[3], vec3[3];
  double area,r1,r2;

  v1[0] = trg->vert[0]->crd[0];
  v1[1] = trg->vert[0]->crd[1];
  v1[2] = trg->vert[0]->crd[2];  
  v2[0] = trg->vert[1]->crd[0];
  v2[1] = trg->vert[1]->crd[1];
  v2[2] = trg->vert[1]->crd[2];  
  v3[0] = trg->vert[2]->crd[0];
  v3[1] = trg->vert[2]->crd[1];
  v3[2] = trg->vert[2]->crd[2];  
  
  subtract_vectors(v2,v1,vec1);
  r1 = random->uniform();
  vec1[0] = vec1[0]*r1;vec1[1] = vec1[1]*r1;vec1[2] = vec1[2]*r1;
  
  add_vectors(v1,vec1,vec2);

  subtract_vectors(v3,vec2,vec1);
  r2 = random->uniform();
  vec1[0] = vec1[0]*r2;vec1[1] = vec1[1]*r2;vec1[2] = vec1[2]*r2;

  add_vectors(vec2,vec1,coord);
}
void Surface::populate_triangle(triangle * trg, double coord[3],double sigma)
{
  double v1[3], v2[3], v3[3];
  double vec1[3], vec2[3], vec3[3];
  double area,r1,r2;
  
  v1[0] = trg->vert[0]->crd[0];
  v1[1] = trg->vert[0]->crd[1];
  v1[2] = trg->vert[0]->crd[2];  
  v2[0] = trg->vert[1]->crd[0];
  v2[1] = trg->vert[1]->crd[1];
  v2[2] = trg->vert[1]->crd[2];  
  v3[0] = trg->vert[2]->crd[0];
  v3[1] = trg->vert[2]->crd[1];
  v3[2] = trg->vert[2]->crd[2];  
  
  subtract_vectors(v2,v1,vec1);
  r1 = random->gaussian(sigma,0.0);
  vec1[0] = vec1[0]*r1;vec1[1] = vec1[1]*r1;vec1[2] = vec1[2]*r1;
  
  add_vectors(v1,vec1,vec2);

  subtract_vectors(v3,vec2,vec1);
  r2 = random->gaussian(sigma,0.0);
  vec1[0] = vec1[0]*r2;vec1[1] = vec1[1]*r2;vec1[2] = vec1[2]*r2;

  add_vectors(vec2,vec1,coord);
}
void Surface::calculate_bounding()
{
  double xmin,ymin,zmin,xmax,ymax,zmax;
  vector<vertex *> ::iterator vertex_i;

  xmin=ymin=zmin=10e-6;
  xmax=ymax=zmax=-10e6;

  for(vertex_i=vertices.begin();vertex_i!=vertices.end();vertex_i++)
    {
      if((*vertex_i)->crd[0]>xmax) xmax = (*vertex_i)->crd[0];
      if((*vertex_i)->crd[1]>ymax) ymax = (*vertex_i)->crd[1];
      if((*vertex_i)->crd[2]>zmax) zmax = (*vertex_i)->crd[2];
      if((*vertex_i)->crd[0]<xmin) xmin = (*vertex_i)->crd[0];
      if((*vertex_i)->crd[1]<ymin) ymin = (*vertex_i)->crd[1];
      if((*vertex_i)->crd[2]<zmin) zmin = (*vertex_i)->crd[2];
    }

  bounding_box[0] = xmax;bounding_box[1] = xmin;
  bounding_box[2] = ymax;bounding_box[3] = ymin;
  bounding_box[4] = zmax;bounding_box[5] = zmin;
  
  //  cout << "xmax: "<<xmax<< " ymax: "<<ymax << " zmax: "<< zmax << endl;
  //  cout << "xmin: "<<xmin<< " ymin: "<<ymin << " zmin: "<< zmin << endl;
  set_xyzbounds();
}

void Surface::set_outside_point()
{
  double epsilon = .1;

  calculate_bounding();

  outside[0] = bounding_box[0]+epsilon;
  outside[1] = bounding_box[2]+epsilon;
  outside[2] = bounding_box[4]+epsilon;

}

bool Surface::inside(double x[3])
{
  
  double normal[3],intersect_point[3];
  double parametric_value;
  int side_of_intersect,intsct_flag;
  double v0[3],v1[3],v2[3];
  double s1[3],s2[3];
  double norm;
  double min = 10000;
  int insout=0;
  vector<triangle *> ::iterator triangle_i;



//   cout << "particle test: " << endl;
//   cout << x[0]<<" " <<x[1]<< " "<<x[2]<<endl;

//   cout << "bounding_box: "<< endl;
//   cout << bounding_box[0] << " "<<
//     bounding_box[1]<<"\t "<<
//     bounding_box[2] << " "<<
//     bounding_box[3]<<"\t "<< 
//     bounding_box[4] << " "<<
//     bounding_box[5]<<"\t "<<endl;

//   cout << "outside point: "<<endl;
//   cout << outside[0]<<" " <<outside[1] << " " << outside[2]<<endl;
  if (!inside(x[0],x[1],x[2])) 
    {
      //  cout << "outside bounding "<< x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
      return false;
    }
  outside[0] = x[0]; outside[1] = x[1]; outside[2] = 10000.0;

  // int skips = 0;
  
  for(triangle_i=triangles.begin();triangle_i!=triangles.end();triangle_i++)
    {

      if (x[0]<(*triangle_i)->xbound[0] || x[0]>(*triangle_i)->xbound[1] ) 
	{
	  //	  skips++;
	  continue;
	}
      if (x[1]<(*triangle_i)->ybound[0] || x[1]>(*triangle_i)->ybound[1] ) 
	{
	  //	  skips ++;
	  continue;
	}
      //if (x[2]>(*triangle_i)->zbound[1]) 
      //	{
	  //	  skips ++;
      //	    continue;
      //    }
      v0[0] = (*triangle_i)->vert[0]->crd[0];
      v0[1] = (*triangle_i)->vert[0]->crd[1];
      v0[2] = (*triangle_i)->vert[0]->crd[2];
      v1[0] = (*triangle_i)->vert[1]->crd[0];
      v1[1] = (*triangle_i)->vert[1]->crd[1];
      v1[2] = (*triangle_i)->vert[1]->crd[2];
      v2[0] = (*triangle_i)->vert[2]->crd[0];
      v2[1] = (*triangle_i)->vert[2]->crd[1];
      v2[2] = (*triangle_i)->vert[2]->crd[2];
      
      subtract_vectors(v1,v0,s1);
      subtract_vectors(v2,v1,s2);
      cross_product(s1,s2,normal);
      
      norm = sqrt(dot_product(normal,normal));
      
      normal[0] = normal[0]/norm;
      normal[1] = normal[1]/norm;
      normal[2] = normal[2]/norm;
      
//       cout <<(*triangle_i)->vert[0]->crd[0]<<" "
// 	   <<(*triangle_i)->vert[0]->crd[1]<<" "
// 	   <<(*triangle_i)->vert[0]->crd[2]<<"\t "
// 	   <<(*triangle_i)->vert[1]->crd[0]<<" "
// 	   <<(*triangle_i)->vert[1]->crd[1]<<" "
// 	   <<(*triangle_i)->vert[1]->crd[2]<<" \t"
// 	   <<(*triangle_i)->vert[2]->crd[0]<<" "
// 	   <<(*triangle_i)->vert[2]->crd[1]<<" "
// 	   <<(*triangle_i)->vert[2]->crd[2]<<endl;
//       cout << "normal: " << normal[0]<<" "<<normal[1] <<" "<<normal[2]<<endl;
      
      if(tri_line_intersect(v0,v1,v2,normal,x,outside,intersect_point,
			    parametric_value,intsct_flag,side_of_intersect) )
	{
// 	  cout << " intersecting triangle: "
// 	       <<parametric_value<<" "<<side_of_intersect <<endl;
	  if(parametric_value<min)
	  {
// 	    cout << " its a min " <<endl;
	    min = parametric_value;
	    insout = side_of_intersect;
	  }
	}
    }
  //  cout<< "skips: "<< skips<<" in if 2 "<<insout<<endl;
  
  if (insout !=2) 
    {
      //      cout <<"not this.."<<insout <<endl;
      return false;
    }
  else  
    {
      //      cout <<"this.."<<insout<<endl;
      return true;
    }

  
}



 bool Surface::inside(double x, double y, double z)
{
  //bounding box screen

  if (x>bounding_box[0] || x<bounding_box[1]) return false;
  if (y>bounding_box[2] || y<bounding_box[3]) return false;
  if (z>bounding_box[4] || z<bounding_box[5]) return false;

  return true;
}

//  bool Surface::inside(double x[3])
// {
//   if(x[0]>bounding_box[0]) return false;
//   if(x[0]<bounding_box[1]) return false;
//   if(x[1]>bounding_box[2]) return false;
//   if(x[1]<bounding_box[3]) return false;
//   if(x[2]>bounding_box[4]) return false;
//   if(x[2]<bounding_box[5]) return false;

// //   cout <<"coord: "<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
// //   cout << "bounding: "<< bounding_box[0]<<" "<<
// //     bounding_box[1]<<"\t "<<
// //     bounding_box[2]<<" "<<
// //     bounding_box[3]<<"\t "<<
// //     bounding_box[4]<<" "<<
// //     bounding_box[5]<<endl;
// //   if(ins) {cout<<"true"<<endl;}
// //   else cout << "false"<<endl;
//   return true;
// }

void Surface::triangle_out(triangle *trg)
{
  cout << "Triangle : "<<endl;
  cout << "vertex 0: "<< trg->vert[0]->crd[0] << "\t";
  cout << trg->vert[0]->crd[1] << "\t";
  cout << trg->vert[0]->crd[2] << "\n";
  cout << "vertex 1: "<< trg->vert[1]->crd[0] << "\t";
  cout << trg->vert[1]->crd[1] << "\t";
  cout << trg->vert[1]->crd[2] << "\n";
  cout << "vertex 2: "<< trg->vert[2]->crd[0] << "\t";
  cout << trg->vert[2]->crd[1] << "\t";
  cout << trg->vert[2]->crd[2] << "\n";
  cout << "Area: " << trg->area << endl;
}

void Surface::cross_product(double a[3],double b[3],double c[3])
{

  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = a[2]*b[0]-a[0]*b[2];
  c[2] = a[0]*b[1]-a[1]*b[0];
}
void Surface::subtract_vectors(double a[3],double b[3],double c[3])
{
  c[0] = a[0]-b[0];c[1] = a[1]-b[1];c[2] = a[2]-b[2];
}

void Surface::add_vectors(double a[3],double b[3],double c[3])
{
  c[0] = a[0]+b[0];c[1] = a[1]+b[1];c[2] = a[2]+b[2];
}
double Surface::dot_product(double a[3],double b[3])
{
  return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

int Surface::pdb_out(ofstream &out, int base)
{
  vector<vertex *>::iterator vertex_i;

  for (vertex_i=vertices.begin() ; vertex_i != vertices.end() ; vertex_i++) {
    base ++;
    out << "ATOM";
    out.width(7);
    out << base;
    out.width(5);
    if (name == "ONE")
      out << "S";
    if (name == "TWO")
      out << "C";
    else
      out << "H";
    out.width(4);
    out << "ALA";
    out.width(6);
    out << 1;
    out << "    ";
    fixed(out);
    out.width(8);
    out << setprecision(3);
    out << (*vertex_i)->crd[0];
    fixed(out);
    out.width(8);
    out << setprecision(3);
    out << (*vertex_i)->crd[1];
    fixed(out);
    out.width(8);
    out << setprecision(3);
    out << (*vertex_i)->crd[2];
    out << "\n";
  }

  return base;
}

int Surface::write_file(ofstream &out)
{
  vector<vertex *>::iterator vertex_i;
  vector<triangle *>::iterator triangle_i;
  vector<edge *>::iterator edge_i;
  int indx;

  if (vertices.size()==0) return 0;

  out << "ChemCell surface" <<endl;
  out << name << "\t"<< vertices.size()<<"\t"<<triangles.size()<<endl;
  out << endl;

  indx = 0;

  for (vertex_i=vertices.begin() ; vertex_i != vertices.end() ; vertex_i++) 
    {
      indx ++;
      out << indx;
      fixed(out);
      out.width(8);
      out << setprecision(3);
      out<< (*vertex_i)->crd[0];
      fixed(out);
      out.width(8);
      out << setprecision(3);
      out<< (*vertex_i)->crd[1];
      fixed(out);
      out.width(8);
      out << setprecision(3);
      out<< (*vertex_i)->crd[2]<<endl;
    }

  out << endl;

  indx = 0;
  for (triangle_i=triangles.begin() ; triangle_i != triangles.end() ; triangle_i++) 
    {
      indx ++;
      out << indx << "\t" << (*triangle_i)->iv[0]<<"\t"
	  <<(*triangle_i)->iv[1]<<"\t"<<(*triangle_i)->iv[2]<<"\t"<<endl;
    }  

  out << endl;
  indx = 0;
  for (triangle_i=triangles.begin() ; triangle_i != triangles.end() ; triangle_i++) 
    {
      indx ++;

      out << indx                 << "\t";
      out << (*triangle_i)->edgtr[0]<< "\t";
      out <<(*triangle_i)->side[0]<< "\t";
      out << (*triangle_i)->edgtr[1]<< "\t";
      out <<(*triangle_i)->side[1]<< "\t";
      out << (*triangle_i)->edgtr[2]<< "\t";
      out <<(*triangle_i)->side[2]<< endl;
    }  
  out << endl;
}
