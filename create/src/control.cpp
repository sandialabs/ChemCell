#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <ios>
#include <iomanip>

#include <control.h>
#include <Surface.h>
#include <system.h>
#include <Distribution.h>
#include <Particle.h>
#include <random.h>
#include <geometry.h>

#define MAXLINE 200

control::control (void)
{

  script_file.clear();
  particles.clear();
  types.clear();
  surfaces.clear();
  surface_vect.clear();

}

control::~control (void)
{
  vector <Particle *>::iterator particle_i;
  for (particle_i=particles.begin() ; particle_i != particles.end() ; particle_i++) 
    {
      (*particle_i)->~Particle();
      delete (*particle_i);
    }
  vector <Surface *>::iterator surface_i;
  for (surface_i=surface_vect.begin() ; surface_i != surface_vect.end() ; surface_i++) 
    {
      // (*surface_i)->~Surface();
      delete (*surface_i);
    } 
  vector <Parttype *>::iterator type_i;
  for (type_i=types.begin() ; type_i != types.end() ; type_i++) 
    {
      // (*type_i)->~Type();
      delete (*type_i);
    } 

  script_file.clear();
  particles.clear();
  types.clear();
  surfaces.clear();
  surface_vect.clear();


}

int control::execute (char *script_file_in)
{
  string token, input_line, *last_line;
  
  list<string *> input_list;
  list<string *>::iterator input_i;
  string::size_type p;

  int seed;

  SurfaceMap::iterator surface_i;  
  Surface *current_surface;
  Surface *current_surface_kill;
  vector <Surface *>::iterator surface_v;
  string surface_name;
  string surface_name_kill;
  string surface_file;
  string output_file;
  string filesurface;
  double translate[3],scale[3],rotate[3];
  bool twod;
  bool ok;

  Parttype *p_type;

  vector<Parttype *>::iterator type_i;

  double sph_center[3],radius[2];
  double exponent, sigma,mu;
  int nclumps;


  int num_particles;
  string particle_name;
  string dist_type;
  bool dist_set = false;
  bool in_surf_set = false;
  bool out_surf_set = false;
  int prt;
  double crd[3];
  Particle *current_particle;

  script_file = script_file_in;  
  ifstream inp (script_file.data());
  int base=0;

  string pdbfile;
  pdbfile = "out.pdb";
  ofstream pdb (pdbfile.data());

  if (!inp) 
    {
      cout << "Unable to open script file: " << script_file << "\n";
      return 1;
    }
  
  last_line = new string();
  
  while (!inp.eof()) 
    {
      getline (inp, input_line);
      
      p = input_line.find_first_of("!");
      if (p <= input_line.size()) 
	{
	  input_line.erase(p,input_line.size());
	}
      if (input_line.size() > 0) 
	{
	  if (*(input_line.data()) == '+') 
	    {
	      *last_line += input_line.substr(1, input_line.size()-1);
	    }
	  else {
	    if ((*last_line).size())
	      input_list.push_back (last_line);
	    last_line = new string();
	    *last_line = input_line;
	  }
	}
    } // read input loop
  inp.close();

  if ((*last_line).size())
    input_list.push_back (last_line);

  for (input_i = input_list.begin() ; input_i != input_list.end() ; input_i++) 
    {
      
      istringstream lin((*input_i)->data());
      lin >> token;    

      if (!strcasecmp(token.data(), "RANDOM")) 
	{
	  lin >> seed;
	  cout << "Seed: "<<seed<<endl;
	  random->set_seed(seed);
	}
      else if (!strcasecmp(token.data(), "create_surf"))
	{
	  lin >> surface_name;
	  lin >> token;       
	  if (!strcasecmp(token.data(), "one")) 
	    {
	      
	      lin >> surface_file;
	      lin >> filesurface;
	      surface_i = surfaces.find(surface_name);    
	      
	      if (surface_i == surfaces.end()) 
		{
		  current_surface = 
		    new Surface(surface_name,surface_file,filesurface);
		  surfaces.insert (SurfacePair(surface_name, current_surface));
		  surface_vect.push_back(current_surface);
		}	  
	      else 
		{
		  current_surface = (*surface_i).second;
		}  
	      cout << "Created surface "<< current_surface->name 
		   <<" from "<<surface_file << endl;
	      lin >> token;
	      if(token == "origin")
		{
		  lin >> current_surface->origin[0];
		  lin >> current_surface->origin[1];
		  lin >> current_surface->origin[2];
		  current_surface->origin_set = true;
		  cout << "origin set: "<<"\t"<< current_surface->origin[0]<<"\t"<<
		    current_surface->origin[1]<<"\t"
		       <<current_surface->origin[2]<<endl;
		  lin >> token;
		}
	      if(token == "translate")
		{
		  translate[0]=0.0;translate[1]=0.0;translate[2]=0.0;
		  lin >> translate[0];
		  lin >> translate[1];
		  lin >> translate[2];
		  current_surface->translate(translate);
		  cout << "translating by: "<<"\t"<< translate[0]<<"\t"<<
		    translate[1]<<"\t"<<translate[2]<<endl;
		  lin >> token;
		}
	      
	      if(token == "scale")
		{
		  lin >> scale[0];
		  lin >> scale[1];
		  lin >> scale[2];
		  
		  cout << "scaling by: "<<"\t"<< scale[0]<<"\t"<<
		    scale[1]<<"\t"<<scale[2]<<endl;
		  current_surface->scale(scale);
		  
		  lin >> token;
		}
	      if(token == "rotate")
		{
		  lin >> rotate[0];
		  lin >> rotate[1];
		  lin >> rotate[2];
		  
		  cout << "rotating by: "<<"\t"<< rotate[0]<<"\t"<<
		    rotate[1]<<"\t"<<rotate[2]<<endl;
		  current_surface->rotate(rotate);
		  
		  lin >> token;
		}
	      current_surface->set_area();
	    }
	  current_surface->set_outside_point();
	} //create surface case
      else if (!strcasecmp(token.data(), "add_surf"))
	{
	  lin >> surface_name;
	  lin >> surface_name_kill;
	  
	  surface_i = surfaces.find(surface_name);    
	  
	  if (surface_i == surfaces.end()) 
	    {
	      current_surface = new Surface(surface_name);
	      surfaces.insert (SurfacePair(surface_name, current_surface));
	      surface_vect.push_back(current_surface);
	    }	  
	  else 
	    {
	      current_surface = (*surface_i).second;
	    }  
	  
	  surface_i = surfaces.find(surface_name_kill);    
	  
	  if (surface_i == surfaces.end()) 
	    {
	      cout << "No surface "<< surface_name_kill<< " to add."<<endl;
	      break;
	    }	  
	  else 
	    {
	      current_surface_kill = (*surface_i).second;
	    }  
	  
	  current_surface->add_surface(current_surface_kill);
	  current_surface->set_area();
	  current_surface->set_outside_point();

	  cout << "Added surface "<<surface_name_kill <<
	    " to surface "<< surface_name << endl;
	}
      else if (!strcasecmp(token.data(), "create_particle")) 
	{
	  lin >> num_particles;
	  lin >> particle_name;
	  p_type=add_particle_type(particle_name);
	  lin >> token;
	  if(token == "2d")
	    {
	      twod = true;
	      p_type->twod = twod;
	      lin >> surface_name;
	      surface_i = surfaces.find(surface_name);      
	      if (surface_i == surfaces.end()) 
		{
		  cout << "Cannot place particles on a non-existent surface "<<
		    surface_name << endl;
		}    
	      else 
		{
		  current_surface = (*surface_i).second;

		  token.clear();
		  lin >> token;
		  if (!strcasecmp(token.data(), "clumps"))
		    {
		      lin >> nclumps;
		      token.clear();
		      lin >> token; 
		      if (!strcasecmp(token.data(), "sigma"))
			lin >> sigma;
		      else 
			sigma = 0.0001;
		      distribution->param_set(current_surface,nclumps,sigma);
		    }
		  else
		    {
		      distribution->param_set(current_surface);
		    }
		  for (prt=0;prt<num_particles;prt++)
		    {
		      distribution->sample_particle(crd);
		      current_particle = new Particle(crd,particle_name,twod);
		      current_particle->surfname = surface_name;
		      particles.push_back(current_particle);
		      p_type->num ++;
		      
		    }
		  cout << "Created "<< num_particles <<" "<<particle_name<<endl;
 		  base = 0;
 		  base = pdb_out(pdb,base);
		}     
	      
	    }
	  else if(token == "sphere")
	    {
	      twod = false;
	      p_type->twod = twod;
	      lin >> sph_center[0];
	      lin >> sph_center[1];
	      lin >> sph_center[2];
	      lin >> radius[0];  // min
	      lin >> radius[1];  // max
	      
	      lin >> token;
	      
	      if (token == "dist") 
		{
		  dist_set = true;
		  lin >> token;
		  if (token == "uniform")
		    distribution->param_set(sph_center,radius);
		  else if (token == "exponential")
		    {
		      lin >> exponent;
		      distribution->param_set(sph_center,radius,exponent);
		    }
		  else if (token == "gaussian")
		    {
		      lin >> sigma;
		      lin>> mu;
		      distribution->param_set(sph_center,radius,mu,sigma);
		    }
		  lin >> token;
		}
	      if (token == "inside")
		{
		  lin >> surface_name;
		  surface_i = surfaces.find(surface_name);      
		  if (surface_i == surfaces.end()) 
		    {cout<< "surface "<<surface_name<<" not defined"<< endl;}
		  else 
		    {
		      current_surface = (*surface_i).second;
		      current_surface->set_xyzbounds();
		      in_surf_set = true;
		    } 
		  lin >> token;
		}
	      if (token == "outside")
		{
		  lin >> surface_name;
		  surface_i = surfaces.find(surface_name);      
		  if (surface_i == surfaces.end()) 
		    {cout<< "surface "<<surface_name<<" not defined"<< endl;}
		  else 
		    {
		      current_surface_kill = (*surface_i).second;
		      current_surface_kill->set_xyzbounds();
		      out_surf_set = true;
		    } 
		  lin >> token;
		}
	      if(!dist_set)
		  distribution->param_set(sph_center,radius);

	      for (prt=0;prt<num_particles;prt++)
		{
		  ok = false;
		  while(!ok)
		    {
		      ok = true;
		      distribution->sample_particle(crd);
		      if(in_surf_set && !current_surface->inside(crd)) ok = false;
		      if(out_surf_set && ok && current_surface_kill->inside(crd)) 
			ok = false;
		    }
		  
		  current_particle = new Particle(crd,particle_name,twod);
		  particles.push_back(current_particle);
		  p_type->num ++;
		}
	      dist_set = false;
	      in_surf_set = false;
	      out_surf_set = false;
	      cout << "Created "<< num_particles <<" "<<particle_name<<endl;

	      base = 0;
	      base = pdb_out(pdb,base);

	    } //sphere case

	} // create particles
	
      
	
      else if (!strcasecmp(token.data(), "write_surf")) 
	{
	  
	  lin >> output_file;
	  ofstream outfile (output_file.data());
	  string pdb_name;
	  int pdbase = 0;

	  cout << "Writing surface to file "<< output_file<<endl;
	  
	  if (!outfile) 
	    {
	      cout << "Unable to open output file: " << output_file << "\n";
	      return -2;
	    }

	  token = "";
	  lin >> token;

	  if(token == "")  
	    {
	      cout << "Writing all surfaces."<<endl;
	      for (surface_v=surface_vect.begin() ; 
		   surface_v != surface_vect.end() ; 
		   surface_v++) 
		{
		  base = (*surface_v)->write_file(outfile);
		  // uncomment the three lines below to generate pdb files of surfaces
// 		  pdb_name = (*surface_v)->name + ".pdb";
// 		  ofstream pdbout(pdb_name.data());
// 		  pdbase = (*surface_v)->pdb_out(pdbout,pdbase);
		}
	    }
	  
	  while(!(token==""))
	    {
	      
	      
	      surface_i = surfaces.find(token);      
	      if (surface_i == surfaces.end() && token!="") 
		{cout<< "surface "<<token<<" not defined"<< endl;}
	      else if(token!="")
		{
		  current_surface = (*surface_i).second;
		  base = current_surface->write_file(outfile);
		} 
	      token = "";
	      lin >> token;
	    }
	}
      else if (!strcasecmp(token.data(), "write_particle")) 
	{
	  
	  lin >> output_file;
	  ofstream outfile (output_file.data());
	  
	  cout << "Writing particles to file "<< output_file<<endl;
	  
	  if (!outfile) 
	    {
	      cout << "Unable to open output file: " << output_file << "\n";
	      return -2;
	    }

	  token = "";
	  lin >> token;

	  if(token == "")  
	    {
	      cout << "Writing all particles."<<endl;
	      base = write_particle_file(outfile,base);
	    }

	  while(!(token==""))
	    {


	      //cycle through types
	      for (type_i=types.begin() ; type_i != types.end() ; type_i++) 
		{
		  if((*type_i)->name == token) break;
		}
	      if(type_i==types.end() && token !="") 
		{cout << "type "<<token<<" not found."<<endl;}
	      
	      else if(token!="")
		{
		  base = write_particle_file(outfile, (*type_i),base);

		} 
	      token = "";
	      lin >> token;
	    }
	}  

	delete (*input_i);
    }//loop over script lines

      return 0;
} //method execute 

Parttype *control::add_particle_type(string &name_in)
{
  int ind;
  Parttype * typ_p;

  vector<Parttype *>::iterator type_i;

  for (type_i=types.begin() ; type_i != types.end() ; type_i++) 
    {
      if(name_in == (*type_i)->name) break;
    }

  if(type_i==types.end()) 
    {
      typ_p = new Parttype;
      typ_p->name = name_in;
      typ_p->num = 0;
      types.push_back(typ_p);

      return typ_p;
    }

  return (*type_i);

}

int control::write_particle_file(ofstream &out, int base)
{
// ChemCell particles              # first line of file
// species-ID N surf-ID            # N = number of particles of species-ID

// 1 x y z                         # list of vertices, numbered 1 to nvert
// ...                             # x,y,z is coord of vertex
// 100 x y z                       # file ends on this line, unless more groups

  vector<Particle *>::iterator particle_i;
  vector<Parttype *>::iterator type_i;
  vector<Surface *>::iterator surface_i;
  vector <Surface::triangle *>::iterator triangle_i;
  int ind,part;

  //cycle through types
  for (type_i=types.begin() ; type_i != types.end() ; type_i++) 
    if((*type_i)->twod)
      {
	//cycle through surfaces
	for (surface_i=surface_vect.begin() ; surface_i != surface_vect.end() ; 
	     surface_i++) 
	  {
	    //cycle through particles
	    part = 0;
	    for (particle_i=particles.begin() ; particle_i != particles.end() ; 
		 particle_i++) 
	      if((*particle_i)->type == (*type_i)->name 
		 && (*particle_i)->surfname == (*surface_i)->name)
		{
		  (*particle_i)->output_flag = true;
		  part++;
		}
	    //write header

	    out << "ChemCell particles" <<endl;
	    //	    out << endl;
	    out << (*type_i)->name << "\t"<<part<<"\t"<<(*surface_i)->name;
	    out << endl << endl;
	    part = 0;
	    //cycle through particles and output
	    for (particle_i=particles.begin() ; particle_i != particles.end() ; 
		 particle_i++) 
	      if((*particle_i)->output_flag)
	      {
		part++;
		out<< part<< "\t";
		out << setprecision(12);
		out<<(*particle_i)->crd[0]<<"  ";
		out<<(*particle_i)->crd[1]<<"  ";    
		out<<(*particle_i)->crd[2]<<endl;
		(*particle_i)->output_flag = false;
	      }

	  }
      }
  
    else
      {
	//cycle through particles
	part = 0;
	for (particle_i=particles.begin() ; particle_i != particles.end() ; 
	     particle_i++) 
	  if((*particle_i)->type == (*type_i)->name)
	    {
	      (*particle_i)->output_flag = true;
	      part++;
	    }
	//write header

	out << "ChemCell particles" <<endl;
	//	out << endl;
	out << (*type_i)->name << "\t"<<part;
	out << endl << endl;
	part = 0;
	//cycle through particles and output
	for (particle_i=particles.begin() ; particle_i != particles.end() ; 
	     particle_i++) 	      
	  if((*particle_i)->output_flag)
	    
	    {
	      part++;
	      out<< part<< "\t"
		 <<(*particle_i)->crd[0]<<" "
		 <<(*particle_i)->crd[1]<<" "
		 <<(*particle_i)->crd[2]<<endl;
	      (*particle_i)->output_flag = false;
	    }	
	out << endl;
      }
  
  return base;

}

int control::write_particle_file(ofstream &out, Parttype * type,int base)
{
// ChemCell particles              # first line of file
// species-ID N surf-ID            # N = number of particles of species-ID

// 1 x y z                         # list of vertices, numbered 1 to nvert
// ...                             # x,y,z is coord of vertex
// 100 x y z                       # file ends on this line, unless more groups

  vector<Particle *>::iterator particle_i;

  vector<Surface *>::iterator surface_i;

  int ind,part;
  
  if(type->twod)
    {
      //cycle through surfaces
      for (surface_i=surface_vect.begin() ; surface_i != surface_vect.end() ; 
	   surface_i++) 
	{
	  //cycle through particles
	  part = 0;
	  for (particle_i=particles.begin() ; particle_i != particles.end() ; 
	       particle_i++) 
	    if((*particle_i)->type == type->name 
	       && (*particle_i)->surfname == (*surface_i)->name)
	      {
		(*particle_i)->output_flag = true;
		part++;
	      }
	  //write header

	  out << "ChemCell particles" <<endl;
	  //	  out << endl;
	  out << type->name << "\t"<<part<<"\t"<<(*surface_i)->name;
	  out << endl << endl;
	  part = 0;
	  //cycle through particles and output
	  for (particle_i=particles.begin() ; particle_i != particles.end() ; 
	       particle_i++) 
	    if((*particle_i)->output_flag)
	      {
		part++;
		out<< part<< "\t";
		fixed(out);
		out.width(16);
		out << setprecision(20);
		out<<(*particle_i)->crd[0]<<"*";
		fixed(out);
		out.width(8);
		out << setprecision(3);
		out<<(*particle_i)->crd[1]<<" ";
		fixed(out);
		out.width(8);    
		out << setprecision(3);
		out<<(*particle_i)->crd[2]<<endl;
		(*particle_i)->output_flag = false;
	      }
	  out << endl;
	}
    }
  
  else
    {
      //cycle through particles
      part = 0;
      for (particle_i=particles.begin() ; particle_i != particles.end() ; 
	   particle_i++) 
	if((*particle_i)->type == type->name)
	  {
	    (*particle_i)->output_flag = true;
	    part++;
	  }
      //write header

      out << "ChemCell particles" <<endl;
      out << endl;
      out << type->name << "\t"<<part;
      out << endl << endl;
      part = 0;
      //cycle through particles and output
      for (particle_i=particles.begin() ; particle_i != particles.end() ; 
	   particle_i++) 	      
	if((*particle_i)->output_flag)
	  
	  {
	    part++;
	    out<< part<< "\t"
	       <<(*particle_i)->crd[0]<<" "
	       <<(*particle_i)->crd[1]<<" "
	       <<(*particle_i)->crd[2]<<endl;
	    (*particle_i)->output_flag = false;
	  }
      out << endl;
    }
  

  
  return base;

}


int control::pdb_out(ofstream &out, int base)
{
  vector<Particle *>::iterator particle_i;

  for (particle_i=particles.begin() ; particle_i != particles.end() ; particle_i++) {
    //    cout << "particle: "<<(*particle_i)->crd[0]<<"\t"<<(*particle_i)->crd[1]
    //	 <<"\t"<<(*particle_i)->crd[2]<<endl;
    base ++;
    out << "ATOM";
    out.width(7);
    out << base;
    out.width(5);
    out << (*particle_i)->type;
    out.width(4);
    out << "ALA";
    out.width(6);
    out << 1;
    out << "    ";
    fixed(out);
    out.width(8);
    out << setprecision(3);
    out << (*particle_i)->crd[0];
    fixed(out);
    out.width(8);
    out << setprecision(3);
    out << (*particle_i)->crd[1];
    fixed(out);
    out.width(8);
    out << setprecision(3);
    out << (*particle_i)->crd[2];
    out << "\n";
  }

  return base;
}
