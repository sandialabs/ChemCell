
#include <iostream>
#include <string>
#include <system.h>

#include <Surface.h>
#include <Distribution.h>
#include <random.h>
#include <control.h>

int main(int argc, char *argv[]) 
{
  //  double part[3];

  System *system;
  system->create();


//   string surfname,filename;
//   double trans[3],scale[3],rotate[3];
//   int base=0;

//   filename = "try.surf";
//   string pdbfile ="dat.pdb";

//   ofstream out (pdbfile.data());
//   if (!out) {
//     cout << "Unable to open pdb file: " << pdbfile << "\n";
//     return -2;
//   }

//   surfname = "ONE";
//   Surface *surf_p = new Surface(surfname,filename);


//   trans[0] = 0.0;trans[1] = 0.0; trans[2] = 0.0;
//   surf_p->translate(trans);
//   scale[0] = 1.0;scale[1] = 1.0; scale[2] = 1.0;
//   surf_p->scale(scale);
//   rotate[0] = 0.0;rotate[1] = 0.0; rotate[2] = 0.0;
//   surf_p->rotate(rotate);

//   surf_p->set_area();  
//   base = surf_p->pdb_out(out,0);


//   surfname = "TWO";
//   Surface *surf_p1 = new Surface(surfname,filename);


//   trans[0] = 0.0;trans[1] = 0.0; trans[2] = 0.0;
//   surf_p1->translate(trans);
//   scale[0] = 3.0;scale[1] = 3.0; scale[2] = 3.0;
//   surf_p1->scale(scale);
//   rotate[0] = 1.0;rotate[1] = 0.0; rotate[2] = 0.0;
//   surf_p1->rotate(rotate);

//   surf_p1->set_area();  
//   base = surf_p1->pdb_out(out,base);



//   out.close();

  control *c;
  if (argc == 2 ) {
    c = new control();
    c->execute(argv[1]);
    
    //c->~control();
       delete c;
  }
  else {
    cout << "usage: main <script file name>\n";
  }


  system->destroy();


}

