#include <Distribution.h>
#include <system.h>
#include <random.h>
#include <Surface.h>

Distribution::Distribution(){}
Distribution::~Distribution(){}

void Distribution::param_set(double orig[3], double reach[2])
{
  dist_type = 0;
  origin[0] = orig[0];origin[1] = orig[1];origin[2] = orig[2];
  rmin = reach[0];
  rmax = reach[1];
}
void Distribution::param_set(double orig[3], double reach[2], double expo)
{
  dist_type = 1;
  origin[0] = orig[0];origin[1] = orig[1];origin[2] = orig[2];
  exponent = expo;
  rmin = reach[0];
  rmax = reach[1];
}
void Distribution::param_set(double orig[3], double reach[2], 
			     double mu_in, double sigma_in)
{
  dist_type = 2;
  origin[0] = orig[0];origin[1] = orig[1];origin[2] = orig[2];
  mu = mu_in;
  rmin = reach[0];
  rmax = reach[1];
}
void Distribution::param_set(Surface * surf_p_in)
{
  dist_type = 3;
  surf_p = surf_p_in;  
  cout << "Setting surface distribution"<<endl;

}
void Distribution::param_set(Surface * surf_p_in, int nclumps,
			     double sigma_in)
{
  dist_type = 4;
  surf_p = surf_p_in;
  sigma = sigma_in;
  centers.clear();
  clumps = nclumps;
  for (int i = 0; i < clumps; i++)
    centers.push_back(random->uniform());
  cout << "Setting surface distribution to "<<nclumps<<
    " clumps with sigma "<<sigma<<endl;

}

void Distribution::sample_particle(double part[3])
{
  
  double sample[3],norm;
  int trg;

  switch ( dist_type )
    {
    case 0:    //uniform
      norm = 1.E15;
      while ((norm > rmax) || (norm < rmin))
	{
	  sample[0] = (random->uniform()-.5)*2.*rmax;
	  sample[1] = (random->uniform()-.5)*2.*rmax;
	  sample[2] = (random->uniform()-.5)*2.*rmax;
	  norm = sqrt(sample[0]*sample[0]+sample[1]*sample[1]+
		      sample[2]*sample[2] );
	}
      
      part[0] = sample[0]+origin[0];
      part[1] = sample[1]+origin[1];
      part[2] = sample[2]+origin[2];
      
      //       cout <<"sample: "<< part[0]<<" "<<part[1]<<" "<<part[2]<<endl;
      //       cout << "norm: "<< norm<<endl;
      //       cout  <<"rmax "<<rmax << "rmin "<<rmin <<endl;
      break;
    case 1:
      norm = 1.E5;
      while  ((norm > rmax) || (norm < rmin))
	{
	  sample[0] = random->exponential(exponent)*rmax;
	  sample[1] = random->exponential(exponent)*rmax;
	  sample[2] = random->exponential(exponent)*rmax;
	  norm = sqrt(sample[0]*sample[0]+sample[1]*sample[1]+
		      sample[2]*sample[2] );
	}
      part[0] = sample[0]+origin[0];
      part[1] = sample[1]+origin[1];
      part[2] = sample[2]+origin[2];
      break;
    case 2:
      
      norm = 1.E15;
      while ((norm > rmax) || (norm < rmin))
	{
	  sample[0] = random->gaussian(sigma,mu)*rmax;
	  sample[1] = random->gaussian(sigma,mu)*rmax;
	  sample[2] = random->gaussian(sigma,mu)*rmax;
	  norm = sqrt(sample[0]*sample[0]+sample[1]*sample[1]+
		      sample[2]*sample[2] );
	}
      part[0] = sample[0]+origin[0];
      part[1] = sample[1]+origin[1];
      part[2] = sample[2]+origin[2];
      break;
    case 3:

      trg = surf_p->find_triangle_bin(random->uniform());
      surf_p->populate_triangle(surf_p->triangles[trg], part);
      break;
    case 4:
 
      trg = (int)((double)clumps*random->uniform());
      trg = surf_p->find_triangle_bin(centers[trg]);
      surf_p->populate_triangle(surf_p->triangles[trg], part,sigma);
      break;
      //    default:
      //cout << "What distribution is this?? " << endl;
    }
}

