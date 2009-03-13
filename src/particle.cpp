/* ----------------------------------------------------------------------
   ChemCell - Cell simulator for particle diffusion and reactions
   Steve Plimpton (sjplimp@sandia.gov)
   Alex Slepoy (alexander.slepoy@nnsa.doe.gov)
   Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/chemcell.html

   Copyright (2004) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level ChemCell directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "particle.h"
#include "simulator.h"
#include "domain.h"
#include "chem.h"
#include "grid.h"
#include "surf.h"
#include "region.h"
#include "move.h"
#include "random.h"
#include "balance.h"
#include "geometry.h"
#include "error.h"
#include "memory.h"

#define AVOGADRO 6.023e23
#define MAXLINE 256
#define BIG 1.0e20
#define FACTOR 1.5
#define DELTA 1

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

enum {COUNT,MOLARITY,UMOLARITY,NMOLARITY};    // same as in output.cpp

/* ---------------------------------------------------------------------- */

Particle::Particle()
{
  nlocal = nghost = ntotal = nghost_last = 0;
  maxpart = 0;
  plist = NULL;

  nspecies = maxspecies = 0;
  name = NULL;

  nalias = maxalias = 0;
  alias = NULL;
  alias2name = NULL;

  diffusivity = NULL;
  dimension = NULL;
  pcount = NULL;
  ccount = NULL;
  scount = NULL;

  size1 = size2 = size3 = 0;
  buf1 = buf2 = buf3 = NULL;
  sizeproc = 0;
  proclist = NULL;

  size_restart = 6;

  MPI_Comm_rank(world,&me);
}

/* ---------------------------------------------------------------------- */

Particle::~Particle()
{
  memory->sfree(plist);

  for (int i = 0; i < nspecies; i++) delete [] name[i];
  memory->sfree(name);

  for (int i = 0; i < nalias; i++) delete [] alias[i];
  memory->sfree(alias);
  memory->sfree(alias2name);

  memory->sfree(diffusivity);
  memory->sfree(dimension);
  memory->sfree(pcount);
  memory->sfree(ccount);
  memory->sfree(scount);

  memory->sfree(buf1);
  memory->sfree(buf2);
  memory->sfree(buf3);
  memory->sfree(proclist);
}

/* ----------------------------------------------------------------------
   init - print particle stats
------------------------------------------------------------------------- */

void Particle::init()
{
  // species stats

  int count2d = 0;
  int count3d = 0;
  double min2d = BIG;
  double min3d = BIG;
  double max2d = -1.0;
  double max3d = -1.0;

  for (int i = 0; i < nspecies; i++) {
    if (dimension[i] == 2) {
      count2d++;
      min2d = MIN(diffusivity[i],min2d);
      max2d = MAX(diffusivity[i],max2d);
    }
    if (dimension[i] == 3) {
      count3d++;
      min3d = MIN(diffusivity[i],min3d);
      max3d = MAX(diffusivity[i],max3d);
    }
  }

  if (me == 0) {
    fprintf(screen,"Particles:\n");
    fprintf(screen,"  number of 3d species = %d\n",count3d);
    if (count3d)
      fprintf(screen,"  min/max 3d diffusivities = %g %g\n",min3d,max3d);
    fprintf(screen,"  number of 2d species = %d\n",count2d);
    if (count2d)
      fprintf(screen,"  min/max 2d diffusivities = %g %g\n",min2d,max2d);

    fprintf(logfile,"Particles:\n");
    fprintf(logfile,"  number of 3d species = %d\n",count3d);
    if (count3d)
      fprintf(logfile,"  min/max 3d diffusivities = %g %g\n",min3d,max3d);
    fprintf(logfile,"  number of 2d species = %d\n",count2d);
    if (count2d)
      fprintf(logfile,"  min/max 2d diffusivities = %g %g\n",min2d,max2d);
  }

  // particle count

  int nall;
  MPI_Allreduce(&nlocal,&nall,1,MPI_INT,MPI_SUM,world);
  if (me == 0) {
    fprintf(screen,"  initial particles = %d\n",nall);
    fprintf(logfile,"  initial particles = %d\n",nall);
  }
}

/* ----------------------------------------------------------------------
   add a particle to end of plist, could be own or ghost particle
   grow particle array if necessary
   set species and xyz coords
   return index of added particle
------------------------------------------------------------------------- */

int Particle::add(int species, double x, double y, double z)
{
  // grow plist by FACTOR if necessary

  if (ntotal == maxpart) {
    if (ntotal == 0) maxpart = 1000;
    else maxpart = static_cast<int> (FACTOR*ntotal);
    plist = (OnePart *) 
      memory->srealloc(plist,maxpart*sizeof(OnePart),"particle:plist");
  }

  plist[ntotal].species = species;
  plist[ntotal].x[0] = x;
  plist[ntotal].x[1] = y;
  plist[ntotal].x[2] = z;
  ntotal++;
  return ntotal-1;
}

/* ----------------------------------------------------------------------
   return species index corresponding to str
   if str is in alias list, return associated species index
   if not, return -1
------------------------------------------------------------------------- */

int Particle::find(char *str)
{
  for (int ialias = 0; ialias < nalias; ialias++)
    if (strcmp(str,alias[ialias]) == 0) return alias2name[ialias];
  return -1;
}

/* ----------------------------------------------------------------------
   create a new species-ID by adding str to name list and alias list
   assumes str is not already a species-ID
   extend all species-based arrays and set default values
   return index of added species
------------------------------------------------------------------------- */

int Particle::add_species(char *str)
{
  // grow name and other species arrays by DELTA if necessary

  if (nspecies == maxspecies) {
    maxspecies += DELTA;
    name = (char **)
      memory->srealloc(name,maxspecies*sizeof(char *),
		       "particle:name");
    dimension = (int *)
      memory->srealloc(dimension,maxspecies*sizeof(int),
		       "particle:dimension");
    diffusivity = (double *) 
      memory->srealloc(diffusivity,maxspecies*sizeof(double),
		       "particle:diffusivity");
    pcount = (int *)
      memory->srealloc(pcount,maxspecies*sizeof(int),
		       "particle:pcount");
    ccount = (double *)
      memory->srealloc(ccount,maxspecies*sizeof(double),
		       "particle:ccount");
    scount = (double *)
      memory->srealloc(scount,maxspecies*sizeof(double),
		       "particle:scount");

    if (simulator->spatial_flag) move->grow_partsurf(maxspecies,surf->nsurf);
  }

  // add str to species root name list and alias list

  int nlen = strlen(str) + 1;
  name[nspecies] = new char[nlen];
  strcpy(name[nspecies],str);
  add_alias(str,nspecies);

  // set species attributes to default values

  dimension[nspecies] = 3;
  diffusivity[nspecies] = 0.0;
  pcount[nspecies] = 0;
  ccount[nspecies] = 0.0;
  if (simulator->spatial_flag) move->default_partsurf(nspecies,-1);

  nspecies++;
  return nspecies-1;
}

/* ----------------------------------------------------------------------
   add a string to alias list and point it at ispecies
   assumes str is not already a species-ID
------------------------------------------------------------------------- */

void Particle::add_alias(char *str, int ispecies)
{
  // grow alias arrays by DELTA if necessary

  if (nalias == maxalias) {
    maxalias += DELTA;
    alias = (char **) memory->srealloc(alias,maxalias*sizeof(char *),
				       "particle:alias");
    alias2name = (int *) memory->srealloc(alias2name,maxalias*sizeof(int),
					  "particle:alias2name");
  }

  // set alias string and point it to root name

  int nlen = strlen(str) + 1;
  alias[nalias] = new char[nlen];
  strcpy(alias[nalias],str);
  alias2name[nalias] = ispecies;
  nalias++;
}

/* ----------------------------------------------------------------------
   define a new species name with optional aliases
------------------------------------------------------------------------- */

void Particle::set_species(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal species command");

  // check if species exists
  // if yes, and no extra aliases, then error
  // if no, create it

  int ispecies = find(arg[0]);
  if (ispecies >= 0 && narg == 1) 
    error->all("Species already exists and no new aliases are defined");
  if (ispecies == -1) ispecies = add_species(arg[0]);

  // loop over alias names and add them to list

  for (int i = 1; i < narg; i++) {
    if (find(arg[i]) >= 0) {
      char *str = new char[128];
      sprintf(str,"Alias %s already exists",arg[i]);
      error->all(str);
    }
    add_alias(arg[i],ispecies);
  }
}

/* ----------------------------------------------------------------------
   set a species diffusion constant
   if species does not exist, it is an error
------------------------------------------------------------------------- */

void Particle::set_diffusion(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal diffusion command");
  double d = atof(arg[1]);
  if (d < 0.0) error->all("Invalid diffusion constant");

  if (strchr(arg[0],'*') == NULL) {
    int ispecies = find(arg[0]);
    if (ispecies == -1) error->all("Unknown species in diffusion command");
    diffusivity[ispecies] = d;
  } else {
    if (strchr(arg[0],'*') != strrchr(arg[0],'*'))
	error->all("Two or more wildcard * in diffusion species ID");
    int ncount = 0;
    for (int ialias = 0; ialias < nalias; ialias++)
      if (match(arg[0],alias[ialias])) {
	diffusivity[alias2name[ialias]] = d;
	ncount++;
      }
    if (me == 0) {
      fprintf(screen,"%d diffusion coeffs set via wildcard\n",ncount);
      fprintf(logfile,"%d diffusion coeffs set via wildcard\n",ncount);
    }
  }
}

/* ----------------------------------------------------------------------
   set dimensionality of a species
   if species does not exist, it is an error
   cannot change dimension of a species with existing particles
   cannot set dimension to 2d for a species with permeability defined
------------------------------------------------------------------------- */

void Particle::set_dimension(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal dimension command");
  int dim = atoi(arg[1]);
  if (dim != 2 && dim != 3) error->all("Invalid dimensionality");

  if (strchr(arg[0],'*') == NULL) {
    int ispecies = find(arg[0]);
    if (ispecies == -1) error->all("Unknown species in dimension command");
    dimension[ispecies] = dim;
    check_dimension(ispecies);
  } else {
    if (strchr(arg[0],'*') != strrchr(arg[0],'*'))
	error->all("Two or more wildcard * in dimension species ID");
    int ncount = 0;
    for (int ialias = 0; ialias < nalias; ialias++)
      if (match(arg[0],alias[ialias])) {
	dimension[alias2name[ialias]] = dim;
	check_dimension(alias2name[ialias]);
	ncount++;
      }
    if (me == 0) {
      fprintf(screen,"%d dimensions set via wildcard\n",ncount);
      fprintf(logfile,"%d dimensions set via wildcard\n",ncount);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Particle::check_dimension(int ispecies)
{
  // check if any particles of these species already exist

  Particle::OnePart *plist = particle->plist;
  int nlocal = particle->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (plist[i].species == ispecies) flag = 1;

  int allflag;
  MPI_Allreduce(&flag,&allflag,1,MPI_INT,MPI_SUM,world);
  if (allflag) error->all("Cannot change dimension of existing particles");

  // for 2d, check if permeability already set to non-defaults

  if (dimension[ispecies] == 2) {
    int flag = 0;
    Move::Part *part = move->part;
    for (int i = 0; i < surf->nsurf; i++) {
      if (part[ispecies].surf[i].in_prob[0] != 1.0) flag = 1;
      if (part[ispecies].surf[i].out_prob[0] != 1.0) flag = 1;
      if (part[ispecies].surf[i].in_species[0] != ispecies) flag = 1;
      if (part[ispecies].surf[i].out_species[0] != ispecies) flag = 1;
    }
    if (flag) error->all("Permeability is set for 2d species");
  }
}

/* ----------------------------------------------------------------------
   set explicit count of a species
   if species does not exist, it is an error
------------------------------------------------------------------------- */

void Particle::set_count(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal count command");

  double count;
  count = atof(arg[1]);

  if (strchr(arg[0],'*') == NULL) {
    int ispecies = find(arg[0]);
    if (ispecies == -1) error->all("Unknown species in count command");
    pcount[ispecies] = 0;
    ccount[ispecies] = 0.0;
    if (simulator->stochastic_flag)
      pcount[ispecies] = static_cast<int> (count);
    else ccount[ispecies] = count;
  } else {
    if (strchr(arg[0],'*') != strrchr(arg[0],'*'))
	error->all("Two or more wildcard * in count species ID");
    int ncount = 0;
    for (int ialias = 0; ialias < nalias; ialias++)
      if (match(arg[0],alias[ialias])) {
	pcount[alias2name[ialias]] = 0;
	ccount[alias2name[ialias]] = 0.0;
	if (simulator->stochastic_flag)
	  pcount[alias2name[ialias]] = static_cast<int> (count);
	else ccount[alias2name[ialias]] = count;
	ncount++;
      }
    if (me == 0) {
      fprintf(screen,"%d counts set via wildcard\n",ncount);
      fprintf(logfile,"%d counts set via wildcard\n",ncount);
    }
  }
}

/* ----------------------------------------------------------------------
   convert pcount or ccount of each species into
     appropriate format/units for stats output
------------------------------------------------------------------------- */

void Particle::compute_count(int flag)
{
  int i;
  double factor;

  // store integer particle count in pcount:
  // if spatial stochastic, sum counts across all procs
  // if ODE, convert ccount molarity to particle count
  // if Gillespie, count already stored in pcount

  if (flag == COUNT) {
    if (simulator->spatial_flag) {
      int *my_count = new int[nspecies];
      for (i = 0; i < nspecies; i++) pcount[i] = my_count[i] = 0;
      for (i = 0; i < nlocal; i++) my_count[plist[i].species]++;
      MPI_Allreduce(my_count,pcount,nspecies,MPI_INT,MPI_SUM,world);
      delete [] my_count;
    } else if (simulator->stochastic_flag == 0) {
      double factor = AVOGADRO * chem->volume;
      for (i = 0; i < nspecies; i++)
	pcount[i] = static_cast<int> (factor*ccount[i]);
    }

  // store molarity in scount in requested units:
  // not possible for spatial since volume is uncertain
  // if ODE, convert ccount
  // if Gillespie, convert pcount to concentration

  } else {
    if (simulator->stochastic_flag == 0) {
      if (flag == MOLARITY) factor = 1.0;
      else if (flag == UMOLARITY) factor = 1.0e6;
      else if (flag == NMOLARITY) factor = 1.0e9;
      for (i = 0; i < nspecies; i++) scount[i] = factor*ccount[i];
    } else {
      if (flag == MOLARITY) factor = 1.0 / (AVOGADRO * chem->volume);
      else if (flag == UMOLARITY) factor = 1.0e6 / (AVOGADRO * chem->volume);
      else if (flag == NMOLARITY) factor = 1.0e9 / (AVOGADRO * chem->volume);
      for (i = 0; i < nspecies; i++) scount[i] = factor*pcount[i];
    }
  }
}

/* ----------------------------------------------------------------------
   read group of particles from input file
   assign each particle to a proc
   if particles diffuse on a 2d surf, assign to a triangle
------------------------------------------------------------------------- */

void Particle::read(int narg, char **arg)
{
  char line[MAXLINE];
  char *err;

  if (me == 0) {
    if (screen) fprintf(screen,"Reading particles ...\n");
    if (logfile) fprintf(logfile,"Reading particles ...\n");
  }

  if (narg != 2 && narg != 3) error->all("Illegal particles command");

  // nprevious = total # of particles before reading new ones

  int nprevious;
  MPI_Allreduce(&nlocal,&nprevious,1,MPI_INT,MPI_SUM,world);

  // check particle species-ID
  // error if doesn't exist

  int ispecies = find(arg[0]);
  if (ispecies == -1) error->all("Unknown species in particles command");
  int n = atoi(arg[1]);
  if (n <= 0) error->all("Illegal particles command");

  // if 2 args, check if existing species is 3d

  if (narg == 2 && dimension[ispecies] == 2)
    error->all("Reading 3d particles for a 2d species");

  // if 3 args, check if existing species is 2d and surf ID exists

  int isurf;
  if (narg == 3) {
    if (dimension[ispecies] == 3)
      error->all("Reading 2d particles for a 3d species");
    isurf = surf->find(arg[2]);
    if (isurf == -1) error->all("Particles surf-ID does not exist");
  }

  // read and bcast list of particle coords

  double **xyz = memory->create_2d_double_array(n,3,"particle:xyz");

  if (me == 0) {
    err = fgets(line,MAXLINE,infile);
    if (err == NULL) error->one("Unexpected end of file");
    int tmp;
    for (int i = 0; i < n; i++) {
      err = fgets(line,MAXLINE,infile);
      sscanf(line,"%d %lg %lg %lg",&tmp,&xyz[i][0],&xyz[i][1],&xyz[i][2]);
    }
    if (err == NULL) error->one("Unexpected end of file");
  }
    
  MPI_Bcast(xyz[0],3*n,MPI_DOUBLE,0,world);

  // scan particle list and keep those whose bin maps to my owned domain
  // add particle to my list (species and coords)
  // initialize other particle quantities (ibin, itri)

  int seed,ix,iy,iz,ibin,ip,itri,iglobal;
  double near,distance;
  double *v0,*v1,*v2;
  Grid::OneBin *blist = grid->blist;
  Surf::OneTri *tlist = surf->tlist;
  double **vlist = surf->vlist;
  Region **rlist = surf->rlist;
  int *region2surf = surf->region2surf;
  Grid::Index *ptr;

  for (int i = 0; i < n; i++) {
    if (!domain->inside(xyz[i])) {
      char str[128];
      sprintf(str,"Particle %d is outside global domain",i+1);
      error->all(str);
    }

    seed = random->read();
    iglobal = grid->whichglobal(xyz[i],&ix,&iy,&iz);
    ibin = grid->global2local(ix,iy,iz);
    if (ibin == -1 || blist[ibin].ghost) continue;

    ip = add(ispecies,xyz[i][0],xyz[i][1],xyz[i][2]);
    nlocal++;
    plist[ip].seed = seed;
    plist[ip].ibin = ibin;
    plist[ip].itri = 0;
    
    // if a 2d particle, find which local triangle it is on
    // loop over all triangles in bin that belong to the correct surf
    // compute nearness of particle to each triangle
    // assign particle to triangle it is nearest to
    
    if (dimension[ispecies] == 2) {
      near = BIG;
      ptr = blist[ibin].tri;
      while (ptr) {
	itri = ptr->index;
	ptr = ptr->ptr;
	if (itri >= 0) {
	  if (tlist[itri].isurf != isurf) continue;
	  v0 = vlist[tlist[itri].vert[0]];
	  v1 = vlist[tlist[itri].vert[1]];
	  v2 = vlist[tlist[itri].vert[2]];
	  if (point_tri_distance(xyz[i],v0,v1,v2,
				 tlist[itri].normal,distance)) {
	    if (distance < near) {
	      near = distance;
	      plist[ip].itri = itri;
	    }
	  }
	} else {
	  if (region2surf[-itri-1] != isurf) continue;
	  distance = rlist[-itri-1]->distance(xyz[i]);
	  if (distance < near) {
	    near = distance;
	    plist[ip].itri = itri;
	  }
	}
      }
      if (near == BIG) {
	char str[128];
	sprintf(str,"Particle %d's triangle does not exist",i+1);
	error->one(str);
      }
    }
  }
    
  // free memory
  
  memory->destroy_2d_double_array(xyz);

  // print stats

  if (me == 0) {
    if (screen) fprintf(screen,"  %d particles\n",n);
    if (logfile) fprintf(logfile,"  %d particles\n",n);
  }

  // check if all particles were assigned to procs

  int nsum;
  MPI_Allreduce(&nlocal,&nsum,1,MPI_INT,MPI_SUM,world);

  if (nprevious + n != nsum) {
    char str[128];
    sprintf(str,"%d particles assigned to procs",nsum-nprevious);
    error->all(str);
  }
}

/* ----------------------------------------------------------------------
   migrate particles after particle motion
   after motion, each particle stores local bin its coords are in
     can be local owned bin or local ghost bin
   final outcome of migration:
     all particles owned by proc that owns particle's bin
     all procs have copies of upwind ghost particles
     particles across PBC have their coords and bin remapped
   both migrating to new procs and ghosting are done in one comm stage
   ghost bin migrate lists must store downwind neighbor info
------------------------------------------------------------------------- */

void Particle::migrate()
{
  int i,ibin;

  // debug test to see if all particles are in correct bins
  // should be unnecessary

  int iflag = 0;
  for (i = 0; i < ntotal; i++) {
    if (plist[i].ibin < 0 || plist[i].ibin >= grid->nbins) {
      iflag++;
      continue;
    }
    ibin = grid->whichlocal(plist[i].x);
    if (ibin != plist[i].ibin) iflag++;
  }
  if (iflag) {
    char str[128];
    sprintf(str,"Pre-migrate: %d particles on proc %d are not in correct bin",
	    iflag,me);
    error->one(str);
  }

  // loop over owned particles
  // each particle has ibin = local bin ID matching its current coords
  // perform any migrate actions this bin stores
  // if ibin is local ghost bin, delete particle from this proc

  Grid::OneBin *blist = grid->blist;
  Grid::Migrate *ptr;

  i = 0;
  int nsend = 0;
  int ncopy = 0;
  while (i < nlocal) {
    ibin = plist[i].ibin;
    ptr = blist[ibin].migrate;
    while (ptr) {
      if (ptr->proc == me) fill_pc(i,ptr,&ncopy);
      else fill_pm(i,ptr,&nsend);
      ptr = ptr->ptr;
    }
    if (blist[ibin].ghost) {
      memcpy(&plist[i],&plist[nlocal-1],sizeof(OnePart));
      nlocal--;
    } else i++;
  }

  ntotal = nlocal;

  // create Zoltan plan, perform unstructured communication, destroy plan
  // insure buf2 is large enough to hold incoming data

  int nrecv;
  Zoltan_Comm_Obj *plan;
  Zoltan_Comm_Create(&plan,nsend,proclist,world,0,&nrecv);

  if (nrecv > size2) {
    memory->sfree(buf2);
    size2 = nrecv;
    buf2 = (Migrate *) memory->smalloc(size2*sizeof(Migrate),"particle:buf2");
  }

  Zoltan_Comm_Do(plan,0,(char *) buf1,sizeof(Migrate),(char *) buf2);
  Zoltan_Comm_Destroy(&plan);

  // add copied and received particles to plist
  // add owned particles 1st, ghost particles last
  // 1st pass = copied owned particles (skips ghost)
  // 2nd pass = received owned particles (skips ghost)
  // 3rd pass = copied ghost particles (skips owned)
  // 4th pass = received ghost particles (skips owned)

  unpack(ncopy,buf3,0);
  unpack(nrecv,buf2,0);
  unpack(ncopy,buf3,1);
  unpack(nrecv,buf2,1);

  nghost_last = nghost;

  // debug test to see if all particles are in correct bins
  // should be unnecessary

  iflag = 0;
  for (i = 0; i < ntotal; i++) {
    if (plist[i].ibin < 0 || plist[i].ibin >= grid->nbins) {
      iflag++;
      continue;
    }
    ibin = grid->whichlocal(plist[i].x);
    if (ibin != plist[i].ibin) iflag++;
  }
  if (iflag) {
    char str[128];
    sprintf(str,"Post-migrate: %d particles on proc %d are not in correct bin",
	    iflag,me);
    error->one(str);
  }
}

/* ----------------------------------------------------------------------
   unpack N copied/received particles from Migrate buffer and add to plist
   ghost = 0/1 if unpacking owned/ghost particles, skip the others
------------------------------------------------------------------------- */

void Particle::unpack(int n, Migrate *buf, int ghostflag)
{
  int m;
  Grid::OneBin *blist = grid->blist;

  for (int i = 0; i < n; i++) {
    if (blist[buf[i].ibin].ghost != ghostflag) continue;
    if (ghostflag) nghost++;
    else nlocal++;
    m = add(buf[i].species,buf[i].x[0],buf[i].x[1],buf[i].x[2]);
    plist[m].seed = buf[i].seed;
    plist[m].ibin = buf[i].ibin;
    plist[m].itri = buf[i].itri;
  }
}

/* ----------------------------------------------------------------------
   add Ith particle from plist to Migrate send buffer
------------------------------------------------------------------------- */

void Particle::fill_pm(int i, Grid::Migrate *ptr, int *pn)
{
  int n = *pn;

  // grow buffer and proclist by 2x if necessary

  if (n >= size1) {
    if (size1 == 0) size1 = 1000;
    else size1 = 2*n;
    buf1 = (Migrate *)
      memory->srealloc(buf1,size1*sizeof(Migrate),"particle:buf1");
  }
  if (n >= sizeproc) {
    if (sizeproc == 0) sizeproc = 1000;
    else sizeproc = 2*n;
    proclist = (int *)
      memory->srealloc(proclist,sizeproc*sizeof(int),"particle:proclist");
  }

  // proc to send to

  proclist[n] = ptr->proc;

  // remap particle across PBC if necessary

  if (ptr->pbc) {
    buf1[n].x[0] = plist[i].x[0] + ptr->xpbc * domain->xsize;
    buf1[n].x[1] = plist[i].x[1] + ptr->ypbc * domain->ysize;
    buf1[n].x[2] = plist[i].x[2] + ptr->zpbc * domain->zsize;
  } else {
    buf1[n].x[0] = plist[i].x[0];
    buf1[n].x[1] = plist[i].x[1];
    buf1[n].x[2] = plist[i].x[2];
  }

  // rest of particle attributes

  buf1[n].species = plist[i].species;
  buf1[n].seed = plist[i].seed;
  buf1[n].ibin = ptr->ibin;
  buf1[n].itri = plist[i].itri;

  *pn = n+1;
}

/* ----------------------------------------------------------------------
   add Ith particle from plist to Migrate on-processor copy buffer
------------------------------------------------------------------------- */

void Particle::fill_pc(int i, Grid::Migrate *ptr, int *pn)
{
  int n = *pn;

  // grow buffer by 2x if necessary

  if (n >= size3) {
    if (size3 == 0) size3 = 1000;
    else size3 = 2*n;
    buf3 = (Migrate *)
      memory->srealloc(buf3,size3*sizeof(Migrate),"particle:buf3");
  }

  // remap particle across PBC if necessary

  if (ptr->pbc) {
    buf3[n].x[0] = plist[i].x[0] + ptr->xpbc * domain->xsize;
    buf3[n].x[1] = plist[i].x[1] + ptr->ypbc * domain->ysize;
    buf3[n].x[2] = plist[i].x[2] + ptr->zpbc * domain->zsize;
  } else {
    buf3[n].x[0] = plist[i].x[0];
    buf3[n].x[1] = plist[i].x[1];
    buf3[n].x[2] = plist[i].x[2];
  }

  // rest of particle attributes

  buf3[n].species = plist[i].species;
  buf3[n].seed = plist[i].seed;
  buf3[n].ibin = ptr->ibin;
  buf3[n].itri = plist[i].itri;

  *pn = n+1;
}

/* ----------------------------------------------------------------------
   swap ghost atoms across PBC remapping their coords
   serial routine that assumes 1 proc owns entire domain
   set seed of ghost atoms to index of original image atom
------------------------------------------------------------------------- */

void Particle::ghost_acquire()
{
  int ibin,ix,iy,iz,xpbc,ypbc,zpbc,species,itri,m;
  double x,y,z;

  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;
  double xsize = domain->xsize;
  double ysize = domain->ysize;
  double zsize = domain->zsize;
  int bx = grid->nbinx - 2;
  int by = grid->nbiny - 2;
  int bz = grid->nbinz - 2;

  for (int i = 0; i < nlocal; i++) {
    ibin = plist[i].ibin;
    grid->local2local(ibin,&ix,&iy,&iz);

    xpbc = ypbc = zpbc = 0;
    if (xperiodic) {
      if (ix == 1) xpbc = 1;
      if (ix == bx) xpbc = -1;
    }
    if (yperiodic) {
      if (iy == 1) ypbc = 1;
      if (iy == by) ypbc = -1;
    }
    if (zperiodic) {
      if (iz == 1) zpbc = 1;
      if (iz == bz) zpbc = -1;
    }

    if (!xpbc && !ypbc && !zpbc) continue;
    
    species = plist[i].species;
    x = plist[i].x[0];
    y = plist[i].x[1];
    z = plist[i].x[2];
    itri = plist[i].itri;

    if (xpbc) {
      m = add(species,x+xpbc*xsize,y,z);
      nghost++;
      plist[m].seed = i;
      plist[m].ibin = ibin + xpbc*bx;
      plist[m].itri = itri;
    }
    if (ypbc) {
      m = add(species,x,y+ypbc*ysize,z);
      nghost++;
      plist[m].seed = i;
      plist[m].ibin = ibin + ypbc*by;
      plist[m].itri = itri;
    }
    if (zpbc) {
      m = add(species,x,y,z+zpbc*zsize);
      nghost++;
      plist[m].seed = i;
      plist[m].ibin = ibin + zpbc*bz;
      plist[m].itri = itri;
    }
    if (xpbc && ypbc) {
      m = add(species,x+xpbc*xsize,y+ypbc*ysize,z);
      nghost++;
      plist[m].seed = i;
      plist[m].ibin = ibin + xpbc*bx * ypbc*by;
      plist[m].itri = itri;
    }
    if (xpbc && zpbc) {
      m = add(species,x+xpbc*xsize,y,z+zpbc*zsize);
      nghost++;
      plist[m].seed = i;
      plist[m].ibin = ibin + xpbc*bx * zpbc*bz;
      plist[m].itri = itri;
    }
    if (ypbc && zpbc) {
      m = add(species,x,y+ypbc*ysize,z+zpbc*zsize);
      nghost++;
      plist[m].seed = i;
      plist[m].ibin = ibin + ypbc*by * zpbc*bz;
      plist[m].itri = itri;
    }
    if (xpbc && ypbc && zpbc) {
      m = add(species,x+xpbc*xsize,y+ypbc*ysize,z+zpbc*zsize);
      nghost++;
      plist[m].seed = i;
      plist[m].ibin = ibin + xpbc*bx + ypbc*by * zpbc*bz;
      plist[m].itri = itri;
    }
  }

  nghost_last = nghost;
}

/* ----------------------------------------------------------------------
   compact the particle list, owned and ghost
   if flag = 0, keep it
   if flag = -1, delete it (spent reactant)
   if flag = 1, delete it (ghost)
   reset nlocal, nghost, ntotal
------------------------------------------------------------------------- */

void Particle::compact()
{
  int i = 0;
  while (i < ntotal) {
    if (plist[i].flag) {
      memcpy(&plist[i],&plist[ntotal-1],sizeof(OnePart));
      ntotal--;
    } else i++;
  }
  nghost = 0;
  nlocal = ntotal;

  // debug test to see if all particles are in my owned domain
  // should be unnecessary

  Grid::OneBin *blist = grid->blist;
  int iflag = 0;
  for (i = 0; i < nlocal; i++)
    if (blist[plist[i].ibin].ghost == 1) iflag++;
  if (iflag) {
    char str[128];
    sprintf(str,
	    "%d local particles on proc %d are in ghost cells after compact",
	    iflag,me);
    error->one(str);
  }
}

/* ----------------------------------------------------------------------
   if pattern matches str, return 1, else return 0
   pattern can have one or zero wildcard char *
   5 cases handled: no *, *, *ab, ab*, ab*cd
------------------------------------------------------------------------- */

int Particle::match(char *pattern, char *str)
{
  if (strchr(pattern,'*') == NULL) {
    if (strcmp(pattern,str) == 0) return 1;
    else return 0;
  }

  if (strchr(pattern,'*') != strrchr(pattern,'*'))
    error->all("Cannot perform match with two or more wildcard *");
  
  // pre and post = strings before and after the *

  char *wild = strchr(pattern,'*');
  char *pre = pattern;
  char *post = wild + 1;

  int flag = 1;
  *wild = '\0';
  if (strstr(str,pre) != str) flag = 0;
   if (strstr(&str[strlen(str)-strlen(post)],post) == NULL) flag = 0;
  *wild = '*';
  return flag;
}

/* ----------------------------------------------------------------------
   match pattern against all species and alias names
   pattern can contain a wildcard character
   set flag to 0/1 for each species depending on match
   return # of species matches
------------------------------------------------------------------------- */

int Particle::allmatch(char *pattern, int *flag)
{
  for (int i = 0; i < nspecies; i++) flag[i] = match(pattern,name[i]);
  for (int i = 0; i < nalias; i++)
    if (match(pattern,alias[i])) flag[alias2name[i]] = 1;

  int count = 0;
  for (int i = 0; i < nspecies; i++) if (flag[i]) count++;
  return count;
}

/* ---------------------------------------------------------------------- */

int Particle::pack_restart(double *buf)
{
  int m = 0;
  for (int i = 0; i < nlocal; i++) {
    buf[m++] = plist[i].x[0];
    buf[m++] = plist[i].x[1];
    buf[m++] = plist[i].x[2];
    buf[m++] = plist[i].species;
    buf[m++] = plist[i].itri;
    buf[m++] = plist[i].seed;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void Particle::unpack_restart(int n, double *buf)
{
  double xyz[3];
  int ispecies,itri,seed,iglobal,ibin,ip,ix,iy,iz;

  Grid::OneBin *blist = grid->blist;
  n /= size_restart;
  int m = 0;

  for (int i = 0; i < n; i++) {
    xyz[0] = buf[m++];
    xyz[1] = buf[m++];
    xyz[2] = buf[m++];
    ispecies = static_cast<int> (buf[m++]);
    itri = static_cast<int> (buf[m++]);
    seed = static_cast<int> (buf[m++]);

    if (!domain->inside(xyz))
      error->one("Particle is outside global domain");

    iglobal = grid->whichglobal(xyz,&ix,&iy,&iz);
    ibin = grid->global2local(ix,iy,iz);
    if (ibin == -1 || blist[ibin].ghost) continue;

    ip = add(ispecies,xyz[0],xyz[1],xyz[2]);
    nlocal++;
    plist[ip].itri = itri;
    plist[ip].ibin = ibin;
    plist[ip].seed = seed;
  }
}

/* ----------------------------------------------------------------------
   return size of plist array
------------------------------------------------------------------------- */

int Particle::memory_usage()
{
  return maxpart*sizeof(OnePart);
}
