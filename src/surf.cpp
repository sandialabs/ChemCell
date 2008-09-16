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

#include "math.h"
#include "string.h"
#include "mpi.h"
#include "surf.h"
#include "simulator.h"
#include "domain.h"
#include "region.h"
#include "grid.h"
#include "move.h"
#include "particle.h"
#include "geometry.h"
#include "memory.h"
#include "error.h"

#define RegionInclude
#include "style.h"
#undef RegionInclude

#define MAXLINE 256

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

Surf::Surf()
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  nsurf = 0;
  name = NULL;
  surf2region = NULL;
  nregion = 0;
  ntri = nvert = 0;

  tlist = NULL;
  vlist = NULL;
  rlist = NULL;
  region2surf = NULL;
}

/* ---------------------------------------------------------------------- */

Surf::~Surf()
{
  for (int i = 0; i < nsurf; i++) delete [] name[i];
  memory->sfree(name);
  memory->sfree(surf2region);

  memory->sfree(tlist);
  memory->destroy_2d_double_array(vlist);
  for (int i = 0; i < nregion; i++) delete rlist[i];
  memory->sfree(rlist);
  memory->sfree(region2surf);
}

/* ---------------------------------------------------------------------- */

void Surf::init()
{
  // print surf stats

  int count = 0;
  for (int i = 0; i < grid->nbins; i++) {
    if (grid->blist[i].id == -1) continue;
    count += grid->blist[i].ntri;
  }
  int all;
  MPI_Allreduce(&count,&all,1,MPI_INT,MPI_SUM,world);
  double tbave = 1.0*all/nprocs/grid->gbins;

  count = 0;
  for (int i = 0; i < grid->nbins; i++) {
    if (grid->blist[i].id == -1) continue;
    count = MAX(count,grid->blist[i].ntri);
  }
  int tbmax;
  MPI_Allreduce(&count,&tbmax,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {
    if (screen) {
      fprintf(screen,"Surfs:\n");
      fprintf(screen,"  surfs regions triangles vertices = %d %d %d %d\n",
	      nsurf,nregion,ntri,nvert);
      fprintf(screen,"  tri/bin: ave max = %g %d\n",tbave,tbmax);
    }
    if (logfile) {
      fprintf(logfile,"Surfs:\n");
      fprintf(logfile,"  surfs regions triangles vertices = %d %d %d %d\n",
	      nsurf,nregion,ntri,nvert);
      fprintf(logfile,"  tri/bin: ave max = %g %d\n",tbave,tbmax);
    }
  }
}

/* ----------------------------------------------------------------------
   return surf index (0 to N-1) for surf-ID str
   return -1 if doesn't exist
------------------------------------------------------------------------- */

int Surf::find(char *str)
{
  for (int isurf = 0; isurf < nsurf; isurf++)
    if (strcmp(str,name[isurf]) == 0) return isurf;
  return -1;
}

/* ----------------------------------------------------------------------
   read triangulated surf and assign triangles to bins
------------------------------------------------------------------------- */

void Surf::read_triangles(int narg, char **arg)
{
  int i,j,k,m;

  char line[MAXLINE];
  char *err;

  if (me == 0) {
    if (screen) fprintf(screen,"Reading surface ...\n");
    if (logfile) fprintf(logfile,"Reading surface ...\n");
  }

  if (narg != 3) error->all("Illegal surface command");

  // check surf-ID
  // if doesn't exist, create new surf and extend/init permeability arrays
  // if exists, triangles are just added to it

  int isurf = find(arg[0]);
  int nv = atoi(arg[1]);
  int nt = atoi(arg[2]);

  if (isurf >= 0 && surf2region[isurf] >= 0)
    error->all("Cannot add triangles to a region surface");

  if (isurf == -1) {
    isurf = add_surf(arg[0]);
    move->grow_partsurf(particle->maxspecies,nsurf);
    move->default_partsurf(-1,isurf);
  }

  // read and bcast list of vertices
  // store coords in vlist
  
  vlist = memory->grow_2d_double_array(vlist,nvert+nv,3,"surf:vlist");

  if (me == 0) {
    err = fgets(line,MAXLINE,infile);
    if (err == NULL) error->one("Unexpected end of file");
    int tmp;
    j = nvert;
    for (i = 0; i < nv; i++) {
      err = fgets(line,MAXLINE,infile);
      sscanf(line,"%d %lg %lg %lg",&tmp,
	     &vlist[j][0],&vlist[j][1],&vlist[j][2]);
      j++;
    }
    if (err == NULL) error->one("Unexpected end of file");
  }
    
  MPI_Bcast(vlist[nvert],3*nv,MPI_DOUBLE,0,world);

  // read and bcast list of triangles

  int **tri = memory->create_2d_int_array(nt,3,"surf:tri");

  if (me == 0) {
    err = fgets(line,MAXLINE,infile);
    if (err == NULL) error->one("Unexpected end of file");
    int tmp;
    for (i = 0; i < nt; i++) {
      err = fgets(line,MAXLINE,infile);
      sscanf(line,"%d %d %d %d",&tmp,&tri[i][0],&tri[i][1],&tri[i][2]);
    }
    if (err == NULL) error->one("Unexpected end of file");
  }
  
  MPI_Bcast(tri[0],3*nt,MPI_INT,0,world);

  // store isurf and vertices in tlist
  // offset by pre-existing vertices
  // -1 to convert to C index

  tlist = (OneTri *) 
    memory->srealloc(tlist,(ntri+nt)*sizeof(OneTri),"surf:tlist");

  j = ntri;
  for (i = 0; i < nt; i++) {
    tlist[j].isurf = isurf;
    tlist[j].vert[0] = nvert + tri[i][0] - 1;
    tlist[j].vert[1] = nvert + tri[i][1] - 1;
    tlist[j].vert[2] = nvert + tri[i][2] - 1;
    j++;
  }

  for (i = ntri; i < nt; i++) tri_normal(i);

  memory->destroy_2d_int_array(tri);

  // read and bcast list of connections
  
  int **connect = memory->create_2d_int_array(nt,6,"surf:connect");
  
  if (me == 0) {
    err = fgets(line,MAXLINE,infile);
    if (err == NULL) error->one("Unexpected end of file");
    int tmp;
    for (i = 0; i < nt; i++) {
      err = fgets(line,MAXLINE,infile);
      sscanf(line,"%d %d %d %d %d %d %d",&tmp,
	     &connect[i][0],&connect[i][1],&connect[i][2],
	     &connect[i][3],&connect[i][4],&connect[i][5]);
    }
    if (err == NULL) error->one("Unexpected end of file");
  }
  
  MPI_Bcast(connect[0],6*nt,MPI_INT,0,world);

  // store connections in tlist
  // offset by pre-existing triangles
  // -1 to convert to C index (0 converted to -1 if no connection)

  j = ntri;
  for (i = 0; i < nt; i++) {
    if (connect[i][0]) tlist[j].connect[0] = ntri + connect[i][0] - 1;
    else tlist[j].connect[0] = -1;
    if (connect[i][2]) tlist[j].connect[1] = ntri + connect[i][2] - 1;
    else tlist[j].connect[1] = -1;
    if (connect[i][4]) tlist[j].connect[2] = ntri + connect[i][4] - 1;
    else tlist[j].connect[2] = -1;
    j++;
  }

  memory->destroy_2d_int_array(connect);

  // check that all vertices are inside global domain

  for (i = 0; i < nv; i++)
    if (!domain->inside(vlist[nvert+i])) {
      char str[128];
      sprintf(str,"Vertex %d in surf %s is outside global domain",i+1,arg[0]);
      error->all(str);
    }

  // scan triangles for those that intersect any of my bins (own or ghost)

  int gbinx = grid->gbinx;
  int gbiny = grid->gbiny;
  int gbinz = grid->gbinz;

  double *v1,*v2,*v3;
  double lo[3],hi[3];
  int imin,jmin,kmin,imax,jmax,kmax,ibin,ix,iy,iz;
  int xflag,yflag,zflag;
  xflag = yflag = zflag = 0;

  for (m = ntri; m < ntri+nt; m++) {

    // lo/hi = bounding box around tri

    v1 = vlist[tlist[m].vert[0]];
    v2 = vlist[tlist[m].vert[1]];
    v3 = vlist[tlist[m].vert[2]];

    lo[0] = MIN(v1[0],v2[0]);
    lo[0] = MIN(lo[0],v3[0]);
    hi[0] = MAX(v1[0],v2[0]);
    hi[0] = MAX(hi[0],v3[0]);

    lo[1] = MIN(v1[1],v2[1]);
    lo[1] = MIN(lo[1],v3[1]);
    hi[1] = MAX(v1[1],v2[1]);
    hi[1] = MAX(hi[1],v3[1]);

    lo[2] = MIN(v1[2],v2[2]);
    lo[2] = MIN(lo[2],v3[2]);
    hi[2] = MAX(v1[2],v2[2]);
    hi[2] = MAX(hi[2],v3[2]);

    // find subset of my bins that are overlapped by tri bounding box

    if (!grid->box_overlap(lo,hi,imin,jmin,kmin,imax,jmax,kmax)) continue;

    // loop over possible intersecting bins

    for (k = kmin; k <= kmax; k++)
      for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++) {
	  grid->coord_corners(i,j,k,lo,hi);
	  if (tri_hex_intersect(v1,v2,v3,tlist[m].normal,lo,hi)) {
	    ibin = grid->local2one(i,j,k);
	    grid->add_tri(ibin,m);
	    grid->local2global(ibin,&ix,&iy,&iz);
	    if (ix <= 0 || ix >= gbinx-1) xflag = 1;
	    if (iy <= 0 || iy >= gbiny-1) yflag = 1;
	    if (iz <= 0 || iz >= gbinz-1) zflag = 1;
	  }
	}
  }

  // error if any tri intersects bin next to periodic boundary

  int flag = 0;
  if (domain->xperiodic && xflag) flag = 1;
  if (domain->yperiodic && yflag) flag = 1;
  if (domain->zperiodic && zflag) flag = 1;

  int all;
  MPI_Allreduce(&flag,&all,1,MPI_INT,MPI_SUM,world);
  if (all) error->all("Triangle intersects bin next to periodic boundary");

  // update counters

  nvert += nv;
  ntri += nt;

  // print stats

  if (me == 0) {
    if (screen) fprintf(screen,"  %d vertices & %d triangles\n",nv,nt);
    if (logfile) fprintf(logfile,"  %d vertices & %d triangles\n",nv,nt);
  }
}

/* ----------------------------------------------------------------------
   add a region to rlist and intersect it with my bins (owned and ghost)
------------------------------------------------------------------------- */

void Surf::add_region(int narg, char **arg)
{
  if (narg < 2) error->all("Illegal region command");

  // check surf-ID to insure doesn't exist
  // create new surf and extend/init permeability arrays

  int isurf = find(arg[0]);
  if (isurf >= 0) error->all("Region surf-ID already exists");
  isurf = add_surf(arg[0]);
  move->grow_partsurf(particle->maxspecies,nsurf);
  move->default_partsurf(-1,isurf);

  // extend rlist and region2surf

  rlist = (Region **) 
    memory->srealloc(rlist,(nregion+1)*sizeof(Region *),"surf:rlist");
  region2surf = (int *) 
    memory->srealloc(region2surf,(nregion+1)*sizeof(int),"surf:region2surf");
  
  // create the new region

  if (strcmp(arg[1],"none") == 0) error->all("Invalid region style");

#define RegionClass
#define RegionStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) \
    rlist[nregion] = new Class(narg-1,&arg[1]);
#include "style.h"
#undef RegionClass

  else error->all("Invalid region style");

  region2surf[nregion] = isurf;
  surf2region[isurf] = nregion;

  // check that region does not extend beyond global domain

  double lo[3],hi[3];
  rlist[nregion]->bbox(lo,hi);

  if (!domain->inside(lo) || !domain->inside(hi))
    error->all("Region extends outside global domain");

  // flag each of my owned or ghost bins that region surface intersects
  // add region to bins that it overlaps

  int i,j,k,ibin,ix,iy,iz;

  int nbinx = grid->nbinx;
  int nbiny = grid->nbiny;
  int nbinz = grid->nbinz;
  int gbinx = grid->gbinx;
  int gbiny = grid->gbiny;
  int gbinz = grid->gbinz;

  int iregion = -nregion - 1;
  int xflag,yflag,zflag;
  xflag = yflag = zflag = 0;
  
  for (k = 0; k < nbinz; k++)
    for (j = 0; j < nbiny; j++)
      for (i = 0; i < nbinx; i++) {
	grid->coord_corners(i,j,k,lo,hi);
	if (rlist[nregion]->hex_intersect(lo,hi)) {
	  ibin = grid->local2one(i,j,k);
	  grid->add_tri(ibin,iregion);
	  grid->local2global(ibin,&ix,&iy,&iz);
	  if (ix <= 0 || ix >= gbinx-1) xflag = 1;
	  if (iy <= 0 || iy >= gbiny-1) yflag = 1;
	  if (iz <= 0 || iz >= gbinz-1) zflag = 1;
	}
      }

  // error if region intersects bin next to periodic boundary

  double pt[3],norm[3];

  int flag = 0;
  if (domain->xperiodic && xflag) {
    if (strcmp(rlist[nregion]->style,"plane") == 0) {
      pt[0] = pt[1] = pt[2] = 0.0;
      rlist[nregion]->compute_normal(pt,norm);
      pt[0] = 1.0;
      if (dot(norm,pt) != 0.0) flag = 1;
    } else if (strcmp(rlist[nregion]->style,"cylinder") == 0) {
      if (((RegionCylinder *) rlist[nregion])->axis != 0) flag = 1;
    } else flag = 1;
  }
  if (domain->yperiodic && yflag) {
    if (strcmp(rlist[nregion]->style,"plane") == 0) {
      pt[0] = pt[1] = pt[2] = 0.0;
      rlist[nregion]->compute_normal(pt,norm);
      pt[1] = 1.0;
      if (dot(norm,pt) != 0.0) flag = 1;
    } else if (strcmp(rlist[nregion]->style,"cylinder") == 0) {
      if (((RegionCylinder *) rlist[nregion])->axis != 1) flag = 1;
    } else flag = 1;
  }
  if (domain->zperiodic && zflag) {
    if (strcmp(rlist[nregion]->style,"plane") == 0) {
      pt[0] = pt[1] = pt[2] = 0.0;
      rlist[nregion]->compute_normal(pt,norm);
      pt[2] = 1.0;
      if (dot(norm,pt) != 0.0) flag = 1;
    } else if (strcmp(rlist[nregion]->style,"cylinder") == 0) {
      if (((RegionCylinder *) rlist[nregion])->axis != 2) flag = 1;
    } else flag = 1;
  }

  int all;
  MPI_Allreduce(&flag,&all,1,MPI_INT,MPI_SUM,world);
  if (all) error->all("Region intersects bin next to periodic boundary");

  nregion++;
}

/* ----------------------------------------------------------------------
   return size of tlist and vlist arrays
------------------------------------------------------------------------- */

int Surf::memory_usage()
{
  int bytes = ntri*sizeof(OneTri) + nvert*3*sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   create a new surf by adding str surf-ID to name
   assume it does not point to a region
   return index of added surf
------------------------------------------------------------------------- */

int Surf::add_surf(char *str)
{
  name = (char **)
    memory->srealloc(name,(nsurf+1)*sizeof(char *),"surf:name");
  surf2region = (int *)
    memory->srealloc(surf2region,(nsurf+1)*sizeof(int),"surf:surf2region");

  int nlen = strlen(str) + 1;
  name[nsurf] = new char[nlen];
  strcpy(name[nsurf],str);

  surf2region[nsurf] = -1;

  nsurf++;
  return nsurf-1;
}

/* ----------------------------------------------------------------------
   compute unit normal of triangle I
   assumes vert has local IDs of 3 vertices
------------------------------------------------------------------------- */

void Surf::tri_normal(int i)
{
  double *v0 = vlist[tlist[i].vert[0]];
  double *v1 = vlist[tlist[i].vert[1]];
  double *v2 = vlist[tlist[i].vert[2]];

  double x01 = v1[0] - v0[0];
  double y01 = v1[1] - v0[1];
  double z01 = v1[2] - v0[2];
  double x02 = v2[0] - v0[0];
  double y02 = v2[1] - v0[1];
  double z02 = v2[2] - v0[2];

  double xn = y01*z02 - z01*y02;
  double yn = z01*x02 - x01*z02;
  double zn = x01*y02 - y01*x02;

  double len = sqrt(xn*xn + yn*yn + zn*zn);

  tlist[i].normal[0] = xn/len;
  tlist[i].normal[1] = yn/len;
  tlist[i].normal[2] = zn/len;
}
