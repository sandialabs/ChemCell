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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "grid.h"
#include "simulator.h"
#include "domain.h"
#include "particle.h"
#include "chem.h"
#include "move.h"
#include "surf.h"
#include "balance.h"
#include "react.h"
#include "memory.h"
#include "error.h"

#define SURF          1
#define PARTICLE      2
#define SURF_PARTICLE 3

#define STATIC_BRICK 1
#define STATIC_BIN   2
#define STATIC_PART  3
#define DYNAMIC_PART 4

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define CHUNK 10000
#define BIG 1000000000

/* ---------------------------------------------------------------------- */

Grid::Grid()
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  setflag = 0;
  balance_style = STATIC_BRICK;
  balance_file = 0;
  fpbal = NULL;

  blist = NULL;

  size1 = size2 = size3 = size4 = 0;
  buf1 = buf2 = buf3 = buf4 = NULL;
  sizeproc1 = sizeproc2 = 0;
  proclist1 = proclist2 = NULL;

  list = (int *) memory->smalloc(nprocs*sizeof(int),"grid:list");
  x2proc = y2proc = z2proc = NULL;
  bininfo = (int **) memory->create_2d_int_array(nprocs,6,"grid:bininfo");

  nichunks = nmchunks = 0;
  ichunks = NULL;
  mchunks = NULL;
  ifreelist = NULL;
  mfreelist = NULL;
}

/* ---------------------------------------------------------------------- */

Grid::~Grid()
{
  if (me == 0 && fpbal) fclose(fpbal);

  memory->sfree(blist);

  memory->sfree(buf1);
  memory->sfree(buf2);
  memory->sfree(buf3);
  memory->sfree(buf4);
  memory->sfree(proclist1);
  memory->sfree(proclist2);

  memory->sfree(list);
  if (x2proc) delete [] x2proc;
  if (y2proc) delete [] y2proc;
  if (z2proc) delete [] z2proc;
  memory->destroy_2d_int_array(bininfo);

  for (int i = 0; i < nichunks; i++) memory->sfree(ichunks[i]);
  memory->sfree(ichunks);
  for (int i = 0; i < nmchunks; i++) memory->sfree(mchunks[i]);
  memory->sfree(mchunks);
}

/* ---------------------------------------------------------------------- */

void Grid::init()
{
  if (setflag == 0) error->all("Must set bins before run");

  // perform initial bin balancing
  // style = STATIC_BRICK, do nothing
  // style = STATIC_BIN, only balance if never have before
  // style = STATIC_PART, balance
  // style = DYNAMIC_PART, balance if never have before or if out-of-balance

  if (balance_style == STATIC_BRICK) {
    next_balance = simulator->laststep + 1;
  } else if (balance_style == STATIC_BIN) {
    weightflag = 0;
    if (decomp_style == STATIC_BRICK) rebalance();
    decomp_style = STATIC_BIN;
    next_balance = simulator->laststep + 1;
  } else if (balance_style == STATIC_PART) {
    weightflag = 1;
    rebalance();
    decomp_style = STATIC_PART;
    next_balance = simulator->laststep + 1;
  } else if (balance_style == DYNAMIC_PART) {
    weightflag = 1;
    if (decomp_style == STATIC_BRICK) rebalance();
    else if (check()) rebalance();
    decomp_style = DYNAMIC_PART;
    int ntimestep = simulator->ntimestep;
    next_balance = ntimestep + nbalance;
    if (nbalance && ntimestep % nbalance != 0)
      next_balance = (ntimestep/nbalance)*nbalance + nbalance;
  }

  // zero out count for calls to balancer

  ncount = 0;

  // print bin stats

  if (me == 0) {
    if (screen) {
      fprintf(screen,"Geometry Bins:\n");
      fprintf(screen,"  xyz size = %g %g %g\n",xbinsize,ybinsize,zbinsize);
      fprintf(screen,"  xyz counts = %d %d %d = %d total\n",
	      gbinx,gbiny,gbinz,gbins);
    }
    if (logfile) {
      fprintf(logfile,"Geometry Bins:\n");
      fprintf(logfile,"  xyz size = %g %g %g\n",xbinsize,ybinsize,zbinsize);
      fprintf(logfile,"  xyz counts = %d %d %d = %d total\n",
	      gbinx,gbiny,gbinz,gbins);
    }
  }
}

/* ----------------------------------------------------------------------
   check if dynamic load-balancing needs to be done
   test current particle imbalance against threshhold
   return 1 if threshhold exceeded, 0 if not
------------------------------------------------------------------------- */

int Grid::check()
{
  if (nprocs == 1) return 0;

  int ntotal,nmax;
  MPI_Allreduce(&particle->nlocal,&ntotal,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&particle->nlocal,&nmax,1,MPI_INT,MPI_MAX,world);
  double ave = ntotal/nprocs;

  if (ntotal == 0) return 0;
  if (nmax/ave <= threshhold) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   perform dynamic balancing of bins
------------------------------------------------------------------------- */

void Grid::dynamic()
{
  // test particle count against threshhold
  // increment count and next timestep to balance on

  next_balance += nbalance;
  if (check() == 0) return;
  rebalance();
  ncount++;
}

/* ----------------------------------------------------------------------
   create bins for 1st time and assign them to procs
   also create reaction bins
------------------------------------------------------------------------- */

void Grid::create(int narg, char **arg)
{
  if (domain->setflag == 0) error->all("Must set global domain before bins");
  if (setflag == 1) error->all("Bins are already setup");

  if (me == 0) {
    if (screen) fprintf(screen,"Setting up bins ...\n");
    if (logfile) fprintf(logfile,"Setting up bins ...\n");
  }

  setflag = 1;
  decomp_style = STATIC_BRICK;

  // compute bin size, counts, gbin for global domain

  global(narg,arg);

  // construct a 3d brick decomposition
  // set nbin & offset info for my local sub-domain
  // create x2proc,y2proc,z2proc mappings that describe 3d brick decomp

  decomp3d();

  // allocate blist and initialize bin settings

  blist = (OneBin *) memory->smalloc(nbins*sizeof(OneBin),"grid:blist");
  topology();

  // allocate reaction bins and initialize reaction bin settings

  chem->create();

  // initialize particle's knowledge of reaction bins

  particle->topology();
}

/* ----------------------------------------------------------------------
   global bin creation
   sets gbin xyz, binsize, origin
------------------------------------------------------------------------- */

void Grid::global(int narg, char **arg)
{
  // calculate grid bin size and reaction bin size
  // diff, gsize, gcount induce grid bin size(s)
  // react, rsize, rcount induce reaction bin size(s)

  double value,xvalue,yvalue,zvalue;
  double gxsize,gysize,gzsize;
  double rxsize,rysize,rzsize;
  int nx,ny,nz,rperx,rpery,rperz;

  int gbinflag = 0;
  int rbinflag = 0;
  gxsize = gysize = gzsize = 0.0;
  rxsize = rysize = rzsize = 0.0;
  rperx = rpery = rperz = 0;

  int iarg = 0;
  while (iarg < narg) {

    if (strcmp(arg[iarg],"diff") == 0) {
      if (narg < iarg+2) error->all("Illegal bin command");
      value = atof(arg[iarg+1]);
      if (value < 0.0) error->all("Illegal bin command");
      value = move->maxbin(value);
      if (me == 0) {
	fprintf(screen,"  bin size induced by diff constraint = %g\n",value);
	fprintf(logfile,"  bin size induced by diff constraint = %g\n",value);
      }
      gxsize = MAX(gxsize,value);
      gysize = MAX(gysize,value);
      gzsize = MAX(gzsize,value);
      gbinflag = 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"gsize") == 0) {
      if (narg < iarg+4) error->all("Illegal bin command");
      xvalue = atof(arg[iarg+1]);
      yvalue = atof(arg[iarg+2]);
      zvalue = atof(arg[iarg+3]);
      if (xvalue <= 0.0 || yvalue <= 0.0 || zvalue <= 0.0)
	error->all("Illegal bin command");
      if (me == 0) {
	fprintf(screen,"  bin size induced by gsize constraint = %g %g %g\n",
		xvalue,yvalue,zvalue);
	fprintf(logfile,"  bin size induced by gsize constraint = %g %g %g\n",
		xvalue,yvalue,zvalue);
      }
      gxsize = MAX(gxsize,xvalue);
      gysize = MAX(gysize,yvalue);
      gzsize = MAX(gzsize,zvalue);
      gbinflag = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg],"gcount") == 0) {
      if (narg < iarg+4) error->all("Illegal bin command");
      nx = atoi(arg[iarg+1]);
      ny = atoi(arg[iarg+2]);
      nz = atoi(arg[iarg+3]);
      if (nx <= 0 || ny <= 0 || nz <= 0) error->all("Illegal bin command");
      xvalue = domain->xsize / nx;
      yvalue = domain->ysize / ny;
      zvalue = domain->zsize / nz;
      if (me == 0) {
	fprintf(screen,"  bin size induced by gcount constraint = %g %g %g\n",
		xvalue,yvalue,zvalue);
	fprintf(logfile,"  bin size induced by gcount constraint = %g %g %g\n",
		xvalue,yvalue,zvalue);
      }
      gxsize = MAX(gxsize,xvalue);
      gysize = MAX(gysize,yvalue);
      gzsize = MAX(gzsize,zvalue);
      gbinflag = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg],"react") == 0) {
      if (narg < iarg+2) error->all("Illegal bin command");
      value = atof(arg[iarg+1]);
      if (value < 0.0) error->all("Illegal bin command");
      if (value == 0.0) value = chem->maxbin();
      if (me == 0) {
	fprintf(screen,"  bin size induced by react constraint = %g\n",value);
	fprintf(logfile,"  bin size induced by react constraint = %g\n",value);
      }
      rxsize = MAX(rxsize,value);
      rysize = MAX(rysize,value);
      rzsize = MAX(rzsize,value);
      rbinflag = 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"rsize") == 0) {
      if (narg < iarg+4) error->all("Illegal bin command");
      xvalue = atof(arg[iarg+1]);
      yvalue = atof(arg[iarg+2]);
      zvalue = atof(arg[iarg+3]);
      if (xvalue <= 0.0 || yvalue <= 0.0 || zvalue <= 0.0)
	error->all("Illegal bin command");
      if (me == 0) {
	fprintf(screen,"  bin size induced by rsize constraint = %g %g %g\n",
		xvalue,yvalue,zvalue);
	fprintf(logfile,"  bin size induced by rsize constraint = %g %g %g\n",
		xvalue,yvalue,zvalue);
      }
      rxsize = MAX(rxsize,xvalue);
      rysize = MAX(rysize,yvalue);
      rzsize = MAX(rzsize,zvalue);
      rbinflag = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg],"rcount") == 0) {
      if (narg < iarg+4) error->all("Illegal bin command");
      nx = atoi(arg[iarg+1]);
      ny = atoi(arg[iarg+2]);
      nz = atoi(arg[iarg+3]);
      if (nx <= 0 || ny <= 0 || nz <= 0) error->all("Illegal bin command");
      if (me == 0) {
	fprintf(screen,"  bin count induced by rcount constraint = %d %d %d\n",
		nx,ny,nz);
	fprintf(logfile,"  bin count induced by rcount constraint = %d %d %d\n",
		nx,ny,nz);
      }
      rperx = nx;
      rpery = ny;
      rperz = nz;
      rbinflag = 1;
      iarg += 4;

    } else error->all("Illegal bin command");
  }

  // set grid bin size = xbinsize,ybinsize,zbinsize
  // calculate default size if grid keywords were not specified

  if (gbinflag == 0)
    xbinsize = ybinsize = zbinsize = move->maxbin(0.0);
  else {
    xbinsize = gxsize;
    ybinsize = gysize;
    zbinsize = gzsize;
  }

  if (xbinsize == 0.0 || ybinsize == 0.0 || zbinsize == 0.0)
    error->all("Geometry bin size of 0.0");

  // set actual grid bin size >= desired size
  // fit integer # of grid bins into global domain
  // must be at least one grid bin in a dimension
  // enforce multiple of 2 grid bins in a periodic dimension

  gbinx = static_cast<int> (domain->xsize / xbinsize);
  gbiny = static_cast<int> (domain->ysize / ybinsize);
  gbinz = static_cast<int> (domain->zsize / zbinsize);

  if (gbinx == 0) gbinx = 1;
  if (gbiny == 0) gbiny = 1;
  if (gbinz == 0) gbinz = 1;

  if (domain->xperiodic && (gbinx % 2) == 1) {
    if (gbinx == 1) 
      error->all("Must be 2 or more bins in periodic dimensions");
    gbinx--;
  }

  if (domain->yperiodic && (gbiny % 2) == 1) {
    if (gbiny == 1) 
      error->all("Must be 2 or more bins in periodic dimensions");
    gbiny--;
  }

  if (domain->zperiodic && (gbinz % 2) == 1) {
    if (gbinz == 1) 
      error->all("Must be 2 or more bins in periodic dimensions");
    gbinz--;
  }

  // gbinx,gbiny,gbinz = actual bin counts in interior of global domain
  
  gbins = gbinx*gbiny*gbinz;

  xbinsize = domain->xsize / gbinx;
  ybinsize = domain->ysize / gbiny;
  zbinsize = domain->zsize / gbinz;

  xbininv = 1.0/xbinsize;
  ybininv = 1.0/ybinsize;
  zbininv = 1.0/zbinsize;

  xbinhalf = 0.5*xbinsize;
  ybinhalf = 0.5*ybinsize;
  zbinhalf = 0.5*zbinsize;

  xorigin = domain->xlo;
  yorigin = domain->ylo;
  zorigin = domain->zlo;

  // set reaction bin counts per grid bin = nperx,npery,nperz
  // default = 1 if reaction keywords were not specified

  if (rbinflag == 0)
    nperx = npery = nperz = 1;
  else {
    nperx = npery = nperz = BIG;
    if (rxsize > 0.0) nperx = static_cast<int> (xbinsize/rxsize);
    if (rysize > 0.0) npery = static_cast<int> (ybinsize/rysize);
    if (rzsize > 0.0) nperz = static_cast<int> (zbinsize/rzsize);
    if (rperx) nperx = MIN(nperx,rperx);
    if (rpery) npery = MIN(npery,rpery);
    if (rperz) nperz = MIN(nperz,rperz);
  }
}

/* ----------------------------------------------------------------------
   initialize grid ID, proc, ghost from current decomp
   zeroes geometry and particle counters
   sets up migrate info from current decomp
------------------------------------------------------------------------- */

void Grid::topology()
{
  int i,j,k,n,n1;
  Grid::Migrate *ptr;

  // set id,ghost,proc
  // zero out geometry and particles

  for (int m = 0; m < nbins; m++) {
    blist[m].id = local2global(m);
    local2local(m,&i,&j,&k);
    if (i == 0 || i == nbinx-1 || j == 0 || j == nbiny-1 ||
	k == 0 || k == nbinz-1) {
      blist[m].ghost = 1;
      blist[m].proc = owner(i+xoffset,j+yoffset,k+zoffset);
    } else {
      blist[m].ghost = 0;
      blist[m].proc = me;
    }

    blist[m].ntri = 0;
    blist[m].tri = NULL;
    blist[m].nparts = 0;
    blist[m].first = -1;
  }

  // setup list of migrate actions for the bin and its particles
  // migrate actions perform 2 tasks:
  //   particles in ghost cells migrate to new owner
  //   procs acquire particles in upwind ghost cells
  // action = send (or copy) particle to another proc (or self)
  //   with particle coords and ibin for where that proc should store it
  //   both coords and ibin are remapped by PBC if necessary
  // if bin ID = -1, it has no actions (particle migrate just deletes it)
  // bins with migrate actions:
  //   all ghost bins
  //   owned bins at lower left of owned domain (i or j or k = 1)
  //   lower left restriction is b/c procs acquire only upwind ghosts
  // each master bin has 8 possible actions, one for each of 8 stencil bins:
  //   stencil bins = 1 for self, 7 for downwind bins
  // each action tells what proc and what bin the master's particles go to
  //   in order to satisfy stencil bin's upwind dependency
  // list of 8 possible actions is culled to insure each recv proc
  //   (including self) gets correct # of copies of particle (one per image)
  // culling rules:
  //   a) if 2 or more actions have same recv proc and ibin, only keep 1
  //   b) no 1st action if master bin is owned since it's a self-copy
  //   test a) then b) so other self-copies get eliminated

  for (int m = 0; m < nbins; m++) {
    local2local(m,&i,&j,&k);
    blist[m].nmigrate = 0;
    blist[m].migrate = NULL;
    if (blist[m].ghost == -1) continue;
    if (blist[m].ghost == 0 && i > 1 && j > 1 && k > 1) continue;

    // 8 stencil bins

    migrate_possible(0,i,j,k,i,j,k);
    migrate_possible(1,i,j,k,i-1,j,k);
    migrate_possible(2,i,j,k,i,j-1,k);
    migrate_possible(3,i,j,k,i,j,k-1);
    migrate_possible(4,i,j,k,i-1,j-1,k);
    migrate_possible(5,i,j,k,i,j-1,k-1);
    migrate_possible(6,i,j,k,i-1,j,k-1);
    migrate_possible(7,i,j,k,i-1,j-1,k-1);

    // cull

    for (n = 0; n < 8; n++) {
      if (possible[n].ptr == NULL) continue;
      for (n1 = n+1; n1 < 8; n1++) {
	if (possible[n1].ptr == NULL) continue;
	if (possible[n].proc != possible[n1].proc) continue;
	if (possible[n].ibin != possible[n1].ibin) continue;
	possible[n1].ptr = NULL;
      }
    }

    if (blist[m].ghost == 0) possible[0].ptr = NULL;

    // store remaining possible actions in bin's migrate list
    // loop backwards so action for self stencil comes first

    for (n = 7; n >= 0; n--) {
      if (possible[n].ptr == NULL) continue;
      ptr = migrate_request();
      ptr->ptr = blist[m].migrate;
      blist[m].migrate = ptr;
      blist[m].nmigrate++;
      ptr->proc = possible[n].proc;
      ptr->ibin = possible[n].ibin;
      ptr->rcheck = possible[n].rcheck;
      ptr->ilo = possible[n].ilo;
      ptr->jlo = possible[n].jlo;
      ptr->klo = possible[n].klo;
      ptr->pbc = possible[n].pbc;
      if (ptr->pbc) {
	ptr->xpbc = possible[n].xpbc;
	ptr->ypbc = possible[n].ypbc;
	ptr->zpbc = possible[n].zpbc;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   use xyz of particle to determine which local bin the particle is in
   assume xyz is in local bin (owned or ghost)
   return local bin ID (0 to nbins-1)
------------------------------------------------------------------------- */

int Grid::whichlocal(double *x)
{
  // ix,iy,iz range from -1 to gbinxyz inclusive

  int ix,iy,iz;

  if (x[0] < xorigin) ix = -1;
  else ix = static_cast<int> ((x[0]-xorigin)*xbininv);
  if (x[1] < yorigin) iy = -1;
  else iy = static_cast<int> ((x[1]-yorigin)*ybininv);
  if (x[2] < zorigin) iz = -1;
  else iz = static_cast<int> ((x[2]-zorigin)*zbininv);

  // offset converts ix,iy,iz to 0 to nbinxyz-1 inclusive
  // ilocal ranges from 0 to nbins-1 inclusive

  int ilocal = (iz-zoffset)*nbinx*nbiny + (iy-yoffset)*nbinx + ix-xoffset;
  return ilocal;
}

/* ----------------------------------------------------------------------
   use xyz of particle to determine which global bin the particle is in
   assume xyz is inside global domain
   set ix,iy,iz to global bin indices (0 to gbinxyz-1)
   return global bin ID (0 to gbins-1)
------------------------------------------------------------------------- */

int Grid::whichglobal(double *x, int *ix, int *iy, int *iz)
{
  // ix,iy,iz range from 0 to gbinxyz-1 inclusive

  *ix = static_cast<int> ((x[0]-xorigin)*xbininv);
  *iy = static_cast<int> ((x[1]-yorigin)*ybininv);
  *iz = static_cast<int> ((x[2]-zorigin)*zbininv);

  // iglobal ranges from 0 to gbins-1 inclusive

  int iglobal = (*iz)*gbinx*gbiny + (*iy)*gbinx + (*ix);
  return iglobal;
}

/* ----------------------------------------------------------------------
   convert global ID to local bin ID
   assume global ID is valid (0 to gbins-1)
   error if global bin does not map to local domain
   return local bin ID (0 to nbins-1)
------------------------------------------------------------------------- */

int Grid::global2local(int id)
{
  // ix,iy,iz range from 0 to nbinxyz-1 inclusive

  int ix = id % gbinx - xoffset;
  int iy = (id/gbinx) % gbiny - yoffset;
  int iz = id / (gbinx*gbiny) - zoffset;

  if (ix < 0 || ix >= nbinx || iy < 0 || iy >= nbiny ||
      iz < 0 || iz >= nbinz) 
    error->one("Global bin does not map to local domain");

  // ilocal ranges from 0 to nbins-1 inclusive

  int ilocal = iz*nbinx*nbiny + iy*nbinx + ix;
  return ilocal;
}

/* ----------------------------------------------------------------------
   convert global bin indices to local bin ID
   assume global indices range from -1 to gbinxyz inclusive
   if global indices do not map to local domain, return -1
   return local bin ID (0 to nbins-1)
------------------------------------------------------------------------- */

int Grid::global2local(int ix, int iy, int iz)
{
  // ixlocal,iylocal,izlocal range from 0 to nbinxyz-1 inclusive

  int ixlocal = ix - xoffset;
  int iylocal = iy - yoffset;
  int izlocal = iz - zoffset;

  if (ixlocal < 0 || ixlocal >= nbinx) return -1;
  if (iylocal < 0 || iylocal >= nbiny) return -1;
  if (izlocal < 0 || izlocal >= nbinz) return -1;

  // ilocal ranges from 0 to nbins-1 inclusive

  int ilocal = izlocal*nbinx*nbiny + iylocal*nbinx + ixlocal;
  return ilocal;
}

/* ----------------------------------------------------------------------
   convert local bin ID to global ID 
   assume local ID is valid (0 to nbins-1)
   if local bin is global ghost, global bin may be periodic or may not exist
   return -1 if it doesn't exist (non-periodic ghost)
   else return global ID (0 to gbins-1) with any needed PBC remap
------------------------------------------------------------------------- */

int Grid::local2global(int ibin)
{
  // ix,iy,iz range from -1 to gbinxyz inclusive

  int ix = ibin % nbinx + xoffset;
  int iy = (ibin/nbinx) % nbiny + yoffset;
  int iz = ibin / (nbinx*nbiny) + zoffset;

  // correct for PBC

  if (ix == -1) {
    if (domain->xperiodic) ix = gbinx - 1;
    else return -1;
  }
  if (ix == gbinx) {
    if (domain->xperiodic) ix = 0;
    else return -1;
  }
  if (iy == -1) {
    if (domain->yperiodic) iy = gbiny - 1;
    else return -1;
  }
  if (iy == gbiny) {
    if (domain->yperiodic) iy = 0;
    else return -1;
  }
  if (iz == -1) {
    if (domain->zperiodic) iz = gbinz - 1;
    else return -1;
  }
  if (iz == gbinz) {
    if (domain->zperiodic) iz = 0;
    else return -1;
  }

  // iglobal ranges from 0 to gbins-1 inclusive

  int iglobal = iz*gbinx*gbiny + iy*gbinx + ix;
  return iglobal;
}

/* ----------------------------------------------------------------------
   convert local bin ID to global bin indices
   assume local ID is valid (0 to nbins-1)
   set ix,iy,iz to global bin indices (-1 to gbinxyz)
------------------------------------------------------------------------- */

void Grid::local2global(int ibin, int *ix, int *iy, int *iz)
{
  // ix,iy,iz range from -1 to gbinxyz inclusive

  *ix = ibin % nbinx + xoffset;
  *iy = (ibin/nbinx) % nbiny + yoffset;
  *iz = ibin / (nbinx*nbiny) + zoffset;
}

/* ----------------------------------------------------------------------
   convert local bin ID to local bin indices
   assume local ID is valid (0 to nbins-1)
   set ix,iy,iz to local bin indices (0 to nbinxyz-1)
------------------------------------------------------------------------- */

void Grid::local2local(int id, int *ix, int *iy, int *iz)
{
  // ix,iy,iz range from 0 to nbinxyz-1 inclusive

  *ix = id % nbinx;
  *iy = (id/nbinx) % nbiny;
  *iz = id / (nbinx*nbiny);
}

/* ----------------------------------------------------------------------
   convert local bin indices to local bin ID
   assume local indices range from 0 to nbinxyz-1 inclusive
   return local bin ID (0 to nbins-1)
------------------------------------------------------------------------- */

int Grid::local2one(int ix, int iy, int iz)
{
  int ilocal = iz*nbinx*nbiny + iy*nbinx + ix;
  return ilocal;
}

/* ----------------------------------------------------------------------
   compute center point of global bin ID
   assume global ID is valid (0 to gbins-1)
   set x,y,z to center point displaced from xyz offset
------------------------------------------------------------------------- */

void Grid::coord_center(int id, double *x, double *y, double *z)
{
  // ix,iy,iz range from 0 to gbinxyz-1 inclusive

  int ix = id % gbinx;
  int iy = (id/gbinx) % gbiny;
  int iz = id / (gbinx*gbiny);

  // xyz = center of bin

  *x = ix*xbinsize + xbinhalf + xorigin;
  *y = iy*ybinsize + ybinhalf + yorigin;
  *z = iz*zbinsize + zbinhalf + zorigin;
}

/* ----------------------------------------------------------------------
   compute center point of global bin with indices ix,iy,iz
   assume global indices range from -1 to gbinxyz inclusive
   set x,y,z to center point displaced from xyz offset
------------------------------------------------------------------------- */

void Grid::coord_center(int ix, int iy, int iz,
			double *x, double *y, double *z)
{
  *x = ix*xbinsize + xbinhalf + xorigin;
  *y = iy*ybinsize + ybinhalf + yorigin;
  *z = iz*zbinsize + zbinhalf + zorigin;
}

/* ----------------------------------------------------------------------
   compute lower edge of global bin index N in dim 012 = xyz
   assume N ranges from -1 to gbinxyz inclusive
   return edge value
------------------------------------------------------------------------- */

double Grid::coord_cut(int dim, int n)
{
  if (dim == 0) return (n*xbinsize + xorigin);
  else if (dim == 1) return (n*ybinsize + yorigin);
  return (n*zbinsize + zorigin);
}

/* ----------------------------------------------------------------------
   compute 2 opposite corner points of global bin ID
   assume global ID is valid (0 to gbins-1)
------------------------------------------------------------------------- */

void Grid::coord_corners(int id, double *x1, double *y1, double *z1,
			 double *x2, double *y2, double *z2)
{
  // ix,iy,iz range from 0 to gbinxyz-1 inclusive

  int ix = id % gbinx;
  int iy = (id/gbinx) % gbiny;
  int iz = id / (gbinx*gbiny);

  // lower-left and upper-right corners of bin

  *x1 = ix*xbinsize + xorigin;
  *y1 = iy*ybinsize + yorigin;
  *z1 = iz*zbinsize + zorigin;
  *x2 = (ix+1)*xbinsize + xorigin;
  *y2 = (iy+1)*ybinsize + yorigin;
  *z2 = (iz+1)*zbinsize + zorigin;
}

/* ----------------------------------------------------------------------
   compute 2 opposite corner points of local bin with indices ix,iy,iz
   assume local indices range from 0 to nbinxyz-1 inclusive
------------------------------------------------------------------------- */

void Grid::coord_corners(int ix, int iy, int iz, double *lo, double *hi)
{
  // lower-left and upper-right corners of bin

  lo[0] = (ix+xoffset)*xbinsize + xorigin;
  lo[1] = (iy+yoffset)*ybinsize + yorigin;
  lo[2] = (iz+zoffset)*zbinsize + zorigin;
  hi[0] = (ix+xoffset+1)*xbinsize + xorigin;
  hi[1] = (iy+yoffset+1)*ybinsize + yorigin;
  hi[2] = (iz+zoffset+1)*zbinsize + zorigin;
}

/* ----------------------------------------------------------------------
   construct a 3d layout of procs and map to global bins
   set all nbin & n2g values for this proc's sub-domain
   set x2proc,y2proc,z2proc for 3d decomp
------------------------------------------------------------------------- */

void Grid::decomp3d()
{
  // map nprocs to global bins
  // mapping is done as 3d grid of procs
  // lo/hi = extent (inclusive) of owned bins by this proc

  int ilo,jlo,klo,ihi,jhi,khi;
  proc3d(me,nprocs,gbinx,gbiny,gbinz,&ilo,&jlo,&klo,&ihi,&jhi,&khi);
  
  // create my 3d array of bins, including ghosts
  // nbins = count of bins on this proc

  nbinx = ihi - ilo + 3;
  nbiny = jhi - jlo + 3;
  nbinz = khi - klo + 3;
  nbins = nbinx*nbiny*nbinz;

  // offset from local indexing to global indexing in each dim
  // global = local + offset
  // local bin #1 is global lo, so offset = lo-1

  xoffset = ilo - 1;
  yoffset = jlo - 1;
  zoffset = klo - 1;

  // allocate and fill x2proc,y2proc,z2proc arrays with 1d mapping
  // decomp1d gives bounds for each proc within actual cells
  // set end cells (0,gbin-1) separately

  x2proc = new int[gbinx];
  y2proc = new int[gbiny];
  z2proc = new int[gbinz];

  int i,j,lo,hi;

  for (i = 0; i < nprocx; i++) {
    decomp1d(gbinx,nprocx,i,&lo,&hi);
    for (j = lo; j <= hi; j++) x2proc[j] = i;
  }
  x2proc[0] = 0;
  x2proc[gbinx-1] = nprocx-1;

  for (i = 0; i < nprocy; i++) {
    decomp1d(gbiny,nprocy,i,&lo,&hi);
    for (j = lo; j <= hi; j++) y2proc[j] = i;
  }
  y2proc[0] = 0;
  y2proc[gbiny-1] = nprocy-1;

  for (i = 0; i < nprocz; i++) {
    decomp1d(gbinz,nprocz,i,&lo,&hi);
    for (j = lo; j <= hi; j++) z2proc[j] = i;
  }
  z2proc[0] = 0;
  z2proc[gbinz-1] = nprocz-1;

  // create copy of all proc's nbin and offset info

  bininfo[me][0] = nbinx;
  bininfo[me][1] = nbiny;
  bininfo[me][2] = nbinz;
  bininfo[me][3] = xoffset;
  bininfo[me][4] = yoffset;
  bininfo[me][5] = zoffset;

  MPI_Allgather(bininfo[me],6,MPI_INT,bininfo[0],6,MPI_INT,world);
}

/* ----------------------------------------------------------------------
   assign nprocs to Nx by Ny by Nz domain as a 3d grid of procs
   minimze surf area of sub-domain
   return lo/hi = ijk extent of each proc's sub-domain (0 to N-1)
------------------------------------------------------------------------- */

void Grid::proc3d(int me, int nprocs, int nx, int ny, int nz,
		 int *ilo, int *jlo, int *klo, int *ihi, int *jhi, int *khi)
{
  // max possible surface area of a sub-domain

  int bestsurf = 2 * (nx*ny + ny*nz + nx*nz);

  // loop thru all possible factorizations of nprocs
  // surf = surface area of a proc sub-domain
  // when done: nprocx,nprocy,nprocz = # of procs in each dimension

  int px,py,pz,nremain,onex,oney,onez,surf;

  px = 1;
  while (px <= nprocs) {
    if (nprocs % px == 0) {
      nremain = nprocs/px;
      py = 1;
      while (py <= nremain) {
	if (nremain % py == 0) {
	  pz = nremain/py;
	  onex = nx/px;
	  oney = ny/py;
	  onez = nz/pz;
	  surf = onex*oney + oney*onez + onex*onez;
	  if (surf < bestsurf) {
	    bestsurf = surf;
	    nprocx = px;
	    nprocy = py;
	    nprocz = pz;
	  }
	}
	py++;
      }
    }
    px++;
  }

  // ipx,ipy,ipz = which proc I am in each dim
  
  int ipx = me % nprocx;
  int ipy = (me/nprocx) % nprocy;
  int ipz = me / (nprocx*nprocy);

  // lo/hi = bins I own at lo/hi end of my sub-domain
  // lo/hi values = 0 to N-1

  decomp1d(nx,nprocx,ipx,ilo,ihi);
  decomp1d(ny,nprocy,ipy,jlo,jhi);
  decomp1d(nz,nprocz,ipz,klo,khi);
}

/* ----------------------------------------------------------------------
   compute a 1d decomp for N values onto NP procs
   return lo-hi range owned by IPth proc
   returned values range from 0 to N-1
------------------------------------------------------------------------- */

void Grid::decomp1d(int n, int np, int ip, int *lo, int *hi)
{
  *lo = static_cast<int> (1.0*ip / np * n + 0.5);
  *hi = static_cast<int> (1.0*(ip+1) / np * n + 0.5) - 1;
}

/* ----------------------------------------------------------------------
   compute proc owner of bin with global indices i,j,k
   decomp_style MUST be current and balance->lb MUST be current
   global indices can be < 0 or > gbinxyz -1
   apply PBC if necessary
   return -1 if global bin doesn't exist (non-periodic ghost)
   else return proc owner of bin
------------------------------------------------------------------------- */

int Grid::owner(int i, int j, int k)
{
  // apply PBC
  // return -1 if global bin is non-periodic ghost

  if (i < 0) {
    if (domain->xperiodic) i += gbinx;
    else return -1;
  }
  if (i >= gbinx) {
    if (domain->xperiodic) i -= gbinx;
    else return -1;
  }
  if (j < 0) {
    if (domain->yperiodic) j += gbiny;
    else return -1;
  }
  if (j >= gbiny) {
    if (domain->yperiodic) j -= gbiny;
    else return -1;
  }
  if (k < 0) {
    if (domain->zperiodic) k += gbinz;
    else return -1;
  }
  if (k >= gbinz) {
    if (domain->zperiodic) k -= gbinz;
    else return -1;
  }

  // bin owner depends on decomposition

  int proc;
  if (decomp_style == STATIC_BRICK) {
    int ipx = x2proc[i];
    int ipy = y2proc[j];
    int ipz = z2proc[k];
    proc = ipz*nprocy*nprocx + ipy*nprocx + ipx;
  } else {
    int id = k*gbinx*gbiny + j*gbinx + i;
    double x[3];
    coord_center(id,&x[0],&x[1],&x[2]);
    Zoltan_LB_Point_Assign(balance->lb,x,&proc);
  }

  return proc;
}

/* ----------------------------------------------------------------------
   check if box defined by lo/hi has any overlap with my 3d array of local bins
   return 1 if yes, 0 if no
   if yes, set ilo,khi to local bin indices of overlap (0 to nbinxyz-1)
------------------------------------------------------------------------- */

int Grid::box_overlap(double *lo, double *hi,
		      int &ilo, int &jlo, int &klo,
		      int &ihi, int &jhi, int &khi)
{
  double xlo = xoffset*xbinsize + xorigin;
  if (xlo > hi[0]) return 0;
  double xhi = (nbinx+xoffset)*xbinsize + xorigin;
  if (xhi < lo[0]) return 0;

  double ylo = yoffset*ybinsize + yorigin;
  if (ylo > hi[1]) return 0;
  double yhi = (nbiny+yoffset)*ybinsize + yorigin;
  if (yhi < lo[1]) return 0;

  double zlo = zoffset*zbinsize + zorigin;
  if (zlo > hi[2]) return 0;
  double zhi = (nbinz+zoffset)*zbinsize + zorigin;
  if (zhi < lo[2]) return 0;

  ilo = static_cast<int> ((lo[0]-xorigin)*xbininv) - xoffset;
  ilo = MAX(ilo,0);
  ihi = static_cast<int> ((hi[0]-xorigin)*xbininv) - xoffset;
  ihi = MIN(ihi,nbinx-1);

  jlo = static_cast<int> ((lo[1]-yorigin)*ybininv) - yoffset;
  jlo = MAX(jlo,0);
  jhi = static_cast<int> ((hi[1]-yorigin)*ybininv) - yoffset;
  jhi = MIN(jhi,nbiny-1);

  klo = static_cast<int> ((lo[2]-zorigin)*zbininv) - zoffset;
  klo = MAX(klo,0);
  khi = static_cast<int> ((hi[2]-zorigin)*zbininv) - zoffset;
  khi = MIN(khi,nbinz-1);

  return 1;
}

/* ----------------------------------------------------------------------
   add geometry itri (triangle or region) to local ibin
------------------------------------------------------------------------- */

void Grid::add_tri(int ibin, int itri)
{
  Index *newptr = index_request();
  newptr->index = itri;
  newptr->ptr = blist[ibin].tri;
  blist[ibin].ntri++;
  blist[ibin].tri = newptr;
}

/* ----------------------------------------------------------------------
   delete linked list of geometry starting with ptr
   return Index objects to free list
------------------------------------------------------------------------- */

void Grid::delete_tri(Index *ptr)
{
  Index *newptr;
  while (ptr) {
    newptr = ptr->ptr;
    index_return(ptr);
    ptr = newptr;
  }
}

/* ----------------------------------------------------------------------
   delete linked list of migrations starting with ptr
   return Migrate objs to free list
------------------------------------------------------------------------- */

void Grid::delete_migrate(Migrate *ptr)
{
  Migrate *newptr;
  while (ptr) {
    newptr = ptr->ptr;
    migrate_return(ptr);
    ptr = newptr;
  }
}

/* ----------------------------------------------------------------------
   set dynamic load-balance parameters
------------------------------------------------------------------------- */

void Grid::set_balance(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal balance command");

  if (strcmp(arg[0],"static") == 0) {
    if (narg < 2 || narg > 3) error->all("Illegal balance command");
    if (strcmp(arg[1],"bin") == 0) balance_style = STATIC_BIN;
    else if (strcmp(arg[1],"particle") == 0) balance_style = STATIC_PART;
    else error->all("Illegal balance command");
    if (narg == 3) {
      balance_file = 1;
      if (me == 0) {
	fpbal = fopen(arg[2],"w");
	if (fpbal == NULL) error->one("Could not open balance file");
      }
    } else balance_file = 0;
  } else if (strcmp(arg[0],"dynamic") == 0) {
    if (narg < 3 || narg > 4) error->all("Illegal balance command");
    balance_style = DYNAMIC_PART;
    nbalance = atoi(arg[1]);
    threshhold = atof(arg[2]);
    if (narg == 4) {
      balance_file = 1;
      if (me == 0) {
	fpbal = fopen(arg[3],"w");
	if (fpbal == NULL) error->one("Could not open balance file");
      }
    } else balance_file = 0;
  } else error->all("Illegal balance command");
}

/* ----------------------------------------------------------------------
   perform dynamic load-balancing
------------------------------------------------------------------------- */

void Grid::rebalance()
{
  int i,m;
  Zoltan_Comm_Obj *plan;
  Particle::OnePart *plist = particle->plist;

  // setup linked list of particles in each grid bin and count them

  link();

  // perform RCB on owned bins with weight = bins or particles
  // RCB returns bounding box for my sub-domain
  // adjust edges of box that are outside global domain

  double xlo,ylo,zlo,xhi,yhi,zhi;
  balance->RCB(weightflag,&xlo,&ylo,&zlo,&xhi,&yhi,&zhi);
  
  xlo = MAX(xlo,domain->xlo);
  ylo = MAX(ylo,domain->ylo);
  zlo = MAX(zlo,domain->zlo);
  xhi = MIN(xhi,domain->xhi);
  yhi = MIN(yhi,domain->yhi);
  zhi = MIN(zhi,domain->zhi);

  if (balance_file) write_balance(xlo,ylo,zlo,xhi,yhi,zhi);

  // save old bin list and decomposition info so can retrieve bin info

  OneBin *bold = blist;

  int decomp_style_old = decomp_style;

  int nbinx_old = nbinx;
  int nbiny_old = nbiny;
  int nbins_old = nbins;

  int xoffset_old = xoffset;
  int yoffset_old = yoffset;
  int zoffset_old = zoffset;

  // new decomposition info for RCB domains
  // convert bounding box to integer bins
  // ilo-khi = lo/hi global bounds (0 to gbinx-1) of owned bins (no ghosts)

  int ilo = static_cast<int> ((xlo+xbinhalf - xorigin) * xbininv);
  int ihi = static_cast<int> ((xhi-xbinhalf - xorigin) * xbininv);
  int jlo = static_cast<int> ((ylo+ybinhalf - yorigin) * ybininv);
  int jhi = static_cast<int> ((yhi-ybinhalf - yorigin) * ybininv);
  int klo = static_cast<int> ((zlo+zbinhalf - zorigin) * zbininv);
  int khi = static_cast<int> ((zhi-zbinhalf - zorigin) * zbininv);

  nbinx = ihi-ilo + 3;
  nbiny = jhi-jlo + 3;
  nbinz = khi-klo + 3;
  nbins = nbinx*nbiny*nbinz;

  xoffset = ilo - 1;
  yoffset = jlo - 1;
  zoffset = klo - 1;

  // create copy of all proc's new nbin and offset info

  bininfo[me][0] = nbinx;
  bininfo[me][1] = nbiny;
  bininfo[me][2] = nbinz;
  bininfo[me][3] = xoffset;
  bininfo[me][4] = yoffset;
  bininfo[me][5] = zoffset;

  MPI_Allgather(bininfo[me],6,MPI_INT,bininfo[0],6,MPI_INT,world);

  // allocate new bins and compute topology connections for new decomposition
  // topology calls owner which requires decomp_style and lb be set

  blist = (OneBin *) memory->smalloc(nbins*sizeof(OneBin),"grid:blist");
  decomp_style = STATIC_BIN;
  balance->swap();
  topology();

  // loop over new own+ghost bins
  // skip if bin is global ghost and non-periodic
  // request SURF and PARTICLE info for my owned bins
  // request just SURF info for my ghost bins
  // fill_ba finds bin owners in old decomposition (3d or RCB lb)
  //   which requires decomp_style and lb be set to old decomposition

  decomp_style = decomp_style_old;
  balance->swap();

  int nsend = 0;
  for (m = 0; m < nbins; m++) {
    if (blist[m].id == -1) continue;
    if (blist[m].ghost == 0) fill_ba(blist[m].id,m,SURF_PARTICLE,&nsend);
    else fill_ba(blist[m].id,m,SURF,&nsend);
  }

  balance->swap();

  // communicate the BinAsk as unstructured communication
  // insure buf2 is large enough to hold incoming data

  int nrecv;
  Zoltan_Comm_Create(&plan,nsend,proclist1,world,0,&nrecv);

  if (nrecv*sizeof(BinAsk) > size2) {
    memory->sfree(buf2);
    size2 = nrecv*sizeof(BinAsk);
    buf2 = (char *) memory->smalloc(size2,"grid:buf2");
  }

  Zoltan_Comm_Do(plan,0,buf1,sizeof(BinAsk),buf2);
  Zoltan_Comm_Destroy(&plan);

  // loop over nrecv from BinAsk communcation
  // convert global bin ID -> old local bin via old decomposition settings
  // always create SurfReply
  // if BinAsk flag = SURF_PARTICLE, also create PartReply
  // pass info from old decomposition to fill_sr and fill_pr

  int iglobal,ix,iy,iz,ilocal;
  BinAsk *bufba = (BinAsk *) buf2;

  int nsends = 0;
  int nsendp = 0;
  for (i = 0; i < nrecv; i++) {
    iglobal = bufba[i].id;
    ix = iglobal % gbinx - xoffset_old;
    iy = (iglobal/gbinx) % gbiny - yoffset_old;
    iz = iglobal / (gbinx*gbiny) - zoffset_old;
    ilocal = iz*nbinx_old*nbiny_old + iy*nbinx_old + ix;
    fill_sr(bold[ilocal].tri,bufba[i].ibin,bufba[i].requester,&nsends);
    if (bufba[i].flag == SURF_PARTICLE)
      fill_pr(bold[ilocal].first,bufba[i].ibin,bufba[i].requester,&nsendp);
  }

  // can now free geometry and migration lists from old decomposition
  // can now free all particles

  for (m = 0; m < nbins_old; m++) {
    delete_tri(bold[m].tri);
    delete_migrate(bold[m].migrate);
  }
  particle->nlocal = particle->nghost = particle->ntotal = 0;

  // communicate the SurfReply as unstructured communication
  // insure buf2 is large enough to hold incoming data

  int nrecvs;
  Zoltan_Comm_Create(&plan,nsends,proclist1,world,0,&nrecvs);

  if (nrecvs*sizeof(SurfReply) > size2) {
    memory->sfree(buf2);
    size2 = nrecvs*sizeof(SurfReply);
    buf2 = (char *) memory->smalloc(size2,"grid:buf2");
  }

  Zoltan_Comm_Do(plan,0,buf1,sizeof(SurfReply),buf2);
  Zoltan_Comm_Destroy(&plan);

  // communicate the PartReply as unstructured communication
  // insure buf4 is large enough to hold incoming data

  int nrecvp;
  Zoltan_Comm_Create(&plan,nsendp,proclist2,world,0,&nrecvp);

  if (nrecvp*sizeof(PartReply) > size4) {
    memory->sfree(buf4);
    size4 = nrecvp*sizeof(PartReply);
    buf4 = (char *) memory->smalloc(size4,"grid:buf4");
  }

  Zoltan_Comm_Do(plan,0,buf3,sizeof(PartReply),buf4);
  Zoltan_Comm_Destroy(&plan);

  // loop over nrecvs from SurfReply communication
  // add each tri to ibin in new decomposition

  SurfReply *bufbr = (SurfReply *) buf2;

  for (i = 0; i < nrecvs; i++)
    add_tri(bufbr[i].ibin,bufbr[i].itri);

  // loop over nrecvp from PartReply communication
  // add each particle to plist

  PartReply *bufpr = (PartReply *) buf4;

  for (i = 0; i < nrecvp; i++) {
    m = particle->add(bufpr[i].species,
		      bufpr[i].x[0],bufpr[i].x[1],bufpr[i].x[2]);
    plist = particle->plist;
    particle->nlocal++;
    plist[m].ibin = bufpr[i].ibin;

    // debug statement
    if (plist[m].ibin >= nbins)
      printf("BAD BIN %d %d %d %d\n",me,m,plist[m].ibin,nbins);

    plist[m].itri = bufpr[i].itri;
    plist[m].seed = bufpr[i].seed;
  }

  // free old blist

  memory->sfree(bold);
}

/* ----------------------------------------------------------------------
   add global bin id to BinAsk buffer
   send to proc that owns it in old decomposition
------------------------------------------------------------------------- */

void Grid::fill_ba(int id, int ibin, int flag, int *pn)
{
  int n = *pn;

  // grow buffer and proclist1 by 2x if necessary

  if (n*sizeof(BinAsk) >= size1) {
    if (size1 == 0) size1 = 1000*sizeof(BinAsk);
    else size1 = 2*n*sizeof(BinAsk);
    buf1 = (char *) memory->srealloc(buf1,size1,"grid:buf1");
  }

  if (n >= sizeproc1) {
    if (sizeproc1 == 0) sizeproc1 = 1000;
    else sizeproc1 = 2*n;
    proclist1 = (int *)
      memory->srealloc(proclist1,sizeproc1*sizeof(int),"grid:proclist1");
  }

  BinAsk *buf = (BinAsk *) buf1;
  buf[n].requester = me;
  buf[n].id = id;
  buf[n].ibin = ibin;
  buf[n].flag = flag;

  // proclist1 = owner of global cell id in old decomposition

  int ix = id % gbinx;
  int iy = (id/gbinx) % gbiny;
  int iz = id / (gbinx*gbiny);
  proclist1[n] = owner(ix,iy,iz);
  *pn = n+1;
}

/* ----------------------------------------------------------------------
   add linked list of geometry to SurfReply buffer to send to obin of proc
------------------------------------------------------------------------- */

void Grid::fill_sr(Index *ptr, int obin, int proc, int *pn)
{
  int n = *pn;
  SurfReply *buf = (SurfReply *) buf1;

  // loop over tris in linked list

  while (ptr) {

    // grow buffer and proclist1 by 2x if necessary

    if (n*sizeof(SurfReply) >= size1) {
      if (size1 == 0) size1 = 1000*sizeof(SurfReply);
      else size1 = 2*n*sizeof(SurfReply);
      buf1 = (char *) memory->srealloc(buf1,size1,"grid:buf1");
      buf = (SurfReply *) buf1;
    }

    if (n >= sizeproc1) {
      if (sizeproc1 == 0) sizeproc1 = 1000;
      else sizeproc1 = 2*n;
      proclist1 = (int *)
	memory->srealloc(proclist1,sizeproc1*sizeof(int),"grid:proclist1");
    }

    buf[n].itri = ptr->index;
    buf[n].ibin = obin;

    proclist1[n] = proc;
    n++;

    // goto next tri in linked list

    ptr = ptr->ptr;
  }

  *pn = n;
}

/* ----------------------------------------------------------------------
   add linked list of particles to PartReply buffer to send to obin of proc
------------------------------------------------------------------------- */

void Grid::fill_pr(int m, int obin, int proc, int *pn)
{
  int n = *pn;
  PartReply *buf = (PartReply *) buf3;
  Particle::OnePart *plist = particle->plist;

  // loop over particles in linked list

  while (m >= 0) {

    // grow buffer and proclist2 by 2x if necessary

    if (n*sizeof(PartReply) >= size3) {
      if (size3 == 0) size3 = 1000*sizeof(PartReply);
      else size3 = 2*n*sizeof(PartReply);
      buf3 = (char *) memory->srealloc(buf3,size3,"grid:buf3");
      buf = (PartReply *) buf3;
    }
    
    if (n >= sizeproc2) {
      if (sizeproc2 == 0) sizeproc2 = 1000;
      else sizeproc2 = 2*n;
      proclist2 = (int *)
	memory->srealloc(proclist2,sizeproc2*sizeof(int),"grid:proclist2");
    }

    // store particle m in PartReply buffer[n]

    buf[n].x[0] = plist[m].x[0];
    buf[n].x[1] = plist[m].x[1];
    buf[n].x[2] = plist[m].x[2];
    buf[n].species = plist[m].species;
    buf[n].seed = plist[m].seed;
    buf[n].itri = plist[m].itri;
    buf[n].ibin = obin;

    proclist2[n] = proc;
    n++;

    // goto next particle in linked list

    m = plist[m].next;
  }

  *pn = n;
}

/* ----------------------------------------------------------------------
   setup particle-to-particle links in plist and counter in blist
   do it for owned and ghost particles
   assumes blist[].first = -1, blist[].nparts = 0 for all bins w/ particles
   set this way from initial topology()
------------------------------------------------------------------------- */

void Grid::link()
{
  int ibin;

  Particle::OnePart *plist = particle->plist;
  int ntotal = particle->ntotal;
  for (int i = 0; i < ntotal; i++) {
    ibin = plist[i].ibin;
    plist[i].next = blist[ibin].first;
    blist[ibin].first = i;
    blist[ibin].nparts++;
  }
}

/* ----------------------------------------------------------------------
   return size of OneBin, Index, Migrate, buffers
------------------------------------------------------------------------- */

int Grid::memory_usage()
{
  int n = nbins*sizeof(OneBin);
  n += nichunks*CHUNK*sizeof(Index);
  n += nmchunks*CHUNK*sizeof(Migrate);
  n += size1 + size2 + size3 + size4;
  return n;
}

/* ----------------------------------------------------------------------
   write load-balance boxes to file as region commands
------------------------------------------------------------------------- */

void Grid::write_balance(double xlo, double ylo, double zlo,
			 double xhi, double yhi, double zhi)
{
  // boxes = sub-domain boundaries from all procs

  double box[6];
  box[0] = xlo;
  box[1] = ylo;
  box[2] = zlo;
  box[3] = xhi;
  box[4] = yhi;
  box[5] = zhi;
  double *boxes = new double[6*nprocs];

  MPI_Gather(box,6,MPI_DOUBLE,boxes,6,MPI_DOUBLE,0,world);

  // write one region command per box with final newline

  if (me == 0) {
    for (int i = 0; i < nprocs; i++)
      fprintf(fpbal,"region %d box %g %g %g %g %g %g\n",
	      i,boxes[6*i+0],boxes[6*i+1],boxes[6*i+2],
	      boxes[6*i+3],boxes[6*i+4],boxes[6*i+5]);
    fprintf(fpbal,"\n");
  }

  // free local memory

  delete [] boxes;
}

/* ----------------------------------------------------------------------
   return one index object from index freelist
------------------------------------------------------------------------- */

Grid::Index *Grid::index_request()
{
  if (!ifreelist) {
    Index *chunk = (Index *)
      memory->smalloc(CHUNK*sizeof(Index),"grid:chunk");
    for (int i = 0; i < CHUNK; i++) chunk[i].ptr = &chunk[i+1];
    chunk[CHUNK-1].ptr = NULL;
    ifreelist = chunk;
    ichunks = (Index **)
      memory->srealloc(ichunks,(nichunks+1)*sizeof(Index *),"grid:ichunks");
    ichunks[nichunks++] = chunk;
  }
  Index *ptr = ifreelist;
  ifreelist = ifreelist->ptr;
  return ptr;
}

/* ----------------------------------------------------------------------
   add one returned index object to index freelist
------------------------------------------------------------------------- */

void Grid::index_return(Index *ptr)
{
  ptr->ptr = ifreelist;
  ifreelist = ptr;
}

/* ----------------------------------------------------------------------
   return one migrate object from migrate freelist
------------------------------------------------------------------------- */

Grid::Migrate *Grid::migrate_request()
{
  if (!mfreelist) {
    Migrate *chunk = (Migrate *) 
      memory->smalloc(CHUNK*sizeof(Migrate),"grid:chunk");
    for (int i = 0; i < CHUNK; i++) chunk[i].ptr = &chunk[i+1];
    chunk[CHUNK-1].ptr = NULL;
    mfreelist = chunk;
    mchunks = (Migrate **)
      memory->srealloc(mchunks,(nmchunks+1)*sizeof(Migrate *),"grid:mchunks");
    mchunks[nmchunks++] = chunk;
  }
  Migrate *ptr = mfreelist;
  mfreelist = mfreelist->ptr;
  return ptr;
}

/* ----------------------------------------------------------------------
   add one returned migrate object to migrate freelist
------------------------------------------------------------------------- */

void Grid::migrate_return(Migrate *ptr)
{
  ptr->ptr = mfreelist;
  mfreelist = ptr;
}

/* ----------------------------------------------------------------------
   compute migrate action for one downwind stencil bin of a local master bin
   iposs = index into possible array for storing action parameters
   im,jm,km = local indices of master bin (0 to nbinxyz-1 inclusive)
   is,js,ks = local indices of stencil bin (-1 to nbinxyz-1 inclusive)
------------------------------------------------------------------------- */

void Grid::migrate_possible(int iposs, int im, int jm, int km,
			    int is, int js, int ks)
{
  // ijk_mglobal = global indices of master bin from -1 to gbinxyz inclusive
  // ijk_sglobal = global indices of stencil bin from -2 to gbinxyz inclusive

  int imglobal = im + xoffset;
  int jmglobal = jm + yoffset;
  int kmglobal = km + zoffset;

  int isglobal = is + xoffset;
  int jsglobal = js + yoffset;
  int ksglobal = ks + zoffset;

  // compute PBC remapping of stencil cell to global interior
  // if stencil cell is inside global domain, no remapping needed
  // if outside and periodic, store the PBC remap for stencil and master bin
  // if outside and non-periodic, no migrate action so return with ptr = NULL

  possible[iposs].ptr = NULL;
  int xpbc,ypbc,zpbc;
  xpbc = ypbc = zpbc = 0;

  if (isglobal < 0) {
    if (!domain->xperiodic) return;
    isglobal += gbinx;
    imglobal += gbinx;
    xpbc = 1;
  }
  if (isglobal >= gbinx) {
    if (!domain->xperiodic) return;
    isglobal -= gbinx;
    imglobal -= gbinx;
    xpbc = -1;
  }
  if (jsglobal < 0) {
    if (!domain->yperiodic) return;
    jsglobal += gbiny;
    jmglobal += gbiny;
    ypbc = 1;
  }
  if (jsglobal >= gbiny) {
    if (!domain->yperiodic) return;
    jsglobal -= gbiny;
    jmglobal -= gbiny;
    ypbc = -1;
  }
  if (ksglobal < 0) {
    if (!domain->zperiodic) return;
    ksglobal += gbinz;
    kmglobal += gbinz;
    zpbc = 1;
  }
  if (ksglobal >= gbinz) {
    if (!domain->zperiodic) return;
    ksglobal -= gbinz;
    kmglobal -= gbinz;
    zpbc = -1;
  }

  // proc = owner of global (remapped) stencil bin
  // if original stencil bin is a owned or ghost cell of this proc
  //   (is,js,ks all >= 0), then proc owner is stored in blist[].proc
  // else call owner() with globally remapped indices

  int proc;
  if (is >= 0 && js >= 0 && ks >= 0) {
    int ibin = ks*nbinx*nbiny + js*nbinx + is;
    proc = blist[ibin].proc;
  } else proc = owner(isglobal,jsglobal,ksglobal);

  // set ibin to local ID of remapped master bin on proc owning stencil bin
  // use ijk_mglobal so master bin has been PBC remapped same as stencil bin
  // bininfo stores nbin and offset values for the all procs

  int ix = imglobal - bininfo[proc][3];
  int iy = jmglobal - bininfo[proc][4];
  int iz = kmglobal - bininfo[proc][5];
  int ibin = iz*bininfo[proc][0]*bininfo[proc][1] + iy*bininfo[proc][0] + ix;

  // if proc is me, check global ID of ibin
  // if global ghost and non-periodic, then return with ptr = NULL
  // this prevents copying a ghost particle that should be deleted

  if (proc == me && blist[ibin].id == -1) return;

  // store values in iposs element of possible vector
  // set ptr to non-NULL to indicate valid migrate action

  possible[iposs].proc = proc;
  possible[iposs].ibin = ibin;
  if (nperx > 1) possible[iposs].ilo = im-is;
  else possible[iposs].ilo = 0;
  if (npery > 1) possible[iposs].jlo = jm-js;
  else possible[iposs].jlo = 0;
  if (nperz > 1) possible[iposs].klo = km-ks;
  else possible[iposs].klo = 0;
  if (possible[iposs].ilo || possible[iposs].jlo || possible[iposs].klo)
    possible[iposs].rcheck = 1;
  else possible[iposs].rcheck = 0;

  if (xpbc || ypbc || zpbc) {
    possible[iposs].pbc = 1;
    possible[iposs].xpbc = xpbc;
    possible[iposs].ypbc = ypbc;
    possible[iposs].zpbc = zpbc;
  } else possible[iposs].pbc = 0;

  possible[iposs].ptr = possible;
}
