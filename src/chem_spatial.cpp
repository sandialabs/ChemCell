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
#include "chem_spatial.h"
#include "simulator.h"
#include "particle.h"
#include "react.h"
#include "domain.h"
#include "grid.h"
#include "random.h"
#include "surf.h"
#include "balance.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

#define REACTANT 1
#define PRODUCT  2

#define AVOGADRO 6.023e23
#define PI 3.1415926
#define BIG 1.0e20

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

ChemSpatial::ChemSpatial() : Chem()
{
  char *str = "spatial";
  int n = strlen(str) + 1;
  style = new char[n];
  strcpy(style,str);

  MPI_Comm_rank(world,&me);

  allocated = 0;

  ncolor = 8;
  firstcolor = new int[ncolor];

  nstencil = 14;
  stencil1 = new int[nstencil];
  stencil2 = new int[nstencil];

  size1 = size2 = size3 = 0;
  buf1 = buf2 = buf3 = NULL;
  sizeproc = 0;
  proclist = NULL;

  monoprob = NULL;
  monoprobsum = NULL;
  pairprob = NULL;
  pairprobsum = NULL;
  distsq = NULL;
}

/* ---------------------------------------------------------------------- */

ChemSpatial::~ChemSpatial()
{
  if (allocated) free_arrays();

  delete [] firstcolor;
  delete [] stencil1;
  delete [] stencil2;

  memory->sfree(buf1);
  memory->sfree(buf2);
  memory->sfree(buf3);
  memory->sfree(proclist);
}

/* ---------------------------------------------------------------------- */

void ChemSpatial::init()
{
  int ireact,ispecies,jspecies,imono,ipair;

  // memory allocation

  if (allocated) free_arrays();
  allocated = 1;

  int nreactions = react->nreactions;
  double *rate = react->rate;

  nreactant = react->nreactant;
  reactants = react->reactants;
  nproduct = react->nproduct;
  products = react->products;
  nmono = react->nmono;
  monolist = react->monolist;
  npair = react->npair;
  pairlist = react->pairlist;

  // error if any zero-reactant reactions

  for (ireact = 0; ireact < nreactions; ireact++)
    if (nreactant[ireact] == 0) {
      char *str = "Cannot have 0-reactant reactions with spatial simulation";
      error->all(str);
    }

  // memory allocation

  int nspecies = particle->nspecies;

  int max = 0;
  for (ispecies = 0; ispecies < nspecies; ispecies++)
    max = MAX(max,nmono[ispecies]);

  monoprob = memory->create_2d_double_array(nspecies,max,"spatial:monoprob");
  monoprobsum = (double *)
    memory->smalloc(nspecies*sizeof(double),"spatial:monoprobsum");

  max = 0;
  for (ispecies = 0; ispecies < nspecies; ispecies++)
    for (jspecies = 0; jspecies < nspecies; jspecies++)
      max = MAX(max,npair[ispecies][jspecies]);

  pairprob = 
    memory->create_3d_double_array(nspecies,nspecies,max,"spatial:pairprob");
  pairprobsum = 
    memory->create_2d_double_array(nspecies,nspecies,"spatial:pairprobsum");
  distsq = memory->create_2d_double_array(nspecies,nspecies,"spatial:distsq");

  // set probability factor for mono reactions
  // prob = reaction rate * timestep

  for (ispecies = 0; ispecies < nspecies; ispecies++)
    for (imono = 0; imono < nmono[ispecies]; imono++) {
      ireact = monolist[ispecies][imono];
      monoprob[ispecies][imono] = rate[ireact] * simulator->dt;
    }

  // set probability factor and cutoff distance for dual reactions
  // 1.0E15 factor is conversion from liters to um^3

  // if prob_style = max style (0):
  // kmax = fastest dual reaction rate
  // max_prob = user-set probability for kmax
  // distance = separation at which kmax reaction occurs with max_prob
  // distance is same for all reactions
  // prob is set for each reaction to give correct k

  if (prob_style == 0) {
    int iwhich = -1;
    double kmax = 0.0;
    for (ireact = 0; ireact < nreactions; ireact++)
      if (nreactant[ireact] == 2 && rate[ireact] > kmax) {
	kmax = rate[ireact];
	iwhich = ireact;
      }

    double distance = 0.0;
    if (iwhich >= 0) {
      double vol_overlap = kmax * simulator->dt / (AVOGADRO * max_prob);
      vol_overlap *= 1.0E15;
      double radius = 0.5 * pow(0.75*vol_overlap/PI,1.0/3.0);
      distance = 2.0*radius;
    }

    for (ispecies = 0; ispecies < nspecies; ispecies++) {
      for (jspecies = 0; jspecies < nspecies; jspecies++) {
	for (ipair = 0; ipair < npair[ispecies][jspecies]; ipair++) {
	  ireact = pairlist[ispecies][jspecies][ipair];

	  distsq[ispecies][jspecies] = distance*distance;

	  double prob;
	  double vol_overlap = 4.0/3.0 * PI * distance*distance*distance;
	  if (distance == 0.0) prob = 0.0;
	  else prob = rate[ireact] * simulator->dt / 
		 (AVOGADRO * vol_overlap/1.0e15);
	  pairprob[ispecies][jspecies][ipair] = prob;
	}
      }
    }
  }

  // if prob_style = diff style (1):
  // distance = diff_factor times sum of 2 ave diff distances for reactants
  // prob is set for each reaction to give correct k
  // if setdist or setprob was user-specified, they override these values

  if (prob_style == 1) {
    double *setdist = react->setdist;
    double *setprob = react->setprob;

    for (ispecies = 0; ispecies < nspecies; ispecies++) {
      for (jspecies = 0; jspecies < nspecies; jspecies++) {
	for (ipair = 0; ipair < npair[ispecies][jspecies]; ipair++) {
	  ireact = pairlist[ispecies][jspecies][ipair];

	  double distance;
	  if (setdist[ireact] == -1.0) {
	    double d1ave = 4.0/sqrt(PI) * 
	      sqrt(1.0e8*particle->diffusivity[ispecies] * simulator->dt);
	    double d2ave = 4.0/sqrt(PI) * 
	      sqrt(1.0e8*particle->diffusivity[jspecies] * simulator->dt);
	    distance = diff_factor * (d1ave + d2ave);
	  } else distance = setdist[ireact];
	  distsq[ispecies][jspecies] = distance*distance;

	  double prob;
	  double vol_overlap = 4.0/3.0 * PI * distance*distance*distance;
	  if (distance == 0.0) prob = 0.0;
	  else if (setprob[ireact] == -1.0)
	    prob = rate[ireact] * simulator->dt / 
	      (AVOGADRO * vol_overlap/1.0E15);
	  else prob = setprob[ireact];
	  pairprob[ispecies][jspecies][ipair] = prob;
	}
      }
    }
  }

  // monoprobsum & pairprobsum = sum of individual reaction probabilities
  // set monoprob & pairprob = cummulative probablity for all reactions
  //   involving same reactants
  // check that no individual or summed probability > 1.0

  for (ispecies = 0; ispecies < nspecies; ispecies++) {
    monoprobsum[ispecies] = 0.0;
    for (imono = 0; imono < nmono[ispecies]; imono++) {
      monoprobsum[ispecies] += monoprob[ispecies][imono];
      if (imono) monoprob[ispecies][imono] += monoprob[ispecies][imono-1];
      if (monoprob[ispecies][imono] > 1.0) {
	char str[128];
	sprintf(str,"Probability %g for reaction %s is too large",
		monoprob[ispecies][imono],
		react->name[monolist[ispecies][imono]]);
	error->all(str);
      }
      if (monoprobsum[ispecies] > 1.0) {
	char str[128];
	sprintf(str,"Summed probability %g for species %s is too large",
		monoprobsum[ispecies],particle->name[ispecies]);
	error->all(str);
      }
    }
  }

  for (ispecies = 0; ispecies < nspecies; ispecies++) {
    for (jspecies = 0; jspecies < nspecies; jspecies++) {
      pairprobsum[ispecies][jspecies] = 0.0;
      for (ipair = 0; ipair < npair[ispecies][jspecies]; ipair++) {
	pairprobsum[ispecies][jspecies] +=
	  pairprob[ispecies][jspecies][ipair];
	if (ipair) pairprob[ispecies][jspecies][ipair] += 
		     pairprob[ispecies][jspecies][ipair-1];
	if (pairprob[ispecies][jspecies][ipair] > 1.0) {
	  char str[128];
	  sprintf(str,"Probability %g for reaction %s is too large",
		  pairprob[ispecies][jspecies][ipair],
		  react->name[pairlist[ispecies][jspecies][ipair]]);
	  error->all(str);
	}
	if (pairprobsum[ispecies][jspecies] > 1.0) {
	  char str[128];
	  sprintf(str,"Summed probability %g for species %s & %s is too large",
		  pairprobsum[ispecies][jspecies],
		  particle->name[ispecies],particle->name[jspecies]);
	  error->all(str);
	}
      }
    }
  }

  // print reaction stats

  double dmin = BIG;
  double pmin = BIG;
  double dmax = -1.0;
  double pmax = -1.0;

  for (ispecies = 0; ispecies < nspecies; ispecies++) {
    for (imono = 0; imono < nmono[ispecies]; imono++) {
      pmin = MIN(pmin,monoprob[ispecies][imono]);
      pmax = MAX(pmax,monoprob[ispecies][imono]);
    }
  }

  for (ispecies = 0; ispecies < nspecies; ispecies++) {
    for (jspecies = 0; jspecies < nspecies; jspecies++) {
      for (ipair = 0; ipair < npair[ispecies][jspecies]; ipair++) {
	dmin = MIN(dmin,sqrt(distsq[ispecies][jspecies]));
	dmax = MAX(dmax,sqrt(distsq[ispecies][jspecies]));
	pmin = MIN(pmin,pairprob[ispecies][jspecies][ipair]);
	pmax = MAX(pmax,pairprob[ispecies][jspecies][ipair]);
      }
    }
  }

  if (me == 0) {
    fprintf(screen,"Chemistry:\n");
    fprintf(screen,"  min/max reaction distance = %g %g\n",dmin,dmax);
    fprintf(screen,"  min/max reaction probabilities = %g %g\n",pmin,pmax);

    fprintf(logfile,"Chemistry:\n");
    fprintf(logfile,"  min/max reaction distance = %g %g\n",dmin,dmax);
    fprintf(logfile,"  min/max reaction probabilities = %g %g\n",pmin,pmax);
  }

  // error if max reaction distance exceeds bin size

  if (dmax > grid->xbinsize || dmax > grid->ybinsize || 
      dmax > grid->zbinsize) error->all("Max reaction distance > bin size");

  // setup reaction stencils for initial decomp
  // setup_stencil also performs bin coloring for decomp

  setup_stencil();

  // initialize reaction counters

  rcount = new int[nreactions];
  for (int ireact = 0; ireact < nreactions; ireact++) rcount[ireact] = 0;
  nbinbin = nbinpair = ndist = noverlap = 0;
}

/* ----------------------------------------------------------------------
   perform all particle reactions
   2-reactant and 1-reactant cases
   do in subsets of colored bins to avoid parallel conflicts
   communicate changed particles to other procs after each color
------------------------------------------------------------------------- */

void ChemSpatial::reactions()
{
  int i,icolor,ip,jp,kp,mp,istencil,ireact,nsend,nrecv,ncopy;
  int ispecies,jspecies,kspecies,ibin,jbin,kbin,ibinnew;
  int npossible,ipossible;
  double delx,dely,delz,rsq;
  double rn;
  Zoltan_Comm_Obj *plan;

  simulator->ctime += simulator->dt;

  // if no reactions defined, erase ghost particles and return

  if (react->nreactions == 0) {
    particle->nghost = 0;
    particle->ntotal = particle->nlocal;
    return;
  }

  // setup linked list of particles in each bin and sort them within each bin

  particle->link();
  particle->sort();

  // flag all particles as not reacted: own = 0, ghost = 1

  Particle::OnePart *plist = particle->plist;
  int nlocal = particle->nlocal;
  int ntotal = particle->ntotal;

  for (i = 0; i < nlocal; i++) plist[i].flag = 0;
  for (i = nlocal; i < ntotal; i++) plist[i].flag = 1;

  // loop over each color

  Grid::OneBin *blist = grid->blist;
  Grid::Migrate *ptr;

  for (icolor = 0; icolor < ncolor; icolor++) {

    // loop over bins of a particular color
   
    nsend = ncopy = 0;
    kbin = firstcolor[icolor];

    while (kbin >= 0) {

      // consider all dual reactions with this stencil origin = kbin
      // loop over stencil pairs ibin,jbin = 2 interacting bins

      nbinbin += nstencil;

      for (istencil = 0; istencil < nstencil; istencil++) {
	ibin = kbin + stencil1[istencil];
	jbin = kbin + stencil2[istencil];

	// if either bin is global non-periodic ghost, go to next stencil pair
	// if either bin has no particles, go to next stencil pair

	if (blist[ibin].id == -1 || blist[jbin].id == -1) continue;
	if (blist[ibin].first == -1 || blist[jbin].first == -1) continue;
	if (ibin != jbin) nbinpair += blist[ibin].nparts * blist[jbin].nparts;
	else nbinpair += (blist[ibin].nparts * (blist[jbin].nparts-1)) / 2;

	// ip,jp = 2 particles
	// double loop over all ip in ibin, all jp in jbin
	// if ibin = jbin, then start jp one particle beyond ip
	// skip any ip or jp that has already reacted (flag = -1)

	ip = blist[ibin].first;
	while (ip >= 0) {

	  if (plist[ip].flag == -1) goto idone;
	  ispecies = plist[ip].species;

	  if (ibin == jbin) jp = plist[ip].next;
	  else jp = blist[jbin].first;
	  while (jp >= 0) {

	    if (plist[jp].flag == -1) goto jdone;
	    jspecies = plist[jp].species;

	    // check if any pair reactions between ispecies & jspecies

	    npossible = npair[ispecies][jspecies];
	    if (npossible == 0) goto jdone;

	    // rsq = distance between 2 particles
	    // check if within reaction cutoff for the pair

	    delx = plist[ip].x[0] - plist[jp].x[0];
	    dely = plist[ip].x[1] - plist[jp].x[1];
	    delz = plist[ip].x[2] - plist[jp].x[2];
	    rsq = delx*delx + dely*dely + delz*delz;
	    ndist++;

	    if (rsq >= distsq[ispecies][jspecies]) goto jdone;
	    noverlap++;

	    // roll MC dice to see if a reaction takes place

	    rn = random->react(plist[ip].seed,plist[jp].seed);
	    if (rn >= pairprobsum[ispecies][jspecies]) goto jdone;

	    // ireact = which reaction is selected
	    // multiple reactions possible for the species pair

	    if (npossible == 1) ireact = pairlist[ispecies][jspecies][0];
	    else
	      for (ipossible = 0; ipossible < npossible; ipossible++) {
		if (rn < pairprob[ispecies][jspecies][ipossible]) {
		  ireact = pairlist[ispecies][jspecies][ipossible];
		  break;
		}
	      }
	    rcount[ireact]++;

	    // flag the 2 reactants as non-reactive

	    plist[ip].flag = -1;
	    plist[jp].flag = -1;

	    // notify procs and image bins that reactants are now non-reactive
	    // this is all procs/bins that store ibin/jbin as owned or
	    //   upwind ghost bin including PBC images on same proc
	    // bin migrate action list stores this info
	    // one exception: if ibin/jbin is ghost it has a migrate action
	    //   to same ibin/jbin on this proc, which is not needed

	    ptr = blist[ibin].migrate;
	    while (ptr) {
	      if (ptr->proc == me && ptr->ibin == ibin) {
		ptr = ptr->ptr;
		continue;
	      }
	      if (ptr->proc == me) fill_rc(&plist[ip],ptr,REACTANT,&ncopy);
	      else fill_rm(&plist[ip],ptr,REACTANT,&nsend);
	      ptr = ptr->ptr;
	    }

	    ptr = blist[jbin].migrate;
	    while (ptr) {
	      if (ptr->proc == me && ptr->ibin == jbin) {
		ptr = ptr->ptr;
		continue;
	      }
	      if (ptr->proc == me) fill_rc(&plist[jp],ptr,REACTANT,&ncopy);
	      else fill_rm(&plist[jp],ptr,REACTANT,&nsend);
	      ptr = ptr->ptr;
	    }

	    // create products as ghost particles
	    // if product is in ghost bin
	    //   set flag = 1, so will delete it
	    //   send to owning proc as stored in bin's 1st migrate action
	    //   no need to tell downwind procs, b/c products are non-reactive
	    // if product is in owned bin
	    //   set flag = 0, so will save it when compacting particles

	    for (int m = 0; m < nproduct[ireact]; m++) {
	      mp = react->create_product(m,ireact,ip,jp);
	      particle->nghost++;
	      plist = particle->plist;
	      ibinnew = plist[mp].ibin;
	      if (blist[ibinnew].ghost) {
		plist[mp].flag = 1;
		if (blist[ibinnew].migrate->proc == me)
		  fill_rc(&plist[mp],blist[ibinnew].migrate,PRODUCT,&ncopy);
		else
		  fill_rm(&plist[mp],blist[ibinnew].migrate,PRODUCT,&nsend);
	      } else plist[mp].flag = 0;
	    }

	    // done with particle ip
	    // go to next particle in ibin

	    goto idone;

	  jdone:
	    jp = plist[jp].next;
	  }
	idone:
	  ip = plist[ip].next;
	}
      }

      // consider mono reactions for particles in this bin

      kp = blist[kbin].first;
      while (kp >= 0) {

	// skip a kp that has already reacted

	if (plist[kp].flag) goto kdone;
	kspecies = plist[kp].species;

	// check if any mono reactions for kspecies

	npossible = nmono[kspecies];
	if (npossible == 0) goto kdone;

	// roll MC dice to see if a reaction takes place
	
	rn = random->react(plist[kp].seed,0);
	if (rn >= monoprobsum[kspecies]) goto kdone;

	// ireact = which reaction is selected
	// multiple reactions possible for the species pair
	
	if (npossible == 1) ireact = monolist[kspecies][0];
	else
	  for (ipossible = 0; ipossible < npossible; ipossible++) {
	    if (rn < monoprob[kspecies][ipossible]) {
	      ireact = monolist[kspecies][ipossible];
	      break;
	    }
	  }
	rcount[ireact]++;

	// flag the reactant as non-reactive

	plist[kp].flag = -1;
	    
	// notify procs and image bins that reactants are now non-reactive
	// this is all procs/bins that store owned kbin as
	//   upwind ghost bin including PBC images on same proc
	// bin migrate action list stores this info

	ptr = blist[kbin].migrate;
	while (ptr) {
	  if (ptr->proc == me) fill_rc(&plist[kp],ptr,REACTANT,&ncopy);
	  else fill_rm(&plist[kp],ptr,REACTANT,&nsend);
	  ptr = ptr->ptr;
	}

	// create reaction products as ghost particles
	// indicate mono reaction via last arg to create_product = -1
	// if product is in ghost bin (due to being placed epsilon away from surf)
	//   set flag = 1, so will delete it
	//   send to owning proc as stored in bin's 1st migrate action
	//   no need to tell downwind procs, b/c products are non-reactive
	// if product is in owned bin
	//   set flag = 0, so will save it when compacting particles

	for (int m = 0; m < nproduct[ireact]; m++) {
	  mp = react->create_product(m,ireact,kp,-1);
	  particle->nghost++;
	  plist = particle->plist;
	  ibinnew = plist[mp].ibin;
	  if (blist[ibinnew].ghost) {
	    plist[mp].flag = 1;
	    if (blist[ibinnew].migrate->proc == me)
	      fill_rc(&plist[mp],blist[ibinnew].migrate,PRODUCT,&ncopy);
	    else
	      fill_rm(&plist[mp],blist[ibinnew].migrate,PRODUCT,&nsend);
	  } else plist[mp].flag = 0;
	}
	
      kdone:
	kp = plist[kp].next;
      }

      // next bin of same color

      kbin = blist[kbin].next;
    }

    // send and receive REACTANT/PRODUCT particles before next color
    // create Zoltan plan, perform unstructured communication, destroy plan
    // insure buf2 is large enough to hold incoming data

    timer->sub_stamp();

    Zoltan_Comm_Create(&plan,nsend,proclist,world,0,&nrecv);

    if (nrecv > size2) {
      memory->sfree(buf2);
      size2 = nrecv;
      buf2 = (Migrate *) memory->smalloc(size2*sizeof(Migrate),"spatial:buf2");
    }

    Zoltan_Comm_Do(plan,0,(char *) buf1,sizeof(Migrate),(char *) buf2);
    Zoltan_Comm_Destroy(&plan);

    timer->sub_stamp(TIME_REACT_COMM);

    // process copied and received REACTANT and PRODUCT particles

    unpack(ncopy,buf3);
    unpack(nrecv,buf2);
    plist = particle->plist;
  }

  // clear bins before compacting particles
  // compact the particle list: no deleted or ghost particles

  particle->unlink(0);
  particle->compact();
}

/* ----------------------------------------------------------------------
   process received or copied particles
   for REACTANT particles, find them in ibin and flag as non-reactive
   for PRODUCT particles, I am new owner, so add them to plist as ghosts
   flag them as 0 so will be saved when compacting particles
------------------------------------------------------------------------- */

void ChemSpatial::unpack(int n, Migrate *buf)
{
  int mbin,mp;
  Grid::OneBin *blist = grid->blist;
  Particle::OnePart *plist = particle->plist;

  for (int i = 0; i < n; i++) {
    if (buf[i].flag == REACTANT) {
      mbin = buf[i].ibin;
      mp = blist[mbin].first;
      while (mp >= 0) {
	if (match(&buf[i],mp)) break;
	mp = plist[mp].next;
      }
      if (mp == -1) error->one("Did not find matching reaction particle");
      plist[mp].flag = -1;
    } else {
      mp = particle->add(buf[i].species,
			 buf[i].x[0],buf[i].x[1],buf[i].x[2]);
      plist = particle->plist;
      particle->nghost++;
      plist[mp].seed = buf[i].seed;
      plist[mp].ibin = buf[i].ibin;
      plist[mp].itri = buf[i].itri;
      plist[mp].flag = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   add Ith particle from plist to Migrate send buffer
   use ibin since it is local bin ID in receiving proc
------------------------------------------------------------------------- */

void ChemSpatial::fill_rm(Particle::OnePart *p, Grid::Migrate *ptr,
			  int flag, int *pn)
{
  int n = *pn;

  // grow buffer and proclist by 2x if necessary

  if (n >= size1) {
    if (size1 == 0) size1 = 1000;
    else size1 = 2*n;
    buf1 = (Migrate *)
      memory->srealloc(buf1,size1*sizeof(Migrate),"spatial:buf1");
  }

  if (n >= sizeproc) {
    if (sizeproc == 0) sizeproc = 1000;
    else sizeproc = 2*n;
    proclist = (int *)
      memory->srealloc(proclist,sizeproc*sizeof(int),"spatial:proclist");
  }

  // proc to send to

  proclist[n] = ptr->proc;

  // remap particle across PBC if necessary

  if (ptr->pbc) {
    buf1[n].x[0] = p->x[0] + ptr->xpbc * domain->xsize;
    buf1[n].x[1] = p->x[1] + ptr->ypbc * domain->ysize;
    buf1[n].x[2] = p->x[2] + ptr->zpbc * domain->zsize;
  } else {
    buf1[n].x[0] = p->x[0];
    buf1[n].x[1] = p->x[1];
    buf1[n].x[2] = p->x[2];
  }

  // rest of particle attributes

  buf1[n].species = p->species;
  buf1[n].seed = p->seed;
  buf1[n].ibin = ptr->ibin;
  buf1[n].itri = p->itri;
  buf1[n].flag = flag;

  *pn = n+1;
}

/* ----------------------------------------------------------------------
   add Ith particle from plist to Migrate on-processor copy buffer
   use ibin since it is local bin ID in self
------------------------------------------------------------------------- */

void ChemSpatial::fill_rc(Particle::OnePart *p, Grid::Migrate *ptr,
			  int flag, int *pn)
{
  int n = *pn;

  // grow buffer by 2x if necessary

  if (n >= size3) {
    if (size3 == 0) size3 = 1000;
    else size3 = 2*n;
    buf3 = (Migrate *)
      memory->srealloc(buf3,size3*sizeof(Migrate),"spatial:buf3");
  }

  // remap particle across PBC if necessary

  if (ptr->pbc) {
    buf3[n].x[0] = p->x[0] + ptr->xpbc * domain->xsize;
    buf3[n].x[1] = p->x[1] + ptr->ypbc * domain->ysize;
    buf3[n].x[2] = p->x[2] + ptr->zpbc * domain->zsize;
  } else {
    buf3[n].x[0] = p->x[0];
    buf3[n].x[1] = p->x[1];
    buf3[n].x[2] = p->x[2];
  }

  // rest of particle attributes
  // convert itri to global IDs

  buf3[n].species = p->species;
  buf3[n].seed = p->seed;
  buf3[n].ibin = ptr->ibin;
  buf3[n].itri = p->itri;
  buf3[n].flag = flag;

  *pn = n+1;
}

/* ----------------------------------------------------------------------
   match Migrate buffer particle to plist particle i
   match is based on seed, species, xyz equality (if not PBC)
   return 1 if a match, 0 if not
------------------------------------------------------------------------- */

int ChemSpatial::match(Migrate *buf, int i)
{
  Particle::OnePart *p = &particle->plist[i];

  if (buf->seed != p->seed) return 0;
  if (buf->species != p->species) return 0;
  if (!domain->xperiodic && buf->x[0] != p->x[0]) return 0;
  if (!domain->yperiodic && buf->x[1] != p->x[1]) return 0;
  if (!domain->zperiodic && buf->x[2] != p->x[2]) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   setup reaction stencils based on current decomposition
------------------------------------------------------------------------- */

void ChemSpatial::setup_stencil()
{
  int nbinx = grid->nbinx;
  int nbiny = grid->nbiny;

  stencil1[ 0] = 0;          stencil2[ 0] = 0;
  stencil1[ 1] = 0;          stencil2[ 1] = 1;
  stencil1[ 2] = 0;          stencil2[ 2] = nbinx;
  stencil1[ 3] = 0;          stencil2[ 3] = nbinx + 1;
  stencil1[ 4] = 1;          stencil2[ 4] = nbinx;
  stencil1[ 5] = 0;          stencil2[ 5] = nbinx*nbiny;
  stencil1[ 6] = 0;          stencil2[ 6] = nbinx*nbiny + 1;
  stencil1[ 7] = 0;          stencil2[ 7] = nbinx*nbiny + nbinx;
  stencil1[ 8] = 0;          stencil2[ 8] = nbinx*nbiny + nbinx + 1;
  stencil1[ 9] = 1;          stencil2[ 9] = nbinx*nbiny;
  stencil1[10] = 1;          stencil2[10] = nbinx*nbiny + nbinx;
  stencil1[11] = nbinx;      stencil2[11] = nbinx*nbiny;
  stencil1[12] = nbinx;      stencil2[12] = nbinx*nbiny + 1;
  stencil1[13] = nbinx + 1;  stencil2[13] = nbinx*nbiny;

  setup_colors();
}

/* ----------------------------------------------------------------------
   setup bin-to-bin ptrs for each reaction color
   only for owned bins
   loop in reverse order so that 1st bins are 1st in color list
------------------------------------------------------------------------- */

void ChemSpatial::setup_colors() 
{
  int i,icolor;

  for (i = 0; i < ncolor; i++) firstcolor[i] = -1;

  Grid::OneBin *blist = grid->blist;
  for (i = grid->nbins-1; i >= 0; i--) {
    if (!blist[i].ghost) {
      icolor = whichcolor(blist[i].id);
      blist[i].next = firstcolor[icolor];
      firstcolor[icolor] = i;
    }
  }
}

/* ----------------------------------------------------------------------
   compute which of 8 colors the global bin ID is
   return color # (0-7)
------------------------------------------------------------------------- */

int ChemSpatial::whichcolor(int id) 
{
  // ix,iy,iz = global bin indices

  int ix = id % grid->gbinx;
  int iy = (id/grid->gbinx) % grid->gbiny;
  int iz = id / (grid->gbinx * grid->gbiny);

  // cx,cy,cz = color indices
  // add 1 so that lowest-left global bin is color 0

  int cx = (ix+1) % 2;
  int cy = (iy+1) % 2;
  int cz = (iz+1) % 2;

  int icolor = cz*4 + cy*2 + cx;
  return icolor;
}

/* ----------------------------------------------------------------------
   free arrays
------------------------------------------------------------------------- */

void ChemSpatial::free_arrays() 
{
  memory->destroy_2d_double_array(monoprob);
  memory->sfree(monoprobsum);
  memory->destroy_3d_double_array(pairprob);
  memory->destroy_2d_double_array(pairprobsum);
  memory->destroy_2d_double_array(distsq);
  if (rcount) delete [] rcount;
}
