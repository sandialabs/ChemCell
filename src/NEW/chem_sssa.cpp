/* ----------------------------------------------------------------------
   ChemCell - Cell simulator for particle diffusion and reactions
   Steve Plimpton (sjplimp@sandia.gov), Alex Slepoy (aslepoy@sandia.gov)
   Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/chemcell.html

   Copyright (2004) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level ChemCell directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "chem_sssa.h"
#include "simulator.h"
#include "move.h"
#include "particle.h"
#include "grid.h"
#include "domain.h"
#include "random.h"
#include "react.h"
#include "output.h"
#include "surf.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#define AVOGADRO 6.023e23
#define PI 3.1415926
#define BIG 1.0e20

#define CHUNK 16

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

ChemSpatialSSA::ChemSpatialSSA() : Chem()
{
  char *str = "spatial/ssa";
  int n = strlen(str) + 1;
  style = new char[n];
  strcpy(style,str);

  allocated = 0;

  nstencil1 = 14;
  stencil1 = new int[nstencil1];
  nstencil2 = 27;
  stencil2 = new int[nstencil2];

  distsq = NULL;
  rate = NULL;

  maxtotal = 0;
  elist = NULL;
  maxelist = NULL;
  nelist = NULL;

  rlist = NULL;
}

/* ---------------------------------------------------------------------- */

ChemSpatialSSA::~ChemSpatialSSA()
{
  if (allocated) free_arrays();

  delete [] stencil1;
  delete [] stencil2;
}

/* ----------------------------------------------------------------------
   init a simulation
------------------------------------------------------------------------- */

void ChemSpatialSSA::init()
{
  int ispecies,jspecies;

  // memory allocation

  if (allocated) free_arrays();
  allocated = 1;

  int nspecies = particle->nspecies;
  nreactions = react->nreactions;
  nreactant = react->nreactant;
  reactants = react->reactants;
  nproduct = react->nproduct;
  products = react->products;
  nmono = react->nmono;
  monolist = react->monolist;
  npair = react->npair;
  pairlist = react->pairlist;

  // n = value such that 2^n >= nreactions

  int n = 0;
  while ((1 << n) < nreactions) n++;

  // create tree of length ntotal
  // init all leaves to 0.0

  int ntotal = 2*(1 << n) - 1;
  tree = new double[ntotal];
  offset = (1 << n) - 1;
  for (int i = offset; i < ntotal; i++) tree[i] = 0.0;

  // initialize rlist

  resize_rlist();

  // init reaction rates for SSA spatial algorithm

  distsq = memory->create_2d_double_array(nspecies,nspecies,"update:distsq");
  for (int ispecies = 0; ispecies < nspecies; ispecies++)
    for (int jspecies = 0; jspecies < nspecies; jspecies++)
      distsq[ispecies][jspecies] = 0.0;

  rate = new double[nreactions];

  if (prob_style == 0)
    error->all("Must use probability diff with run style spatial/ssa");

  if (prob_style == 1) {
    double *setdist = react->setdist;

    for (int ireact = 0; ireact < nreactions; ireact++) {
      if (nreactant[ireact] == 0) {
	char *str = "Cannot have 0-reactant reactions with "
	  "run style spatial/ssa";
	error->all(str);
      }
      if (nreactant[ireact] == 1) rate[ireact] = react->rate[ireact];
      if (nreactant[ireact] == 2) {
	ispecies = reactants[ireact][0];
	jspecies = reactants[ireact][1];

	double distance;
	if (setdist[ireact] == -1.0) {
	  double d1ave = 4.0/sqrt(PI) * 
	    sqrt(1.0e8*particle->diffusivity[ispecies] * simulator->dt);
	  double d2ave = 4.0/sqrt(PI) * 
	    sqrt(1.0e8*particle->diffusivity[jspecies] * simulator->dt);
	  distance = diff_factor * (d1ave + d2ave);
	} else distance = setdist[ireact];
	distsq[ispecies][jspecies] = distance*distance;
	distsq[jspecies][ispecies] = distance*distance;

	double prob;
	double vol_overlap = 4.0/3.0 * PI * distance*distance*distance;
	if (distance == 0.0) prob = 0.0;
	else prob = react->rate[ireact] / (AVOGADRO * vol_overlap/1.0E15);
	rate[ireact] = prob;
      }
    }
  }

  // print reaction stats

  double dmin = BIG;
  double dmax = -1.0;

  for (ispecies = 0; ispecies < nspecies; ispecies++) {
    for (jspecies = 0; jspecies < nspecies; jspecies++) {
      for (int ipair = 0; ipair < npair[ispecies][jspecies]; ipair++) {
	dmin = MIN(dmin,sqrt(distsq[ispecies][jspecies]));
	dmax = MAX(dmax,sqrt(distsq[ispecies][jspecies]));
      }
    }
  }

  int me;
  MPI_Comm_rank(world,&me);

  if (me == 0) {
    fprintf(screen,"Chemistry:\n");
    fprintf(screen,"  min/max reaction distance = %g %g\n",dmin,dmax);

    fprintf(logfile,"Chemistry:\n");
    fprintf(logfile,"  min/max reaction distance = %g %g\n",dmin,dmax);
  }

  // error if max reaction distance exceeds bin size

  if (dmax > grid->xbinsize || dmax > grid->ybinsize || 
      dmax > grid->zbinsize) error->all("Max reaction distance > bin size");

  // error if max reaction distance >= 1/2 domain size in a periodic dimension
  // because A-B interaction will be put in edge list for B and its ghost

  int iflag = 0;
  if (domain->xperiodic && dmax >= 0.5*domain->xsize) iflag = 1;
  if (domain->yperiodic && dmax >= 0.5*domain->ysize) iflag = 1;
  if (domain->zperiodic && dmax >= 0.5*domain->zsize) iflag = 1;
  if (iflag) error->all("Max reaction distance > 1/2 periodic domain");

  // setup stencils for neighbor finding
 
  setup_stencil();

  // initialize reaction counters

  rcount = new int[nreactions];
  for (int ireact = 0; ireact < nreactions; ireact++) rcount[ireact] = 0;
  nbinbin = nbinpair = ndist = noverlap = 0;
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::reactions()
{
  int i,m,n,ireact,nr,np,mp,part1,part2,ibin;
  double r1,r2,r3,delta;

  // if no reactions defined, erase ghost particles and return

  if (react->nreactions == 0) {
    particle->nghost = 0;
    particle->ntotal = particle->nlocal;
    return;
  }

  // acquire all ghosts across PBC, not just upwind ones provided by migrate
  // sets ghost atom seed to index of owned particle it is image of
  //   so can always use owned particles in elist and rlist

  timer->sub_stamp();

  particle->nghost = 0;
  particle->ntotal = particle->nlocal;
  particle->ghost_acquire();

  timer->sub_stamp(TIME_REACT_COMM);

  // flag all particles as not reacted: own = 0, ghost = 1

  Grid::OneBin *blist = grid->blist;
  Particle::OnePart *plist = particle->plist;
  int nlocal = particle->nlocal;
  int ntotal = particle->ntotal;

  for (i = 0; i < nlocal; i++) plist[i].flag = 0;
  for (i = nlocal; i < ntotal; i++) plist[i].flag = 1;

  // setup linked list of particles in each bin including ghost bins
  // clear all reaction lists
  // resize edge list to current # of owned particles
  // build edge list for all particle pairs and associated reaction list

  particle->link();
  for (m = 0; m < nreactions; m++) nrlist[m] = 0;
  resize_elist();
  build_elist_all();

  // compute initial propensity for each reaction using rlist counts
  // sum the tree initially

  propensity = &tree[offset];
  for (m = 0; m < nreactions; m++) propensity[m] = compute_propensity(m);
  sum();

  // loop until reach diffusive time

  double target = simulator->ntimestep * simulator->dt;
  double ctime = simulator->ctime;

  // int loopcount = 0;

  while (ctime < target) {

    // loopcount++;
    // consistency();

    if (tree[0] > 0.0) {
      r1 = random->gillespie();
      delta = 1.0/tree[0] * log(1.0/r1);
    } else break;
    ctime += delta;

    r2 = random->gillespie();
    ireact = find(r2*tree[0]);
    r3 = random->gillespie();
    n = static_cast<int> (r3 * nrlist[ireact]);
    rcount[ireact]++;

    // delete reactants of chosen reaction, their edges, and edge's reactions
    // must set part2 before deleting edges since reactions are rearranged

    nr = nreactant[ireact];

    if (nr) {
      part1 = rlist[ireact][n].part1;
      if (nr == 2) part2 = rlist[ireact][n].part2;
      delete_edges(part1);
      plist[part1].flag = -1;
      if (nr == 2) {
	delete_edges(part2);
	plist[part2].flag = -1;
      }
    }

    // create reaction products as owned particles
    // set flag = 0 so will participate in subsequent reactions
    // add new product to front of its bin list
    // create edges for new product and their associated reactions

    np = nproduct[ireact];
    for (m = 0; m < np; m++) {
      if (nr == 2) mp = react->create_product(m,ireact,part1,part2);
      else mp = react->create_product(m,ireact,part1,-1);
      particle->nlocal++;
      plist = particle->plist;
      
      plist[mp].flag = 0;
      ibin = plist[mp].ibin;
      plist[mp].next = blist[ibin].first;
      blist[ibin].first = mp;
      blist[ibin].nparts++;

      resize_elist();
      build_elist_one(mp);
    }

    // recompute all propensities
    
    for (int m = 0; m < nreactions; m++)
      propensity[m] = compute_propensity(m);
    sum();
  }

  // consistency();

  // set simulator time to ctime
  // ctime is > target which will be carried over to next set of reactions

  simulator->ctime = ctime;

  // clear edge list before compacting particles since nlocal will change
  // clear bins with current particle list
  // compact the particle list: no deleted or ghost particles

  for (i = 0; i < nreactions; i++)
    if (nreactant[i] == 2) noverlap += nrlist[i];
  
  ntotal = particle->ntotal;
  for (i = 0; i < ntotal; i++) nelist[i] = 0;

  particle->unlink(1);
  particle->compact();
}

/* ----------------------------------------------------------------------
   check self-consistency of elist,rlist data structs
------------------------------------------------------------------------- */

void ChemSpatialSSA::consistency()
{
  int i,j,part1,edge1,part2,edge2;
  int type,ireact,active,type1,ireact1,active1,type2,ireact2,active2;
  Particle::OnePart *plist = particle->plist;

  printf("Consistency checking ...\n");

  // examine rlist

  for (i = 0; i < nreactions; i++) {
    if (nrlist[i] < 0 || nrlist[i] > maxrlist[i]) {
      printf("BAD: %d %d %d\n",i,nrlist[i],maxrlist[i]);
      error->one("Invalid nrlist length");
    }
    for (j = 0; j < nrlist[i]; j++) {
      part1 = rlist[i][j].part1;
      edge1 = rlist[i][j].edge1;
      if (nreactant[i] == 2) {
	part2 = rlist[i][j].part2;
	edge2 = rlist[i][j].edge2;
      } else part2 = edge2 = 0;
      if (part1 < 0 || part1 > maxtotal ||
	  part2 < 0 || part2 > maxtotal) {
	printf("BAD: %d %d %d\n",part1,part2,maxtotal);
	error->one("Invalid rlist particle params");
      }
      if (edge1 < 0 || edge1 > nelist[part1] || 
	  edge2 < 0 || edge2 > nelist[part2]) {
	printf("BAD: %d %d: %d %d\n",edge1,nelist[part1],edge2,nelist[part2]);
	error->one("Invalid rlist edge params");
      }
      if (plist[part1].flag || (nreactant[i] == 2 && plist[part2].flag)) {
	printf("BAD: %d %d: %d %d\n",part1,part2,
	       plist[part1].flag,plist[part2].flag);
	error->one("Invalid rlist particle existence");
      }
      type1 = elist[part1][edge1].type;
      ireact1 = elist[part1][edge1].ireact;
      active1 = elist[part1][edge1].active;
      if (nreactant[i] == 2) {
	type2 = elist[part2][edge2].type;
	ireact2 = elist[part2][edge2].ireact;
	active2 = elist[part2][edge2].active;
      } else {
	type2 = i;
	ireact2 = j;
	active2 = 1;
      }
      if (type1 != i || type2 != i || ireact1 != j || ireact2 != j ||
	  active1 != 1 || active2 != 1) {
	printf("BAD: %d %d %d: %d %d %d: %d %d\n",
	       type1,type2,i,ireact1,ireact2,j,active1,active2);
	error->one("Invalid rlist->elist reaction params");
      }
    }
  }

  // examine elist

  if (particle->ntotal > maxtotal) {
    error->one("Invalid elist length");
  }

  for (i = 0; i < particle->ntotal; i++) {
    if (nelist[i] < 0 || nelist[i] > maxelist[i]) {
      printf("BAD: %d %d %d\n",i,nelist[i],maxelist[i]);
      error->one("Invalid nelist length");
    }
    for (j = 0; j < nelist[i]; j++) {
      type = elist[i][j].type;
      ireact = elist[i][j].ireact;
      active = elist[i][j].active;
      if (!active) continue;
      if (type < 0 || type > nreactions) {
	printf("BAD: %d %d\n",type,nreactions);
	error->one("Invalid elist type params");
      }
      if (ireact < 0 || ireact > nrlist[type]) {
	printf("BAD: %d %d\n",ireact,nrlist[type]);
	error->one("Invalid elist ireact params");
      }
      part1 = rlist[type][ireact].part1;
      edge1 = rlist[type][ireact].edge1;
      if (nreactant[type] == 2) {
	part2 = rlist[type][ireact].part2;
	edge2 = rlist[type][ireact].edge2;
      } else part2 = edge2 = -1;
      if (part1 != i && part2 != i) {
	printf("BAD: %d %d %d %d: %d %d\n",i,j,type,ireact,part1,part2);
	error->one("Invalid elist->rlist particle params");
      }
      if (edge1 != j && edge2 != j) {
	printf("BAD: %d %d %d %d: %d %d\n",i,j,type,ireact,edge1,edge2);
	error->one("Invalid elist->rlist edge params");
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::build_elist_all()
{
  int j,k,m,istencil,npossible;
  int ibin,jbin,ispecies,jspecies;
  double delx,dely,delz,rsq;

  // loop over half neighbors of each particle via half stencil1
  // this will examine each pair of particles once

  Grid::OneBin *blist = grid->blist;
  Particle::OnePart *plist = particle->plist;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++) {
    ibin = plist[i].ibin;
    ispecies = plist[i].species;

    // add all mono reactions this particle is involved in to rlist
    // npossible = # of mono reactions this species is involved in

    npossible = nmono[ispecies];
    for (k = 0; k < npossible; k++) {
      m = monolist[ispecies][k];
      add_reaction(m,i);
    }
    
    // add all dual reactions this particle is involved in to rlist
    // loop over neighbor bins in half stencil
    // loop over particles in each bin
    // 1st stencil is 0 so loop is over remaining particles in ibin
    // if stencil jbin is a non-periodic ghost, skip it
    // npossible = # of dual reactions these species are involved in

    for (istencil = 0; istencil < nstencil1; istencil++) {
      jbin = ibin + stencil1[istencil];
      if (blist[jbin].id == -1) continue;
      if (ibin != jbin) nbinpair += blist[ibin].nparts * blist[jbin].nparts;
      else nbinpair += (blist[ibin].nparts * (blist[jbin].nparts-1)) / 2;
      if (istencil == 0) j = plist[i].next;
      else j = blist[jbin].first;
      while (j >= 0) {
	jspecies = plist[j].species;

	npossible = npair[ispecies][jspecies];
	if (npossible == 0) {
	  j = plist[j].next;
	  continue;
	}
	
	// no need to apply PBC because ghost particles are already remapped

	delx = plist[i].x[0] - plist[j].x[0];
	dely = plist[i].x[1] - plist[j].x[1];
	delz = plist[i].x[2] - plist[j].x[2];
	rsq = delx*delx + dely*dely + delz*delz;
	ndist++;

	if (rsq >= distsq[ispecies][jspecies]) {
	  j = plist[j].next;
	  continue;
	}

	// if j particle is ghost, add reaction to original image particle

	for (k = 0; k < npossible; k++) {
	  m = pairlist[ispecies][jspecies][k];
	  if (blist[jbin].ghost) add_reaction(m,i,plist[j].seed);
	  else add_reaction(m,i,j);
	}
	
	j = plist[j].next;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::build_elist_one(int i)
{
  int j,k,m,istencil,npossible,flag;
  int ibin,jbin,ispecies,jspecies;
  double delx,dely,delz,rsq;

  // loop over all neighbors of particle I via full stencil2
  // this will examine all neighbors of I once

  Grid::OneBin *blist = grid->blist;
  Particle::OnePart *plist = particle->plist;
  
  ibin = plist[i].ibin;
  ispecies = plist[i].species;
    
  // add all mono reactions this particle is involved in to rlist
  // npossible = # of mono reactions this species is involved in

  npossible = nmono[ispecies];
  for (k = 0; k < npossible; k++) {
    m = monolist[ispecies][k];
    add_reaction(m,i);
  }
    
  // add all dual reactions this particle is involved in to rlist
  // loop over neighbor bins in full stencil
  // loop over particles in each bin
  // if stencil jbin is a non-periodic ghost, skip it
  
  for (istencil = 0; istencil < nstencil2; istencil++) {
    jbin = ibin + stencil2[istencil];
    if (blist[jbin].id == -1) continue;
    j = blist[jbin].first;
    while (j >= 0) {

      // skip if j is already reacted particle
      // if j is ghost, check original image particle's flag

      if (blist[jbin].ghost) flag = plist[plist[j].seed].flag;
      else flag = plist[j].flag;
      if (flag == -1) {
	j = plist[j].next;
	continue;
      }

      // skip if i = j

      if (i == j) {
	j = plist[j].next;
	continue;
      }

      jspecies = plist[j].species;

      // npossible = # of dual reactions these species are involved in

      npossible = npair[ispecies][jspecies];
      if (npossible == 0) {
	j = plist[j].next;
	continue;
      }
	
      // no need to apply PBC because ghost particles are already remapped

      delx = plist[i].x[0] - plist[j].x[0];
      dely = plist[i].x[1] - plist[j].x[1];
      delz = plist[i].x[2] - plist[j].x[2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= distsq[ispecies][jspecies]) {
	j = plist[j].next;
	continue;
      }

      // if j particle is ghost, add reaction to original image particle
	
      for (k = 0; k < npossible; k++) {
	m = pairlist[ispecies][jspecies][k];
	if (blist[jbin].ghost) add_reaction(m,i,plist[j].seed);
	else add_reaction(m,i,j);
      }
	
      j = plist[j].next;
    }
  }
}

/* ----------------------------------------------------------------------
   compute propensity of a single reaction
   for mono or dual reaction: propensity = count * rate
------------------------------------------------------------------------- */

double ChemSpatialSSA::compute_propensity(int m)
{
  double p;
  if (nreactant[m] == 0) p = rate[m];
  else p = nrlist[m] * rate[m];
  return p;
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::add_reaction(int type, int part1)
{
  int n = nrlist[type]++;
  if (n == maxrlist[type]) resize_rlist_one(type);

  rlist[type][n].part1 = part1;
  rlist[type][n].edge1 = add_edge(part1,type,n);
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::add_reaction(int type, int part1, int part2)
{
  int n = nrlist[type]++;
  if (n == maxrlist[type]) resize_rlist_one(type);

  rlist[type][n].part1 = part1;
  rlist[type][n].part2 = part2;
  add_edge(part1,part2,type,n,rlist[type][n].edge1,rlist[type][n].edge2);
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::delete_reaction(int type, int ireact)
{
  /*
  if (type < 0 || type >= nreactions)
    error->all("Illegal type of reaction to delete");
  if (ireact < 0 || ireact >= nrlist[type]) {
    printf("Illegal reaction %d %d\n",ireact,nrlist[type]);
    error->all("Illegal reaction to delete");
  }
  */

  int last = nrlist[type] - 1;
  nrlist[type] = last;
  Reaction *r1 = &rlist[type][ireact];
  Reaction *r2 = &rlist[type][last];

  elist[r1->part1][r1->edge1].active = 0;
  r1->part1 = r2->part1;
  r1->edge1 = r2->edge1;
  elist[r1->part1][r1->edge1].ireact = ireact;
  if (nreactant[type] == 1) return;

  elist[r1->part2][r1->edge2].active = 0;
  r1->part2 = r2->part2;
  r1->edge2 = r2->edge2;
  elist[r1->part2][r1->edge2].ireact = ireact;
}

/* ---------------------------------------------------------------------- */

int ChemSpatialSSA::add_edge(int part1, int type, int ireact)
{
  int edge1 = nelist[part1]++;
  if (edge1 == maxelist[part1]) resize_elist_one(part1);

  elist[part1][edge1].type = type;
  elist[part1][edge1].ireact = ireact;
  elist[part1][edge1].active = 1;

  return edge1;
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::add_edge(int part1, int part2, int type, int ireact,
			  int &edge1, int &edge2)
{
  edge1 = nelist[part1]++;
  if (edge1 == maxelist[part1]) resize_elist_one(part1);
  edge2 = nelist[part2]++;
  if (edge2 == maxelist[part2]) resize_elist_one(part2);

  elist[part1][edge1].type = type;
  elist[part1][edge1].ireact = ireact;
  elist[part1][edge1].active = 1;

  elist[part2][edge2].type = type;
  elist[part2][edge2].ireact = ireact;
  elist[part2][edge2].active = 1;
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::delete_edges(int part1)
{
  Edge *list = elist[part1];
  int n = nelist[part1];
  for (int i = 0; i < n; i++) {
    if (!list[i].active) continue;
    delete_reaction(list[i].type,list[i].ireact);
  }
}

/* ----------------------------------------------------------------------
   sum entire tree, all nodes are computed
------------------------------------------------------------------------- */

void ChemSpatialSSA::sum()
{
  int child1,child2;

  for (int parent = offset-1; parent >= 0; parent--) {
    child1 = 2*parent + 1;
    child2 = 2*parent + 2;
    tree[parent] = tree[child1] + tree[child2];
  }
}

/* ----------------------------------------------------------------------
   value = uniform RN from 0 to tree[0]
   return index (0 to M-1) of propensity bin it falls in
------------------------------------------------------------------------- */

int ChemSpatialSSA::find(double value)
{
  int i,leftchild;

  // i walks tree from root to appropriate leaf
  // value is modified when right branch of tree is traversed

  i = 0;
  while (i < offset) {
    leftchild = 2*i + 1;
    if (value <= tree[leftchild]) i = leftchild;
    else {
      value -= tree[leftchild];
      i = leftchild + 1;
    }
  }
  return i - offset;
}

/* ----------------------------------------------------------------------
   setup reaction stencils based on current decomposition
   stencil1 is 1/2 a full stencil (self + 13)
   stencil2 is a full stencil (self + 26)
   self bin is always stencil[0]
------------------------------------------------------------------------- */

void ChemSpatialSSA::setup_stencil()
{
  int nbinx = grid->nbinx;
  int nbiny = grid->nbiny;

  stencil1[ 0] = 0;
  stencil1[ 1] = 1;
  stencil1[ 2] = nbinx - 1;
  stencil1[ 3] = nbinx;
  stencil1[ 4] = nbinx + 1;
  stencil1[ 5] = nbiny*nbinx - nbinx - 1;
  stencil1[ 6] = nbiny*nbinx - nbinx;
  stencil1[ 7] = nbiny*nbinx - nbinx + 1;
  stencil1[ 8] = nbiny*nbinx - 1;
  stencil1[ 9] = nbiny*nbinx;
  stencil1[10] = nbiny*nbinx + 1;
  stencil1[11] = nbiny*nbinx + nbinx - 1;
  stencil1[12] = nbiny*nbinx + nbinx;
  stencil1[13] = nbiny*nbinx + nbinx + 1;

  stencil2[ 0] = 0;
  stencil2[ 1] = -nbiny*nbinx - nbinx - 1;
  stencil2[ 2] = -nbiny*nbinx - nbinx;
  stencil2[ 3] = -nbiny*nbinx - nbinx + 1;
  stencil2[ 4] = -nbiny*nbinx - 1;
  stencil2[ 5] = -nbiny*nbinx;
  stencil2[ 6] = -nbiny*nbinx + 1;
  stencil2[ 7] = -nbiny*nbinx + nbinx - 1;
  stencil2[ 8] = -nbiny*nbinx + nbinx;
  stencil2[ 9] = -nbiny*nbinx + nbinx + 1;
  stencil2[10] = -nbinx - 1;
  stencil2[11] = -nbinx;
  stencil2[12] = -nbinx + 1;
  stencil2[13] = -1;
  stencil2[14] = 1;
  stencil2[15] = nbinx - 1;
  stencil2[16] = nbinx;
  stencil2[17] = nbinx + 1;
  stencil2[18] = nbiny*nbinx - nbinx - 1;
  stencil2[19] = nbiny*nbinx - nbinx;
  stencil2[20] = nbiny*nbinx - nbinx + 1;
  stencil2[21] = nbiny*nbinx - 1;
  stencil2[22] = nbiny*nbinx;
  stencil2[23] = nbiny*nbinx + 1;
  stencil2[24] = nbiny*nbinx + nbinx - 1;
  stencil2[25] = nbiny*nbinx + nbinx;
  stencil2[26] = nbiny*nbinx + nbinx + 1;
}

/* ----------------------------------------------------------------------
   free arrays used by Gillspie solver
------------------------------------------------------------------------- */

void ChemSpatialSSA::free_arrays()
{
  delete [] tree;
  delete [] rcount;

  memory->destroy_2d_double_array(distsq);

  delete [] rate;

  for (int i = 0; i < nreactions; i++) memory->sfree(rlist[i]);
  memory->sfree(maxrlist);
  memory->sfree(nrlist);
  memory->sfree(rlist);

  for (int i = 0; i < maxtotal; i++) memory->sfree(elist[i]);
  memory->sfree(maxelist);
  memory->sfree(nelist);
  memory->sfree(elist);

  maxtotal = 0;
  elist = NULL;
  maxelist = NULL;
  nelist = NULL;
}

/* ----------------------------------------------------------------------
   allocate elist to current length of ntotal
------------------------------------------------------------------------- */

void ChemSpatialSSA::resize_elist()
{
  if (particle->ntotal > maxtotal) {
    int nmax = particle->ntotal + CHUNK;
    elist = (Edge **) memory->srealloc(elist,nmax*sizeof(Edge *),"sssa:elist");
    maxelist = (int *) 
      memory->srealloc(maxelist,nmax*sizeof(int),"sssa:maxelist");
    nelist = (int *) memory->srealloc(nelist,nmax*sizeof(int),"sssa:nelist");

    for (int i = maxtotal; i < nmax; i++) {
      elist[i] = NULL;
      maxelist[i] = nelist[i] = 0;
    }
    maxtotal = nmax;
  }
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::resize_elist_one(int i)
{
  maxelist[i] += CHUNK;
  elist[i] = (Edge *) 
    memory->srealloc(elist[i],maxelist[i]*sizeof(Edge),"sssa:elist");
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::resize_rlist()
{
  rlist = (Reaction **) 
    memory->smalloc(nreactions*sizeof(Reaction *),"sssa:rlist");
  maxrlist = (int *) memory->smalloc(nreactions*sizeof(int),"sssa:maxrlist");
  nrlist = (int *) memory->smalloc(nreactions*sizeof(int),"sssa:nrlist");

  for (int i = 0; i < nreactions; i++) {
    rlist[i] = NULL;
    maxrlist[i] = nrlist[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */

void ChemSpatialSSA::resize_rlist_one(int i)
{
  maxrlist[i] += CHUNK;
  rlist[i] = (Reaction *) 
    memory->srealloc(rlist[i],maxrlist[i]*sizeof(Reaction),"sssa:rlist");
}
