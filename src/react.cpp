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

#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "react.h"
#include "particle.h"
#include "surf.h"
#include "region.h"
#include "random.h"
#include "memory.h"
#include "error.h"

#define DEFAULT 1
#define AT      2
#define NEAR    3
#define ONE     4
#define TWO     5
#define ONETWO  6
#define INSIDE  7
#define OUTSIDE 8
#define INOUT   9

#define MAX_PRODUCT  5
#define EPSILON 1.0e-4
#define BIG 1.0e20

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

React::React()
{
  allocated = 0;

  nreactions = 0;
  name = NULL;
  nreactant = NULL;
  reactants = NULL;
  nproduct = NULL;
  products = NULL;
  rate = NULL;
  setwhich = NULL;
  setstyle = NULL;
  setdir = NULL;
  setdist = NULL;
  setprob = NULL;

  wreactant = NULL;
  wproduct = NULL;

  npair = NULL;
  pairlist = NULL;
  nmono = NULL;
  monolist = NULL;

  locwhich = NULL;
  locstyle = NULL;
  locdir = NULL;
}

/* ---------------------------------------------------------------------- */

React::~React()
{
  for (int i = 0; i < nreactions; i++) delete [] name[i];
  memory->sfree(name);

  memory->sfree(nreactant);
  memory->destroy_2d_int_array(reactants);
  memory->sfree(nproduct);
  memory->destroy_2d_int_array(products);
  memory->sfree(rate);
  memory->destroy_2d_int_array(setwhich);
  memory->destroy_2d_int_array(setstyle);
  memory->destroy_2d_int_array(setdir);
  memory->sfree(setdist);
  memory->sfree(setprob);
  memory->destroy_2d_double_array(wreactant);
  memory->destroy_2d_double_array(wproduct);

  if (allocated) free_arrays();
}

/* ----------------------------------------------------------------------
   add a new reaction
------------------------------------------------------------------------- */

void React::add(int narg, char **arg)
{
  if (narg < 3) error->all("Illegal reaction command");

  // store ID

  if (find(arg[0]) >= 0) {
    char *str = new char[128];
    sprintf(str,"Reaction ID %s already exists",arg[0]);
    error->all(str);
  }

  int n = nreactions + 1;
  name = (char **) memory->srealloc(name,n*sizeof(char *),"react:name");
  int nlen = strlen(arg[0]) + 1;
  name[nreactions] = new char[nlen];
  strcpy(name[nreactions],arg[0]);

  // grow all reaction arrays

  nreactant = (int *) 
    memory->srealloc(nreactant,n*sizeof(int),"react:nreactnant");
  reactants = memory->grow_2d_int_array(reactants,n,2,"react:reactants");
  nproduct = (int *) memory->srealloc(nproduct,n*sizeof(int),"react:nproduct");
  products = 
    memory->grow_2d_int_array(products,n,MAX_PRODUCT,"react:products");
  rate = (double *) memory->srealloc(rate,n*sizeof(double),"react:rate");
  setwhich = 
    memory->grow_2d_int_array(setwhich,n,MAX_PRODUCT,"react:setwhich");
  setstyle = 
    memory->grow_2d_int_array(setstyle,n,MAX_PRODUCT,"react:setstyle");
  setdir = 
    memory->grow_2d_int_array(setdir,n,MAX_PRODUCT,"react:setdir");
  setdist = (double *)
    memory->srealloc(setdist,n*sizeof(double),"react:setdist");
  setprob = (double *)
    memory->srealloc(setprob,n*sizeof(double),"react:setprob");
  wreactant = memory->grow_2d_double_array(wreactant,n,2,"react:wreactant");
  wproduct =
    memory->grow_2d_double_array(wproduct,n,MAX_PRODUCT,"react:wproduct");

  // find which arg is numeric reaction rate

  char c;
  int iarg = 1;
  while (iarg < narg) {
    c = arg[iarg][0];
    if ((c >= '0' && c <= '9') || c == '+' || c == '-' || c == '.') break;
    iarg++;
  }

  // error checks

  if (iarg == narg) error->all("Reaction has no numeric rate");
  if (iarg < 1 || iarg > 3) 
    error->all("Reaction must have 0,1,2 reactants");
  if (narg-1 - iarg > MAX_PRODUCT) 
    error->all("Reaction cannot have more than MAX_PRODUCT products");

  // extract reactant and product species names
  // if any species does not exist, create it

  nreactant[nreactions] = 0;
  for (int i = 1; i < iarg; i++) {
    int ispecies = particle->find(arg[i]);
    if (ispecies == -1) error->all("Unknown species in reaction command");
    reactants[nreactions][i-1] = ispecies;
    nreactant[nreactions]++;
  }

  rate[nreactions] = atof(arg[iarg]);

  nproduct[nreactions] = 0;
  for (int i = iarg+1; i < narg; i++) {
    int ispecies = particle->find(arg[i]);
    if (ispecies == -1) error->all("Unknown species in reaction command");
    products[nreactions][i - (iarg+1)] = ispecies;
    nproduct[nreactions]++;
  }

  // set defaults for product locations and dist/prob and volume weights

  for (int i = 0; i < nproduct[nreactions]; i++)
    setwhich[nreactions][i] = setstyle[nreactions][i] = 
      setdir[nreactions][i] = DEFAULT;
  setdist[nreactions] = -1.0;
  setprob[nreactions] = -1.0;

  for (int i = 0; i < nreactant[nreactions]; i++)
    wreactant[nreactions][i] = 1.0;
  for (int i = 0; i < nproduct[nreactions]; i++)
    wproduct[nreactions][i] = 1.0;
  
  nreactions++;
}

/* ----------------------------------------------------------------------
   delete a reaction
------------------------------------------------------------------------- */

void React::delete_reaction(char *str)
{
  // find which reaction it is and delete it

  int ireact = find(str);
  if (ireact == -1) error->all("Could not find unreact ID");
  delete [] name[ireact];

  // move other reactions down in list one slot

  for (int i = ireact+1; i < nreactions; i++) {
    name[i-1] = name[i];
    nreactant[i-1] = nreactant[i];
    nproduct[i-1] = nproduct[i];

    for (int j = 0; j < nreactant[i-1]; j++) {
      reactants[i-1][j] = reactants[i][j];
      wreactant[i-1][j] = wreactant[i][j];
    }
    for (int j = 0; j < nproduct[i-1]; j++) {
      products[i-1][j] = products[i][j];
      wproduct[i-1][j] = wproduct[i][j];
      setwhich[i-1][j] = setwhich[i][j];
      setstyle[i-1][j] = setstyle[i][j];
      setdir[i-1][j] = setdir[i][j];
    }

    rate[i-1] = rate[i];
    setdist[i-1] = setdist[i];
    setprob[i-1] = setprob[i];
  }
  nreactions--;
}

/* ----------------------------------------------------------------------
   modify the settings of a reaction
------------------------------------------------------------------------- */

void React::modify(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal react_modify command");

  int ireact = find(arg[0]);
  if (ireact == -1) error->all("Could not find reaction ID");

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"rate") == 0) {
      if (narg < iarg+2) error->all("Illegal react_modify command");
      rate[ireact] = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"loc") == 0) {
      if (narg < iarg+5) error->all("Illegal react_modify command");

      int iproduct;
      if (strcmp(arg[iarg+1],"all") == 0) iproduct = -1;
      else iproduct = atoi(arg[iarg+1]);
      if (iproduct == 0 || iproduct > nproduct[ireact])
	error->all("Illegal react_modify command");

      int value1,value2,value3;
      if (strcmp(arg[iarg+2],"def") == 0) value1 = DEFAULT;
      else if (strcmp(arg[iarg+2],"at") == 0) value1 = AT;
      else if (strcmp(arg[iarg+2],"near") == 0) value1 = NEAR;
      else error->all("Illegal react_modify command");

      if (strcmp(arg[iarg+3],"def") == 0) value2 = DEFAULT;
      else if (strcmp(arg[iarg+3],"1") == 0) value2 = ONE;
      else if (strcmp(arg[iarg+3],"2") == 0) value2 = TWO;
      else if (strcmp(arg[iarg+3],"1/2") == 0) value2 = ONETWO;
      else error->all("Illegal react_modify command");

      if (strcmp(arg[iarg+4],"def") == 0) value3 = DEFAULT;
      else if (strcmp(arg[iarg+4],"in") == 0) value3 = INSIDE;
      else if (strcmp(arg[iarg+4],"out") == 0) value3 = OUTSIDE;
      else if (strcmp(arg[iarg+4],"in/out") == 0) value3 = INOUT;
      else error->all("Illegal react_modify command");

      if (iproduct == -1) {
	for (int i = 0; i < nproduct[ireact]; i++) {
	  setstyle[ireact][i] = value1;
	  setwhich[ireact][i] = value2;
	  setdir[ireact][i] = value3;
	}
      } else {
	setstyle[ireact][iproduct-1] = value1;
	setwhich[ireact][iproduct-1] = value2;
	setdir[ireact][iproduct-1] = value3;
      }

      iarg += 5;
    } else if (strcmp(arg[iarg],"dist") == 0) {
      if (narg < iarg+2) error->all("Illegal react_modify command");
      setdist[ireact] = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"prob") == 0) {
      if (narg < iarg+2) error->all("Illegal react_modify command");
      setprob[ireact] = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"weight") == 0) {
      if (narg < iarg+4) error->all("Illegal react_modify command");
      int which;
      if (strcmp(arg[iarg+1],"reactant") == 0) which = 0;
      else if (strcmp(arg[iarg+1],"product") == 0) which = 1;
      else error->all("Illegal react_modify command");
      int index = atoi(arg[iarg+2]);
      double value = atof(arg[iarg+3]);
      if (which == 0 && (index <= 0 || index > nreactant[ireact]))
	error->all("Illegal react_modify command");
      if (which == 1 && (index <= 0 || index > nproduct[ireact]))
	error->all("Illegal react_modify command");
      if (value < 0.0) error->all("Illegal react_modify command");
      if (which == 0) wreactant[ireact][index-1] = value;
      else wproduct[ireact][index-1] = value;
      iarg += 4;
    } else error->all("Illegal react_modify command");
  }
}

/* ----------------------------------------------------------------------
   initialize various reaction counters, lists, and settings
------------------------------------------------------------------------- */

void React::init()
{
  int ispecies,jspecies,ireact;

  if (allocated) free_arrays();
  allocated = 1;

  // mono reaction arrays
  
  int nspecies = particle->nspecies;
  nmono = new int[nspecies];
  
  for (ispecies = 0; ispecies < nspecies; ispecies++)
    nmono[ispecies] = 0;
  for (ireact = 0; ireact < nreactions; ireact++)
    if (nreactant[ireact] == 1) {
      ispecies = reactants[ireact][0];
      nmono[ispecies]++;
    }
  
  int max = 0;
  for (ispecies = 0; ispecies < nspecies; ispecies++)
    max = MAX(max,nmono[ispecies]);
  
  monolist = memory->create_2d_int_array(nspecies,max,"react:monolist");
  
  for (ispecies = 0; ispecies < nspecies; ispecies++) nmono[ispecies] = 0;
  for (ireact = 0; ireact < nreactions; ireact++)
    if (nreactant[ireact] == 1) {
      ispecies = reactants[ireact][0];
      monolist[ispecies][nmono[ispecies]++] = ireact;
    }
  
  // dual (pair) reaction arrays
  // max = max # of reactions for any pair of species
			 
  npair = memory->create_2d_int_array(nspecies,nspecies,"react:npair");
  
  for (ispecies = 0; ispecies < nspecies; ispecies++)
    for (jspecies = 0; jspecies < nspecies; jspecies++)
      npair[ispecies][jspecies] = 0;
  
  for (ireact = 0; ireact < nreactions; ireact++)
    if (nreactant[ireact] == 2) {
      ispecies = reactants[ireact][0];
      jspecies = reactants[ireact][1];
      npair[ispecies][jspecies]++;
      if (ispecies != jspecies) npair[jspecies][ispecies]++;
    }

  max = 0;
  for (ispecies = 0; ispecies < nspecies; ispecies++)
    for (jspecies = 0; jspecies < nspecies; jspecies++)
      max = MAX(max,npair[ispecies][jspecies]);

  pairlist = 
    memory->create_3d_int_array(nspecies,nspecies,max,"react:pairlist");

  for (ispecies = 0; ispecies < nspecies; ispecies++)
    for (jspecies = 0; jspecies < nspecies; jspecies++)
      npair[ispecies][jspecies] = 0;
  
  for (ireact = 0; ireact < nreactions; ireact++)
    if (nreactant[ireact] == 2) {
      ispecies = reactants[ireact][0];
      jspecies = reactants[ireact][1];
      pairlist[ispecies][jspecies][npair[ispecies][jspecies]++] = ireact;
      if (ispecies != jspecies)
	pairlist[jspecies][ispecies][npair[jspecies][ispecies]++] = ireact;
    }

  // reaction arrays

  locwhich = 
    memory->create_2d_int_array(nreactions,MAX_PRODUCT,"react:locwhich");
  locstyle = 
    memory->create_2d_int_array(nreactions,MAX_PRODUCT,"react:locstyle");
  locdir = 
    memory->create_2d_int_array(nreactions,MAX_PRODUCT,"react:locdir");
  
  // set location flags for placement of each reaction product
  // idim,jdim = dimensionality of reactant species
  // kdim = dimensionality of product species
  // idiff,jdiff = diffusion coeffs of 2 reactants
  // check for errors in user requests

  for (ireact = 0; ireact < nreactions; ireact++)

    // mono reactions

    if (nreactant[ireact] == 1) {
      int idim = particle->dimension[reactants[ireact][0]];
      for (int k = 0; k < nproduct[ireact]; k++) {
	int kdim = particle->dimension[products[ireact][k]];

	// check all possible which/style/dir settings

	if (setwhich[ireact][k] == DEFAULT) locwhich[ireact][k] = ONE;
	if (setwhich[ireact][k] == ONE) locwhich[ireact][k] = ONE;
	if (setwhich[ireact][k] == TWO) invalid(ireact);
	if (setwhich[ireact][k] == ONETWO) invalid(ireact);

	if (setstyle[ireact][k] == DEFAULT) {
	  if (idim == 3 && kdim == 3) locstyle[ireact][k] = AT;
	  else if (idim == 2 && kdim == 2) locstyle[ireact][k] = AT;
	  else if (idim == 3 && kdim == 2) invalid(ireact);
	  else if (idim == 2 && kdim == 3) locstyle[ireact][k] = NEAR;
	}
	if (setstyle[ireact][k] == AT) {
	  if (idim == 3 && kdim == 3) locstyle[ireact][k] = AT;
	  else if (idim == 2 && kdim == 2) locstyle[ireact][k] = AT;
	  else if (idim == 3 && kdim == 2) invalid(ireact);
	  else if (idim == 3 && kdim == 2) invalid(ireact);
	}
	if (setstyle[ireact][k] == NEAR) {
	  if (idim == 3 && kdim == 3) invalid(ireact);
	  else if (idim == 2 && kdim == 2) invalid(ireact);
	  else if (idim == 3 && kdim == 2) invalid(ireact);
	  else if (idim == 3 && kdim == 2) locstyle[ireact][k] = NEAR;
	}

	if (setdir[ireact][k] == DEFAULT) {
	  if (idim == 2 && kdim == 3) locdir[ireact][k] = INSIDE;
	}
	if (setdir[ireact][k] == INSIDE) {
	  if (idim == 3 && kdim == 3) invalid(ireact);
	  else if (idim == 2 && kdim == 2) invalid(ireact);
	  else if (idim == 3 && kdim == 2) invalid(ireact);
	  else if (idim == 2 && kdim == 3) locdir[ireact][k] = INSIDE;
	}
	if (setdir[ireact][k] == OUTSIDE) {
	  if (idim == 3 && kdim == 3) invalid(ireact);
	  else if (idim == 2 && kdim == 2) invalid(ireact);
	  else if (idim == 3 && kdim == 2) invalid(ireact);
	  else if (idim == 2 && kdim == 3) locdir[ireact][k] = OUTSIDE;
	}
	if (setdir[ireact][k] == INOUT) {
	  if (idim == 3 && kdim == 3) invalid(ireact);
	  else if (idim == 2 && kdim == 2) invalid(ireact);
	  else if (idim == 3 && kdim == 2) invalid(ireact);
	  else if (idim == 2 && kdim == 3) locdir[ireact][k] = INOUT;
	}
      }

      // dual reactions
	
    } else if (nreactant[ireact] == 2) {
      int idim = particle->dimension[reactants[ireact][0]];
      int jdim = particle->dimension[reactants[ireact][1]];
      double idiff = particle->diffusivity[reactants[ireact][0]];
      double jdiff = particle->diffusivity[reactants[ireact][1]];

      for (int k = 0; k < nproduct[ireact]; k++) {
	int kdim = particle->dimension[products[ireact][k]];

	// check all possible which/style/dir settings

	if (setwhich[ireact][k] == DEFAULT) {
	  if (idim == 3 && jdim == 3 && kdim == 3) {
	    if (idiff <= jdiff) locwhich[ireact][k] = ONE;
	    else locwhich[ireact][k] = TWO;
	  } else if (idim == 2 && jdim == 2 && kdim == 2) {
	    if (idiff <= jdiff) locwhich[ireact][k] = ONE;
	    else locwhich[ireact][k] = TWO;
	  } else if (idim == 3 && jdim == 2 && kdim == 3)
	    locwhich[ireact][k] = ONE;
	  else if (idim == 2 && jdim == 3 && kdim == 3)
	    locwhich[ireact][k] = TWO;
	  else if (idim == 3 && jdim == 2 && kdim == 2)
	    locwhich[ireact][k] = TWO;
	  else if (idim == 2 && jdim == 3 && kdim == 2)
	    locwhich[ireact][k] = ONE;
	  else if (idim == 3 && jdim == 3 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 2 && kdim == 3) {
	    if (idiff <= jdiff) locwhich[ireact][k] = ONE;
	    else locwhich[ireact][k] = TWO;
	  }
	}
	if (setwhich[ireact][k] == ONE) {
	  if (idim == 3 && jdim == 3 && kdim == 3)
	    locwhich[ireact][k] = ONE;
	  else if (idim == 2 && jdim == 2 && kdim == 2)
	    locwhich[ireact][k] = ONE;
	  else if (idim == 3 && jdim == 2 && kdim == 3)
	    locwhich[ireact][k] = ONE;
	  else if (idim == 2 && jdim == 3 && kdim == 3)
	    locwhich[ireact][k] = ONE;
	  else if (idim == 3 && jdim == 2 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 3 && kdim == 2)
	    locwhich[ireact][k] = ONE;
	  else if (idim == 3 && jdim == 3 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 2 && kdim == 3)
	    locwhich[ireact][k] = ONE;
	}
	if (setwhich[ireact][k] == TWO) {
	  if (idim == 3 && jdim == 3 && kdim == 3)
	    locwhich[ireact][k] = TWO;
	  else if (idim == 2 && jdim == 2 && kdim == 2)
	    locwhich[ireact][k] = TWO;
	  else if (idim == 3 && jdim == 2 && kdim == 3)
	    locwhich[ireact][k] = TWO;
	  else if (idim == 2 && jdim == 3 && kdim == 3)
	    locwhich[ireact][k] = TWO;
	  else if (idim == 3 && jdim == 2 && kdim == 2)
	    locwhich[ireact][k] = TWO;
	  else if (idim == 2 && jdim == 3 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 3 && jdim == 3 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 2 && kdim == 3)
	    locwhich[ireact][k] = TWO;
	}
	if (setwhich[ireact][k] == ONETWO) {
	  if (idim == 3 && jdim == 3 && kdim == 3)
	    locwhich[ireact][k] = ONETWO;
	  else if (idim == 2 && jdim == 2 && kdim == 2)
	    locwhich[ireact][k] = ONETWO;
	  else if (idim == 3 && jdim == 2 && kdim == 3)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 3 && kdim == 3)
	    invalid(ireact);
	  else if (idim == 3 && jdim == 2 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 3 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 3 && jdim == 3 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 2 && kdim == 3)
	    locwhich[ireact][k] = ONETWO;
	}

	if (setstyle[ireact][k] == DEFAULT) {
	  if (idim == 3 && jdim == 3 && kdim == 3)
	    locstyle[ireact][k] = AT;
	  else if (idim == 2 && jdim == 2 && kdim == 2) 
	    locstyle[ireact][k] = AT;
	  else if (idim == 3 && jdim == 2 && kdim == 3) {
	    if (locwhich[ireact][k] == TWO) invalid(ireact);
	    else locstyle[ireact][k] = AT;
	  } else if (idim == 2 && jdim == 3 && kdim == 3) {
	    if (locwhich[ireact][k] == ONE) invalid(ireact);
	    else locstyle[ireact][k] = AT;
	  } else if (idim == 3 && jdim == 2 && kdim == 2) {
	    if (locwhich[ireact][k] == ONE) invalid(ireact);
	    else locstyle[ireact][k] = AT;
	  } else if (idim == 2 && jdim == 3 && kdim == 2) {
	    if (locwhich[ireact][k] == TWO) invalid(ireact);
	    else locstyle[ireact][k] = AT;
	  } else if (idim == 3 && jdim == 3 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 2 && kdim == 3)
	    locstyle[ireact][k] = NEAR;
	}
	if (setstyle[ireact][k] == AT) {
	  if (idim == 3 && jdim == 3 && kdim == 3)
	    locstyle[ireact][k] = AT;
	  else if (idim == 2 && jdim == 2 && kdim == 2) 
	    locstyle[ireact][k] = AT;
	  else if (idim == 3 && jdim == 2 && kdim == 3) {
	    if (locwhich[ireact][k] == TWO) invalid(ireact);
	    else locstyle[ireact][k] = AT;
	  } else if (idim == 2 && jdim == 3 && kdim == 3) {
	    if (locwhich[ireact][k] == ONE) invalid(ireact);
	    else locstyle[ireact][k] = AT;
	  } else if (idim == 3 && jdim == 2 && kdim == 2) {
	    if (locwhich[ireact][k] == ONE) invalid(ireact);
	    else locstyle[ireact][k] = AT;
	  } else if (idim == 2 && jdim == 3 && kdim == 2) {
	    if (locwhich[ireact][k] == TWO) invalid(ireact);
	    else locstyle[ireact][k] = AT;
	  } else if (idim == 3 && jdim == 3 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 2 && kdim == 3)
	    invalid(ireact);
	}
	if (setstyle[ireact][k] == NEAR) {
	  if (idim == 3 && jdim == 3 && kdim == 3)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 2 && kdim == 2) 
	    invalid(ireact);
	  else if (idim == 3 && jdim == 2 && kdim == 3) {
	    if (locwhich[ireact][k] == ONE) invalid(ireact);
	    else locstyle[ireact][k] = NEAR;
	  } else if (idim == 2 && jdim == 3 && kdim == 3) {
	    if (locwhich[ireact][k] == TWO) invalid(ireact);
	    else locstyle[ireact][k] = NEAR;
	  } else if (idim == 3 && jdim == 2 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 3 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 3 && jdim == 3 && kdim == 2)
	    invalid(ireact);
	  else if (idim == 2 && jdim == 2 && kdim == 3)
	    locstyle[ireact][k] = NEAR;
	}

	if (setdir[ireact][k] == DEFAULT) locdir[ireact][k] = INSIDE;
	else locdir[ireact][k] = setdir[ireact][k];
      }
    }

  // print reaction stats

  int count0 = 0;
  int count1 = 0;
  int count2 = 0;
  double min0 = BIG;
  double min1 = BIG;
  double min2 = BIG;
  double max0 = -1.0;
  double max1 = -1.0;
  double max2 = -1.0;

  for (ireact = 0; ireact < nreactions; ireact++) {
    if (nreactant[ireact] == 0) {
      count0++;
      min0 = MIN(rate[ireact],min0);
      max0 = MAX(rate[ireact],max0);
    }
    if (nreactant[ireact] == 1) {
      count1++;
      min1 = MIN(rate[ireact],min1);
      max1 = MAX(rate[ireact],max1);
    }
    if (nreactant[ireact] == 2) {
      count2++;
      min2 = MIN(rate[ireact],min2);
      max2 = MAX(rate[ireact],max2);
    }
  }

  int me;
  MPI_Comm_rank(world,&me);

  if (me == 0) {
    fprintf(screen,"Reactions:\n");
    fprintf(screen,"  number of 0-reactant reactions = %d\n",count0);
    if (count0)
      fprintf(screen,"  min/max 0-reactant rates = %g %g\n",min0,max0);
    fprintf(screen,"  number of 1-reactant reactions = %d\n",count1);
    if (count1)
      fprintf(screen,"  min/max 1-reactant rates = %g %g\n",min1,max1);
    fprintf(screen,"  number of 2-reactant reactions = %d\n",count2);
    if (count1)
      fprintf(screen,"  min/max 2-reactant rates = %g %g\n",min2,max2);

    fprintf(logfile,"Reactions:\n");
    fprintf(logfile,"  number of 0-reactant reactions = %d\n",count0);
    if (count0)
      fprintf(logfile,"  min/max 0-reactant rates = %g %g\n",min0,max0);
    fprintf(logfile,"  number of 1-reactant reactions = %d\n",count1);
    if (count1)
      fprintf(logfile,"  min/max 1-reactant rates = %g %g\n",min1,max1);
    fprintf(logfile,"  number of 2-reactant reactions = %d\n",count2);
    if (count1)
      fprintf(logfile,"  min/max 2-reactant rates = %g %g\n",min2,max2);
  }
}

/* ----------------------------------------------------------------------
   return reaction index (0 to N-1) for a reaction ID name
   return -1 if doesn't exist
------------------------------------------------------------------------- */

int React::find(char *str)
{
  for (int i = 0; i < nreactions; i++)
    if (strcmp(str,name[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   print an error message for an invalid reaction
------------------------------------------------------------------------- */

void React::invalid(int ireact)
{
  char str[128];
  sprintf(str,"Invalid settings for reaction %s",name[ireact]);
  error->all(str);
}

/* ----------------------------------------------------------------------
   create Mth product particle of IREACT reaction
------------------------------------------------------------------------- */

int React::create_product(int m, int ireact, int ip, int jp)
{
  int kp,mp;
  double *normal;
  double regnormal[3];
  Particle::OnePart *plist = particle->plist;

  if (jp >= 0) {

    // dual reaction
    // locwhich determines which reactant location to use
    //   ONE or TWO or ONETWO (random)
    // locstyle determines whether creation is AT or NEAR a reactant
    // locdir determines (for NEAR) which direction from surface
    //   the loc is at: INSIDE or OUTSIDE or INOUT (random)
    // use paticle->add() to add the particle
    // update plist in case it was re-allocated
    // set seed & bin & triangle for new products

    if (locwhich[ireact][m] == ONE) {
      if (plist[ip].species == reactants[ireact][0]) kp = ip;
      else kp = jp;
    } else if (locwhich[ireact][m] == TWO) {
      if (plist[ip].species == reactants[ireact][1]) kp = ip;
      else kp = jp;
    } else {
      if (random->move(&plist[ip].seed) < 0.5) kp = ip;
      else kp = jp;
    }

    if (locstyle[ireact][m] == AT) {
      mp = particle->add(products[ireact][m],
			 plist[kp].x[0],plist[kp].x[1],plist[kp].x[2]);
      plist = particle->plist;

      plist[mp].seed = random->product(&plist[ip].seed);
      if (kp == ip) plist[mp].ibin = plist[ip].ibin;
      else plist[mp].ibin = plist[jp].ibin;
      if (particle->dimension[products[ireact][m]] == 3)
	plist[mp].itri = 0;
      else plist[mp].itri = plist[kp].itri;
    
    } else {
      double sign;
      if (locdir[ireact][m] == INSIDE) sign = -1.0;
      else if (locdir[ireact][m] == OUTSIDE) sign = 1.0;
      else if (random->move(&plist[ip].seed) < 0.5) sign = -1.0;
      else sign = 1.0;
      
      int itri = plist[kp].itri;
      if (itri >= 0) normal = surf->tlist[itri].normal;
      else {
	normal = regnormal;
	surf->rlist[-itri-1]->compute_normal(plist[kp].x,normal);
      }
      
      double xnew[3];
      xnew[0] = plist[kp].x[0] + sign*EPSILON*normal[0];
      xnew[1] = plist[kp].x[1] + sign*EPSILON*normal[1];
      xnew[2] = plist[kp].x[2] + sign*EPSILON*normal[2];
      
      mp = particle->add(products[ireact][m],xnew[0],xnew[1],xnew[2]);
      plist = particle->plist;
      
      plist[mp].seed = random->product(&plist[ip].seed);
      plist[mp].ibin = grid->whichlocal(xnew);
      plist[mp].itri = 0;
    }
    
  } else {

    // mono reaction
    // locstyle determines whether creation is AT or NEAR a reactant
    // locdir determines (for NEAR) which direction from surface
    //   the loc is at: INSIDE or OUTSIDE or INOUT (random)
    // use paticle->add() to add the particle
    // update plist in case it was re-allocated
    // set seed & bin & triangle for new products

    if (locstyle[ireact][m] == AT) {
      mp = particle->add(products[ireact][m],
			 plist[ip].x[0],plist[ip].x[1],plist[ip].x[2]);
      plist = particle->plist;
      
      plist[mp].seed = random->product(&plist[ip].seed);
      plist[mp].ibin = plist[ip].ibin;
      if (particle->dimension[products[ireact][m]] == 3)
	plist[mp].itri = 0;
      else plist[mp].itri = plist[ip].itri;
    
    } else {
      double sign;
      if (locdir[ireact][m] == INSIDE) sign = -1.0;
      else if (locdir[ireact][m] == OUTSIDE) sign = 1.0;
      else if (random->move(&plist[ip].seed) < 0.5) sign = -1.0;
      else sign = 1.0;
      
      int itri = plist[ip].itri;
      if (itri >= 0) normal = surf->tlist[itri].normal;
      else {
	normal = regnormal;
	surf->rlist[-itri-1]->compute_normal(plist[kp].x,normal);
      }
      
      double xnew[3];
      xnew[0] = plist[ip].x[0] + sign*EPSILON*normal[0];
      xnew[1] = plist[ip].x[1] + sign*EPSILON*normal[1];
      xnew[2] = plist[ip].x[2] + sign*EPSILON*normal[2];

      mp = particle->add(products[ireact][m],xnew[0],xnew[1],xnew[2]);
      plist = particle->plist;
      
      plist[mp].seed = random->product(&plist[ip].seed);
      plist[mp].ibin = grid->whichlocal(xnew);
      plist[mp].itri = 0;
    }
  }

  return mp;
}

/* ----------------------------------------------------------------------
   free arrays
------------------------------------------------------------------------- */

void React::free_arrays()
{
  if (nmono) delete [] nmono;
  memory->destroy_2d_int_array(monolist);

  memory->destroy_2d_int_array(npair);
  memory->destroy_3d_int_array(pairlist);

  memory->destroy_2d_int_array(locwhich);
  memory->destroy_2d_int_array(locstyle);
  memory->destroy_2d_int_array(locdir);
}
