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

#include "stdio.h"
#include "string.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

#define FixInclude
#include "style.h"
#undef FixInclude

#define DELTA 1

// mask settings - same as mask settings in fix.cpp

#define INITIAL  1
#define FINAL    2
#define CLEANUP  4

/* ---------------------------------------------------------------------- */

Modify::Modify()
{
  nfix = maxfix = 0;
  n_initial = n_final = n_cleanup = 0;

  fix = NULL;
  fmask = NULL;
  list_initial = NULL;
  list_final = NULL;
  list_cleanup = NULL;
}

/* ---------------------------------------------------------------------- */

Modify::~Modify()
{
  // delete all fixes

  for (int i = 0; i < nfix; i++) delete fix[i];
  memory->sfree(fix);
  memory->sfree(fmask);

  delete [] list_initial;
  delete [] list_final;
  delete [] list_cleanup;
}

/* ----------------------------------------------------------------------
   initialize all fixes and lists of fixes
------------------------------------------------------------------------- */

void Modify::init()
{
  int i;

  // init each fix

  for (i = 0; i < nfix; i++) fix[i]->init();

  // create lists of fixes to call at each stage of run

  list_init(INITIAL,n_initial,list_initial);
  list_init(INITIAL,n_final,list_final);
  list_init(CLEANUP,n_cleanup,list_cleanup);
}

/* ----------------------------------------------------------------------
   beginning of timestep
------------------------------------------------------------------------- */

void Modify::initial()
{
  for (int i = 0; i < n_initial; i++)
    fix[list_initial[i]]->initial();
}

/* ----------------------------------------------------------------------
   end of timestep
------------------------------------------------------------------------- */

void Modify::final()
{
  for (int i = 0; i < n_final; i++)
    fix[list_final[i]]->final();
}

/* ----------------------------------------------------------------------
   end of run
------------------------------------------------------------------------- */

void Modify::cleanup()
{
  for (int i = 0; i < n_cleanup; i++)
    fix[list_cleanup[i]]->cleanup();
}

/* ----------------------------------------------------------------------
   add a new fix or replace one with same ID
------------------------------------------------------------------------- */

void Modify::add_fix(int narg, char **arg)
{
  if (narg < 2) error->all("Illegal fix command");

  // if fix ID exists:
  //   set newflag = 0 so create new fix in same location in fix list
  //   error if new style does not match old style
  //     since can't replace it (all when-to-invoke ptrs would be invalid)
  //   delete old fix
  //   set ptr to NULL in case new fix scans list of fixes
  // if fix ID does not exist:
  //   set newflag = 1 so create new fix
  //   extend fix and fmask lists as necessary

  int ifix,newflag;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(arg[0],fix[ifix]->id) == 0) break;

  if (ifix < nfix) {
    newflag = 0;
    if (strcmp(arg[1],fix[ifix]->style) != 0)
      error->all("Replacing a fix, but new style != old style");
    delete fix[ifix];
    fix[ifix] = NULL;
  } else {
    newflag = 1;
    if (nfix == maxfix) {
      maxfix += DELTA;
      fix = (Fix **) memory->srealloc(fix,maxfix*sizeof(Fix *),"modify:fix");
      fmask = (int *) 
	memory->srealloc(fmask,maxfix*sizeof(int),"modify:fmask");
    }
  }

  // create the Fix

  if (0) return;         // dummy line to enable else-if macro expansion

#define FixClass
#define FixStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) fix[ifix] = new Class(narg,arg);
#include "style.h"
#undef FixClass

  else error->all("Invalid fix style");

  // if fix is new, set it's mask values and increment nfix

  if (newflag) {
    fmask[ifix] = fix[ifix]->setmask();
    nfix++;
  }
}

/* ----------------------------------------------------------------------
   delete a Fix from list of Fixes
------------------------------------------------------------------------- */

void Modify::delete_fix(char *id)
{
  int ifix;

  // find which fix it is and delete it

  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,fix[ifix]->id) == 0) break;
  if (ifix == nfix) error->all("Could not find unfix ID");

  delete fix[ifix];

  // move other Fixes and fmask down in list one slot

  for (int i = ifix+1; i < nfix; i++) fix[i-1] = fix[i];
  for (int i = ifix+1; i < nfix; i++) fmask[i-1] = fmask[i];
  nfix--;
}

/* ----------------------------------------------------------------------
   create list of fix indices for fixes which match mask
------------------------------------------------------------------------- */

void Modify::list_init(int mask, int &n, int *&list)
{
  delete [] list;
  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];
  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) list[n++] = i;
}
