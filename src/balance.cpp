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

#include "mpi.h"
#include "balance.h"
#include "grid.h"

#define BIN_WEIGHT 0.0001

/* ---------------------------------------------------------------------- */

Balance::Balance()
{
  float version;
  Zoltan_Initialize(0,NULL,&version);
  RCB_setup();
  MPI_Comm_rank(world,&me);
}

/* ---------------------------------------------------------------------- */

Balance::~Balance()
{
  if (lb) Zoltan_Destroy(&lb);
  if (lb_new) Zoltan_Destroy(&lb_new);
}

/* ----------------------------------------------------------------------
   setup the bin balancer - called once
------------------------------------------------------------------------- */

void Balance::RCB_setup()
{
  // create the bin balancers

  lb = Zoltan_Create(world);
  lb_new = Zoltan_Create(world);

  // use RCB to balance bins
  // RCB_RECTILINEAR_BLOCKS and AVERAGE_CUTS insure cuts are between bins
  // REMAP = 0 turns off processor remapping at end of RCB
  // define balancer callbacks
  // callback for weight computation is set at time of RCB

  Zoltan_Struct *lbptr;
  for (int i = 0; i < 2; i++) {
    if (i == 0) lbptr = lb;
    else lbptr = lb_new;

    Zoltan_Set_Param(lbptr,"DEBUG_LEVEL","0");
    Zoltan_Set_Param(lbptr,"DEBUG_MEMORY","0");
    Zoltan_Set_Param(lbptr,"CHECK_GEOM","0");
    Zoltan_Set_Param(lbptr,"OBJ_WEIGHT_DIM","1");
    Zoltan_Set_Param(lbptr,"LB_METHOD","RCB");
    Zoltan_Set_Param(lbptr,"RETURN_LISTS","NONE");
    Zoltan_Set_Param(lbptr,"KEEP_CUTS","1");
    Zoltan_Set_Param(lbptr,"RCB_RECTILINEAR_BLOCKS","1");
    Zoltan_Set_Param(lbptr,"AVERAGE_CUTS","1");
    // Zoltan_Set_Param(lbptr,"REMAP","0");

    Zoltan_Set_Fn(lbptr,ZOLTAN_NUM_GEOM_FN_TYPE,
		  (void (*)()) geometry_dimension,NULL); 
    Zoltan_Set_Fn(lbptr,ZOLTAN_NUM_OBJ_FN_TYPE,
		  (void (*)()) bin_number,NULL); 
    Zoltan_Set_Fn(lbptr,ZOLTAN_GEOM_MULTI_FN_TYPE,
		  (void (*)()) bin_coords,NULL); 
  }
}

/* ----------------------------------------------------------------------
   balance bins across procs via RCB on bin centers
------------------------------------------------------------------------- */

void Balance::RCB(int weightflag, 
		 double *xlo, double *ylo, double *zlo, 
		 double *xhi, double *yhi, double *zhi)
{
  int change,num_gid,num_lid,num_import,num_export,dim;
  ZOLTAN_ID_PTR import_global_ids,import_local_ids;
  ZOLTAN_ID_PTR export_global_ids,export_local_ids;
  int *import_procs,*export_procs;

  // set correct weight-computation function based on weightflag

  if (weightflag == 0)
    Zoltan_Set_Fn(lb_new,ZOLTAN_OBJ_LIST_FN_TYPE,
		  (void (*)()) bin_wts_bin,NULL); 
  else
    Zoltan_Set_Fn(lb_new,ZOLTAN_OBJ_LIST_FN_TYPE,
		  (void (*)()) bin_wts_particle,NULL); 

  // generate debug files for Zoltan

  // Zoltan_Generate_Files(lb_new,"tmp.bad",1,1,0,0);

  // perform balance using lb_new

  Zoltan_LB_Balance(lb_new,&change,&num_gid,&num_lid, 
		    &num_import,&import_global_ids,&import_local_ids,
		    &import_procs,
		    &num_export,&export_global_ids,&export_local_ids,
		    &export_procs);

  // convert RCB cuts to lo/hi bounds

  Zoltan_RCB_Box(lb_new,me,&dim,xlo,ylo,zlo,xhi,yhi,zhi);

  // free LB memory

  Zoltan_LB_Free_Data(&import_global_ids,&import_local_ids,&import_procs, 
		      &export_global_ids,&export_local_ids,&export_procs);
}

/* ----------------------------------------------------------------------
   swap lb ptrs so new decomp becomes current and vice versa
------------------------------------------------------------------------- */

void Balance::swap()
{
  Zoltan_Struct *tmp = lb;
  lb = lb_new;
  lb_new = tmp;
}

/* ----------------------------------------------------------------------
   static callback functions
   cannot reference anything in Balance class
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   geometry_dimension: return number of dimensions for global problem
------------------------------------------------------------------------- */

int Balance::geometry_dimension(void *null, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}

/* ----------------------------------------------------------------------
   bin_number: return number of bins owned by this proc (no ghosts)
------------------------------------------------------------------------- */

int Balance::bin_number(void *null, int *ierr)
{
  int nbins = (grid->nbinx - 2) * (grid->nbiny - 2) * (grid->nbinz - 2);
  *ierr = ZOLTAN_OK;
  return nbins;
}

/* ----------------------------------------------------------------------
   bin_coords: return x,y,z centers for a list of requested bins
------------------------------------------------------------------------- */

void Balance::bin_coords(void *null, int num_global, int num_local, 
			int nbins, ZOLTAN_ID_PTR globals, ZOLTAN_ID_PTR locals,
			int ndim, double *coords, int *ierr)
{
  int n = 0;
  for (int i = 0; i < nbins; i++) {
    grid->coord_center(globals[i],&coords[n],&coords[n+1],&coords[n+2]);
    n += 3;
  }
  *ierr = ZOLTAN_OK;
}

/* ----------------------------------------------------------------------
   bin_wts_bin: fill lists of global & local IDs & weights for owned bins
   weight of each bin = BIN_WEIGHT
------------------------------------------------------------------------- */

void Balance::bin_wts_bin(void *null, int num_global, int num_local,
			 ZOLTAN_ID_PTR globals, ZOLTAN_ID_PTR locals,
			 int wgt_dim, float *weights, int *ierr)
{
  Grid::OneBin *blist = grid->blist;
  int nbins = grid->nbins;

  int n = 0;
  for (int i = 0; i < nbins; i++)
    if (blist[i].ghost == 0) {
      globals[n] = blist[i].id;
      locals[n] = i;
      weights[n] = BIN_WEIGHT;
      n++;
    }
  *ierr = ZOLTAN_OK;
}

/* ----------------------------------------------------------------------
   bin_wts_particle: fill lists of global & local IDs & weights for owned bins
   weight of each bin = BIN_WEIGHT + # of particles in bin
------------------------------------------------------------------------- */

void Balance::bin_wts_particle(void *null, int num_global, int num_local,
			      ZOLTAN_ID_PTR globals, ZOLTAN_ID_PTR locals,
			      int wgt_dim, float *weights, int *ierr)
{
  Grid::OneBin *blist = grid->blist;
  int nbins = grid->nbins;

  int n = 0;
  for (int i = 0; i < nbins; i++)
    if (blist[i].ghost == 0) {
      globals[n] = blist[i].id;
      locals[n] = i;
      weights[n] = BIN_WEIGHT + blist[i].nparts;
      n++;
    }
  *ierr = ZOLTAN_OK;
}
