#include "string.h"
#include "stdlib.h"
#include "mpi.h"
#include "zoltan.h"

#define BIG 1.0e20

struct Zoltan_Comm_Obj {
  int n;
};

/* ------------------------------------------------------------------------ */

struct Zoltan_Struct *Zoltan_Create(MPI_Comm communicator)
{
  return NULL;
}

void Zoltan_Destroy(struct Zoltan_Struct **zz) {}

int Zoltan_Initialize(int argc, char **argv, float *ver)
{
  return 0;
}

int Zoltan_Set_Param(struct Zoltan_Struct *zz, const char *name,
		     const char *val)
{
  return 0;
}

int Zoltan_Set_Fn(struct Zoltan_Struct *zz, ZOLTAN_FN_TYPE fn_type,
		  void (*fn_ptr)(), void *data_ptr)
{
  return 0;
}

int Zoltan_LB_Balance(
  struct Zoltan_Struct *zz,
  int *changes,
  int *num_gid_entries,
  int *num_lid_entries,
  int *num_import,
  ZOLTAN_ID_PTR *import_global_ids,
  ZOLTAN_ID_PTR *import_local_ids,
  int **import_procs,
  int *num_export,
  ZOLTAN_ID_PTR *export_global_ids,
  ZOLTAN_ID_PTR *export_local_ids,
  int **export_procs)
{
  return 0;
}

int Zoltan_LB_Free_Data (
  ZOLTAN_ID_PTR *import_global_ids,
  ZOLTAN_ID_PTR *import_local_ids,
  int **import_procs, 
  ZOLTAN_ID_PTR *export_global_ids,
  ZOLTAN_ID_PTR *export_local_ids,
  int **export_procs)
{
  return 0;
}

int Zoltan_LB_Point_Assign(struct Zoltan_Struct *zz, double *coords, int *proc)
{ 
  *proc = 0;

  return 0;
}

int Zoltan_LB_Box_Assign(
  struct Zoltan_Struct *zz,
  double xmin,
  double ymin,
  double zmin,
  double xmax,
  double ymax,
  double zmax,
  int *procs,
  int *numprocs)
{
  procs[0] = 0;
  *numprocs = 1;

  return 0;
}

int Zoltan_RCB_Box(
  struct Zoltan_Struct *zz,
  int     part,
  int     *dim,
  double *xmin,
  double *ymin,
  double *zmin,
  double *xmax,
  double *ymax,
  double *zmax)
{
  *dim = 3;
  *xmin = -BIG;
  *ymin = -BIG;
  *zmin = -BIG;
  *xmax = BIG;
  *ymax = BIG;
  *zmax = BIG;

  return 0;
}

int Zoltan_Comm_Do(Zoltan_Comm_Obj *plan, int tag,
		   char *sendbuf, int nsize, char *recvbuf)
{
  memcpy(recvbuf,sendbuf,nsize*plan->n);
  return 0;
}

int Zoltan_Comm_Destroy(Zoltan_Comm_Obj **pplan)
{
  free(*pplan);
  *pplan = NULL;
  return 0;
}

int Zoltan_Comm_Create(Zoltan_Comm_Obj **pplan, int nsend, int *proclist,
		       MPI_Comm world, int tag, int *nrecv)
{
  Zoltan_Comm_Obj *plan;
  plan = (Zoltan_Comm_Obj *) malloc(sizeof(Zoltan_Comm_Obj));
  plan->n = nsend;

  *pplan = plan;
  *nrecv = nsend;
  return 0;
}

