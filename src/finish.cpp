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
#include "mpi.h"
#include "finish.h"
#include "simulator.h"
#include "particle.h"
#include "chem.h"
#include "move.h"
#include "grid.h"
#include "surf.h"
#include "react.h"
#include "timer.h"
#include "memory.h"

/* ---------------------------------------------------------------------- */

void Finish::end()
{
  int i;
  int histo[10];
  double time,tmp,ave,max,min;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // adjust TIME_REACT for TIME_REACT_COMM

  timer->array[TIME_REACT] -= timer->array[TIME_REACT_COMM];

  // total loop time

  double time_other = timer->array[TIME_LOOP] -
    (timer->array[TIME_MOVE] + timer->array[TIME_MIGRATE] + 
     timer->array[TIME_REACT] + timer->array[TIME_REACT_COMM] + 
     timer->array[TIME_OUTPUT] + timer->array[TIME_BALANCE]);

  double time_loop = timer->array[TIME_LOOP];
  MPI_Allreduce(&time_loop,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time_loop = tmp/nprocs;

  if (me == 0) {
    fprintf(screen,
	    "Loop time of %g on %d procs for %d steps\n\n",
	    time_loop,nprocs,simulator->nsteps);
    fprintf(logfile,
	    "Loop time of %g on %d procs for %d steps\n\n",
	    time_loop,nprocs,simulator->nsteps);
  }

  // time for each portion of simulation

  time = timer->array[TIME_MOVE];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  if (me == 0) {
    fprintf(screen,"Move  time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    fprintf(logfile,"Move  time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  time = timer->array[TIME_MIGRATE];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  if (me == 0) {
    fprintf(screen,"Migrt time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    fprintf(logfile,"Migrt time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  time = timer->array[TIME_REACT];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  if (me == 0) {
    fprintf(screen,"React time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    fprintf(logfile,"React time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  time = timer->array[TIME_REACT_COMM];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  if (me == 0) {
    fprintf(screen,"RComm time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    fprintf(logfile,"RComm time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  time = timer->array[TIME_OUTPUT];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  if (me == 0) {
    fprintf(screen,"Outpt time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    fprintf(logfile,"Outpt time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  time = timer->array[TIME_BALANCE];
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  if (me == 0) {
    fprintf(screen,"Balnc time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    fprintf(logfile,"Balnc time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }


  time = time_other;
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  if (me == 0) {
    fprintf(screen,"Other time (%%) = %g (%g)\n",time,time/time_loop*100.0);
    fprintf(logfile,"Other time (%%) = %g (%g)\n",time,time/time_loop*100.0);
  }

  if (me == 0) {
    fprintf(screen,"\n");
    fprintf(logfile,"\n");
  }

  if (simulator->spatial_flag) {

    // per processor counts

    tmp = particle->nlocal;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      fprintf(screen,"Nlocal:    %g ave %g max %g min\n",ave,max,min);
      fprintf(screen,"Histogram:");
      for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
      fprintf(screen,"\n");
      fprintf(logfile,"Nlocal:    %g ave %g max %g min\n",ave,max,min);
      fprintf(logfile,"Histogram:");
      for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
      fprintf(logfile,"\n");
    }
    
    tmp = particle->nghost_last;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      fprintf(screen,"Nghost:    %g ave %g max %g min\n",ave,max,min);
      fprintf(screen,"Histogram:");
      for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
      fprintf(screen,"\n");
      fprintf(logfile,"Nghost:    %g ave %g max %g min\n",ave,max,min);
      fprintf(logfile,"Histogram:");
      for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
      fprintf(logfile,"\n");
    }
    
    tmp = grid->nbinx * grid->nbiny * grid->nbinz;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      fprintf(screen,"Nbin:    %g ave %g max %g min\n",ave,max,min);
      fprintf(screen,"Histogram:");
      for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
      fprintf(screen,"\n");
      fprintf(logfile,"Nbin:    %g ave %g max %g min\n",ave,max,min);
      fprintf(logfile,"Histogram:");
      for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
      fprintf(logfile,"\n");
    }
  
    if (me == 0) {
      fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }
  }

  // total counts for moves, reactions, balancing

  int nsteps = simulator->nsteps;
  if (nsteps == 0) nsteps = 1;

  if (simulator->spatial_flag) {
    int nmove,ncheck,nreflect,nnear,nstick,nfar,nthru;
    MPI_Allreduce(&move->nmove,&nmove,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&move->ncheck,&ncheck,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&move->nreflect,&nreflect,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&move->nnear,&nnear,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&move->nstick,&nstick,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&move->nfar,&nfar,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&move->nthru,&nthru,1,MPI_INT,MPI_SUM,world);
    
    if (me == 0) {
      if (screen) {
	fprintf(screen,"Move statistics (total & per-step):\n");
	fprintf(screen,"  moves      = %d %g\n",nmove,1.0*nmove/nsteps);
	fprintf(screen,"  tri checks = %d %g\n",ncheck,1.0*ncheck/nsteps);
	fprintf(screen,"  refl  hits = %d %g\n",nreflect,1.0*nreflect/nsteps);
	fprintf(screen,"  near  hits = %d %g\n",nnear,1.0*nnear/nsteps);
	fprintf(screen,"  stick hits = %d %g\n",nstick,1.0*nstick/nsteps);
	fprintf(screen,"  far   hits = %d %g\n",nfar,1.0*nfar/nsteps);
	fprintf(screen,"  thru  hits = %d %g\n",nthru,1.0*nthru/nsteps);
      }
      if (logfile) {
	fprintf(logfile,"Move statistics (total & per-step):\n");
	fprintf(logfile,"  moves      = %d %g\n",nmove,1.0*nmove/nsteps);
	fprintf(logfile,"  tri checks = %d %g\n",ncheck,1.0*ncheck/nsteps);
	fprintf(logfile,"  refl  hits = %d %g\n",nreflect,1.0*nreflect/nsteps);
	fprintf(logfile,"  near  hits = %d %g\n",nnear,1.0*nnear/nsteps);
	fprintf(logfile,"  stick hits = %d %g\n",nstick,1.0*nstick/nsteps);
	fprintf(logfile,"  far   hits = %d %g\n",nfar,1.0*nfar/nsteps);
	fprintf(logfile,"  thru  hits = %d %g\n",nthru,1.0*nthru/nsteps);
      }
    }
  }

  int nbinbin;
  MPI_Allreduce(&chem->nbinbin,&nbinbin,1,MPI_INT,MPI_SUM,world);
  int nbinpair;
  MPI_Allreduce(&chem->nbinpair,&nbinpair,1,MPI_INT,MPI_SUM,world);
  int ndist;
  MPI_Allreduce(&chem->ndist,&ndist,1,MPI_INT,MPI_SUM,world);
  int noverlap;
  MPI_Allreduce(&chem->noverlap,&noverlap,1,MPI_INT,MPI_SUM,world);
  
  int nreactions = react->nreactions;
  int *rcount = new int[nreactions];
  MPI_Allreduce(chem->rcount,rcount,nreactions,MPI_INT,MPI_SUM,world);

  int nreact = 0;
  for (int i = 0; i < nreactions; i++) nreact += rcount[i];

  if (me == 0) {
    if (screen) {
      fprintf(screen,"Reaction statistics (total & per-step):\n");
      fprintf(screen,"  bin-bin     = %d %g\n",nbinbin,1.0*nbinbin/nsteps);
      fprintf(screen,"  bin pairs   = %d %g\n",nbinpair,1.0*nbinpair/nsteps);
      fprintf(screen,"  dist checks = %d %g\n",ndist,1.0*ndist/nsteps);
      fprintf(screen,"  overlaps    = %d %g\n",noverlap,1.0*noverlap/nsteps);
      fprintf(screen,"  reactions   = %d %g\n",nreact,1.0*nreact/nsteps);
      fprintf(screen,"  count of each reaction: ");
      for (int i = 0; i < nreactions; i++)
	fprintf(screen,"(%s %d) ",react->name[i],rcount[i]);
      fprintf(screen,"\n");
    }
    if (logfile) {
      fprintf(logfile,"Reaction statistics (total & per-step):\n");
      fprintf(logfile,"  bin-bin     = %d %g\n",nbinbin,1.0*nbinbin/nsteps);
      fprintf(logfile,"  bin pairs   = %d %g\n",nbinpair,1.0*nbinpair/nsteps);
      fprintf(logfile,"  dist checks = %d %g\n",ndist,1.0*ndist/nsteps);
      fprintf(logfile,"  overlaps    = %d %g\n",noverlap,1.0*noverlap/nsteps);
      fprintf(logfile,"  reactions   = %d %g\n",nreact,1.0*nreact/nsteps);
      fprintf(logfile,"  count of each reaction: ");
      for (int i = 0; i < nreactions; i++)
	fprintf(logfile,"(%s %d) ",react->name[i],rcount[i]);
      fprintf(logfile,"\n");
    }
  }

  delete [] rcount;

  if (simulator->spatial_flag && me == 0) {
    if (screen) fprintf(screen,"Number of balancing calls = %d\n",grid->ncount);
    if (logfile) fprintf(logfile,"Number of balancing calls = %d\n",grid->ncount);
  }

  memory_usage();

  // print mapping of particle species name and type just to log file
  // useful for visualizer output

  if (me == 0) {
    fprintf(logfile,"Equivalance map of species & type\n");
    for (int i = 0; i < particle->nspecies; i++)
      fprintf(logfile,"  map %s %d\n",particle->name[i],i+1);
  }
}

/* ----------------------------------------------------------------------
   print out memory usage
------------------------------------------------------------------------- */

void Finish::memory_usage()
{
  int pbyte = particle->memory_usage();
  int bbyte = 0;
  int sbyte = 0;
  if (simulator->spatial_flag) {
    bbyte = grid->memory_usage();
    sbyte = surf->memory_usage();
  }

  int byte = pbyte + bbyte + sbyte;

  int pbyte_ave,pbyte_max,pbyte_min;
  MPI_Allreduce(&pbyte,&pbyte_ave,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&pbyte,&pbyte_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&pbyte,&pbyte_min,1,MPI_INT,MPI_MIN,world);

  int bbyte_ave,bbyte_max,bbyte_min;
  MPI_Allreduce(&bbyte,&bbyte_ave,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&bbyte,&bbyte_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&bbyte,&bbyte_min,1,MPI_INT,MPI_MIN,world);

  int sbyte_ave,sbyte_max,sbyte_min;
  MPI_Allreduce(&sbyte,&sbyte_ave,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&sbyte,&sbyte_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&sbyte,&sbyte_min,1,MPI_INT,MPI_MIN,world);

  int byte_ave,byte_max,byte_min;
  MPI_Allreduce(&byte,&byte_ave,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&byte,&byte_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&byte,&byte_min,1,MPI_INT,MPI_MIN,world);

  pbyte_ave /= nprocs;
  bbyte_ave /= nprocs;
  sbyte_ave /= nprocs;
  byte_ave /= nprocs;

  double s = 1.0/1024/1024;

  if (me == 0) {
    if (screen) {
      fprintf(screen,"Memory usage in Mbyte/proc (ave/max/min)\n");
      fprintf(screen,"  parts = %g %g %g\n",
	      s*pbyte_ave,s*pbyte_max,s*pbyte_min);
      fprintf(screen,"  bins  = %g %g %g\n",
	      s*bbyte_ave,s*bbyte_max,s*bbyte_min);
      fprintf(screen,"  surfs = %g %g %g\n",
	      s*sbyte_ave,s*sbyte_max,s*sbyte_min);
      fprintf(screen,"  total = %g %g %g\n",
	      s*byte_ave,s*byte_max,s*byte_min);
    }
    if (logfile) {
      fprintf(logfile,"Memory usage in Mbyte/proc (ave/max/min)\n");
      fprintf(logfile,"  parts = %g %g %g\n",
	      s*pbyte_ave,s*pbyte_max,s*pbyte_min);
      fprintf(logfile,"  bins  = %g %g %g\n",
	      s*bbyte_ave,s*bbyte_max,s*bbyte_min);
      fprintf(logfile,"  surfs = %g %g %g\n",
	      s*sbyte_ave,s*sbyte_max,s*sbyte_min);
      fprintf(logfile,"  total = %g %g %g\n",
	      s*byte_ave,s*byte_max,s*byte_min);
    }
  }
}

/* ----------------------------------------------------------------------
   histogram stats
------------------------------------------------------------------------- */

void Finish::stats(int n, double *data, 
		   double *pave, double *pmax, double *pmin,
		   int nhisto, int *histo)
{
  int i,m;
  int *histotmp;

  double min = 1.0e20;
  double max = -1.0e20;
  double ave = 0.0;
  for (i = 0; i < n; i++) {
    ave += data[i];
    if (data[i] < min) min = data[i];
    if (data[i] > max) max = data[i];
  }

  int ntotal;
  MPI_Allreduce(&n,&ntotal,1,MPI_INT,MPI_SUM,world);
  double tmp;
  MPI_Allreduce(&ave,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  ave = tmp/ntotal;
  MPI_Allreduce(&min,&tmp,1,MPI_DOUBLE,MPI_MIN,world);
  min = tmp;
  MPI_Allreduce(&max,&tmp,1,MPI_DOUBLE,MPI_MAX,world);
  max = tmp;

  for (i = 0; i < nhisto; i++) histo[i] = 0;

  double del = max - min;
  for (i = 0; i < n; i++) {
    if (del == 0.0) m = 0;
    else m = static_cast<int> ((data[i]-min)/del * nhisto);
    if (m > nhisto-1) m = nhisto-1;
    histo[m]++;
  }

  histotmp = (int *) memory->smalloc(nhisto*sizeof(int),"finish:histotmp");
  MPI_Allreduce(histo,histotmp,nhisto,MPI_INT,MPI_SUM,world);
  for (i = 0; i < nhisto; i++) histo[i] = histotmp[i];
  memory->sfree(histotmp);

  *pave = ave;
  *pmax = max;
  *pmin = min;
}
