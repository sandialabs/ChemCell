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

// universal classes for ChemCell styles

#ifdef CommandInclude
#include "move_test.h"
#include "read_restart.h"
#include "run.h"
#include "write_restart.h"
#endif

#ifdef CommandClass
CommandStyle(move_test,MoveTest)
CommandStyle(read_restart,ReadRestart)
CommandStyle(run,Run)
CommandStyle(write_restart,WriteRestart)
#endif

#ifdef FixInclude
#include "fix_conc_random.h"
#include "fix_conc_set.h"
#endif

#ifdef FixClass
FixStyle(conc/random,FixConcRandom)
FixStyle(conc/set,FixConcSet)
#endif

#ifdef RegionInclude
#include "region_box.h"
#include "region_cylinder.h"
#include "region_plane.h"
#include "region_sphere.h"
#endif

#ifdef RegionClass
RegionStyle(box,RegionBox)
RegionStyle(cylinder,RegionCylinder)
RegionStyle(plane,RegionPlane)
RegionStyle(sphere,RegionSphere)
#endif

#ifdef SimulatorInclude
#include "sim_gillespie.h"
#include "sim_ode.h"
#include "sim_ode_rk.h"
#include "sim_spatial.h"
#endif

#ifdef SimulatorClass
SimulatorStyle(gillespie,SimGillespie)
SimulatorStyle(ode,SimODE)
SimulatorStyle(ode/rk,SimODERK)
SimulatorStyle(spatial,SimSpatial)
#endif
