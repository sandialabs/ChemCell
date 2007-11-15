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

#ifndef CHEM_SPATIAL_SSA_H
#define CHEM_SPATIAL_SSA_H

#include "chem.h"

class ChemSpatialSSA : public Chem {
 public:
  ChemSpatialSSA();
  ~ChemSpatialSSA();
  void init();
  void reactions();
  void setup_stencil();

 private:
  struct Edge {               // single edge = reaction for a particle
    int type,ireact;          // type and index of associated reaction
    int active;               // 1 if active reaction, 0 if deleted
  };

  struct Reaction {           // single reaction
    int part1,edge1;          // particle and edge involved in reaction
    int part2,edge2;          // 2nd particle and edge involved in reaction
                              // part2,edge2 are unset if 1 reactant
  };

  // data structures for storing reactions used by spatial/ssa algorithm
  // a reaction is between a pair of nearby particles
  // each reaction stored once in rlist
  // dual reactions are stored twice in elist (with each atom)

  Edge **elist;               // list of Edges for each owned atom
  int *nelist,*maxelist;      // # of Edges stored for each ntotal atoms
  int maxtotal;               // max # of total atoms in elist

  Reaction **rlist;           // list of Reactions for each reaction type
  int *nrlist,*maxrlist;      // # of Reactions stored for each reaction type

  int nreactions;
  int *nreactant;
  int **reactants;
  int *nproduct;
  int **products;
  int *nmono;
  int **monolist;
  int **npair;
  int ***pairlist;

  int allocated;
  double *tree;              // tree structure of summed propensities
  double *propensity;        // propensity of each reaction
  int offset;                // loc in tree where 1st leaf starts
  double *rate;              // reaction rates for spatial SSA

  double **distsq;

  int nstencil1,nstencil2;
  int *stencil1,*stencil2;

  void consistency();
  void build_elist_all();
  void build_elist_one(int);
  double compute_propensity(int);
  void add_reaction(int, int);
  void add_reaction(int, int, int);
  void delete_reaction(int, int);
  int add_edge(int, int, int);
  void add_edge(int, int, int, int, int &, int &);
  void delete_edges(int);
  void sum();
  int find(double);
  void setup_stencil(int, int, int);
  void free_arrays();
  void resize_elist();
  void resize_elist_one(int);
  void resize_rlist();
  void resize_rlist_one(int);
};

#endif
