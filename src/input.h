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

#ifndef INPUT_H
#define INPUT_H

#include "stdio.h"
#include "system.h"

class Variable;

class Input : public System {
 public:
  int narg;                      // # of command args
  char **arg;                    // parsed args for command
  Variable *variable;            // defined variables

  Input(int, char **);
  ~Input();
  char *next();                  // return next executable command from script
  char *one(char *);             // execute a single passed-in command
  void substitute(char *, int);  // substitute for variables in a string

 private:
  int me;                      // proc ID
  char *command;               // ptr to current command
  int maxarg;                  // max # of args in arg
  char *line,*copy,*work;      // input line & copy of it
  int echo_screen;             // 0 = no, 1 = yes
  int echo_log;                // 0 = no, 1 = yes
  int nfile,maxfile;           // current # and max # of open input files
  int label_active;            // 0 = no label, 1 = looking for label
  char *labelstr;              // label string being looked for
  int jump_skip;               // 1 if skipping next jump, 0 otherwise

  FILE **infiles;              // list of open input files

  void parse();                // parse an input text line
  int execute_command();       // execute a single command

  void balance();              // individual commands
  void bin();
  void boundary();
  void cd();
  void check();
  void clear();
  void count();
  void debug();
  void diffusion();
  void dimension();
  void dump();
  void dump_modify();
  void echo();
  void fix();
  void global();
  void include();
  void jump();
  void label();
  void log();
  void move_style();
  void move_test();
  void next_command();
  void particles();
  void permeable();
  void print();
  void probability();
  void react_modify();
  void reaction();
  void region();
  void restart();
  void run();
  void run_style();
  void seed();
  void species();
  void stats();
  void stats_modify();
  void timestep();
  void triangles();
  void undump();
  void unfix();
  void unreact();
  void variable_command();
  void volume();
};

#endif