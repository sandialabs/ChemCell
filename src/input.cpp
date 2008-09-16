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

#include "unistd.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "input.h"
#include "variable.h"
#include "universe.h"
#include "simulator.h"
#include "particle.h"
#include "surf.h"
#include "domain.h"
#include "chem.h"
#include "move.h"
#include "modify.h"
#include "random.h"
#include "output.h"
#include "react.h"
#include "grid.h"
#include "memory.h"
#include "error.h"

#define CommandInclude
#define SimulatorInclude
#include "style.h"
#undef CommandInclude
#undef SimulatorInclude

#define MAXLINE 20000
#define DELTA 4

/* ---------------------------------------------------------------------- */

Input::Input(int argc, char **argv)
{
  MPI_Comm_rank(world,&me);

  line = new char[MAXLINE];
  copy = new char[MAXLINE];
  work = new char[MAXLINE];
  narg = maxarg = 0;
  arg = NULL;

  echo_screen = 0;
  echo_log = 1;

  label_active = 0;
  labelstr = NULL;
  jump_skip = 0;

  if (me == 0) {
    nfile = maxfile = 1;
    infiles = (FILE **) memory->smalloc(sizeof(FILE *),"input:infiles");
    infiles[0] = infile;
  } else infiles = NULL;

  variable = new Variable;

  // process command-line args
  // check for args "-var" and "-echo"
  // caller has already checked that sufficient arguments exist

  int iarg = 0;
  while (iarg < argc) {
    if (strcmp(argv[iarg],"-var") == 0) {
      variable->set(argv[iarg+1],argv[iarg+2]);
      iarg += 3;
    } else if (strcmp(argv[iarg],"-echo") == 0) {
      narg = 1;
      char **tmp = arg;        // trick echo() into using argv instead of arg
      arg = &argv[iarg+1];
      echo();
      arg = tmp;
      iarg += 2;
    } else iarg++;
  }
}

/* ---------------------------------------------------------------------- */

Input::~Input()
{
  // don't free command and arg strings
  // they just point to other allocated memory

  delete variable;
  delete [] line;
  delete [] copy;
  delete [] work;
  if (labelstr) delete [] labelstr;
  if (arg) memory->sfree(arg);
  if (infiles) memory->sfree(infiles);
}

/* ----------------------------------------------------------------------
   process all input from infile
   infile = stdin or file if command-line arg "-in" was used
------------------------------------------------------------------------- */

void Input::file()
{
  int n;

  while (1) {
    
    // read one line from input script
    // if line ends in continuation char '&', concatenate next line(s)
    // n = str length of line
    
    if (me == 0) {
      if (fgets(line,MAXLINE,infile) == NULL) n = 0;
      else n = strlen(line) + 1;
      while (n >= 3 && line[n-3] == '&') {
	if (fgets(&line[n-3],MAXLINE-n+3,infile) == NULL) n = 0;
	else n = strlen(line) + 1;
      }
    }

    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      if (label_active) error->all("Label wasn't found in input script");
      if (me == 0) {
	if (infile != stdin) fclose(infile);
	nfile--;
      }
      MPI_Bcast(&nfile,1,MPI_INT,0,world);
      if (nfile == 0) break;
      if (me == 0) infile = infiles[nfile-1];
      continue;
    }

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // if n = MAXLINE, line is too long

    if (n == MAXLINE) {
      char str[MAXLINE+32];
      sprintf(str,"Input line too long: %s",line);
      error->all(str);
    }

    // echo the command unless scanning for label

    if (me == 0 && label_active == 0) {
      if (echo_screen && screen) fprintf(screen,"%s",line); 
      if (echo_log && logfile) fprintf(logfile,"%s",line);
    }

    // parse the line
    // if no command, skip to next line in input script

    parse();
    if (command == NULL) continue;

    // if scanning for label, skip command unless it's a label command

    if (label_active && strcmp(command,"label") != 0) continue;

    // execute the command

    if (execute_command()) {
      char str[MAXLINE];
      sprintf(str,"Unknown command: %s",line);
      error->all(str);
    }
  }
}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void Input::file(const char *filename)
{
  // error if another nested file still open
  // if single open file is not stdin, close it
  // open new filename and set infile, infiles[0]

  if (me == 0) {
    if (nfile > 1)
      error->one("Another input script is already being processed");
    if (infile != stdin) fclose(infile);
    infile = fopen(filename,"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",filename);
      error->one(str);
    }
    infiles[0] = infile;
  } else infile = NULL;

  file();
}

/* ----------------------------------------------------------------------
   parse the command in single and execute it
   return command name to caller
------------------------------------------------------------------------- */

char *Input::one(const char *single)
{
  strcpy(line,single);

  // echo the command unless scanning for label
  
  if (me == 0 && label_active == 0) {
    if (echo_screen && screen) fprintf(screen,"%s",line); 
    if (echo_log && logfile) fprintf(logfile,"%s",line);
  }

  // parse the line
  // if no command, just return NULL

  parse();
  if (command == NULL) return NULL;

  // if scanning for label, skip command unless it's a label command

  if (label_active && strcmp(command,"label") != 0) return NULL;

  // execute the command and return its name

  if (execute_command()) {
    char str[MAXLINE];
    sprintf(str,"Unknown command: %s",line);
    error->all(str);
  }

  return command;
}

/* ----------------------------------------------------------------------
   parse copy of command line
   strip comment = all chars from # on
   replace all $ via variable substitution
   command = first word
   narg = # of args
   arg[] = individual args
   treat text between double quotes as one arg
------------------------------------------------------------------------- */

void Input::parse()
{
  // make a copy to work on

  strcpy(copy,line);

  // strip any # comment by resetting string terminator
  // do not strip # inside double quotes

  int level = 0;
  char *ptr = copy;
  while (*ptr) {
    if (*ptr == '#' && level == 0) {
      *ptr = '\0';
      break;
    }
    if (*ptr == '"') {
      if (level == 0) level = 1;
      else level = 0;
    }
    ptr++;
  }

  // perform $ variable substitution (print changes)
  // except if searching for a label since earlier variable may not be defined

  if (!label_active) substitute(copy,1);

  // command = 1st arg

  command = strtok(copy," \t\n\r\f");
  if (command == NULL) return;

  // point arg[] at each subsequent arg
  // treat text between double quotes as one arg
  // insert string terminators in copy to delimit args

  narg = 0;
  while (1) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) memory->srealloc(arg,maxarg*sizeof(char *),"input:arg");
    }
    arg[narg] = strtok(NULL," \t\n\r\f");
    if (arg[narg] && arg[narg][0] == '\"') {
      arg[narg] = &arg[narg][1];
      if (arg[narg][strlen(arg[narg])-1] == '\"')
	arg[narg][strlen(arg[narg])-1] = '\0';
      else {
	arg[narg][strlen(arg[narg])] = ' ';
	ptr = strtok(arg[narg],"\"");
	if (ptr == NULL) error->all("Unbalanced quotes in input line");
      }
    }
    if (arg[narg]) narg++;
    else break;
  }
}

/* ----------------------------------------------------------------------
   substitute for $ variables in str
   print updated string if flag is set and not searching for label
------------------------------------------------------------------------- */

void Input::substitute(char *str, int flag)
{
  // use work[] as scratch space to expand str
  // do not replace $ inside double quotes as flagged by level
  // var = pts at variable name, ended by NULL
  //   if $ is followed by '{', trailing '}' becomes NULL
  //   else $x becomes x followed by NULL
  // beyond = pts at text following variable

  char *var,*value,*beyond;
  int level = 0;
  char *ptr = str;

  while (*ptr) {
    if (*ptr == '$' && level == 0) {
      if (*(ptr+1) == '{') {
	var = ptr+2;
	int i = 0;
	while (var[i] != '\0' && var[i] != '}') i++;
	if (var[i] == '\0') error->one("Invalid variable name");
	var[i] = '\0';
	beyond = ptr + strlen(var) + 3;
      } else {
	var = ptr;
	var[0] = var[1];
	var[1] = '\0';
	beyond = ptr + strlen(var) + 1;
      }
      value = variable->retrieve(var);
      if (value == NULL) error->one("Substitution for illegal variable");

      *ptr = '\0';
      strcpy(work,str);
      if (strlen(work)+strlen(value) >= MAXLINE)
	error->one("Input line too long after variable substitution");
      strcat(work,value);
      if (strlen(work)+strlen(beyond) >= MAXLINE)
	error->one("Input line too long after variable substitution");
      strcat(work,beyond);
      strcpy(str,work);
      ptr += strlen(value);
      if (flag && me == 0 && label_active == 0) {
	if (echo_screen && screen) fprintf(screen,"%s",str); 
	if (echo_log && logfile) fprintf(logfile,"%s",str);
      }
      continue;
    }
    if (*ptr == '"') {
      if (level == 0) level = 1;
      else level = 0;
    }
    ptr++;
  }
}

/* ----------------------------------------------------------------------
   process a single parsed command
   return 0 if successful, -1 if did not recognize command
------------------------------------------------------------------------- */

int Input::execute_command()
{
  int flag = 1;

  if (!strcmp(command,"clear")) clear();
  else if (!strcmp(command,"echo")) echo();
  else if (!strcmp(command,"if")) ifthenelse();
  else if (!strcmp(command,"include")) include();
  else if (!strcmp(command,"jump")) jump();
  else if (!strcmp(command,"label")) label();
  else if (!strcmp(command,"log")) log();
  else if (!strcmp(command,"next")) next_command();
  else if (!strcmp(command,"print")) print();
  else if (!strcmp(command,"variable")) variable_command();

  else if (!strcmp(command,"balance")) balance();
  else if (!strcmp(command,"bin")) bin();
  else if (!strcmp(command,"boundary")) boundary();
  else if (!strcmp(command,"check")) check();
  else if (!strcmp(command,"count")) count();
  else if (!strcmp(command,"debug")) debug();
  else if (!strcmp(command,"diffusion")) diffusion();
  else if (!strcmp(command,"dimension")) dimension();
  else if (!strcmp(command,"dump")) dump();
  else if (!strcmp(command,"dump_modify")) dump_modify();
  else if (!strcmp(command,"fix")) fix();
  else if (!strcmp(command,"global")) global();
  else if (!strcmp(command,"move_style")) move_style();
  else if (!strcmp(command,"particles")) particles();
  else if (!strcmp(command,"permeable")) permeable();
  else if (!strcmp(command,"probability")) probability();
  else if (!strcmp(command,"react_modify")) react_modify();
  else if (!strcmp(command,"reaction")) reaction();
  else if (!strcmp(command,"region")) region();
  else if (!strcmp(command,"restart")) restart();
  else if (!strcmp(command,"run_style")) run_style();
  else if (!strcmp(command,"seed")) seed();
  else if (!strcmp(command,"species")) species();
  else if (!strcmp(command,"stats")) stats();
  else if (!strcmp(command,"stats_modify")) stats_modify();
  else if (!strcmp(command,"timestep")) timestep();
  else if (!strcmp(command,"triangles")) triangles();
  else if (!strcmp(command,"undump")) undump();
  else if (!strcmp(command,"unfix")) unfix();
  else if (!strcmp(command,"unreact")) unreact();
  else if (!strcmp(command,"volume")) volume();
  
  else flag = 0;

  // return if command was listed above

  if (flag) return 0;

  // check if command is added via style.h

  if (0) return 0;      // dummy line to enable else-if macro expansion

#define CommandClass
#define CommandStyle(key,Class)         \
  else if (strcmp(command,#key) == 0) { \
    Class key;			\
    key.command(narg,arg);              \
    return 0;                           \
  }
#include "style.h"
#undef CommandClass

  // unrecognized command

  return -1;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Input::clear()
{
  if (narg > 0) error->all("Illegal clear command");
  sleep(2);              // pause to insure RNG does not start with same RN
  System::destroy();
  System::create();
}

/* ---------------------------------------------------------------------- */

void Input::echo()
{
  if (narg != 1) error->all("Illegal echo command");

  if (strcmp(arg[0],"none") == 0) {
    echo_screen = 0;
    echo_log = 0;
  } else if (strcmp(arg[0],"screen") == 0) {
    echo_screen = 1;
    echo_log = 0;
  } else if (strcmp(arg[0],"log") == 0) {
    echo_screen = 0;
    echo_log = 1;
  } else if (strcmp(arg[0],"both") == 0) {
    echo_screen = 1;
    echo_log = 1;
  } else error->all("Illegal echo command");
}

/* ---------------------------------------------------------------------- */

void Input::ifthenelse()
{
  if (narg != 5 && narg != 7) error->all("Illegal if command");

  int flag = 0;
  if (strcmp(arg[1],"==") == 0) {
    if (atof(arg[0]) == atof(arg[2])) flag = 1;
  } else if (strcmp(arg[1],"!=") == 0) {
    if (atof(arg[0]) != atof(arg[2])) flag = 1;
  } else if (strcmp(arg[1],"<") == 0) {
    if (atof(arg[0]) < atof(arg[2])) flag = 1;
  } else if (strcmp(arg[1],"<=") == 0) {
    if (atof(arg[0]) <= atof(arg[2])) flag = 1;
  } else if (strcmp(arg[1],">") == 0) {
    if (atof(arg[0]) > atof(arg[2])) flag = 1;
  } else if (strcmp(arg[1],">=") == 0) {
    if (atof(arg[0]) >= atof(arg[2])) flag = 1;
  } else error->all("Illegal if command");

  if (strcmp(arg[3],"then") != 0) error->all("Illegal if command");
  if (narg == 7 && strcmp(arg[5],"else") != 0) 
    error->all("Illegal if command");

  char str[128] = "\0";
  if (flag) strcpy(str,arg[4]);
  else if (narg == 7) strcpy(str,arg[6]);
  strcat(str,"\n");

  if (strlen(str) > 1) char *tmp = one(str);
}

/* ---------------------------------------------------------------------- */

void Input::include()
{
  if (narg != 1) error->all("Illegal include command");

  if (me == 0) {
    if (nfile == maxfile) {
      maxfile++;
      infiles = (FILE **) 
        memory->srealloc(infiles,maxfile*sizeof(FILE *),"input:infiles");
    }
    infile = fopen(arg[0],"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",arg[0]);
      error->one(str);
    }
    infiles[nfile++] = infile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::jump()
{
  if (narg < 1 || narg > 2) error->all("Illegal jump command");

  if (jump_skip) {
    jump_skip = 0;
    return;
  }

  if (me == 0) {
    if (infile != stdin) fclose(infile);
    infile = fopen(arg[0],"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",arg[0]);
      error->one(str);
    }
    infiles[nfile-1] = infile;
  }

  if (narg == 2) {
    label_active = 1;
    if (labelstr) delete [] labelstr;
    int n = strlen(arg[1]) + 1;
    labelstr = new char[n];
    strcpy(labelstr,arg[1]);
  }
}

/* ---------------------------------------------------------------------- */

void Input::label()
{
  if (narg != 1) error->all("Illegal label command");
  if (label_active && strcmp(labelstr,arg[0]) == 0) label_active = 0;
}

/* ---------------------------------------------------------------------- */

void Input::log()
{
  if (narg != 1) error->all("Illegal log command");

  if (me == 0) {
    if (logfile) fclose(logfile);
    if (strcmp(arg[0],"none") == 0) logfile = NULL;
    else {
      logfile = fopen(arg[0],"w");
      if (logfile == NULL) {
	char str[128];
	sprintf(str,"Cannot open logfile %s",arg[0]);
	error->one(str);
      }
    }
    if (universe->nworlds == 1) universe->ulogfile = logfile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::next_command()
{
  if (variable->next(narg,arg)) jump_skip = 1;
}

/* ---------------------------------------------------------------------- */

void Input::print()
{
  if (narg != 1) error->all("Illegal print command");

  // substitute for $ variables (no printing)

  substitute(arg[0],0);

  if (me == 0) {
    if (screen) fprintf(screen,"%s\n",arg[0]);
    if (logfile) fprintf(logfile,"%s\n",arg[0]);
  }
}

/* ---------------------------------------------------------------------- */

void Input::variable_command()
{
  variable->set(narg,arg);
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   one function for each ChemCell-specific input script command
------------------------------------------------------------------------- */

void Input::balance()
{
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use balance command with non-spatial simulation");
  grid->set_balance(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::bin()
{
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use bin command with non-spatial simulation");
  grid->create(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::boundary()
{
  domain->boundary(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::check()
{
  if (!simulator) error->all("Must set run_style first");
  if (!move) error->all("Run style does not support checking");
  move->set_check(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::count()
{
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 1)
    error->all("Cannot use count command with spatial simulation");
  particle->set_count(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::debug()
{
  if (!simulator) error->all("Must set run_style first");
  if (narg != 3) error->all("Illegal debug command");
  move->debug_proc = atoi(arg[0]);
  move->debug_step = atoi(arg[1]);
  move->debug_index = atoi(arg[2]);
}

/* ---------------------------------------------------------------------- */

void Input::diffusion()
{
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use diffusion command with non-spatial simulation");
  particle->set_diffusion(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dimension()
{
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use dimension command with non-spatial simulation");
  particle->set_dimension(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump()
{
  if (simulator->spatial_flag == 0)
    error->all("Cannot use dump command with non-spatial simulation");
  output->add_dump(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump_modify()
{
  output->dump_modify(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::fix()
{
  if (!simulator) error->all("Must set run_style first");
  modify->add_fix(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::global()
{
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use global command with non-spatial simulation");
  domain->global(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::move_style()
{
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use move_style command with non-spatial simulation");
  move->set_style(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::particles() {
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use particle command with non-spatial simulation");
  if (grid->setflag == 0) error->all("Must set bins before reading particles");
  particle->read(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::permeable()
{
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use permeable command with non-spatial simulation");
  move->set_permeable(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::probability()
{
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use probability command with non-spatial simulation");
  chem->set_prob(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::react_modify()
{
  if (!simulator) error->all("Must set run_style first");
  react->modify(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::reaction()
{
  if (!simulator) error->all("Must set run_style first");
  react->add(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::region() {
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use region command with non-spatial simulation");
  if (grid->setflag == 0) error->all("Must set bins before defining region");
  surf->add_region(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::restart()
{
  output->create_restart(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::run_style() {
  if (simulator) delete simulator;

  if (strcmp(arg[0],"none") == 0) error->all("Invalid run style");

#define SimulatorClass
#define SimulatorStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) simulator = new Class(narg,arg);
#include "style.h"
#undef SimulatorClass

  else error->all("Invalid run style");
}

/* ---------------------------------------------------------------------- */

void Input::seed()
{
  if (narg != 1) error->all("Illegal seed command");
  random->seed = atoi(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::species()
{
  if (!simulator) error->all("Must set run_style first");
  particle->set_species(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::stats()
{
  output->set_stats(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::stats_modify()
{
  output->stats_modify(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::timestep()
{
  if (!simulator) error->all("Must set run_style first");
  if (strcmp(simulator->style,"gillespie") == 0)
    error->all("Cannot use timestep command with gillespie simulation");
  if (narg != 1) error->all("Illegal timestep command");
  simulator->dt = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::triangles() {
  if (!simulator) error->all("Must set run_style first");
  if (simulator->spatial_flag == 0)
    error->all("Cannot use triangles command with non-spatial simulation");
  if (grid->setflag == 0) error->all("Must set bins before reading surface");
  surf->read_triangles(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::undump()
{
  if (narg != 1) error->all("Illegal undump command");
  output->delete_dump(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::unfix()
{
  if (!simulator) error->all("Must set run_style first");
  if (narg != 1) error->all("Illegal unfix command");
  modify->delete_fix(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::unreact()
{
  if (!simulator) error->all("Must set run_style first");
  if (narg != 1) error->all("Illegal unreact command");
  react->delete_reaction(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::volume()
{
  if (narg != 1) error->all("Illegal volume command");
  if (simulator->spatial_flag == 1)
    error->all("Cannot use volume command with spatial simulation");
  chem->volume = atof(arg[0]);
}
