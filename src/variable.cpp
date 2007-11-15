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
#include "stdlib.h"
#include "string.h"
#include "unistd.h"
#include "variable.h"
#include "universe.h"
#include "simulator.h"
#include "particle.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

#define VARDELTA 4

enum{INDEX,LOOP,EQUAL,WORLD,UNIVERSE,ULOOP};

/* ---------------------------------------------------------------------- */

Variable::Variable()
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  nvar = maxvar = 0;
  names = NULL;
  style = NULL;
  num = NULL;
  index = NULL;
  data = NULL;
}

/* ---------------------------------------------------------------------- */

Variable::~Variable()
{
  for (int i = 0; i < nvar; i++) {
    delete [] names[i];
    for (int j = 0; j < num[i]; j++) delete [] data[i][j];
    delete [] data[i];
  }
  memory->sfree(names);
  memory->sfree(style);
  memory->sfree(num);
  memory->sfree(index);
  memory->sfree(data);
}

/* ----------------------------------------------------------------------
   called by variable command in input script
------------------------------------------------------------------------- */

void Variable::set(int narg, char **arg)
{
  if (narg < 3) error->all("Illegal variable command");

  // if var already exists, just skip (except EQUAL vars)

  if (find(arg[0]) >= 0 && strcmp(arg[1],"equal") != 0) return;
      
  // make space for new variable

  if (nvar == maxvar) {
    maxvar += VARDELTA;
    names = (char **)
      memory->srealloc(names,maxvar*sizeof(char *),"var:names");
    style = (int *) memory->srealloc(style,maxvar*sizeof(int),"var:style");
    num = (int *) memory->srealloc(num,maxvar*sizeof(int),"var:num");
    index = (int *) memory->srealloc(index,maxvar*sizeof(int),"var:index");
    data = (char ***) 
      memory->srealloc(data,maxvar*sizeof(char **),"var:data");
  }

  // INDEX
  // num = listed args, index = 1st value, data = copied args

  if (strcmp(arg[1],"index") == 0) {
    style[nvar] = INDEX;
    num[nvar] = narg - 2;
    index[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

  // LOOP
  // num = N, index = 1st value, data = list of NULLS since never used

  } else if (strcmp(arg[1],"loop") == 0) {
    style[nvar] = LOOP;
    num[nvar] = atoi(arg[2]);
    index[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    for (int i = 0; i < num[nvar]; i++) data[nvar][i] = NULL;
    
  // EQUAL
  // remove pre-existing var if also style EQUAL (allows it to be reset)
  // num = 2, index = 1st value
  // data = 2 values, 1st is string to eval, 2nd is filled on retrieval

  } else if (strcmp(arg[1],"equal") == 0) {
    if (find(arg[0]) >= 0) {
      if (style[find(arg[0])] != EQUAL)
	error->all("Cannot redefine variable as a different style");
      remove(find(arg[0]));
    }
    style[nvar] = EQUAL;
    num[nvar] = 2;
    index[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(1,&arg[2],data[nvar]);
    data[nvar][1] = NULL;
    
  // WORLD
  // num = listed args, index = partition this proc is in, data = copied args
  // error check that num = # of worlds in universe

  } else if (strcmp(arg[1],"world") == 0) {
    style[nvar] = WORLD;
    num[nvar] = narg - 2;
    if (num[nvar] != universe->nworlds)
      error->all("World variable count doesn't match # of partitions");
    index[nvar] = universe->iworld;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

  // UNIVERSE and ULOOP
  // for UNIVERSE: num = listed args, data = copied args
  // for ULOOP: num = N, data = list of NULLS since never used
  // index = partition this proc is in
  // universe proc 0 creates lock file
  // error check that all other universe/uloop variables are same length

  } else if (strcmp(arg[1],"universe") == 0 || strcmp(arg[1],"uloop") == 0) {
    if (strcmp(arg[1],"universe") == 0) {
      style[nvar] = UNIVERSE;
      num[nvar] = narg - 2;
      data[nvar] = new char*[num[nvar]];
      copy(num[nvar],&arg[2],data[nvar]);
    } else {
      style[nvar] = ULOOP;
      num[nvar] = atoi(arg[2]);
      data[nvar] = new char*[num[nvar]];
      for (int i = 0; i < num[nvar]; i++) data[nvar][i] = NULL;
    }

    if (num[nvar] < universe->nworlds)
      error->all("Universe/uloop variable count < # of partitions");
    index[nvar] = universe->iworld;

    if (universe->me == 0) {
      FILE *fp = fopen("tmp.ccell.variable","w");
      fprintf(fp,"%d\n",universe->nworlds);
      fclose(fp);
    }

    for (int jvar = 0; jvar < nvar; jvar++)
      if (num[jvar] && (style[jvar] == UNIVERSE || style[jvar] == ULOOP) && 
	  num[nvar] != num[jvar])
	error->all("All universe/uloop variables must have same # of values");

    if (me == 0) {
      if (universe->uscreen)
	fprintf(universe->uscreen,
		"Initial ${%s} setting: value %d on partition %d\n",
		arg[0],index[nvar]+1,universe->iworld);
      if (universe->ulogfile)
	fprintf(universe->ulogfile,
		"Initial ${%s} setting: value %d on partition %d\n",
		arg[0],index[nvar]+1,universe->iworld);
    }

  } else error->all("Illegal variable command");

  // set variable name (after possible EQUAL remove)

  names[nvar] = new char[strlen(arg[0]) + 1];
  strcpy(names[nvar],arg[0]);
  nvar++;
}

/* ----------------------------------------------------------------------
   single-value INDEX variable created by command-line argument
------------------------------------------------------------------------- */

void Variable::set(char *name, char *value)
{
  int ivar = find(name);
  if (ivar >= 0) error->all("Command-line variable already exists");

  if (nvar == maxvar) {
    maxvar += VARDELTA;
    names = (char **)
      memory->srealloc(names,maxvar*sizeof(char *),"var:names");
    style = (int *) memory->srealloc(style,maxvar*sizeof(int),"var:style");
    num = (int *) memory->srealloc(num,maxvar*sizeof(int),"var:num");
    index = (int *) memory->srealloc(index,maxvar*sizeof(int),"var:index");
    data = (char ***) 
      memory->srealloc(data,maxvar*sizeof(char **),"var:data");
  }

  names[nvar] = new char[strlen(name) + 1];
  strcpy(names[nvar],name);
  style[nvar] = INDEX;
  num[nvar] = 1;
  index[nvar] = 0;
  data[nvar] = new char*[num[nvar]];
  copy(1,&value,data[nvar]);
  nvar++;
}

/* ----------------------------------------------------------------------
   increment variable(s)
   return 0 if OK if successfully incremented
   return 1 if any variable is exhausted, free the variable to allow re-use
------------------------------------------------------------------------- */

int Variable::next(int narg, char **arg)
{
  int ivar;

  if (narg == 0) error->all("Illegal next command");

  // check that variables exist and are all the same style
  // exception: UNIVERSE and ULOOP variables can be mixed in same next command

  for (int iarg = 0; iarg < narg; iarg++) {
    ivar = find(arg[iarg]);
    if (ivar == -1) error->all("Invalid variable in next command");
    if (style[ivar] == ULOOP && style[find(arg[0])] == UNIVERSE) continue;
    else if (style[ivar] == UNIVERSE && style[find(arg[0])] == ULOOP) continue;
    else if (style[ivar] != style[find(arg[0])])
      error->all("All variables in next command must be same style");
  }

  // check for invalid styles EQUAL or WORLD

  int istyle = style[find(arg[0])];
  if (istyle == EQUAL || istyle == WORLD)
    error->all("Invalid variable style with next command");

  // increment all variables in list
  // if any variable is exhausted, set flag = 1 and remove var to allow re-use

  int flag = 0;

  if (istyle == INDEX || istyle == LOOP) {
    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      index[ivar]++;
      if (index[ivar] >= num[ivar]) {
	flag = 1;
	remove(ivar);
      }
    }

  } else if (istyle == UNIVERSE || istyle == ULOOP) {

    // wait until lock file can be created and owned by proc 0 of this world
    // read next available index and Bcast it within my world
    // set all variables in list to nextindex

    int nextindex;
    if (me == 0) {
      while (1) {
	if (!rename("tmp.ccell.variable","tmp.ccell.variable.lock")) break;
	usleep(100000);
      }
      FILE *fp = fopen("tmp.ccell.variable.lock","r");
      fscanf(fp,"%d",&nextindex);
      fclose(fp);
      fp = fopen("tmp.ccell.variable.lock","w");
      fprintf(fp,"%d\n",nextindex+1);
      fclose(fp);
      rename("tmp.ccell.variable.lock","tmp.ccell.variable");
      if (universe->uscreen)
	fprintf(universe->uscreen,
		"Increment via next: value %d on partition %d\n",
		nextindex+1,universe->iworld);
      if (universe->ulogfile)
      	fprintf(universe->ulogfile,
      		"Increment via next: value %d on partition %d\n",
      		nextindex+1,universe->iworld);
    }
    MPI_Bcast(&nextindex,1,MPI_INT,0,world);

    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      index[ivar] = nextindex;
      if (index[ivar] >= num[ivar]) {
	flag = 1;
	remove(ivar);
      }
    }
  }

  return flag;
}

/* ----------------------------------------------------------------------
   return ptr to the data text associated with a variable
   return NULL if no variable or index is bad, caller must respond
------------------------------------------------------------------------- */

char *Variable::retrieve(char *name)
{
  int ivar = find(name);
  if (ivar == -1) return NULL;
  if (index[ivar] >= num[ivar]) return NULL;

  char *str;
  if (style[ivar] == INDEX || style[ivar] == WORLD || 
      style[ivar] == UNIVERSE) {
    str = data[ivar][index[ivar]];
  } else if (style[ivar] == LOOP || style[ivar] == ULOOP) {
    char *value = new char[16];
    sprintf(value,"%d",index[ivar]+1);
    int n = strlen(value) + 1;
    if (data[ivar][0]) delete [] data[ivar][0];
    data[ivar][0] = new char[n];
    strcpy(data[ivar][0],value);
    delete [] value;
    str = data[ivar][0];
  } else if (style[ivar] == EQUAL) {
    char *value = evaluate(data[ivar][0]);
    int n = strlen(value) + 1;
    if (data[ivar][1]) delete [] data[ivar][1];
    data[ivar][1] = new char[n];
    strcpy(data[ivar][1],value);
    delete [] value;
    str = data[ivar][1];
  }
  return str;
}

/* ----------------------------------------------------------------------
   search for name in list of variables names
   return index or -1 if not found
------------------------------------------------------------------------- */
  
int Variable::find(char *name)
{
  for (int i = 0; i < nvar; i++)
    if (strcmp(name,names[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   remove Nth variable from list and compact list
------------------------------------------------------------------------- */
  
void Variable::remove(int n)
{
  delete [] names[n];
  for (int i = 0; i < num[n]; i++) delete [] data[n][i];
  delete [] data[n];

  for (int i = n+1; i < nvar; i++) {
    names[i-1] = names[i];
    style[i-1] = style[i];
    num[i-1] = num[i];
    index[i-1] = index[i];
    data[i-1] = data[i];
  }
  nvar--;
}

/* ----------------------------------------------------------------------
   copy narg strings from **from to **to 
------------------------------------------------------------------------- */
  
void Variable::copy(int narg, char **from, char **to)
{
  int n;

  for (int i = 0; i < narg; i++) {
    n = strlen(from[i]) + 1;
    to[i] = new char[n];
    strcpy(to[i],from[i]);
  }
}

/* ----------------------------------------------------------------------
   recursive method to evaluate a string
   string can be a number: 0.0, -5.45, etc
   string can be a keyword: lx, ly, lz, vol, etc
   string can be a vector: x[123], y[3], vx[34], etc
     value inside brackets must be global ID from 1 to Natoms
   string can be a function: div(x,y), mult(x,y), add(x,y), neg(x), etc
     each function arg can be a number, keyword, vector, or function
   see lists of valid keywords, vectors, and functions below
   when string is evaluated, put result in a newly allocated string
   return address of result string (will be freed by caller)
------------------------------------------------------------------------- */

char *Variable::evaluate(char *str)
{
  // allocate a new string for the eventual result

  char *result = new char[32];
  double answer;

  // if string is "function", grab one or two args, evaulate args recursively

  if (strchr(str,'(')) {
    if (str[strlen(str)-1] != ')')
      error->all("Cannot evaluate variable equal command");

    char *ptr = strchr(str,'(');
    int n = ptr - str;
    char *func = new char[n+1];
    strncpy(func,str,n);
    func[n] = '\0';

    char *comma = ++ptr;
    int level = 0;
    while (1) {
      if (*comma == '\0')
	error->all("Cannot evaluate variable equal command");
      else if (*comma == ',' && level == 0) break;
      else if (*comma == ')' && level == 0) break;
      else if (*comma == '(') level++;
      else if (*comma == ')') level--;
      comma++;
    }

    char *arg1,*arg2;
    n = comma - ptr;
    arg1 = new char[n+1];
    strncpy(arg1,ptr,n);
    arg1[n] = '\0';

    if (*comma == ',') {
      ptr = comma + 1;
      comma = &str[strlen(str)-1];
      n = comma - ptr;
      arg2 = new char[n+1];
      strncpy(arg2,ptr,n);
      arg2[n] = '\0';
    } else arg2 = NULL;
    
    double value1,value2;
    char *strarg1 = NULL;
    char *strarg2 = NULL;

    // customize by adding function to this list and to if statement
    // math functions: add(x,y),sub(x,y),mult(x,y),div(x,y),
    //                 neg(x),pow(x,y),exp(x),ln(x),sqrt(x)

    if (strcmp(func,"add") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      strarg2 = evaluate(arg2);
      value1 = atof(strarg1);
      value2 = atof(strarg2);
      answer = value1 + value2;

    } else if (strcmp(func,"sub") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      strarg2 = evaluate(arg2);
      value1 = atof(strarg1);
      value2 = atof(strarg2);
      answer = value1 - value2;

    } else if (strcmp(func,"mult") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      strarg2 = evaluate(arg2);
      value1 = atof(strarg1);
      value2 = atof(strarg2);
      answer = value1 * value2;

    } else if (strcmp(func,"div") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      strarg2 = evaluate(arg2);
      value1 = atof(strarg1);
      value2 = atof(strarg2);
      if (value2 == 0.0)
	error->all("Cannot evaluate variable equal command");
      answer = value1 / value2;

    } else if (strcmp(func,"neg") == 0) {
      if (arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      value1 = atof(strarg1);
      answer = -value1;

    } else if (strcmp(func,"pow") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      strarg2 = evaluate(arg2);
      value1 = atof(strarg1);
      value2 = atof(strarg2);
      if (value2 == 0.0)
	error->all("Cannot evaluate variable equal command");
      answer = pow(value1,value2);

    } else if (strcmp(func,"exp") == 0) {
      if (arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      value1 = atof(strarg1);
      answer = exp(value1);

    } else if (strcmp(func,"ln") == 0) {
      if (arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      value1 = atof(strarg1);
      if (value1 == 0.0)
	error->all("Cannot evaluate variable equal command");
      answer = log(value1);

    } else if (strcmp(func,"sqrt") == 0) {
      if (arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      value1 = atof(strarg1);
      if (value1 == 0.0)
	error->all("Cannot evaluate variable equal command");
      answer = sqrt(value1);

    } else error->all("Cannot evaluate variable equal command");

    delete [] func;
    delete [] arg1;
    delete [] strarg1;
    if (arg2) {
      delete [] arg2;
      delete [] strarg2;
    }
    sprintf(result,"%.10g",answer);

  // if string is "vector", find which proc owns atom, grab vector value

  } else if (strchr(str,'[')) {
    if (str[strlen(str)-1] != ']')
      error->all("Cannot evaluate variable equal command");
    if (nprocs > 1) error->all("Cannot use vectors in parallel");

    char *ptr = strchr(str,'[');
    int n = ptr - str;
    char *vector = new char[n+1];
    strncpy(vector,str,n);
    vector[n] = '\0';

    char *arg;
    ptr++;
    char *ptr2 = &str[strlen(str)-1];
    n = ptr2 - ptr;
    arg = new char[n+1];
    strncpy(arg,ptr,n);
    arg[n] = '\0';
    int i = atoi(arg);

    double mine;

    // customize by adding vector to this list and to if statement
    // x,y,z,sp,bin,tri,seed

    if (strcmp(vector,"x") == 0) mine = particle->plist[i].x[0];
    else if (strcmp(vector,"y") == 0) mine = particle->plist[i].x[1];
    else if (strcmp(vector,"z") == 0) mine = particle->plist[i].x[2];
    else if (strcmp(vector,"sp") == 0) mine = particle->plist[i].species;
    else if (strcmp(vector,"bin") == 0) mine = particle->plist[i].ibin;
    else if (strcmp(vector,"tri") == 0) mine = particle->plist[i].itri;
    else if (strcmp(vector,"seed") == 0) mine = particle->plist[i].seed;
    
    else error->one("Invalid vector in variable equal command");

    delete [] vector;
    delete [] arg;
    sprintf(result,"%.10g",answer);

  // if string is "keyword", compute appropriate value

  } else if (str[0] - 'a' >= 0 && str[0] - 'a' < 26) {

    // customize by adding keyword to this list and to if statement
    // step,npart,lx,ly,lz,vol

    if (!simulator)
      error->all("Using variable equal keyword before sim box is defined");

    // process these keywords here
    // they are allowed anytime after simulation box is defined

    if (strcmp(str,"step") == 0)
      answer = simulator->ntimestep;
    else if (strcmp(str,"npart") == 0)
      answer = particle->nlocal;
    else if (strcmp(str,"lx") == 0)
      answer = domain->xsize;
    else if (strcmp(str,"ly") == 0)
      answer = domain->ysize;
    else if (strcmp(str,"lz") == 0)
      answer = domain->zsize;
    else if (strcmp(str,"vol") == 0)
      answer = domain->xsize * domain->ysize * domain->zsize;

    else error->all("Unknown keyword for variable equal command");
    
    sprintf(result,"%.10g",answer);

  // string is a number, just copy to result

  } else strcpy(result,str);

  // return newly allocated string

  return result;
}
