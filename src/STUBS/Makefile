# Makefile for MPI stubs - edit this for your platform

SHELL = /bin/sh
.IGNORE:

# Files

SRC =		mpi.c zoltan.c
INC =		mpi.h zoltan.h

# Definitions

EXE =		libstub.a
OBJ = 		$(SRC:.c=.o)

# System-specific settings

CC =		g++
CCFLAGS =	-g
ARCHIVE =	ar
ARCHFLAG =	r

# Target

$(EXE):	$(OBJ)
	$(ARCHIVE) $(ARCHFLAG) $(EXE) $(OBJ)
	ranlib $(EXE)

# Clean

clean:
	rm *.o libstub.a

# Compilation rules

.cpp.o:
	$(CC) $(CCFLAGS) -c $<

.c.o:
	$(CC) $(CCFLAGS) -c $<

# Individual dependencies

$(OBJ):	$(INC)

