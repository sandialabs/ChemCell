# ChemCell multiple-machine Makefile

SHELL = /bin/sh

# Definitions

ROOT =	ccell
EXE =	$(ROOT)_$@
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# Targets

help:
	@echo 'Type "make target" where target is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	for file in $$files; do head -1 $$file; done

clean:
	rm -r Obj_*

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp -p $(SRC) $(INC) Obj_$@
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE)  "OBJ = $(OBJ)" "INC = $(INC)" "EXE = ../$(EXE)" ../$(EXE)
	@if [ -d Obj_$@ ]; then cd Obj_$@; rm $(SRC) $(INC) Makefile*; fi
