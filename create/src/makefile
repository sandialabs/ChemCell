CCFILES = Surface.cpp Distribution.cpp random.cpp system.cpp control.cpp Particle.cpp geometry.cpp
HFILES = $(CCFILES:.cpp=.h)
CCOBJS = $(CCFILES:.cpp=.o)
CCFLAGS = -g -LANG:std
LIBS = -lm
CC = g++

main:	main.o $(CCOBJS)
	g++ $(CCFLAGS) -o ../run/chemcreate main.o $(CCOBJS) $(LIBS)

main.o: main.cpp
system.o: system.cpp $(HFILES)
Surface.o: Surface.cpp $(HFILES)
Distribution.o: Distribution.cpp $(HFILES)
random.o: random.cpp  $(HFILES)
particle.o: Particle.cpp  $(HFILES)
control.o: control.cpp $(HFILES)
geometry.o: geometry.cpp $(HFILES)



.cpp.o:
	g++ $(CCFLAGS) -c -I. $*.cpp

clean:
	rm -rf core *.o ../run/chemcreate ii_files

clobber: clean
	rm -f ../run/chemcreate
