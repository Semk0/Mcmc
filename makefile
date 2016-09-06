CFLAGS = -c -Wall

all: hello

hello: mcmc.o parameters.o datastr.o spec.o
	g++ -o mcmc mcmc.o parameters.o datastr.o spec.o
	subl out/mcmc.txt
	./mcmc

mcmc.o: mcmc.cpp
	g++ $(CFLAGS) mcmc.cpp

parameters.o: parameters.cpp
	g++ $(CFLAGS) parameters.cpp

datastr.o: datastr.cpp
	g++ $(CFLAGS) datastr.cpp

spec.o: spec.cpp
	g++ $(CFLAGS) spec.cpp