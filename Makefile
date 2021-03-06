# 
# This relies on installation of the Boost library. In particular, we use the
# Fibonacci Heap data structure. Set BOOST_ROOT to point to the root of your
# Boost library install as an environment variable, otherwise it won't compile.
#

CC = g++-4.9
OPENMP = -fopenmp
#OPENMP =
CFLAGS = -O3
LIBS =
GPROF = -pg
CPP_FLAG = -std=c++11

TARGETS = apsp serial

all:	$(TARGETS)

apsp: apsp.o common.o 
#	$(CC) $(GPROF) -I$(BOOST_ROOT) -o $@ $(LIBS) $(OPENMP) apsp.o common.o
	$(CC) $(GPROF) $(CPP_FLAG) -o $@ $(LIBS) $(OPENMP) apsp.o common.o

serial: serial.o common.o
	$(CC) $(GPROF) $(CPP_FLAG) -o $@ $(LIBS) $(OPENMP) serial.o common.o

serial.o: serial.cpp common.h
#	$(CC) $(GPROF) -I$(BOOST_ROOT) -c $(OPENMP) $(CFLAGS) apsp.cpp
	$(CC) $(GPROF) $(CPP_FLAG) -c $(OPENMP) $(CFLAGS) serial.cpp

apsp.o: apsp.cpp common.h
#	$(CC) $(GPROF) -I$(BOOST_ROOT) -c $(OPENMP) $(CFLAGS) apsp.cpp
	$(CC) $(GPROF) $(CPP_FLAG) -c $(OPENMP) $(CFLAGS) apsp.cpp

common.o: common.cpp common.h
#	$(CC) $(GPROF) -I$(BOOST_ROOT) -c $(OPENMP) $(CFLAGS) common.cpp
	$(CC) $(GPROF) $(CPP_FLAG) -c $(OPENMP) $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS)
