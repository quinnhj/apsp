
CC = g++
#OPENMP = -fopenmp
OPENMP =
CFLAGS = -O3
LIBS =


TARGETS = apsp

all:	$(TARGETS)

apsp: apsp.o common.o 
	$(CC) -o $@ $(LIBS) $(OPENMP) apsp.o common.o

apsp.o: apsp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) apsp.cpp

common.o: common.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
