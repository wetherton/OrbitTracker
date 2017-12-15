# Warnings
WFLAGS	:= -Wall -Wextra

# Optimization and architecture
OPT		:= -O2

# Linker Options
LDOPT   := $(OPT)
LDFLAGS := -fopenmp -lgsl -lgslcblas -lm


# Language standard
CCSTD	:= -std=c99
CXXSTD	:= -std=c++14

.DEFAULT_GOAL := all

EXEC := orbitsolve

all : Makefile $(EXEC)

debug:
	gcc -g -std=c99 -o orbitsolve solveorbit_noMPI.c gdareader.c interp2.c -fopenmp -lgsl -lgslcblas

orbitsolve: 
	mpicc -std=c99 -o orbitsolve -O3 solveorbit.c gdareader.c interp2.c -fopenmp -lgsl -lgslcblas

noMPI:
	gcc -O3 -std=c99 -o orbitsolve solveorbit_noMPI.c gdareader.c interp2.c -fopenmp -lgsl -lgslcblas

.PHONY: clean
clean:
	@ rm -f $(EXEC) *.o
