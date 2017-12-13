# Warnings
WFLAGS	:= -Wall -Wextra

# Optimization and architecture
OPT		:= -O3

# Language standard
CCSTD	:= -std=c99
CXXSTD	:= -std=c++14

.DEFAULT_GOAL := all

EXEC := orbitsolve

all : Makefile $(EXEC)

orbitsolve: solveorbit.c gdareader.c interp2.c
	mpicc $(CXXSTD) $(WFLAGS) $(OPT) -fopenmp -o $@ $<

.PHONY: clean
clean:
	@ rm -f $(EXEC) *.o
