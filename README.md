# OrbitTracker
Final Project for CS759. Tracks particle orbits in parallel. Inputs based on VPIC simulations.
The Makefile will make the executable if you have loaded OpenMPI. 
sample of usage:
mpirun -np 4 ./orbitsolve
Sample files fit to orbit.cfg in repo, so it will run. If you want to speed it up, Nvx, Nvy, and Nvz may be decreased. 
fulltrack is currently 1, so the traces will be output to their own files. Set to zero for the mapping alone. 
omp:nthreads can also be adjusted in orbit.cfg. 
