#ifndef SOLVEORBIT_H
#define SOLVEORBIT_H

typedef struct{
  double Ex, Ey, Ez, Bx, By, Bz;
} field; //6 doubles to categorize E&M fields at a point

typedef struct{
  double **Ex, **Ey, **Ez, **Bx, **By, **Bz;
  double *x, *z;
} fieldgrid; //2D arrays for E&M fields, 1D arrays to specify the grid


typedef struct{
  double x, y, z, vx, vy, vz;
}posvel; //Simple structure to hold phase space coordinate

typedef struct{
  double x, z;
}pos; //Interpolation position

int dxdt(double t, const double y[],double dydt[],void * params);
int dxdtB(double t, const double y[],double dydt[],void * params);
int jacobian(double t, const double y[],double *, double dfdt[], void * params);
posvel solveorbit(posvel IC, fieldgrid MF);
posvel solveorbitB(posvel IC, fieldgrid MF);
int initialize();
int outcheck(double, double);

#endif
