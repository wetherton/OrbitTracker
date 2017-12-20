#ifndef SOLVEORBIT_H
#define SOLVEORBIT_H

typedef struct{
  double Ex, Ey, Ez, Bx, By, Bz;
} field;

typedef struct{
  double **Ex, **Ey, **Ez, **Bx, **By, **Bz;
  double *x, *z;
} fieldgrid;


typedef struct{
  double x, y, z, vx, vy, vz;
}posvel;

typedef struct{
  double x, z;
}pos;

int dxdt(double t, const double y[],double dydt[],void * params);
int dxdtB(double t, const double y[],double dydt[],void * params);
int jacobian(double t, const double y[],double *, double dfdt[], void * params);
posvel solveorbit(posvel IC, fieldgrid MF);
posvel solveorbitB(posvel IC, fieldgrid MF);
int initialize();
int outcheck(double, double);

#endif
