#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gslmatrix.h>
#include <gsl/gsl_odeiv2.h>

#define USEJAC 0;
#define DIM 6;
#define TILEX 20;
#define TILEZ 20;

typedef struct{
  double Ex, Ey, Ez, Bx, By, Bz;
} field;


typedef struct{
  double x, y, z, vx, vy, vz;
}posvel;

int dxdt(double t, const double y[], double dydt[], void *params){
  field masterfield[NX][NZ] = (field *) params;
  field fields = interpfields(y[0],y[3],tilefield);
  (void) (t);
  dydt[0] = y[3];
  dydt[1] = y[4];
  dydt[2] = y[5];
  dydt[3] = -Ex + y[5]*fields.By - y[4]*fields.Bz;
  dydt[4] = -Ey + y[3]*fields.Bz - y[5]*fields.Bx;
  dydt[5] = -Ez + y[4]*fields.Bx - y[3]*fields.By;
  return GSL_SUCCESS;
}

int jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  (void) (t);
  (void) (params);
  (void) (y);
  dfdy = NULL;
  for(int i = 0; i< DIM; i++){
    dfdt[i] = 0.0;
  }
  
  return GSL_SUCCESS;
}

int solveorbit(posvel IC){
  field masterfield[NX][NZ] = readfields(NX,NZ);
  gsl_odeiv2_system sys = {dxdt, jacobian, DIM, (void *) masterfield}
}
