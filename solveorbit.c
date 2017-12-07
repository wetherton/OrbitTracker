#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gslmatrix.h>
#include <gsl/gsl_odeiv2.h>
#include "solveorbit.h"

#define USEJAC 0;
#define DIM 6;
#define TILEX 20;
#define TILEZ 20;

int dxdt(double t, const double y[], double dydt[], void *params){
  fieldgrid masterfield = (masterfield *) params;
  field fields = interpfields(y[0],y[3],masterfield);
  (void) (t);
  dydt[0] = y[3];
  dydt[1] = y[4];
  dydt[2] = y[5];
  dydt[3] = -fields.Ex + y[5]*fields.By - y[4]*fields.Bz;
  dydt[4] = -fields.Ey + y[3]*fields.Bz - y[5]*fields.Bx;
  dydt[5] = -fields.Ez + y[4]*fields.Bx - y[3]*fields.By;
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
  int nts = 10000;
  char filestring[];
  filestring = sprintf(filestring, "x%04.0fz%04.0fvx%03.0fvy%03.0fvz%03.0f.txt",IC.x,IC.z,IC.vx,IC.vy,IC.vz);
  FILE *fp = fopen(filestring,"w");
  fieldgrid masterfield = readfields(NX,NZ);
  gsl_odeiv2_system sys = {dxdt, jacobian, DIM, (void *) masterfield};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,1e-6, 1e-6, 0.0);
  double t = 0.0, tend = 10000.0;
  double y[6] = {IC.x, IC.y, IC.z, IC.vx, IC.vy, IC.vz};
  for(i = 1; i<nts; i++){
    double ti = i*tend/nts;
    int status = gsl_odeiv2_driver_apply(d,&t,ti,y);

    if (status!= GSL_SUCCESS){
      fprintf(fp,"Error, return value = %d\n",status);
      break;
    }

    fprintf("%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",y[0],y[1],y[2],y[3],y[4],y[5],ti);
    
    if(outcheck(y[0],y[2])){
      break;
    }
  }
  gsl_odeiv2_driver_free(d);
}
