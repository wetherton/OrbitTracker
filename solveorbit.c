#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gslmatrix.h>
#include <gsl/gsl_odeiv2.h>
#include "solveorbit.h"
#include "gdareader.h"

#define USEJAC 0;
#define DIM 6;
#define TILEX 20;
#define TILEZ 20;
#define OMPNUM 8;
#define MAGIC 4;

char jobname[128], outdir[256],gdadir[256];
int Nvx, Nvy, Nvz;
double vxmax, vymax, vzmax, Lx, Lz, minx, minz, maxx, maxz;
int Npoints, nts, nx, nz;
double tend;
double *x0, *z0;

int main(int argc,char *argv[]){
  MPI_Init(&argc, &argv);
  int myid, numprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  initialize();
  double *vx, *vy, *vz;
  int nICs = Nvx*Nvy*Nvz;
  int size = nICs*sizeof(double);
  vx = (double *) malloc(size);
  vy = (double *) malloc(size);
  vz = (double *) malloc(size);
  if(Nvx == 1) dvx =  0;
  else dvx = 2*vxmax/(Nvx-1);
  if(Nvy == 1) dvy =  0;
  else dvy = 2*vymax/(Nvy-1);
  if(Nvz == 1) dvz =  0;
  else dvz = 2*vzmax/(Nvz-1);
  
  for(int i = 0; i< Nvx; i++){
    for(int j = 0; j< Nvy; j++){
      for(int k = 0; k<Nvz; k++){
	int ind = k+Nvz*j+Nvz*Nvy*i;
	vx[ind] = vxmax-i*dvx;
	vy[ind] = vymax-j*dvy;
	vz[ind] = vzmax-k*dvz;
      }
    }
  }
  for(int npos = 0; npos<Npoints; npos++){
    for (int i = rank; i<nICs; i+=numprocs*OMPNUM*MAGIC){
#pragma omp_parallel_for()
      for(int j =0; j<OMPNUM*MAGIC; j++){
	int index = j + i;
	if (index>nICs) break;
	posvel IC;
	IC.x = x0[npos];
	IC.y = 0.0;
	IC.z = z0[npos];
	IC.vx = vx[index];
	IC.vy = vy[index];
	IC.vz = vz[index];
	solveorbit(IC);
      }
    }
  }
  MPI_Finalize();
}

int dxdt(double t, const double y[], double dydt[], void *params){
  fieldgrid masterfield = (masterfield *) params;
  field fields = interpfields(y[0],y[3],masterfield);
  (void) (t);
  double gamma = sqrt(1+y[3]*y[3]+y[4]*y[4]+y[5]*y[5]);
  double vx = y[3]/gamma;
  double vy = y[4]/gamma;
  double vz = y[5]/gamma;
  dydt[0] = vx;
  dydt[1] = vy;
  dydt[2] = vz;
  dydt[3] = -fields.Ex + vz*fields.By - vy*fields.Bz;
  dydt[4] = -fields.Ey + vx*fields.Bz - vz*fields.Bx;
  dydt[5] = -fields.Ez + vy*fields.Bx - vx*fields.By;
  return GSL_SUCCESS;
}

int dxdtB(double t, const double y[], double dydt[], void *params){
  fieldgrid masterfield = (masterfield *) params;
  field fields = interpfields(y[0],y[3],masterfield);
  double gamma = sqrt(1+y[3]*y[3]+y[4]*y[4]+y[5]*y[5]);
  double vx = y[3]/gamma;
  double vy = y[4]/gamma;
  double vz = y[5]/gamma;
  (void) (t);
  dydt[0] = vx;
  dydt[1] = vy;
  dydt[2] = vz;
  dydt[3] = -fields.Ex - vz*fields.By + vy*fields.Bz;
  dydt[4] = -fields.Ey - vx*fields.Bz + vz*fields.Bx;
  dydt[5] = -fields.Ez - vy*fields.Bx + vx*fields.By;
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
  char filestring[256];
  filestring = sprintf(filestring, "%sFx%04.0fz%04.0fvx%03.0fvy%03.0fvz%03.0f.txt",jobname,IC.x,IC.z,IC.vx,IC.vy,IC.vz);
  FILE *fp = fopen(filestring,"w");
  fieldgrid masterfield = readfields(NX,NZ);
  gsl_odeiv2_system sys = {dxdt, jacobian, DIM, (void *) masterfield};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,1e-6, 1e-6, 0.0);
  double t = 0.0;
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

int solveorbitB(posvel IC){
  char filestring[256];
  filestring = sprintf(filestring, "%sBx%04.0fz%04.0fvx%03.0fvy%03.0fvz%03.0f.txt",jobname,IC.x,IC.z,IC.vx,IC.vy,IC.vz);
  FILE *fp = fopen(filestring,"w");
  fieldgrid masterfield = readfields(NX,NZ);
  gsl_odeiv2_system sys = {dxdtB, jacobian, DIM, (void *) masterfield};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,1e-6, 1e-6, 0.0);
  double t = 0.0;
  double y[6] = {IC.x, IC.y, IC.z, -IC.vx, -IC.vy, -IC.vz};
  for(i = 1; i<nts; i++){
    double ti = i*tend/nts;
    int status = gsl_odeiv2_driver_apply(d,&t,ti,y);

    if (status!= GSL_SUCCESS){
      fprintf(fp,"Error, return value = %d\n",status);
      break;
    }

    fprintf("%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",y[0],y[1],y[2],-y[3],-y[4],-y[5],ti);

    if(outcheck(y[0],y[2])){
      break;
    }
  }
  gsl_odeiv2_driver_free(d);
}


int initialize(){
  FILE * fp = fopen("orbit.cfg","r");
  fscanf(fp,"job-name\t%s",jobname);
  fscanf(fp,"directory:outputs\t%s",outdir);
  fscanf(fp,"directory:gdas\t%s",gdadir);
  fscanf(fp, "vxmax\t%lf",vxmax);
  fscanf(fp, "vymax\t%lf",vymax);
  fscanf(fp, "vzmax\t%lf",vzmax);
  fscanf(fp, "Nvx\t%d",Nvx);
  fscanf(fp, "Nvy\t%d",Nvy);
  fscanf(fp, "Nvz\t%d",Nvz);
  fscanf(fp, "Npoints\t%d",Npoints);
  x0 = (double *)malloc(Npoints*sizeof(double));
  z0 = (double *)malloc(Npoints*sizeof(double));
  fscanf(fp, "x0\t%lf",x0[0]);
  for(int i = 1; i<Npoints; i++) fscanf("%lf",x0[i]);
  fscanf(fp, "z0\t%lf",z0[0]);
  for(int i = 1; i<Npoints; i++) fscanf("%lf",z0[i]);
  fscanf(fp, "tend\t%lf",tend);
  fscanf(fp, "nts\t%d",nts);
  fscanf(fp, "minx\t%lf",minx);
  fscanf(fp, "maxx\t%lf",maxx);
  fscanf(fp, "minz\t%lf",minz);
  fscanf(fp, "maxz\t%lf",maxz);
  fscanf(fp, "Lx\t%lf",Lx);
  fscanf(fp, "Lz\t%lf",Lz);
  fscanf(fp, "nx\t%d",nx);
  fscanf(fp, "nz\t%d",nz); 
}
