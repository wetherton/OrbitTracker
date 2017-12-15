#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "solveorbit.h"
#include "gdareader.h"
#include "interp2.h"
#include <string.h>

#define USEJAC 0
#define DIM 6
#define TILEX 20
#define TILEZ 20
#define MAGIC 4

const char jobname[128], outdir[256],gdadir[256];
int Nvx, Nvy, Nvz;
double vxmax, vymax, vzmax, Lx, Lz, minx, minz, maxx, maxz;
int Npoints, nts, nx, nz, slice, nthreads;
double tend;
double *x0, *z0;

int main(int argc,char *argv[]){
  //MPI_Init(&argc, &argv);
  int rank, numprocs;
  rank = 0; numprocs = 1;
  //MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  //MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  initialize();
  omp_set_num_threads(nthreads);
  fieldgrid masterfield;
  masterfield = loadfields(nx,nz,Lx,Lz,slice,gdadir);
  double *vx, *vy, *vz;
  int nICs = Nvx*Nvy*Nvz;
  int size = nICs*sizeof(double);
  vx = (double *) malloc(size);
  vy = (double *) malloc(size);
  vz = (double *) malloc(size);
  double dvx,dvy,dvz;
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
    for (int i = rank; i<nICs; i+=numprocs*nthreads*MAGIC){
#pragma omp_parallel_for()
      for(int j =0; j<nthreads*MAGIC; j++){
	int index = j + i;
	if (index>nICs) break;
	posvel IC;
	IC.x = x0[npos];
	IC.y = 0.0;
	IC.z = z0[npos];
	IC.vx = vx[index];
	IC.vy = vy[index];
	IC.vz = vz[index];
	solveorbit(IC,masterfield);
	solveorbitB(IC,masterfield);
      }
    }
  }
  //MPI_Finalize();
  return 0;
}

int dxdt(double t, const double y[], double dydt[], void *params){
  fieldgrid masterfield = *(fieldgrid *) params;
  pos xyz;
  xyz.x = y[0]; xyz.z = y[2];
  field fields = interpfield(masterfield,xyz,nx/2,nz/2,Lx,Lz);
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
  fieldgrid masterfield = *(fieldgrid *) params;
  pos xyz;
  xyz.x = y[0]; xyz.z = y[2];
  field fields = interpfield(masterfield,xyz,nx/2,nz/2,Lx,Lz);
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

int solveorbit(posvel IC, fieldgrid masterfield){
  char filestring[256];
  char buffer[1024];
  char *outstring = (char *) malloc(sizeof(char)*1024*nts);
  *outstring = '\0';
  sprintf(filestring, "%sFx%04.0fz%04.0fvx%03.0fvy%03.0fvz%03.0f.txt\0",jobname,IC.x,IC.z,IC.vx*1000,IC.vy*1000,IC.vz*1000);
  char *out1 = concat(outdir,"/");
  char *out = concat(out1,filestring);
  FILE *fp = fopen(out,"w");
  gsl_odeiv2_system sys = {dxdt, jacobian, DIM, (void *)&masterfield};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,1e-6, 1e-6, 0.0);
  double t = 0.0;
  double y[6] = {IC.x, IC.y, IC.z, IC.vx, IC.vy, IC.vz};
  for(int i = 1; i<nts; i++){
    double ti = i*tend/nts;
    int status = gsl_odeiv2_driver_apply(d,&t,ti,y);

    if (status!= GSL_SUCCESS){
      fprintf(fp,"Error, return value = %d\n",status);
      break;
    }

    sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n\0",y[0],y[1],y[2],y[3],y[4],y[5],ti);
    strcat(outstring,buffer);
    if(outcheck(y[0],y[2])){
      break;
    }
  }
  gsl_odeiv2_driver_free(d);
  fprintf(fp,"%s",outstring);
  fclose(fp);
  free(outstring);
  return 0;
}

int solveorbitB(posvel IC, fieldgrid masterfield){
  char filestring[256];
  char buffer[1024];
  char *outstring = (char *) malloc(sizeof(char)*1024*nts);
  *outstring = '\0';
  sprintf(filestring, "%sBx%04.0fz%04.0fvx%03.0fvy%03.0fvz%03.0f.txt\0",jobname,IC.x,IC.z,1000*IC.vx,1000*IC.vy,1000*IC.vz);
  char *out1 = concat(outdir,"/");
  char *out = concat(out1,filestring);
  FILE *fp = fopen(out,"w");
  gsl_odeiv2_system sys = {dxdtB, jacobian, DIM, (void *)&masterfield};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,1e-6, 1e-6, 0.0);
  double t = 0.0;
  double y[6] = {IC.x, IC.y, IC.z, -IC.vx, -IC.vy, -IC.vz};
  for(int i = 1; i<nts; i++){
    double ti = i*tend/nts;
    int status = gsl_odeiv2_driver_apply(d,&t,ti,y);

    if (status!= GSL_SUCCESS){
      fprintf(fp,"Error, return value = %d\n",status);
      break;
    }

    sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n\0",y[0],y[1],y[2],-y[3],-y[4],-y[5],-ti);
    strcat(outstring,buffer);
    if(outcheck(y[0],y[2])){
      break;
    }
  }
  gsl_odeiv2_driver_free(d);
  fprintf(fp,"%s",outstring);
  fclose(fp);
  free(outstring);
  return 0;
}


int initialize(){
  char buff[256];
  FILE * fp = fopen("orbit.cfg","r");
  fscanf(fp,"job-name %s\n",jobname);
  fscanf(fp,"directory:outputs %s\n",outdir);
  fscanf(fp,"directory:gdas %s\n",gdadir);
  fscanf(fp, "vxmax %lf\n",&vxmax);
  fscanf(fp, "vymax %lf\n",&vymax);
  fscanf(fp, "vzmax %lf\n",&vzmax);
  fscanf(fp, "Nvx %d\n",&Nvx);
  fscanf(fp, "Nvy %d\n",&Nvy);
  fscanf(fp, "Nvz %d\n",&Nvz);
  fscanf(fp, "Npoints %d\n",&Npoints);
  x0 = (double *)malloc(Npoints*sizeof(double));
  z0 = (double *)malloc(Npoints*sizeof(double));
  if(Npoints == 1){
    fscanf(fp, "x0 %lf\n",&x0[0]);
    fscanf(fp, "z0 %lf\n",&z0[0]);
  }
  else{
    fscanf(fp, "x0 %lf",&x0[0]);
    for(int i = 1; i<(Npoints-1); i++) fscanf(fp, " %lf",&x0[i]);
    fscanf(fp,"%lf\n",&x0[Npoints-1]);
    
    fscanf(fp, "z0 %lf",&z0[0]);
    for(int i = 1; i<(Npoints-1); i++) fscanf(fp, " %lf",&z0[i]);
    fscanf(fp,"%lf\n",&z0[Npoints-1]);
  }
  fscanf(fp, "tend %lf\n",&tend);
  fscanf(fp, "nts %d\n",&nts);
  fscanf(fp, "minx %lf\n",&minx);
  fscanf(fp, "maxx %lf\n",&maxx);
  fscanf(fp, "minz %lf\n",&minz);
  fscanf(fp, "maxz %lf\n",&maxz);
  fscanf(fp, "Lx %lf\n",&Lx);
  fscanf(fp, "Lz %lf\n",&Lz);
  fscanf(fp, "nx %d\n",&nx);
  fscanf(fp, "nz %d\n",&nz);
  fscanf(fp, "slice %d\n",&slice);
  fscanf(fp,"omp:nthreads %d\n",&nthreads);
  fclose(fp);
  return 0;
}

int outcheck(double x,double z){
  if((x>maxx)||(x<minx)||(z>maxz)||(z<minz)) return 1;
  return 0;
}
