#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solveorbit.h"
#include "gdareader.h"

char* concat(const char *s1, const char *s2){
  char *result = malloc(strlen(s1)+strlen(s2)+1);
  strcpy(result, s1);
  strcat(result, s2);
  return result;
}

double ** loadgda(const char q[], int slice, int nx, int nz, const char datadir[]){
  char * im1 = concat(datadir,"/");
  char * im2 = concat(im1, q);
  char * fname = concat(im2,".gda");
  double **outdata;
  outdata = setup2Darray(nx,nz);
  float *buffer = (float *) malloc(sizeof(float)*nx*nz);
  FILE *fp = fopen(fname, "rb");
  for(int nslice = 0; nslice<slice; nslice++){
    fread(buffer, sizeof(float), nx*nz, fp);
  }

  fclose(fp);
  for(int i = 0; i<nz; i++){
    for(int j = 0; j<nx; j++){
      outdata[i][j] = (double) buffer[j+nx*i];
    }
  }
  free(buffer);
  // Decimation process begins here

  double ** outdata1;
  outdata1 = setup2Darray(nx/2,nz);

  for(int i = 0; i < nz; i++){
    for(int j = 0; j< nx/2; j++){
      outdata1[i][j] = (outdata[i][2*j] + outdata[i][2*j+1])/2;
    }
  }

  free(outdata);

  double **outdata2 = setup2Darray(nx/2,nz/2);
  

  for(int i = 0; i < nz/2; i++){
    for(int j = 0; j< nx/2; j++){
      outdata2[i][j] = (outdata1[2*i][j] + outdata1[2*i+1][j])/2;
    }
  }

  free(outdata1);
  return outdata2;
}

fieldgrid loadfields(int NX,int NZ, double LX, double LZ, int slice, const char datadir[]){
  fieldgrid masterfield = masterfieldmalloc(NX,NZ);
  masterfield.Ex = loadgda("ex",slice, NX, NZ, datadir);
  masterfield.Ey = loadgda("ey",slice, NX, NZ, datadir);
  masterfield.Ez = loadgda("ez",slice, NX, NZ, datadir);
  masterfield.Bx = loadgda("bx",slice, NX, NZ, datadir);
  masterfield.By = loadgda("by",slice, NX, NZ, datadir);
  masterfield.Bz = loadgda("bz",slice, NX, NZ, datadir);
  masterfield.x =loadx(LX,NX/2);
  masterfield.z =loadz(LZ,NZ/2);
  return masterfield;
}

fieldgrid masterfieldmalloc(int nx,int nz){
  fieldgrid * temp = (fieldgrid *)malloc(sizeof(double **)*6+sizeof(double*)*2);
  fieldgrid masterfield = *temp;
  masterfield.Ex = (double **)malloc(sizeof(double*)*nz/2);
  masterfield.Ex[0] = (double *)malloc(sizeof(double)*nx*nz/4);
  masterfield.Ey = (double **)malloc(sizeof(double*)*nz/2);
  masterfield.Ey[0] = (double *)malloc(sizeof(double)*nx*nz/4);
  masterfield.Ez = (double **)malloc(sizeof(double*)*nz/2);
  masterfield.Ez[0] = (double *)malloc(sizeof(double)*nx*nz/4);
  masterfield.Bx = (double **)malloc(sizeof(double*)*nz/2);
  masterfield.Bx[0] = (double *)malloc(sizeof(double)*nx*nz/4);
  masterfield.By = (double **)malloc(sizeof(double*)*nz/2);
  masterfield.By[0] = (double *)malloc(sizeof(double)*nx*nz/4);
  masterfield.Bz = (double **)malloc(sizeof(double*)*nz/2);
  masterfield.Bz[0] = (double *)malloc(sizeof(double)*nx*nz/4);
  masterfield.x = (double *)malloc(sizeof(double)*nx/2);
  masterfield.z = (double *)malloc(sizeof(double)*nz/2);
  return masterfield;
}

double ** setup2Darray(int nx, int nz){
  double ** array = (double **) malloc(sizeof(double*)*nz);
  array[0] = (double *)malloc(sizeof(double)*nx*nz);
  for(int i = 0; i< nz; i++){
    array[i] = (*array + nx*i);
  }
  return array;
}

double *loadx(double Lx, int nx){
  double *x = (double *) malloc(sizeof(double)*nx);
  double dx = Lx/(nx);
  for(int i = 0; i< nx; i++){
    x[i] = i*dx;
  }
  return x;
}

double *loadz(double Lz, int nz){
  double *z = (double *) malloc(nz*sizeof(double));
  double dz = Lz/(nz);
  for(int i = 0; i< nz; i++){
    z[i] = -(Lz+dz)/2+i*dz;
  }
  return z;
}
