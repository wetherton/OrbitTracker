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

void loadgda(double ** outdata2, const char q[], int slice, int nx, int nz, const char datadir[]){
  char * im1 = concat(datadir,"/");
  char * im2 = concat(im1, q);
  char * fname = concat(im2,".gda");
  float **outdata;
  outdata = (float **)malloc(sizeof(float *)*nz);
  outdata[0] = (float *)malloc(sizeof(float)*nx*nz);
  for(int i = 0; i< nz; i++){
    outdata[i] = (*outdata + nx*i);
  }
  float *buffer = (float *) malloc(sizeof(float)*nx*nz);
  FILE *fp = fopen(fname, "rb");
  for(int nslice = 0; nslice<slice; nslice++){
    fread(buffer, sizeof(float), nx*nz, fp);
  }

  fclose(fp);
  for(int i = 0; i<nz; i++){
    for(int j = 0; j<nx; j++){
      outdata[i][j] = buffer[j+nx*i];
    }
  }
  free(buffer);
  // Decimation process begins here

  float ** outdata1;
  outdata1 = (float **)malloc(sizeof(float *)*nz);
  outdata1[0] = (float *)malloc(sizeof(float)*nx*nz/2);
  for(int i = 0; i< nz; i++){
    outdata1[i] = (*outdata1 + nx*i/2);
  }

  for(int i = 0; i < nz; i++){
    for(int j = 0; j< nx/2; j++){
      outdata1[i][j] = (outdata[i][2*j] + outdata[i][2*j+1])/2;
    }
  }

  free(outdata);
  
  outdata2 = (double **)malloc(sizeof(double *)*nz/2);
  outdata2[0] = (double *)malloc(sizeof(double)*nx*nz/4);
  for(int i = 0; i< nz/2; i++){
    outdata2[i] = (*outdata2 + nx*i/2);
  }

  for(int i = 0; i < nz/2; i++){
    for(int j = 0; j< nx/2; j++){
      outdata2[i][j] = (double)(outdata1[2*i][j] + outdata1[2*i+1][j])/2;
    }
  }

  free(outdata1);
}

void loadfields(fieldgrid masterfield,int NX,int NZ, double LX, double LZ, int slice, const char datadir[]){
  loadgda(masterfield.Ex, "ex",slice, NX, NZ, datadir);
  loadgda(masterfield.Ey, "ey",slice, NX, NZ, datadir);
  loadgda(masterfield.Ez, "ez",slice, NX, NZ, datadir);
  loadgda(masterfield.Bx, "bx",slice, NX, NZ, datadir);
  loadgda(masterfield.By, "by",slice, NX, NZ, datadir);
  loadgda(masterfield.Bz, "bz",slice, NX, NZ, datadir);
  loadx(masterfield.x,LX,NX);
  loadz(masterfield.z,LZ,NZ);
}

void masterfieldmalloc(fieldgrid masterfield,int nx,int nz){
  fieldgrid * temp = (fieldgrid *)malloc(sizeof(double **)*6+sizeof(double*)*2);
  masterfield = *temp;
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
}

void loadx(double *x,double Lx, int nx){
  double dx = Lx/(nx-1);
  for(int i = 0; i< nx; i++){
    x[i] = i*dx;
  }
}

void loadz(double *z, double Lz, int nz){
  double dz = Lz/(nz-1);
  for(int i = 0; i< nz; i++){
    z[i] = -Lz/2+i*dz;
  }
}
