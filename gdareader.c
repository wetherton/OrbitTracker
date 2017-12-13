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

float ** loadgda(const char q[], int slice, int nx, int nz,const char datadir[]){
  char * im1 = concat(datadir,"/");
  char * im2 = concat(im1, q);
  char * fname = concat(im2,".gda");
  float ** outdata;
  outdata = (float **)malloc(sizeof(float *)*nz);
  outdata[0] = (float *)malloc(sizeof(float)*nx*nz);
  for(int i = 0; i< nz; i++){
    outdata[i] = (*outdata + nx*i);
  }
  FILE *fp = fopen(fname, "rb");
  free(im1); free(im2); free(fname);
  for(int nslice = 0; nsclice<slice; slice++){
    fread(outdata, sizeof(float), nx*nz, fp);
  }
  // Decimation process begins here

  float ** outdata1;
  outdata = (float **)malloc(sizeof(float *)*nz);
  outdata[0] = (float *)malloc(sizeof(float)*nx*nz/2);
  for(int i = 0; i< nz; i++){
    outdata[i] = (*outdata + nx*i/2);
  }

  for(int i = 0; i < nz; i++){
    for(int j = 0; j< nx/2; j++){
      outdata1[i][j] = (outdata[i][2*j] + outdata[i][2*j+1])/2;
    }
  }

  free(outdata);
  
  float ** outdata2;
  outdata = (float **)malloc(sizeof(float *)*nz/2);
  outdata[0] = (float *)malloc(sizeof(float)*nx*nz/4);
  for(int i = 0; i< nz/2; i++){
    outdata[i] = (*outdata + nx*i/2);
  }

  for(int i = 0; i < nz/2; i++){
    for(int j = 0; j< nx/2; j++){
      outdata2[i][j] = (outdata1[2*i][j] + outdata1[2*i+1][j])/2;
    }
  }

  free(outdata1);
  return outdata2;
}

void loadfields(fieldgrid masterfield,int NX,int NZ, double LX, double LZ, int slice, const char datadir[]){
  masterfield.Ex = loadgda("ex",slice, NX, NZ, datadir);
  masterfield.Ey = loadgda("ey",slice, NX, NZ, datadir);
  masterfield.Ez = loadgda("ez",slice, NX, NZ, datadir);
  masterfield.Bx = loadgda("bx",slice, NX, NZ, datadir);
  masterfield.By = loadgda("by",slice, NX, NZ, datadir);
  masterfield.Bz = loadgda("bz",slice, NX, NZ, datadir);
  masterfield.x  = loadx(LX,NX);
  masterfield.z  = loadz(LZ,NZ);
}
fieldgrid masterfield;

void masterfieldmalloc(fieldgrid masterfield,int nx,int nz){
  masterfield = (fieldgrid *)malloc(sizeof(float **)*6+sizeof(float*)*2);
  masterfield.Ex = (float **)malloc(sizeof(float*)*nz/2);
  masterfield.Ex[0] = (float *)malloc(sizeof(float)*nx*nz/4);
  masterfield.Ey = (float **)malloc(sizeof(float*)*nz/2);
  masterfield.Ey[0] = (float *)malloc(sizeof(float)*nx*nz/4);
  masterfield.Ez = (float **)malloc(sizeof(float*)*nz/2);
  masterfield.Ez[0] = (float *)malloc(sizeof(float)*nx*nz/4);
  masterfield.Bx = (float **)malloc(sizeof(float*)*nz/2);
  masterfield.Bx[0] = (float *)malloc(sizeof(float)*nx*nz/4);
  masterfield.By = (float **)malloc(sizeof(float*)*nz/2);
  masterfield.By[0] = (float *)malloc(sizeof(float)*nx*nz/4);
  masterfield.Bz = (float **)malloc(sizeof(float*)*nz/2);
  masterfield.Bz[0] = (float *)malloc(sizeof(float)*nx*nz/4);
  masterfield.x = (float *)malloc(sizeof(float)*nx/2);
  masterfield.z = (float *)malloc(sizeof(float)*nz/2);
}

double *loadx(double Lx, int nx){
  double * out = (double *) malloc (sizeof(double)*nx);
  double dx = Lx/(nx-1);
  for(int i = 0; i< nx; i++){
    out[i] = i*dx;
  }
  return out;
}

double *loadz(double Lz, int nz){
  double * out = (double *) malloc (sizeof(double)*nz);
  double dz = Lz/(nz-1);
  for(int i = 0; i< nz; i++){
    out[i] = -Lz/2+i*dz;
  }
  return out;
}
