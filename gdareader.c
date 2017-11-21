#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
  return outdata;
}
