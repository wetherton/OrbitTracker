#ifndef GDAREADER_H
#define GDAREADER_H
char* concat(const char *s1, const char *s2);
float ** loadgda(const char q[], int slice, int nx, int nz,const char datadir[]);
void loadfields(fieldgrid masterfield,int NX,int NZ, double LX, double LZ, int slice, const char datadir[]);
void masterfieldmalloc(fieldgrid , int, int);
double * loadx(double, int);
double * loadz(double, int);
#endif
