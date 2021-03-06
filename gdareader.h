
#ifndef GDAREADER_H
#define GDAREADER_H
char* concat(const char *s1, const char *s2);
double ** loadgda(const char q[], int slice, int nx, int nz,const char datadir[]);
fieldgrid loadfields(int NX,int NZ, double LX, double LZ, int slice, const char datadir[]);
fieldgrid masterfieldmalloc( int, int);
double * loadx(double, int);
double * loadz(double, int);
double ** setup2Darray( int, int);
double ** loadgdaPeter(const char q[], int slice, int nx, int nz,const char datadir[]);
fieldgrid loadfieldsPeter(int NX,int NZ, double LX, double LZ, int slice, const char datadir[]);
double ** smoothfield(double **, int, int, int);
#endif
