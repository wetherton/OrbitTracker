char* concat(const char *s1, const char *s2);
float ** loadgda(const char q[], int slice, int nx, int nz,const char datadir[]);
void loadfields(fieldgrid* masterfield, int NX, int NZ);
void masterfieldmalloc(fieldgrid *, int, int);
double * loadx(double, int);
double * loadz(double, int);
