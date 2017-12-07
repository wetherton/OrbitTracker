typedef struct{
  double Ex, Ey, Ez, Bx, By, Bz;
} field;

typedef struct{
  double **Ex, **Ey, **Ez, **Bx, **By, **Bz;
  double *x, *z;
} fieldgrid;


typedef struct{
  double x, y, z, vx, vy, vz;
}posvel;

typedef struct{
  double x, z;
}pos;

masterfieldmalloc(fieldgrid *, int, int);
