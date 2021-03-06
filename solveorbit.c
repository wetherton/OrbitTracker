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
#include <time.h>

#define USEJAC 0
#define DIM 6
#define TILEX 20 //Eventually planned to be used to manage memory for interp grid
#define TILEZ 20

//globals allow fewer arguments to be passed.

const char jobname[128], outdir[256],gdadir[256];
int Nvx, Nvy, Nvz;
double vxmax, vymax, vzmax, Lx, Lz, minx, minz, maxx, maxz;
int Npoints, nts, nx, nz, slice, nthreads, DoTrack, fileext;
double tend;
double *x0, *z0;
double wpewce;

int main(int argc,char *argv[]){
  clock_t start = clock(), diff;
  MPI_Init(&argc, &argv);
  int rank, numprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  initialize(); //reads "orbit.cfg" and sets the globals
  omp_set_num_threads(nthreads); //for scaling characterization
  fieldgrid masterfield;
  if(slice==0) //Slice = 0 implies gdas are made for a specific timestep.
    masterfield = loadfieldsPeter(nx,nz,Lx,Lz,fileext,gdadir); 
  else
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
  //set up IC arrays
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
  for(int npos = 0; npos<Npoints; npos++){ //go by position
    char endfilename[512];
    *endfilename = '\0';
    strcat(endfilename,outdir);
    strcat(endfilename,"/");
    char temp[256];
    sprintf(temp, "%sx%04.0fz%04.0frank%d.tsv",jobname,x0[npos],z0[npos],rank);
    strcat(endfilename,temp);
    FILE *endpointfile = fopen(endfilename,"w");
#pragma omp parallel for schedule(dynamic) //allow omp to dynamically parallelize most of the ICs.
    for(int j = rank; j<(nICs); j+=numprocs){ //Assign ranks ICs on an alternating basis (hope close in velocity space implies close in load)
	posvel IC;
	IC.x = x0[npos];
	IC.y = 0.0;
	IC.z = z0[npos];
	IC.vx = vx[j];
	IC.vy = vy[j];
	IC.vz = vz[j];
	if(DoTrack)
	  //Only forward track for traces
	  solveorbit(IC,masterfield);
	posvel endpoint = solveorbitB(IC,masterfield); //creates mapping
        #pragma omp critical
	{
	  //Only allows one thread to write to the mapping file at a time.
	  fprintf(endpointfile,  "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n"\
		  ,IC.x,IC.z,IC.vx,IC.vy,IC.vz,endpoint.x,endpoint.z,endpoint.vx,endpoint.vy,endpoint.vz);
	}
      }
    fclose(endpointfile);
  }
  diff = clock() - start;
  int sec = diff/CLOCKS_PER_SEC; //A poor measure of time
  printf("Time: %d s\n", sec); 
  MPI_Finalize();
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
  //straight up Lorentz force
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
  //Lorentz force under time reversal.
}


int jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  //could be used with implicit method
  (void) (t);
  (void) (params);
  (void) (y);
  dfdy = NULL;
  for(int i = 0; i< DIM; i++){
    dfdt[i] = 0.0;
  }
  
  return GSL_SUCCESS;
}

/*posvel solveorbit(posvel IC, fieldgrid masterfield){
  //Forward track
  char filestring[256];
  char buffer[1024];
  char *outstring = (char *) malloc(sizeof(char)*1024*nts);
  *outstring = '\0';
  sprintf(filestring, "%sFx%04.0fz%04.0fvx%03.0fvy%03.0fvz%03.0f.tsv",jobname,IC.x,IC.z,IC.vx*1000,IC.vy*1000,IC.vz*1000);
  char *out1 = concat(outdir,"/");
  char *out = concat(out1,filestring);
  FILE *fp = fopen(out,"w");
  gsl_odeiv2_system sys2 = {dxdt, jacobian, DIM, (void *)&masterfield};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,1e-6, 1e-6, 0.0); //Runge-Kutta 8 now, could be an issue
  double t = 0.0;
  double y[6] = {IC.x, IC.y, IC.z, IC.vx, IC.vy, IC.vz};
  sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",y[0],y[1],y[2],y[3],y[4],y[5],t);
  strcat(outstring,buffer);
  for(int i = 1; i<nts; i++){
    double ti = i*tend/nts;
    int status = gsl_odeiv2_driver_apply(d,&t,ti,y);

    if (status!= GSL_SUCCESS){
      fprintf(fp,"Error, return value = %d\n",status);
      return IC;
    }

    sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",y[0],y[1],y[2],y[3],y[4],y[5],ti);
    if(DoTrack)
      strcat(outstring,buffer);
    if(outcheck(y[0],y[2])){
      break;
    }
  }
  posvel FC;
  FC.x = y[0]; FC.y = y[1]; FC.z = y[2]; FC.vx = y[3]; FC.vy = y[4]; FC.vz = y[5];
  gsl_odeiv2_driver_free(d);
  fprintf(fp,"%s",outstring);
  fclose(fp);
  free(outstring);
  return FC;
  }*/

posvel solveorbit(posvel IC, fieldgrid masterfield){
  char filestring[256];
  char buffer[1024];
  char *outstring = (char *) malloc(sizeof(char)*1024*nts);
  *outstring = '\0';
  sprintf(filestring, "%sFx%04.0fz%04.0fvx%03.0fvy%03.0fvz%03.0f.tsv",jobname,IC.x,IC.z,1000*IC.vx,1000*IC.vy,1000*IC.vz);
  char *out1 = concat(outdir,"/");
  char *out = concat(out1,filestring);
  FILE *fp = fopen(out,"w");
  //gsl_odeiv2_system sys = {dxdtB, jacobian, DIM, (void *)&masterfield};
  //gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk2,1e-6, 1e-6, 0.0);
  double t = 0.0;
  double dt = wpewce/20;
  //double y[6] = {IC.x, IC.y, IC.z, -IC.vx, -IC.vy, -IC.vz}; //Time reversal
  posvel part;
  part.x = IC.x; part.y = IC.y; part.z = IC.z; part.vx = IC.vx; part.vy = IC.vy; part.vz = IC.vz;
  if(DoTrack){
    //sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",y[0],y[1],y[2],y[3],y[4],y[5],t);
    sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",part.x,part.y,part.z,part.vx,part.vy,part.vz,t);
    strcat(outstring,buffer);
  }
  pos xyz;
  xyz.x = IC.x; xyz.z = IC.z;
  field f = interpfield(masterfield,xyz,nx/2,nz/2,Lx,Lz);
  part = BorisUpdate(part,f,-0.5*dt);
  for(int i = 1; i<nts; i++){
    double ti = i*tend/nts;
    for(double tp = (i-1)*tend/nts; tp<ti; tp+=dt){
      xyz.x = part.x; xyz.z = part.z;
      f = interpfield(masterfield,xyz,nx/2,nz/2,Lx,Lz);
      part = BorisUpdate(part,f,dt);
      part = particlepush(part,dt);
    } 
    /*int status = gsl_odeiv2_driver_apply(d,&t,ti,y);

    if (status!= GSL_SUCCESS){
      if(DoTrack)
	fprintf(fp,"Error, return value = %d\n",status);
      return IC;
      }*/
    if(DoTrack){
      //sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",y[0],y[1],y[2],-y[3],-y[4],-y[5],-ti);
      sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",part.x,part.y,part.z,part.vx,part.vy,part.vz,ti);
      strcat(outstring,buffer);
    }
    if(outcheck(part.x,part.z)){
      break;
    }
  }
  posvel FC;
  FC.x = part.x; FC.y = part.y; FC.z = part.z; FC.vx = part.vx; FC.vy = part.vy; FC.vz = part.vz; 
  //gsl_odeiv2_driver_free(d);
  if(DoTrack)
    //Since this can be run without the full track, needs a conditional
    fprintf(fp,"%s",outstring);
  fclose(fp);
  free(outstring);
  return FC;
}

posvel solveorbitB(posvel IC, fieldgrid masterfield){
  //Back track
  char filestring[256];
  char buffer[1024];
  char *outstring = (char *) malloc(sizeof(char)*1024*nts);
  *outstring = '\0';
  sprintf(filestring, "%sBx%04.0fz%04.0fvx%03.0fvy%03.0fvz%03.0f.tsv",jobname,IC.x,IC.z,1000*IC.vx,1000*IC.vy,1000*IC.vz);
  char *out1 = concat(outdir,"/");
  char *out = concat(out1,filestring);
  FILE *fp;
  if(DoTrack)
    fp = fopen(out,"w");
  //gsl_odeiv2_system sys = {dxdtB, jacobian, DIM, (void *)&masterfield};
  //gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk2,1e-6, 1e-6, 0.0);
  double t = 0.0;
  double dt = wpewce/(20);
  //double y[6] = {IC.x, IC.y, IC.z, -IC.vx, -IC.vy, -IC.vz}; //Time reversal
  posvel part;
  part.x = IC.x; part.y = IC.y; part.z = IC.z; part.vx = IC.vx; part.vy = IC.vy; part.vz = IC.vz;
  if(DoTrack){
    //sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",y[0],y[1],y[2],-y[3],-y[4],-y[5],-t);
    sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",part.x,part.y,part.z,-part.vx,-part.vy,-part.vz,-t);
    strcat(outstring,buffer);
  }
  pos xyz;
  xyz.x = IC.x; xyz.z = IC.z;
  field f = interpfield(masterfield,xyz,nx/2,nz/2,Lx,Lz);
;
  part = BorisUpdate(part,f,0.5*dt);
  for(int i = 1; i<nts; i++){
    double ti = i*tend/nts;
    for(double tp = (i-1)*tend/nts; tp<ti; tp+=dt){
      xyz.x = part.x; xyz.z = part.z;
      f = interpfield(masterfield,xyz,nx/2,nz/2,Lx,Lz);
      part = BorisUpdate(part,f,-dt);
      part = particlepush(part,-dt);
    } 
    /*int status = gsl_odeiv2_driver_apply(d,&t,ti,y);

    if (status!= GSL_SUCCESS){
      if(DoTrack)
	fprintf(fp,"Error, return value = %d\n",status);
      return IC;
      }*/
    if(DoTrack){
      //sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",y[0],y[1],y[2],-y[3],-y[4],-y[5],-ti);
      sprintf(buffer, "%.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\n",part.x,part.y,part.z,-part.vx,-part.vy,-part.vz,-ti);
      strcat(outstring,buffer);
    }
    if(outcheck(part.x,part.z)){
      break;
    }
  }
  posvel FC;
  FC.x = part.x; FC.y = part.y; FC.z = part.z; FC.vx = -part.vx; FC.vy = -part.vy; FC.vz = -part.vz; //time reversal
  //gsl_odeiv2_driver_free(d);
  if(DoTrack){
    //Since this can be run without the full track, needs a conditional
    fprintf(fp,"%s",outstring);
    fclose(fp);
  }
  free(outstring);
  return FC;
}


int initialize(){
  //reads in the parameters desired. Not powerful now, as the file needs to follow this exact format.
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
  fscanf(fp, "wpewce %lf\n",&wpewce);
  fscanf(fp, "slice %d\n",&slice);
  if(slice==0)
    fscanf(fp, "timestep %d\n", &fileext);
  fscanf(fp,"omp:nthreads %d\n",&nthreads);
  fscanf(fp, "fulltrack %d\n", &DoTrack);
  fclose(fp);
  return 0;
}

int outcheck(double x,double z){
  //Return 1 if you've left the bounding box
  if((x>maxx)||(x<minx)||(z>maxz)||(z<minz)) return 1;
  return 0;
}

posvel BorisUpdate(posvel part, field f, double dt) {
  //stolen mostly from VPIC Boris method (best to match, right?)
  double v_minus[3];
  double v_plus[3];
  double v0, v1, v2, v3, v4;
  const double one_third  = 1.0/3.0;
  const double two_fifteenths = 2.0/15.0;
  /*v minus, electric field half update*/
   v_minus[0] = part.vx - f.Ex*0.5*dt; //relativistic momentum
   v_minus[1] = part.vy - f.Ey*0.5*dt;
   v_minus[2] = part.vz - f.Ez*0.5*dt;


   double gamma = sqrt(1.0+(v_minus[0]*v_minus[0]+(v_minus[1]*v_minus[1]+v_minus[2]*v_minus[2]))); //Middle of timestep

   v0   = -dt*gamma/2.0;
   /**/                                      // Boris - scalars
   v1   = f.Bx*f.Bx + (f.By*f.By + f.Bz*f.Bz);
   v2   = (v0*v0)*v1;
   v3   = v0*(1.0+v2*(one_third+v2*two_fifteenths));
   v4   = v3/(1.0+v1*(v3*v3));
   v4  += v4;
   v0   = v_minus[0] + v3*( v_minus[1]*f.Bz - v_minus[2]*f.By );       // Boris - uprime
   v1   = v_minus[1] + v3*( v_minus[2]*f.Bx - v_minus[0]*f.Bz );
   v2   = v_minus[2] + v3*( v_minus[0]*f.By - v_minus[1]*f.Bx );
   v_plus[0] = v_minus[0] + v4*( v1*f.Bz - v2*f.By );            // Boris - rotation
   v_plus[1] = v_minus[1] + v4*( v2*f.Bx - v0*f.Bz );
   v_plus[2] = v_minus[2] + v4*( v0*f.By - v1*f.Bx );
   
   part.vx = v_plus[0] - f.Ex*0.5*dt;
   part.vy = v_plus[1] - f.Ey*0.5*dt;
   part.vz = v_plus[2] - f.Ez*0.5*dt;
   
  return part;
}

posvel particlepush(posvel part, double dt){
  double gamma = sqrt(1+part.vx*part.vx+part.vy*part.vy+part.vz*part.vz);
  double vx = part.vx/gamma;
  double vy = part.vy/gamma;
  double vz = part.vz/gamma;

  part.x += vx*dt;
  part.y += vy*dt;
  part.z += vz*dt;
  return part;
}
