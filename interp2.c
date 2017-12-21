#include "interp2.h"
#include "solveorbit.h"

field interpfield(fieldgrid mfield, pos xinterp, int NX, int NZ, double Lx, double Lz){
  //Interpolates from field defined by mfield to position defined by xinterp.
  field outfield;
  int xhighindex, xlowindex = 0, zhighindex,zlowindex = 0;
  double dx = Lx/NX; //Centers of boxes defined by x and z, thus domain actually extend to length Lx + dx. 
  double dz = Lz/NZ; //Same for z
  if(xinterp.x>0) xlowindex = (int) xinterp.x/dx; //takes advantage of uniform grid
  if(xlowindex<(NX-1)) xhighindex = xlowindex+1; //takes care of points not on the right edge
  else xhighindex = xlowindex; // Right edge
  if(xinterp.x<0) xhighindex = 0;

  if(xinterp.z>-Lz/2) zlowindex = (int) (xinterp.z+Lz/2)/dz; //interior grid
  if(zlowindex<(NZ-1)) zhighindex = zlowindex+1; //interior grid
  else zhighindex = zlowindex; //top edge
  if(xinterp.z<-Lz/2) zhighindex = 0; //bottom edge
  
  double wll, whl, wlh, whh; //Weights for linear interpolation
  whh = (xinterp.x - mfield.x[xlowindex])/dx*(mfield.z[zlowindex]-xinterp.z)/dz;
  whl = (xinterp.x - mfield.x[xlowindex])/dx*(-mfield.z[zhighindex]+xinterp.z)/dz;
  wlh = (-xinterp.x + mfield.x[xhighindex])/dx*(mfield.z[zlowindex]-xinterp.z)/dz;
  wll = (xinterp.x - mfield.x[xhighindex])/dx*(mfield.z[zhighindex]-xinterp.z)/dz;

  if((zhighindex == 0)||(zlowindex == (NZ-1))){
    //top and bottom edge
    whh = 0;
    wlh = (mfield.z[xlowindex]-xinterp.z)/dz;
    wll = (-mfield.z[xhighindex]+xinterp.z)/dz;
    whl = 0;
  }

  if((xhighindex == 0)||(xlowindex == (NX-1))){
    if((zhighindex == 0)||(zlowindex == (NZ-1))){
      //corner
      whh = 0; whl = 0; wlh = 0; wll = 1;
    }
    else{
      //left and right edge
      whh = 0;
      wlh = 0;
      wll = (mfield.x[xhighindex]-xinterp.x)/dx;
      whl = (xinterp.x-mfield.x[xlowindex])/dx;
    }
  }

  //interpolate all of the fields
  outfield.Ex = wll*mfield.Ex[zlowindex][xlowindex] +wlh*mfield.Ex[zlowindex][xhighindex] + whl*mfield.Ex[zhighindex][xlowindex] + whh*mfield.Ex[zhighindex][xhighindex];
  outfield.Ey = wll*mfield.Ey[zlowindex][xlowindex] +wlh*mfield.Ey[zlowindex][xhighindex] + whl*mfield.Ey[zhighindex][xlowindex] + whh*mfield.Ey[zhighindex][xhighindex];
  outfield.Ez = wll*mfield.Ez[zlowindex][xlowindex] +wlh*mfield.Ez[zlowindex][xhighindex] + whl*mfield.Ez[zhighindex][xlowindex] + whh*mfield.Ez[zhighindex][xhighindex];
  outfield.Bx = wll*mfield.Bx[zlowindex][xlowindex] +wlh*mfield.Bx[zlowindex][xhighindex] + whl*mfield.Bx[zhighindex][xlowindex] + whh*mfield.Bx[zhighindex][xhighindex];
  outfield.By = wll*mfield.By[zlowindex][xlowindex] +wlh*mfield.By[zlowindex][xhighindex] + whl*mfield.By[zhighindex][xlowindex] + whh*mfield.By[zhighindex][xhighindex];
  outfield.Bz = wll*mfield.Bz[zlowindex][xlowindex] +wlh*mfield.Bz[zlowindex][xhighindex] + whl*mfield.Bz[zhighindex][xlowindex] + whh*mfield.Bz[zhighindex][xhighindex];
  return outfield;
}
