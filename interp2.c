#include "interp2.h"
#include "solveorbit.h"

field interpfield(fieldgrid mfield, pos xinterp, int NX, int NZ, double Lx, double Lz){
  field outfield;
  int xhighindex, xlowindex, zhighindex,zlowindex;
  double dx = Lx/(NX-1);
  double dz = Lz/(NZ-1);
  if(xinterp.x>0) xlowindex = (int) xinterp.x/dx;
  if(xlowindex<(NX-1)) xhighindex = xlowindex+1;
  else xhighindex = xlowindex;

  if(xinterp.z>-Lz/2) zlowindex = (int) (xinterp.z+Lz/2)/dz;
  if(zlowindex<(NZ-1)) zhighindex = zlowindex+1;
  else zhighindex = zlowindex;
  double wll, whl, wlh, whh;
  whh = (xinterp.x - mfield.x[xlowindex])/dx*(mfield.z[zlowindex]-xinterp.z)/dz;
  whl = (xinterp.x - mfield.x[xlowindex])/dx*(-mfield.z[zhighindex]+xinterp.z)/dz;
  wlh = (-xinterp.x + mfield.x[xhighindex])/dx*(mfield.z[zlowindex]-xinterp.z)/dz;
  wll = (xinterp.x - mfield.x[xhighindex])/dx*(mfield.z[zhighindex]-xinterp.z)/dz;

  if((zhighindex == 0)||(zlowindex == (NZ-1))){
    whh = 0;
    wlh = (mfield.z[xlowindex]-xinterp.z)/dz;
    wll = (-mfield.z[xhighindex]+xinterp.z)/dz;
    whl = 0;
  }

  if((xhighindex == 0)||(xlowindex == (NX-1))){
    if((zhighindex == 0)||(zlowindex == (NZ-1))){
      whh = 0; whl = 0; wlh = 0; wll = 1;
    }
    else{
      whh = 0;
      wlh = 0;
      wll = (mfield.x[xhighindex]-xinterp.x)/dx;
      whl = (xinterp.x-mfield.x[xlowindex])/dx;
    }
  }


  outfield.Ex = wll*mfield.Ex[xlowindex][zlowindex] +wlh*mfield.Ex[xlowindex][zhighindex] + whl*mfield.Ex[xhighindex][zlowindex] + whh*mfield.Ex[xhighindex][zhighindex];
  outfield.Ey = wll*mfield.Ey[xlowindex][zlowindex] +wlh*mfield.Ey[xlowindex][zhighindex] + whl*mfield.Ey[xhighindex][zlowindex] + whh*mfield.Ey[xhighindex][zhighindex];
  outfield.Ez = wll*mfield.Ez[xlowindex][zlowindex] +wlh*mfield.Ez[xlowindex][zhighindex] + whl*mfield.Ez[xhighindex][zlowindex] + whh*mfield.Ez[xhighindex][zhighindex];
  outfield.Bx = wll*mfield.Bx[xlowindex][zlowindex] +wlh*mfield.Bx[xlowindex][zhighindex] + whl*mfield.Bx[xhighindex][zlowindex] + whh*mfield.Bx[xhighindex][zhighindex];
  outfield.By = wll*mfield.By[xlowindex][zlowindex] +wlh*mfield.By[xlowindex][zhighindex] + whl*mfield.By[xhighindex][zlowindex] + whh*mfield.By[xhighindex][zhighindex];
  outfield.Bz = wll*mfield.Bz[xlowindex][zlowindex] +wlh*mfield.Bz[xlowindex][zhighindex] + whl*mfield.Bz[xhighindex][zlowindex] + whh*mfield.Bz[xhighindex][zhighindex];
  return outfield;
}
