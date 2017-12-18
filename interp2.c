#include "interp2.h"
#include "solveorbit.h"

field interpfield(fieldgrid mfield, pos xinterp, int NX, int NZ, double Lx, double Lz){
  field outfield;
  int xhighindex, xlowindex = 0, zhighindex,zlowindex = 0;
  double dx = Lx/NX;
  double dz = Lz/NZ;
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


  outfield.Ex = wll*mfield.Ex[zlowindex][xlowindex] +wlh*mfield.Ex[zlowindex][xhighindex] + whl*mfield.Ex[zhighindex][xlowindex] + whh*mfield.Ex[zhighindex][xhighindex];
  outfield.Ey = wll*mfield.Ey[zlowindex][xlowindex] +wlh*mfield.Ey[zlowindex][xhighindex] + whl*mfield.Ey[zhighindex][xlowindex] + whh*mfield.Ey[zhighindex][xhighindex];
  outfield.Ez = wll*mfield.Ez[zlowindex][xlowindex] +wlh*mfield.Ez[zlowindex][xhighindex] + whl*mfield.Ez[zhighindex][xlowindex] + whh*mfield.Ez[zhighindex][xhighindex];
  outfield.Bx = wll*mfield.Bx[zlowindex][xlowindex] +wlh*mfield.Bx[zlowindex][xhighindex] + whl*mfield.Bx[zhighindex][xlowindex] + whh*mfield.Bx[zhighindex][xhighindex];
  outfield.By = wll*mfield.By[zlowindex][xlowindex] +wlh*mfield.By[zlowindex][xhighindex] + whl*mfield.By[zhighindex][xlowindex] + whh*mfield.By[zhighindex][xhighindex];
  outfield.Bz = wll*mfield.Bz[zlowindex][xlowindex] +wlh*mfield.Bz[zlowindex][xhighindex] + whl*mfield.Bz[zhighindex][xlowindex] + whh*mfield.Bz[zhighindex][xhighindex];
  return outfield;
}
