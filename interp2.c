fields interpfield(fieldgrid field, pos xinterp){
  fields outfield;
  outfield.x = xinterp.x;
  outfield.z = xinterp.z;
  int xhighindex, xlowindex, zhighindex,zlowindex;
  for(int i = 0; i< NX; i++){
    if (xinterp.x>field.x[i])
      {
	xhighindex = i;
	if(i>0) xlowindex = i-1;
	break;
      }
    if(i==NX-1){
      xhighindex = NX-1;
      xlowindex = NX -1;
    }
  }
  
  for(int i = 0; i< NZ; i++){
    if (xinterp.z<field.z[i])
      {
        zhighindex = i;
        if(i>0) zlowindex = i-1;
        break;
      }
    if(i==NZ-1){
      zhighindex = NZ-1;
      zlowindex = NZ -1;
    }
  }
  double wll, whl, wlh, whh;
  double dx = field.x[1]-field.x[0];
  double dz = field.z[0]-field.z[1];
  whh = (xinterp.x - field.x[xlowindex])/dx*(field.z[zlowindex]-xinterp.z)/dz;
  whl = (xinterp.x - field.x[xlowindex])/dx*(-field.z[zhighindex]+xinterp.z)/dz;
  wlh = (-xinterp.x + field.x[xhighindex])/dx*(field.z[zlowindex]-xinterp.z)/dz;
  wll = (xinterp.x - field.x[xhighindex])/dx*(field.z[zhighindex]-xinterp.z)/dz;

  if((zhighindex == 0)||(zlowindex == (NZ-1))){
    whh = 0;
    wlh = (field.z[xlowindex]-xinterp.z)/dz;
    wll = (-field.z[xhighindex]+xinterp.z)/dz;
    whl = 0;
  }

  if((xhighindex == 0)||(xlowindex == (NX-1))){
    if((zhighindex == 0)||(zlowindex == (NZ-1))){
      whh = 0; whl = 0; wlh = o; wll = 1;
    }
    else{
      whh = 0;
      wlh = 0;
      wll = (field.x[xhighindex]-xinterp.x)/dx;
      whl = (xinterp.x-field.x[lowindex])/dx;
    }
  }


  outfield.Ex = wll*field.Ex[xlowindex][zlowindex] +wlh*field.Ex[xlowindex][zhighindex] + whl*field.Ex[xhighindex][zlowindex] + whh*field.Ex[xhighindex][zhighindex];
  outfield.Ey = wll*field.Ey[xlowindex][zlowindex] +wlh*field.Ey[xlowindex][zhighindex] + whl*field.Ey[xhighindex][zlowindex] + whh*field.Ey[xhighindex][zhighindex];
  outfield.Ez = wll*field.Ez[xlowindex][zlowindex] +wlh*field.Ez[xlowindex][zhighindex] + whl*field.Ez[xhighindex][zlowindex] + whh*field.Ez[xhighindex][zhighindex];
  outfield.Bx = wll*field.Bx[xlowindex][zlowindex] +wlh*field.Bx[xlowindex][zhighindex] + whl*field.Bx[xhighindex][zlowindex] + whh*field.Bx[xhighindex][zhighindex];
  outfield.By = wll*field.By[xlowindex][zlowindex] +wlh*field.By[xlowindex][zhighindex] + whl*field.By[xhighindex][zlowindex] + whh*field.By[xhighindex][zhighindex];
  outfield.Bz = wll*field.Bz[xlowindex][zlowindex] +wlh*field.Bz[xlowindex][zhighindex] + whl*field.Bz[xhighindex][zlowindex] + whh*field.Bz[xhighindex][zhighindex];
  return outfield;
}
