#include <stdio.h>
#include <math.h>
#include <DataTypes.h>

//Real absolute(Real);
//const Real pi = 4.0*atan(1.0);

const Real pi = 3.14159265359;

ErrorType SetVelocityField(DM da, int ic, Vec u, Vec v, Vec w, void *g)
{
  StructGrid	*grid = (StructGrid*) g;
  DMDALocalInfo	info;
  VarArray	***arru,***arrv,***arrw;
  Integer 	i,istart,iend;
  Integer 	j,jstart,jend;
  Integer 	k,kstart,kend;
  Real 		dcsi,deta,dzeta;

/*
  DMDAGetLocalInfo(da,&info);
  istart = info.xs;
  iend = info.xs + info.xm;
  jstart = info.ys;
  jend = info.ys + info.ym;
  kstart = info.zs;
  kend = info.zs + info.zm;
*/

  Integer xs,xm,ys,ym,zs,zm;
  DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  istart=xs;	iend=xs+xm;
  jstart=ys;	jend=ys+ym;
  kstart=zs;	kend=zs+zm;

  DMDAVecGetArray(da,u,&arru);
  DMDAVecGetArray(da,v,&arrv);
  DMDAVecGetArray(da,w,&arrw);


  dcsi 		= grid->dcsi;
  deta 		= grid->deta;
  dzeta		= grid->dzeta;

/*
  if (ic == 1) {
    for (i=istart; i<iend; i++) {
    	for (j=jstart; j<jend; j++) {
    	    for (k=kstart; k<kend; k++) {

      		Real csi  = (grid->csimin  + i*dcsi);
      		Real eta  = (grid->etamin  + j*deta);
      		Real zeta = (grid->zetamin + k*dzeta);

		// rotational velocity , divergence free
   		arru[k][j][i] = (pi/3.14)*(-eta + 0.5)	;
     		arrv[k][j][i] = (pi/3.14)*( csi - 0.5)	;
     		arrw[k][j][i] = 0.0			;

    	  }
       }
     }
  }


  else{
    printf("Error (SetVelocityField): Invalid choice of initial conditions.\n");
    return(1);
  }
*/

    for (i=istart; i<iend; i++) {
    	for (j=jstart; j<jend; j++) {
    	    for (k=kstart; k<kend; k++) {

		// rotational velocity , divergence free
   		arru[k][j][i] = 0.0	;
     		arrv[k][j][i] = 0.0	;
     		arrw[k][j][i] = 0.0	;

    	  }
       }
     }

  DMDAVecRestoreArray(da,u,&arru);
  DMDAVecRestoreArray(da,v,&arrv);
  DMDAVecRestoreArray(da,w,&arrw);
  return(0);
}
