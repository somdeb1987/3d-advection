#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <DataTypes.h>

//Real absolute(Real);
//const Real pi = 4.0*atan(1.0);

//const Real pi = 3.14159265359;

ErrorType SetScalarField(DM da, int ic, Vec Scalar, void *g)
{
  StructGrid	*grid = (StructGrid*) g;
  DMDALocalInfo	info;
  VarArray	***phi;
  Integer 	i,istart,iend;
  Integer 	j,jstart,jend;
  Integer 	k,kstart,kend;
  Real 		dcsi,deta,dzeta;
  printf("SetScalarfunction\n");

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

  //printf("\nlocal imax = %d\n",info.xm);
  //printf("\nlocal jmax = %d\n",info.ym);
  //printf("\nlocal kmax = %d\n",info.zm);

  DMDAVecGetArray(da,Scalar,&phi);
  dcsi 		= grid->dcsi;
  deta 		= grid->deta;
  dzeta 	= grid->dzeta;

  double r0=0.35,xc=0.5,yc=0.5;
  double xL=0.45,xR=0.55;
//  double yD=0.750,syc=0.35;

  if (ic == 1) {
    for (i=istart; i<iend; i++) {
    	for (j=jstart; j<jend; j++) {
	    for (k=kstart; k<kend; k++) {

      		Real csi  = (grid->csimin  + i*dcsi);
      		Real eta  = (grid->etamin  + j*deta);
      		Real zeta = (grid->zetamin + k*dzeta);

		/* deformation field case (follows the velocity field)*/
/*
		if ((r0*r0-(csi-xc)*(csi-xc)-(eta-yc)*(eta-yc))>=0) {
        		phi[k][j][i] = 1.0;
      		} 
		else {
        		phi[k][j][i] = 0.0;
      		}
*/
		/*----------------------------------------------------*/
		/* zalesak disk case */

		if (csi>=xL && csi <= xR && eta <=yc ) {
        		phi[k][j][i] =0.0;
      		} 		
		else if ((r0*r0-(csi-xc)*(csi-xc)-(eta-yc)*(eta-yc))>=0) {
        		phi[k][j][i] = 1.0;
      		} 
		else {
        		phi[k][j][i] = 0.0;
      		}

		/*----------------------------------------------------*/
		/*          box case        */
/*
      		if (csi< 0.6 && csi>0.4 && eta< 0.35 && eta>0.15 ) {
        		phi[k][j][i] = 1.0;
      		} 
		else {
        		phi[k][j][i] = 0.0;
      		}
*/
		/*----------------------------------------------------*/
    	  }
       }
     }
  }



/* need to make 2d 
  else if (params->ic == -1) {
    Integer i;
    for (i=istart; i<iend; i++) {
      Real x = (xmin + i*dx) - params->ax*t;
      if (x <= -0.8) {
        u[i] = 0;
      } else if (x <= -0.6) {
        u[i] = exp(-log(2.0)*(x+0.7)*(x+0.7)/0.0009);
      } else if (x <= -0.4) {
        u[i] = 0;
      } else if (x <= -0.2) {
        u[i] = 1.0;
      } else if (x <= 0) {
        u[i] = 0;
      } else if (x <= 0.2) {
        u[i] = 1-absolute(10*(x-0.1));
      } else if (x <= 0.4) {
        u[i] = 0;
      } else if (x <= 0.6) {
        u[i] = sqrt(1-100*(x-0.5)*(x-0.5));
      } else {
        u[i] = 0;
      }
    }
  }
  else if (params->ic == -2) {
    Integer i;
    for (i=istart; i<iend; i++) {
      Real x = (xmin + i*dx) - params->ax*t;
      if (x <= 0.4) {
        u[i] = 0;
      } else if (x <= 0.6) {
        u[i] = 1.0;
      } else {
        u[i] = 0;
      }
    }
  }
*/
  else{
    printf("Error (SetScalarFunction): Invalid choice of initial conditions.\n");
    return(1);
  }

  DMDAVecRestoreArray(da,Scalar,&phi);
  return(0);
}
