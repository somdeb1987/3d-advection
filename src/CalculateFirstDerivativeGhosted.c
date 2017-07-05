#include <stdio.h>
#include <DataTypes.h>

ErrorType CalculateFirstDerivativeGhosted(DM da,Vec pLoc, Vec DpDcsi, Vec DpDeta, Vec DpDzeta,void *g)
{
  StructGrid    	*grid	= (StructGrid*)	g;
  VarArray 		***parr;
  VarArray		***ddcsi, ***ddeta,***ddzeta;
  double		dx,dy,dz;
  double  		dxinv,dyinv,dzinv;
  Integer 		i, is, ie;
  Integer 		j, js, je;
  Integer 		k, ks, ke;
/*
  Integer 		LeftB  ,RightB;
  Integer 		BottomB,TopB  ;
  Integer 		FrontB ,BackB ;


  LeftB		=	grid->LeftB	;
  RightB 	=	grid->RightB	;
  BottomB 	=	grid->BottomB	;
  TopB		=	grid->TopB	;
  FrontB 	=	grid->FrontB	;
  BackB		=	grid->BackB	;
*/

  dx	=	grid->dcsi;
  dy	=	grid->deta; 
  dz	=	grid->dzeta;
  dxinv = 1.0/dx;
  dyinv = 1.0/dy;
  dzinv = 1.0/dz;


  Integer xs,xm,ys,ym,zs,zm;
  DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  is=xs;	ie=xs+xm;
  js=ys;	je=ys+ym;
  ks=zs;	ke=zs+zm;


  DMDAVecGetArray(da,pLoc,&parr);
  DMDAVecGetArray(da,DpDcsi,&ddcsi);
  DMDAVecGetArray(da,DpDeta,&ddeta);
  DMDAVecGetArray(da,DpDzeta,&ddzeta);


  for (i=is-1; i<ie+1; i++) {
     for (j=js-1; j<je+1; j++) {
	for (k=ks-1; k<ke+1; k++) {
		//Periodic case
		{
		ddcsi[k][j][i] 	= 0.5*dxinv*(parr[k][j][i+1]-parr[k][j][i-1]);
		ddeta[k][j][i] 	= 0.5*dyinv*(parr[k][j+1][i]-parr[k][j-1][i]);
		ddzeta[k][j][i] = 0.5*dzinv*(parr[k+1][j][i]-parr[k-1][j][i]);
		}	

	}
     }
  }

  DMDAVecRestoreArray(da,pLoc,&parr);
  DMDAVecRestoreArray(da,DpDcsi,&ddcsi);
  DMDAVecRestoreArray(da,DpDeta,&ddeta);
  DMDAVecRestoreArray(da,DpDzeta,&ddzeta);

  return(0);
}
