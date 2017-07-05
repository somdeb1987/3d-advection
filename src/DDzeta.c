#include <stdio.h>
#include <DataTypes.h>

ErrorType DDzeta(DM da, Vec p , Vec Dp, void *g)
{
  StructGrid    	*grid	= (StructGrid*)	g;
  VarArray 		***parr,***dparr;
  double		dz;
  double  		dzinv;
  Integer 		i, is, ie;
  Integer 		j, js, je;
  Integer 		k, ks, ke;


  dz	= grid->dzeta;
  dzinv = 1.0/dz;



  Integer xs,xm,ys,ym,zs,zm;
  DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  is=xs;	ie=xs+xm;
  js=ys;	je=ys+ym;
  ks=zs;	ke=zs+zm;


  DMDAVecGetArray(da,p,&parr);
  DMDAVecGetArray(da,Dp,&dparr);



  for (i=is; i<ie; i++) {
     for (j=js; j<je; j++) {
	for (k=ks; k<ke; k++) {
		//Periodic case
		{
		dparr[k][j][i] 	+=  0.5*dzinv*(parr[k+1][j][i]-parr[k-1][j][i]);
		}	

	}
     }
  }

  DMDAVecRestoreArray(da,p,&parr);
  DMDAVecRestoreArray(da,Dp,&dparr);

  return(0);
}
