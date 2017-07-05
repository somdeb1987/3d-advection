#include <stdio.h>
#include <DataTypes.h>

ErrorType CalculateCsiUpWind(DM da,Vec ucLoc, Vec upw)
{
  VarArray 		***ucarr, ***upwarr;

  Integer 		i, is, ie;
  Integer 		j, js, je;
  Integer 		k, ks, ke;

/*
  DMDALocalInfo		info;
  DMDAGetLocalInfo(da,&info);
  is = info.xs;
  ie = info.xs + info.xm;
  js = info.ys;
  je = info.ys + info.ym;
  ks = info.zs;
  ke = info.zs + info.zm;
*/

  Integer xs,xm,ys,ym,zs,zm;
  DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  is=xs;	ie=xs+xm;
  js=ys;	je=ys+ym;
  ks=zs;	ke=zs+zm;



  DMDAVecGetArray(da,upw,&upwarr);
  DMDAVecGetArray(da,ucLoc,&ucarr);



  for (i=is; i<ie+1; i++) {
	for (j=js; j<je; j++) {
		for (k=ks; k<ke; k++) {
			//Periodic case
        		double uL = ucarr[k][j][i-1];
        		double uR = ucarr[k][j][i];
        		if      ((uL > 0) && (uR > 0)) upwarr[k][j][i] =  1.0;
        		else if ((uL < 0) && (uR < 0)) upwarr[k][j][i] = -1.0;
        		else                           upwarr[k][j][i] =  0.0;

		}
	}
  }

  DMDAVecRestoreArray(da,upw,&upwarr);
  DMDAVecRestoreArray(da,ucLoc,&ucarr);

  return(0);
}
