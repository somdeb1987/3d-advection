#include <stdio.h>
#include <DataTypes.h>

ErrorType CalculateZetaUpWind(DM da,Vec wcLoc, Vec upw)
{
  VarArray 		***wcarr, ***upwarr;
  DMDALocalInfo		info;
  Integer 		i, is, ie;
  Integer 		j, js, je;
  Integer 		k, ks, ke;

/*
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
  DMDAVecGetArray(da,wcLoc,&wcarr);



  for (i=is; i<ie; i++) {
	for (j=js; j<je; j++) {
		for (k=ks; k<ke+1; k++) {
			//Periodic case
        		double wL = wcarr[k][j-1][i];
        		double wR = wcarr[k][j][i];
        		if      ((wL > 0) && (wR > 0)) upwarr[k][j][i] =  1.0;
        		else if ((wL < 0) && (wR < 0)) upwarr[k][j][i] = -1.0;
        		else                           upwarr[k][j][i] =  0.0;
		}
	}
  }

  DMDAVecRestoreArray(da,upw,&upwarr);
  DMDAVecRestoreArray(da,wcLoc,&wcarr);

  return(0);
}
