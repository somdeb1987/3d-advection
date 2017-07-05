#include <stdio.h>
#include <DataTypes.h>

ErrorType CalculateEtaUpWind(DM da,Vec vcLoc, Vec upw)
{
  VarArray 		***vcarr, ***upwarr;
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
  DMDAVecGetArray(da,vcLoc,&vcarr);



  for (i=is; i<ie; i++) {
	for (j=js; j<je+1; j++) {
		for (k=ks; k<ke; k++) {
			//Periodic case
        		double vL = vcarr[k][j-1][i];
        		double vR = vcarr[k][j][i];
        		if      ((vL > 0) && (vR > 0)) upwarr[k][j][i] =  1.0;
        		else if ((vL < 0) && (vR < 0)) upwarr[k][j][i] = -1.0;
        		else                           upwarr[k][j][i] =  0.0;
		}
	}
  }

  DMDAVecRestoreArray(da,upw,&upwarr);
  DMDAVecRestoreArray(da,vcLoc,&vcarr);
  return(0);
}
