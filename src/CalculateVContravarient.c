#include <stdio.h>
#include <DataTypes.h>

ErrorType CalculateVContravarient(DM da,Vec v, Vec vc,void* g)
{
  StructGrid  		*grid 	= (StructGrid*)	g;
  VarArray 		***arrv, ***arrvc;
  DMDALocalInfo		info;
  Integer 		i, is, ie;
  Integer 		j, js, je;
  Integer 		k, ks, ke;

  DMDAVecGetArray(da,v,&arrv);
  DMDAVecGetArray(da,vc,&arrvc);

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
  DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  is=xs;	ie=xs+xm;
  js=ys;	je=ys+ym;
  ks=zs;	ke=zs+zm;



  for (i=is; i<ie; i++) {
	for (j=js; j<je; j++) {
		for (k=ks; k<ke; k++) {
    			arrvc[k][j][i] = arrv[k][j][i];
		}
	}
  }

  DMDAVecRestoreArray(da,v,&arrv);
  DMDAVecRestoreArray(da,vc,&arrvc);

  return(0);
}
