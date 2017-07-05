#include <stdio.h>
#include <DataTypes.h>

ErrorType CalculateUContravarient(DM da,Vec u, Vec uc,void* g)
{
  StructGrid  		*grid 	= (StructGrid*)	g;
  VarArray 		***arru, ***arruc;
  DMDALocalInfo		info;
  Integer 		i, is, ie;
  Integer 		j, js, je;
  Integer 		k, ks, ke;

  DMDAVecGetArray(da,u,&arru);
  DMDAVecGetArray(da,uc,&arruc);

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
    		arruc[k][j][i] = arru[k][j][i];
	}
    }
  }

  DMDAVecRestoreArray(da,u,&arru);
  DMDAVecRestoreArray(da,uc,&arruc);

  return(0);
}
