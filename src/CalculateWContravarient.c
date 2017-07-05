#include <stdio.h>
#include <DataTypes.h>

ErrorType CalculateWContravarient(DM da,Vec w, Vec wc,void* g)
{
  StructGrid  		*grid 	= (StructGrid*)	g;
  VarArray 		***arrw, ***arrwc;
  DMDALocalInfo		info;
  Integer 		i, is, ie;
  Integer 		j, js, je;
  Integer 		k, ks, ke;

  DMDAVecGetArray(da,w,&arrw);
  DMDAVecGetArray(da,wc,&arrwc);

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
    		arrwc[k][j][i] = arrw[k][j][i];
	}
    }
  }

  DMDAVecRestoreArray(da,w,&arrw);
  DMDAVecRestoreArray(da,wc,&arrwc);

  return(0);
}
