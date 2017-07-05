#include <stdio.h>
#include <math.h>
#include <DataTypes.h>

ErrorType SetZeroVectors(DM da,Vec q)
{
  DMDALocalInfo	info;
  VarArray	***arrq;
  Integer i,istart,iend;
  Integer j,jstart,jend;
  Integer k,kstart,kend;


  DMDAGetLocalInfo(da,&info);
  istart = info.xs;
  iend = info.xs + info.xm;
  jstart = info.ys;
  jend = info.ys + info.ym;
  kstart = info.zs;
  kend = info.zs + info.zm;

  DMDAVecGetArray(da,q,&arrq);

  for (i=istart; i<iend; i++) {
	for (j=jstart; j<jend; j++) {
	    	for (k=kstart; k<kend; k++) {

   			arrq[k][j][i] = 0.0;

		}
	}
  }

  DMDAVecRestoreArray(da,q,&arrq);

  return(0);
}
