#include <stdio.h>
#include <DataTypes.h>

ErrorType CalculateEtaFlux(DM da,Vec pLoc, Vec vcLoc, Vec feta,void *g)
{
  StructGrid    	*grid	= (StructGrid*)	g;
  VarArray 		***parr, ***vcarr, ***farr;
  DMDALocalInfo		info;
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

/*
  DMDAGetLocalInfo(da,&info);
  is = info.gxs;
  ie = info.gxs + info.gxm;
  js = info.gys;
  je = info.gys + info.gym;
  ks = info.gzs;
  ke = info.gzs + info.gzm;
*/
  Integer xs,xm,ys,ym,zs,zm;
  DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  is=xs;	ie=xs+xm;
  js=ys;	je=ys+ym;
  ks=zs;	ke=zs+zm;


  DMDAVecGetArray(da,pLoc,&parr);
  DMDAVecGetArray(da,vcLoc,&vcarr);
  DMDAVecGetArray(da,feta   ,&farr);





  for (i=is; i<ie; i++) {
     for (j=js; j<je; j++) {
	for (k=ks; k<ke; k++) {
		//Periodic case
		{farr[k][j][i] = vcarr[k][j][i]*parr[k][j][i];}	

	}
     }
  }

  DMDAVecRestoreArray(da,pLoc,&parr);
  DMDAVecRestoreArray(da,vcLoc,&vcarr);
  DMDAVecRestoreArray(da,feta,&farr);
  return(0);
}
