#include <stdio.h>
#include <DataTypes.h>

ErrorType CalculateCsiFlux(DM da,Vec pLoc, Vec ucLoc, Vec fcsi,void *g)
{
  StructGrid    	*grid	= (StructGrid*)	g;
  VarArray 		***parr, ***ucarr, ***farr;

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

  //printf("\nRightB= %d\n",RightB);
/*
  DMDALocalInfo		info;
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
  DMDAVecGetArray(da,ucLoc,&ucarr);
  DMDAVecGetArray(da,fcsi   ,&farr);





  for (i=is; i<ie; i++) {
     for (j=js; j<je; j++) {
	for (k=ks; k<ke; k++) {
		//Periodic case
		{farr[k][j][i] = ucarr[k][j][i]*parr[k][j][i];}	
	}
     }
  }

  DMDAVecRestoreArray(da,pLoc,&parr);
  DMDAVecRestoreArray(da,ucLoc,&ucarr);
  DMDAVecRestoreArray(da,fcsi,&farr);
  return(0);
}
