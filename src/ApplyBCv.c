#include <stdio.h>
#include <DataTypes.h>

ErrorType ApplyBCv(DM da,Vec v, Vec vl,void *g)
{
  StructGrid    	*grid	= (StructGrid*)	g;
  VarArray 		***varr, ***vlarr;

  Vec 			pLoc;
  Integer 		i, is, ie;
  Integer 		j, js, je;
  Integer 		k, ks, ke;

  Integer 		LeftB  ,RightB;
  Integer 		BottomB,TopB  ;
  Integer 		FrontB ,BackB ;


  LeftB		=	grid->LeftB	;
  RightB 	=	grid->RightB	;
  BottomB 	=	grid->BottomB	;
  TopB		=	grid->TopB	;
  FrontB 	=	grid->FrontB	;
  BackB		=	grid->BackB	;

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
  
  DMGetLocalVector(da,&pLoc);
  VecZeroEntries(pLoc);
  DMGlobalToLocalBegin(da,v,INSERT_VALUES,pLoc);
  DMGlobalToLocalEnd  (da,v,INSERT_VALUES,pLoc);  


  DMDAVecGetArray(da,pLoc,&varr);
  DMDAVecGetArray(da,vl,&vlarr);


  for (i=is; i<ie; i++) {
     for (j=js; j<je; j++) {
	for (k=ks; k<ke; k++) {
		//Periodic case
//		{farr[k][j][i] = ucarr[k][j][i]*parr[k][j][i];}	
		//Dirichlet BC	
    		if(i>=LeftB 	&& 	i<=RightB 
		&& j>=BottomB  	&& 	j<=TopB
		&& k>=FrontB   	&& 	k<=BackB   ){vlarr[k][j][i] = varr[k][j][i];}

    		if(i==LeftB-1		 	   ){vlarr[k][j][i] = - varr[k][j][LeftB+0];}
    		if(i==LeftB-2		 	   ){vlarr[k][j][i] = - varr[k][j][LeftB+1];}
    		if(i==LeftB-3		 	   ){vlarr[k][j][i] = - varr[k][j][LeftB+2];}

    		if(i==RightB+1		 	   ){vlarr[k][j][i] = - varr[k][j][RightB-0];}
    		if(i==RightB+2		 	   ){vlarr[k][j][i] = - varr[k][j][RightB-1];}
    		if(i==RightB+3		 	   ){vlarr[k][j][i] = - varr[k][j][RightB-2];}


    		if(j==BottomB-1		 	   ){vlarr[k][j][i] = - varr[k][BottomB+0][i];}
    		if(j==BottomB-2		 	   ){vlarr[k][j][i] = - varr[k][BottomB+1][i];}
    		if(j==BottomB-3		 	   ){vlarr[k][j][i] = - varr[k][BottomB+2][i];}

    		if(j==TopB+1 		 	   ){vlarr[k][j][i] = - varr[k][TopB-0][i];}
    		if(j==TopB+2		 	   ){vlarr[k][j][i] = - varr[k][TopB-1][i];}
    		if(j==TopB+3		 	   ){vlarr[k][j][i] = - varr[k][TopB-2][i];}

    		if(k==FrontB-1		 	   ){vlarr[k][j][i] = - varr[FrontB+0][j][i];}
    		if(k==FrontB-2		 	   ){vlarr[k][j][i] = - varr[FrontB+1][j][i];}
    		if(k==FrontB-3		 	   ){vlarr[k][j][i] = - varr[FrontB+2][j][i];}

    		if(k==BackB+1		 	   ){vlarr[k][j][i] = - varr[BackB-0][j][i];}
    		if(k==BackB+2		 	   ){vlarr[k][j][i] = - varr[BackB-1][j][i];}
    		if(k==BackB+3		 	   ){vlarr[k][j][i] = - varr[BackB-2][j][i];}

	}
     }
  }

  DMDAVecRestoreArray(da,vl,&vlarr);
  DMDAVecRestoreArray(da,pLoc,&varr);
  DMRestoreLocalVector(da,&pLoc);


  return(0);
}
