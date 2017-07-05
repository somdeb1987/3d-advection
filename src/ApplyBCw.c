#include <stdio.h>
#include <DataTypes.h>

ErrorType ApplyBCw(DM da,Vec w, Vec wl,void *g)
{
  StructGrid    	*grid	= (StructGrid*)	g;
  VarArray 		***warr, ***wlarr;

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
  DMGlobalToLocalBegin(da,w,INSERT_VALUES,pLoc);
  DMGlobalToLocalEnd  (da,w,INSERT_VALUES,pLoc);  

  DMDAVecGetArray(da,pLoc,&warr);
  DMDAVecGetArray(da,wl,&wlarr);


  for (i=is; i<ie; i++) {
     for (j=js; j<je; j++) {
	for (k=ks; k<ke; k++) {
		//Periodic case
//		{farr[k][j][i] = ucarr[k][j][i]*parr[k][j][i];}	
		//Dirichlet BC	
    		if(i>=LeftB 	&& 	i<=RightB 
		&& j>=BottomB  	&& 	j<=TopB
		&& k>=FrontB   	&& 	k<=BackB   ){wlarr[k][j][i] = warr[k][j][i];}

    		if(i==LeftB-1		 	   ){wlarr[k][j][i] = - warr[k][j][LeftB+0];}
    		if(i==LeftB-2		 	   ){wlarr[k][j][i] = - warr[k][j][LeftB+1];}
    		if(i==LeftB-3		 	   ){wlarr[k][j][i] = - warr[k][j][LeftB+2];}

    		if(i==RightB+1		 	   ){wlarr[k][j][i] = - warr[k][j][RightB-0];}
    		if(i==RightB+2		 	   ){wlarr[k][j][i] = - warr[k][j][RightB-1];}
    		if(i==RightB+3		 	   ){wlarr[k][j][i] = - warr[k][j][RightB-2];}


    		if(j==BottomB-1		 	   ){wlarr[k][j][i] = - warr[k][BottomB+0][i];}
    		if(j==BottomB-2		 	   ){wlarr[k][j][i] = - warr[k][BottomB+1][i];}
    		if(j==BottomB-3		 	   ){wlarr[k][j][i] = - warr[k][BottomB+2][i];}

    		if(j==TopB+1 		 	   ){wlarr[k][j][i] = - warr[k][TopB-0][i];}
    		if(j==TopB+2		 	   ){wlarr[k][j][i] = - warr[k][TopB-1][i];}
    		if(j==TopB+3		 	   ){wlarr[k][j][i] = - warr[k][TopB-2][i];}

    		if(k==FrontB-1		 	   ){wlarr[k][j][i] = - warr[FrontB+0][j][i];}
    		if(k==FrontB-2		 	   ){wlarr[k][j][i] = - warr[FrontB+1][j][i];}
    		if(k==FrontB-3		 	   ){wlarr[k][j][i] = - warr[FrontB+2][j][i];}

    		if(k==BackB+1		 	   ){wlarr[k][j][i] = - warr[BackB-0][j][i];}
    		if(k==BackB+2		 	   ){wlarr[k][j][i] = - warr[BackB-1][j][i];}
    		if(k==BackB+3		 	   ){wlarr[k][j][i] = - warr[BackB-2][j][i];}

	}
     }
  }

  DMDAVecRestoreArray(da,wl,&wlarr);
  DMDAVecRestoreArray(da,pLoc,&warr);
  DMRestoreLocalVector(da,&pLoc);


  return(0);
}
