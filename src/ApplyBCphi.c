#include <stdio.h>
#include <DataTypes.h>

ErrorType ApplyBCphi(DM da,Vec phi, Vec phil,void *g)
{
  StructGrid    	*grid	= (StructGrid*)	g;
  VarArray 		***phiarr, ***philarr;

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
  DMGlobalToLocalBegin(da,phi,INSERT_VALUES,pLoc);
  DMGlobalToLocalEnd  (da,phi,INSERT_VALUES,pLoc);





  DMDAVecGetArray(da,pLoc,&phiarr);
  DMDAVecGetArray(da,phil,&philarr);


  for (i=is; i<ie; i++) {
     for (j=js; j<je; j++) {
	for (k=ks; k<ke; k++) {
		//Periodic case
//		{farr[k][j][i] = ucarr[k][j][i]*parr[k][j][i];}	
		//Dirichlet BC	
    		if(i>=LeftB 	&& 	i<=RightB 
		&& j>=BottomB  	&& 	j<=TopB
		&& k>=FrontB   	&& 	k<=BackB   ){philarr[k][j][i] = phiarr[k][j][i];}

    		if(i==LeftB-1		 	   ){philarr[k][j][i] =  phiarr[k][j][LeftB+0];}
    		if(i==LeftB-2		 	   ){philarr[k][j][i] =  phiarr[k][j][LeftB+1];}
    		if(i==LeftB-3		 	   ){philarr[k][j][i] =  phiarr[k][j][LeftB+2];}

    		if(i==RightB+1		 	   ){philarr[k][j][i] =  phiarr[k][j][RightB-0];}
    		if(i==RightB+2		 	   ){philarr[k][j][i] =  phiarr[k][j][RightB-1];}
    		if(i==RightB+3		 	   ){philarr[k][j][i] =  phiarr[k][j][RightB-2];}


    		if(j==BottomB-1		 	   ){philarr[k][j][i] =  phiarr[k][BottomB+0][i];}
    		if(j==BottomB-2		 	   ){philarr[k][j][i] =  phiarr[k][BottomB+1][i];}
    		if(j==BottomB-3		 	   ){philarr[k][j][i] =  phiarr[k][BottomB+2][i];}

    		if(j==TopB+1 		 	   ){philarr[k][j][i] =  phiarr[k][TopB-0][i];}
    		if(j==TopB+2		 	   ){philarr[k][j][i] =  phiarr[k][TopB-1][i];}
    		if(j==TopB+3		 	   ){philarr[k][j][i] =  phiarr[k][TopB-2][i];}

    		if(k==FrontB-1		 	   ){philarr[k][j][i] =  phiarr[FrontB+0][j][i];}
    		if(k==FrontB-2		 	   ){philarr[k][j][i] =  phiarr[FrontB+1][j][i];}
    		if(k==FrontB-3		 	   ){philarr[k][j][i] =  phiarr[FrontB+2][j][i];}

    		if(k==BackB+1		 	   ){philarr[k][j][i] =  phiarr[BackB-0][j][i];}
    		if(k==BackB+2		 	   ){philarr[k][j][i] =  phiarr[BackB-1][j][i];}
    		if(k==BackB+3		 	   ){philarr[k][j][i] =  phiarr[BackB-2][j][i];}

	}
     }
  }

  DMDAVecRestoreArray(da,phil,&philarr);
  DMDAVecRestoreArray(da,pLoc,&phiarr);
  DMRestoreLocalVector(da,&pLoc);


  return(0);
}
