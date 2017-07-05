#include <math.h>
#include <DataTypes.h>

double CalculateErrorL2(DM da, Vec X0, Vec Xf )
{
  int     		i, is, ie;
  int     		j, js, je;
  PetscReal  		**x0,**xf;
  DMDALocalInfo 	info;
  double		err,l2;
  int 			num;

  DMDAGetLocalInfo(da,&info);
  is = info.xs;
  ie = info.xs + info.xm;
  js = info.ys;
  je = info.ys + info.ym;

  DMDAVecGetArray(da,X0,&x0);
  DMDAVecGetArray(da,Xf,&xf);

  err=0;
  num =0;
  
  for (i=is; i<ie+1; i++) {
	  for (j=js; j<je; j++) {
		err= err+(xf[j][i]-x0[j][i])*(xf[j][i]-x0[j][i]);
		num++;

	  }
  }
  printf("\nL2 Error = % lf\n",sqrt(err));
  l2=sqrt(err/num);
  DMDAVecRestoreArray(da,X0,&x0);
  DMDAVecRestoreArray(da,Xf,&xf);

  return(l2);
}
