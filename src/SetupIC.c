#include <stdio.h>
#include <math.h>
#include <DataTypes.h>
#include <Functions.h>

//ErrorType SetScalarFunction(DM,Vec,void*);

ErrorType SetupIC(DM da,int ic, void *g, void *f)
{

  StructGrid		*grid		= (StructGrid*) 	g;
  StructField    	*field	   	= (StructField*) 	f;

  ErrorType ierr;

  printf("SetupIC\n");


  ierr=SetScalarField(da,ic,field->phi,grid);	if (ierr) return (ierr);
  ierr=SetVelocityField(da,ic,field->u,field->v,field->w,grid);	if (ierr) return (ierr);

  ierr=SetZeroVectors(da,field->Au);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Av);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Aphi);	if (ierr) return (ierr);

  ierr=SetZeroVectors(da,field->Du);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Dv);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Dphi);	if (ierr) return (ierr);

  ierr=SetZeroVectors(da,field->Su);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Sv);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Sphi);	if (ierr) return (ierr);

  ierr=SetZeroVectors(da,field->RHSu);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->RHSv);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->RHSphi);if (ierr) return (ierr);

  printf("exit SetupIC\n");




  return(0);
}
