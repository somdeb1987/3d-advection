#include <stdio.h>
#include <DataTypes.h>
#include <Functions.h>

ErrorType CalculateDiagonalDiffusionTerm(DM da,	Vec Local_p, Vec globalDp, void *g, PetscScalar alpha)
{
  StructGrid    	*grid	= (StructGrid*)	g;
  PetscErrorCode  	ierr;
  Vec			DpDcsi,DpDeta,DpDzeta;



  /* Contravarient velocity field */

  DMCreateLocalVector(da,&DpDcsi);
  VecZeroEntries(DpDcsi);
  DMCreateLocalVector(da,&DpDeta);
  VecZeroEntries(DpDcsi);
  DMCreateLocalVector(da,&DpDzeta);
  VecZeroEntries(DpDcsi);


  /* Calculate first derivative */
  ierr = CalculateFirstDerivativeGhosted(da,Local_p,DpDcsi,DpDeta,DpDzeta,grid);
  if (ierr) return(ierr);


  VecScale(DpDcsi,alpha);
  VecScale(DpDeta,alpha);
  VecScale(DpDzeta,alpha);


  ierr = DDcsi(da,DpDcsi,globalDp,grid);  	if (ierr) return(ierr);
  ierr = DDeta(da,DpDeta,globalDp,grid);  	if (ierr) return(ierr);
  ierr = DDzeta(da,DpDzeta,globalDp,grid);  	if (ierr) return(ierr);



 /* Finished computing advective term */

  VecDestroy(&DpDcsi);
  VecDestroy(&DpDeta);
  VecDestroy(&DpDzeta);




  return(0);
}
