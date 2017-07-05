#include <stdio.h>
#include <Functions.h>
#include <DataTypes.h>

ErrorType SetWENOParams(void *p)
{
  WENOParameters    	*WENO	   	= (WENOParameters*) 	p;

  /* WENO options */
  WENO->borges		= PETSC_TRUE;
  WENO->yc		= PETSC_FALSE;
  WENO->no_limiting	= PETSC_FALSE;
  WENO->mapped		= PETSC_FALSE;
  WENO->eps		= 1e-6;
  WENO->p		= 2;

  return(0);
}
