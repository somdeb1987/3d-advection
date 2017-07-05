#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Functions.h>
#include <Definitions.h>
#include <DataTypes.h>

ErrorType Initialize(void *s,void *io,void *TSP,void *g,void *param)
{

  StructINS3D		*Solver		= (StructINS3D*)	s;
  IOParameters    	*IO	   	= (IOParameters*) 	io;
  TSParameters		*TSParam	= (TSParameters*)	TSP;
  StructGrid		*Grid		= (StructGrid*) 	g;
  WENOParameters    	*WENO	   	= (WENOParameters*) 	param;

  /* Basic options */
  //solver.ax	= 1.0;
  //solver.ay	= 1.0;
  Solver->spatial_scheme	= 5;
  Solver->ic			= 1;
/*
  if (Solver->spatial_scheme == _FIRST_ORDER_UPWIND_) {
    Solver->ReconstructionRHS = FirstOrderUpwindRHS;
  } 
  else if (Solver->spatial_scheme == _SECOND_ORDER_CENTRAL_) {
    Solver->ReconstructionRHS = SecondOrderCentralRHS;
  } 
  else if (Solver->spatial_scheme == _FIFTH_ORDER_WENO_) {
    Solver->ReconstructionRHS = FifthOrderWENORHS;
  } 

  else {
    printf("Error: in InitialSetup - invalid choice of spatial_scheme.\n");
    return(1);
  }
*/

  /* IO options */
  strcpy(IO->initial_fname,"initial_solution");
  strcpy(IO->exact_fname,"exact_solution");
  strcpy(IO->final_fname,"final_solution");
  IO->output_interval 	=	100;
  IO->counter 		=	0;

  /* Time-iterations */
  TSParam->T		= 1800;
  TSParam->tf		= 18.84;
  TSParam->dt 	      	= TSParam->tf / TSParam->T;
//  TSParam->cfl 	      = p->ax*p->dt / p->dx;

  /* Grid options */
  Grid->Nx	= 50;
  Grid->Ny	= 50;
  Grid->Nz	= 50;
  Grid->csimin	= 0.0;
  Grid->csimax	= 1.0;
  Grid->etamin	= 0.0;
  Grid->etamax	= 1.0;
  Grid->zetamin	= 0.0;
  Grid->zetamax	= 1.0;

  Grid->LeftB  	= 0;
  Grid->RightB 	= Grid->Nx-1;
  Grid->BottomB	= 0;
  Grid->TopB 	= Grid->Ny-1;
  Grid->FrontB	= 0;
  Grid->BackB	= Grid->Nz-1;

  /* WENO options */
  WENO->borges		= PETSC_TRUE;
  WENO->yc		= PETSC_FALSE;
  WENO->no_limiting	= PETSC_FALSE;
  WENO->mapped		= PETSC_FALSE;
  WENO->eps		= 1e-6;
  WENO->p		= 2;

  return(0);
}
