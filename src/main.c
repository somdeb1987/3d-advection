#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Functions.h>
#include <Definitions.h>
#include <DataTypes.h>
#include <PETSCFunctions.h>

static const char help[] = "Linear advection equation in 1D over a periodic domain\n";


ErrorType Initialize(void*,void*,void*,void*,void*);

int main(int argc, char **argv)
{
  StructINS3D		solver;
  ErrorType		ierr;	  	    	/* Error flag */
  DM			da;

  PetscInitialize(&argc,&argv,(char *)0,help);

  MPI_Comm_rank(PETSC_COMM_WORLD,&solver.MPIVars.rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&solver.MPIVars.nproc);



  /* initialize and allocate arrays */
  ierr= Initialize(&solver,&solver.IOParams,&solver.TSParams,&solver.Grid,&solver.WENOParams);	
  if (ierr) return (ierr);
  //ierr= SetupGrid(&solver);		if (ierr) return (ierr);

  int Nx=solver.Grid.Nx;
  int Ny=solver.Grid.Ny;
  int Nz=solver.Grid.Nz;

  ierr= DMDACreate3d(	PETSC_COMM_WORLD,
			DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,
               		DMDA_STENCIL_STAR,              // no diagonal differencing
			Nx,Ny,Nz,
			PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
			_NDOF_,
			_GHOSTSIZE_,
			PETSC_NULL,PETSC_NULL,PETSC_NULL,
			&da);
  if (ierr) return (ierr);

  ierr=DMSetFromOptions(da);						if (ierr) return (ierr);
  ierr=DMSetUp(da);							if (ierr) return (ierr);


  ierr=CreateSVectors(da,&solver.Field);				if (ierr) return (ierr);

  ierr=SetupGrid(da,&solver.Grid);					if (ierr) return (ierr);

  ierr=SetupIC(da,solver.ic,&solver.Grid,&solver.Field);  		if (ierr) return (ierr);



  ierr=SolvePETSc(da,&solver.IOParams,&solver.TSParams, &solver.Grid, &solver.Field,&solver.MPIVars); 
  if (ierr) return (ierr);
 



  ierr=DestroySVectors(da,&solver.Field);				if (ierr) return (ierr);


  ierr=DMDestroy(&da);							if (ierr) return (ierr);
  PetscFinalize();
  return(0);
}
