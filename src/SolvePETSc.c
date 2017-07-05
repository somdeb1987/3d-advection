#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <DataTypes.h>
#include <PETSCFunctions.h>
#include <Functions.h>

ErrorType SolvePETSc(DM da,void *io, void *TSP, void *g, void *f, void *m)
{
  IOParameters    	*IO	= (IOParameters*) 	io;
  TSParameters		*params	= (TSParameters*)	TSP;
  StructGrid    	*grid	= (StructGrid*)		g;
  StructField		*field	= (StructField*) 	f;
  MPIVariables    	*mpi   	= (MPIVariables*) 	m;

  DMDALocalInfo		info;
  PetscErrorCode  	ierr    = 0;
  TS              	ts; // time integration object  
  Vec             	Y;  // PETSc solution vector    
  Mat             	A;  // Jacobian matrix          
  TSType          	time_scheme;
  PETScTSContext 	context;
  Vec			init,final;

  ierr = DMCreateGlobalVector(da,&init); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&final); CHKERRQ(ierr);
  ierr = VecCopy(field->phi,init); CHKERRQ(ierr);

  // Write an initial solution file 
  ierr = SaveSolution(grid,&IO->counter,0, da,field->u,field->v,field->w,field->phi); CHKERRQ(ierr);
  // set PETSc context 
  context.io		= IO;
  context.params	= params;
  context.grid 		= grid;
  context.field 	= field;
  context.mpi    	= mpi;


  // Register custom time-integration methods, if specified 
  //ierr = PetscRegisterTIMethods(params->ti_filename,mpi->rank); if (ierr) return(ierr);


//  printf("\nIn context, eps = %d",context.weno->eps);

  // create and initialize PETSc TS solution vector and other parameters 
  //DMDAGetLocalInfo(da,&info);
  //int size=info.xm*info.ym*info.zm;

  Integer xm,ym,zm;
  DMDAGetCorners(da,0,0,0,&xm,&ym,&zm);
  int size=xm*ym*zm;




  context.size_u 	= size;
  context.size_v 	= size;
  context.size_w 	= size;
  context.size_phi 	= size;

  context.offset_u 	= 0;
  context.offset_v	= context.offset_u + context.size_u;
  context.offset_w 	= context.offset_v + context.size_v;
  context.offset_phi 	= context.offset_w + context.size_w;
  int total_size 	= context.size_u + context.size_v
			+ context.size_w + context.size_phi;


  ierr = VecCreate(PETSC_COMM_WORLD,&Y);                                  CHKERRQ(ierr);
  ierr = VecSetSizes(Y,total_size,PETSC_DECIDE);                          CHKERRQ(ierr);
  ierr = VecSetUp(Y);                                                     CHKERRQ(ierr);

  // copy initial solution to PETSc's vector 


 ierr = TransferVecToTS(da,field->u,field->v,field->w,field->phi,Y,&context);  if (ierr) return(ierr);


  // Define and initialize the time-integration object 
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);                                CHKERRQ(ierr);
  ierr = TSSetDM(ts,da);
  ierr = TSSetDuration(ts,params->T,params->dt*params->T);      	CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,0.0,params->dt);                       CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);           CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);                                          CHKERRQ(ierr);

  // Set field->dt = 1.0 so that convection and diffusion subroutines are not dt-scaled 
  //field->dt = 1.0;

  // Define the right and left -hand side functions for each time-integration scheme 
  ierr = TSGetType(ts,&time_scheme);                                      CHKERRQ(ierr);

  if ( (!strcmp(time_scheme,TSEULER)) || (!strcmp(time_scheme,TSSSP)) || (!strcmp(time_scheme,TSRK)) ){
    
    // Convection and diffusion are both explicit 
    ierr = TSSetRHSFunction(ts,PETSC_NULL,PetscRHSFunctionExpl,&context); CHKERRQ(ierr);

  } 
  else if (!strcmp(time_scheme,TSARKIMEX)) {

	printf("\nIMEX NOT IMPLEMENTED\n");
	return(1);

    // Convection is explicit, diffusion is implicit 
    //ierr = TSSetRHSFunction(ts,PETSC_NULL,PetscRHSFunctionIMEX,&context); CHKERRQ(ierr);
    //ierr = TSSetIFunction  (ts,PETSC_NULL,PetscIFunctionIMEX,  &context); CHKERRQ(ierr);
    //ierr = MatCreateShell(PETSC_COMM_WORLD,total_size,total_size,PETSC_DETERMINE,PETSC_DETERMINE,&context,&A);CHKERRQ(ierr);
    //ierr = MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX);CHKERRQ(ierr);
    //ierr = MatSetUp(A);                                                   CHKERRQ(ierr);
    //ierr = TSSetIJacobian(ts,A,A,PetscIJacobianIMEX,&context);            CHKERRQ(ierr);

    // set pre-conditioner to none for MatShell 
    //SNES snes;
    //KSP  ksp;
    //PC   pc;
    //TSGetSNES(ts,&snes);
    //SNESGetKSP(snes,&ksp);
    //KSPGetPC(ksp,&pc);
    //PCSetType(pc,PCNONE);

  } else {
    fprintf(stderr,"Error in SolvePETSc: TSType %s not supported.\n",time_scheme);
    return(1);
  }

  // Set pre/post-stage and post-timestep function 
  ierr = TSSetPreStage(ts,PetscPreStage    );                             CHKERRQ(ierr);
  ierr = TSSetPostStage(ts,PetscPostStage  );                             CHKERRQ(ierr);
  ierr = TSSetPostStep(ts,PetscPostTimeStep);                             CHKERRQ(ierr);
  // Set solution vector for TS 
  ierr = TSSetSolution(ts,Y);                                             CHKERRQ(ierr);
  // Set it all up 
  ierr = TSSetUp(ts);                                                     CHKERRQ(ierr);
  // Set application context 
  ierr = TSSetApplicationContext(ts,&context);                            CHKERRQ(ierr);
 

  if (!mpi->rank) printf("** Starting PETSc time integration **\n");
  ierr = TSSolve(ts,Y);                                                   CHKERRQ(ierr);
  if (!mpi->rank) printf("** Completed PETSc time integration **\n");

  // copy final solution from PETSc's vector 
  //ierr = TransferTSToVec(da,field->u,field->v,field->w,field->phi,Y,&context);

  ierr = TransferTSToVec(da,field->u,field->v,field->w,field->phi,Y,&context);
  if (ierr) return(ierr);
/*
  ierr = VecCopy(field->phi,final); CHKERRQ(ierr);

  double error=CalculateErrorL2(da, init, final );
  printf("\nrms Error = % lf\n",error);


  //ierr = SaveSolution(grid,&IO->counter, params->tf, da,field->u,field->v,field->w,field->phi); CHKERRQ(ierr);
  // clean up 
  if (!strcmp(time_scheme,TSARKIMEX)) {
    //ierr = MatDestroy(&A);                                                CHKERRQ(ierr);
  }
*/
  ierr = TSDestroy(&ts);CHKERRQ(ierr);

  ierr = VecDestroy(&Y);CHKERRQ(ierr);
  ierr = VecDestroy(&init);CHKERRQ(ierr);
  ierr = VecDestroy(&final);CHKERRQ(ierr);

  return(0);
}
