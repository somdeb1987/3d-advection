#include <DataTypes.h>
#include <PETSCFunctions.h>
#include <Functions.h>

//PetscErrorCode SaveSolution		(void*,int*,double,DM,Vec);

PetscErrorCode PetscPostTimeStep(TS ts)
{
  PETScTSContext    	*context 	= NULL;
  IOParameters		*io		= NULL;
  TSParameters		*params		= NULL;
  StructGrid    	*grid		= NULL;
  StructField   	*field		= NULL;

  DM			da;
  Vec             	Y;
  PetscErrorCode  	ierr;




  ierr   = TSGetApplicationContext(ts, &context); CHKERRQ(ierr);
  ierr   = TSGetDM(ts, &da); CHKERRQ(ierr);

  io		= context->io;
  params 	= context->params;
  grid		= context->grid;
  field  	= context->field;


  ierr = TSGetSolution(ts,&Y); CHKERRQ(ierr);
  ierr = TransferTSToVec(da,field->u,field->v,field->w,field->phi,Y,context);  CHKERRQ(ierr);

  /* ----------------------------solve pressure poisson equation---------------------------------- */

  Mat			Amat_PPE;
  Vec			rhsPPE,phi_press;
  KSP 			ppe_solver;
  PC 			pc;
  MatNullSpace 		nullspace;
  PetscInt 		its;

  //create coefficient matrix
  ierr = DMSetMatType(da,MATMPIAIJ);						CHKERRQ(ierr);
  ierr = DMCreateMatrix(da,&Amat_PPE); 						CHKERRQ(ierr);
  ierr = MatSetFromOptions(Amat_PPE); 						CHKERRQ(ierr);
  ierr = MatZeroEntries(Amat_PPE);						CHKERRQ(ierr);

  //create solution vector and rhs
  ierr = DMCreateGlobalVector(da,&rhsPPE); 					CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&phi_press); 					CHKERRQ(ierr);

  //Compute the coefficient matrix
  ierr = ComputeMatrixPPE(da,Amat_PPE,grid);					CHKERRQ(ierr);

  //compute rhs
  ierr = ComputeRHSPPE(da,field->u,field->v,field->w,rhsPPE,params->dt, grid);	CHKERRQ(ierr);
  
  //create ksp to solve linear system and set associated A matrix
  ierr = KSPCreate(PETSC_COMM_WORLD, &ppe_solver);				CHKERRQ(ierr);
  ierr = KSPSetOperators(ppe_solver, Amat_PPE, Amat_PPE); 			CHKERRQ(ierr);
  //setup preconditioner
  ierr = KSPGetPC(ppe_solver,&pc);
  ierr = PCSetType(pc,PCBJACOBI);
  //initial guess
  ierr = KSPSetInitialGuessNonzero(ppe_solver, PETSC_TRUE);			CHKERRQ(ierr);


  //setup ksp solver
  ierr = KSPSetTolerances(ppe_solver,1.e-9,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); 	CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ppe_solver); 						CHKERRQ(ierr);
  ierr = KSPSetUp(ppe_solver); 								CHKERRQ(ierr);
  ierr = KSPSetType(ppe_solver, KSPPGMRES);						CHKERRQ(ierr);

  //solve the linear system
  ierr = KSPSolve(ppe_solver,rhsPPE,phi_press); 					CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ppe_solver,&its); 					CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"KSP: %d\n",its); 				CHKERRQ(ierr);

  //detailed analysis of the solver
  ierr = KSPView(ppe_solver,PETSC_VIEWER_STDOUT_WORLD);					CHKERRQ(ierr);




  ierr = KSPDestroy(&ppe_solver); 		CHKERRQ(ierr);
  ierr = MatDestroy(&Amat_PPE); 		CHKERRQ(ierr);
  ierr = VecDestroy(&rhsPPE); 			CHKERRQ(ierr);
  ierr = VecDestroy(&phi_press); 		CHKERRQ(ierr);



  /* -----------------------------------------  done   ------------------------------------------- */
  int iter; 
  ierr = TSGetTimeStepNumber(ts,&iter);

  // Write intermediate solution to file 
  if ((iter+1)%io->output_interval == 0) {
    ierr = SaveSolution(grid,&io->counter, iter*params->dt, da,field->u,field->v,field->w,field->phi); CHKERRQ(ierr);
	    }

  /* copy solution back to PETSc vector */
  ierr = TransferVecToTS(da,field->u,field->v,field->w,field->phi,Y,context);
  ierr = TSSetSolution(ts,Y); CHKERRQ(ierr);
  if (ierr) return(ierr);




  return(0);
}

