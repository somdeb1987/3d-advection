#include <DataTypes.h>
#include <PETSCFunctions.h>
#include <Functions.h>

ErrorType PetscRHSFunctionExpl(TS ts, PetscReal t, Vec Y, Vec F, void *ctxt)
{
  PETScTSContext    	*context 	= (PETScTSContext*) 	ctxt;
  StructGrid    	*grid		= (StructGrid*)		context->grid;
  StructField   	*field		= (StructField*) 	context->field;
  IOParameters		*io		= (IOParameters*) 	context->io;

  DM 			da;
  PetscErrorCode  	ierr;

  Vec			Local_u,Local_v,Local_w,Local_phi;


  TSGetDM(ts,&da);

  /* Boundary condition for local */

  DMCreateLocalVector(da,&Local_u);
  DMCreateLocalVector(da,&Local_v);
  DMCreateLocalVector(da,&Local_w);
  DMCreateLocalVector(da,&Local_phi);



  /* copy solution from PETSc vector */
  ierr = TransferTSToVec(da,field->u,field->v,field->w,field->phi,Y,context);
  if (ierr) return(ierr);


  /* initialize all right-hand side to zero */

  ierr=SetZeroVectors(da,field->Au);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Av);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Aw);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Aphi);	if (ierr) return (ierr);

  ierr=SetZeroVectors(da,field->Du);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Dv);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Dw);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Dphi);	if (ierr) return (ierr);

  ierr=SetZeroVectors(da,field->Su);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Sv);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Sw);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->Sphi);	if (ierr) return (ierr);

  ierr=SetZeroVectors(da,field->RHSu);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->RHSv);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->RHSw);	if (ierr) return (ierr);
  ierr=SetZeroVectors(da,field->RHSphi);if (ierr) return (ierr);



  // Apply Boundary Conditions //
  ierr = ApplyBCu(da,field->u,Local_u,grid); 		if (ierr) return (ierr);
  ierr = ApplyBCv(da,field->v,Local_v,grid); 		if (ierr) return (ierr);
  ierr = ApplyBCw(da,field->w,Local_w,grid);		if (ierr) return (ierr);
  ierr = ApplyBCphi(da,field->phi,Local_phi,grid);	if (ierr) return (ierr);

  
  // Calculate Advection Terms //
  ierr = CalculateAdvectionTerm(da,Local_u,Local_v,Local_w,Local_u  ,field->Au  ,grid);
  if (ierr) return (ierr); 
  ierr = CalculateAdvectionTerm(da,Local_u,Local_v,Local_w,Local_v  ,field->Av  ,grid);
  if (ierr) return (ierr);
  ierr = CalculateAdvectionTerm(da,Local_u,Local_v,Local_w,Local_w  ,field->Aw  ,grid);
  if (ierr) return (ierr);
  ierr = CalculateAdvectionTerm(da,Local_u,Local_v,Local_w,Local_phi,field->Aphi,grid);
  if (ierr) return (ierr);



  // Calculate Diffusion Terms //
  
  PetscScalar coef_nu	= 0.001;
  PetscScalar coef_k	= 1.0;

  ierr = CalculateDiagonalDiffusionTerm(da,Local_u,field->Du,grid, coef_nu);
  if (ierr) return (ierr);
  ierr = CalculateDiagonalDiffusionTerm(da,Local_v,field->Dv,grid, coef_nu);
  if (ierr) return (ierr);
  ierr = CalculateDiagonalDiffusionTerm(da,Local_w,field->Dw,grid, coef_nu);
  if (ierr) return (ierr);
  ierr = CalculateDiagonalDiffusionTerm(da,Local_phi,field->Dphi,grid, coef_k);
  if (ierr) return (ierr);





  //ierr = CalculateDensityProfile(da,Local_phi,rho,grid);
  //ierr = CalculateViscosityProfile(da,Local_phi,mu,grid);





  ierr = VecAYPX(field->RHSu,1.0,field->Au);
  ierr = VecAYPX(field->RHSu,1.0,field->Du);

  ierr = VecAYPX(field->RHSv,1.0,field->Av);
  ierr = VecAYPX(field->RHSv,1.0,field->Dv);

  ierr = VecAYPX(field->RHSw,1.0,field->Aw);
  ierr = VecAYPX(field->RHSw,1.0,field->Dw);




 

  /* Transfer RHS to TS */
  ierr = TransferVecToTS(da,field->RHSu,field->RHSv,field->RHSw,field->RHSphi,F,context);
  if (ierr) return(ierr);


  /* Clean up */

  VecDestroy(&Local_u);
  VecDestroy(&Local_v);
  VecDestroy(&Local_w);
  VecDestroy(&Local_phi);

  return(0);
}

