#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscdmda.h>

PetscErrorCode FirstOrderUpwindRHS		(Vec,Vec, Vec, void*, double,double);
PetscErrorCode SecondOrderCentralRHS		(Vec,Vec, Vec, double,double);
PetscErrorCode FifthOrderWENORHS		(Vec,Vec,Vec,Vec,Vec,Vec,Vec,void*);
PetscErrorCode SetWENOParams			(void*);
PetscErrorCode WENOWeights			(double,double,double,double,double,double,double,double);

PetscErrorCode CreateSVectors			(DM , void*);
PetscErrorCode DestroySVectors			(DM , void*);

PetscErrorCode SetupGrid			(DM ,void*);
PetscErrorCode SetupIC				(DM ,int,void*,void*);
PetscErrorCode SetScalarField			(DM ,int,Vec,void*);
PetscErrorCode SetVelocityField			(DM ,int,Vec,Vec,Vec,void*);
PetscErrorCode SetZeroVectors			(DM ,Vec);
PetscErrorCode SolvePETSc			(DM ,void*,void*,void*,void*,void*);
PetscErrorCode SaveSolution			(void*,int*,double,DM,Vec,Vec,Vec,Vec);
double	       CalculateErrorL2			(DM,Vec,Vec);

PetscErrorCode ApplyBCu				(DM,Vec,Vec,void*);
PetscErrorCode ApplyBCv				(DM,Vec,Vec,void*);
PetscErrorCode ApplyBCw				(DM,Vec,Vec,void*);
PetscErrorCode ApplyBCphi			(DM,Vec,Vec,void*);

PetscErrorCode CalculateAdvectionTerm		(DM,Vec,Vec,Vec,Vec,Vec,void*);

PetscErrorCode CalculateDiagonalDiffusionTerm	(DM,Vec,Vec,void*,PetscScalar);
PetscErrorCode CalculateFirstDerivativeGhosted	(DM,Vec,Vec,Vec,Vec,void*);
PetscErrorCode DDcsi				(DM,Vec,Vec,void*);
PetscErrorCode DDeta				(DM,Vec,Vec,void*);
PetscErrorCode DDzeta				(DM,Vec,Vec,void*);


PetscErrorCode CalculateUContravarient		(DM,Vec,Vec,void*);
PetscErrorCode CalculateVContravarient		(DM,Vec,Vec,void*);
PetscErrorCode CalculateWContravarient		(DM,Vec,Vec,void*);

PetscErrorCode CalculateCsiUpWind		(DM,Vec,Vec);
PetscErrorCode CalculateEtaUpWind		(DM,Vec,Vec);
PetscErrorCode CalculateZetaUpWind		(DM,Vec,Vec);

PetscErrorCode CalculateCsiFlux			(DM,Vec,Vec,Vec,void*);
PetscErrorCode CalculateEtaFlux			(DM,Vec,Vec,Vec,void*);
PetscErrorCode CalculateZetaFlux		(DM,Vec,Vec,Vec,void*);


PetscErrorCode ComputeMatrixPPE			(DM, Mat, void*);
PetscErrorCode ComputeUxRHSPPE			(DM, Vec, Vec, void*);
PetscErrorCode ComputeVyRHSPPE			(DM, Vec, Vec, void*);
PetscErrorCode ComputeWzRHSPPE			(DM, Vec, Vec, void*);
PetscErrorCode ComputeRHSPPE			(DM, Vec, Vec, Vec, Vec, double, void*);



PetscErrorCode FifthOrderCRWENORHS		(Vec,Vec,void*,double);
PetscErrorCode FifthOrderCRWENOLHS		(Mat,void*,int,int);
