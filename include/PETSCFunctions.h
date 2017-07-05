/* include PETSc header files */
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscts.h>


/**************************/
/* petsc time integration */
/**************************/
typedef struct _petsctscontext_ {
  void *io;
  void *params;
  void *grid;
  void *field;
  void *mpi;
  void *weno;
  void *solveLS;

  int size_u, size_v, size_w, size_phi; 
  int offset_u, offset_v, offset_w, offset_phi;

  PetscReal shift;

} PETScTSContext;



/* Copy Functions */
int TransferVecToTS	(DM,Vec,Vec,Vec,Vec,Vec,void*);
int TransferTSToVec	(DM,Vec,Vec,Vec,Vec,Vec,void*);




/* Register custom time-integration RK/ARKIMEX method */
//int PetscRegisterTIMethods(char*,int);

/* Right and left -hand side functions */
PetscErrorCode PetscRHSFunctionExpl(TS,PetscReal,Vec,Vec,void*);
//PetscErrorCode PetscRHSFunctionIMEX(TS,PetscReal,Vec,Vec,void*);
//PetscErrorCode PetscIFunctionIMEX(TS,PetscReal,Vec,Vec,Vec,void*);

/* Jacobian functions for left-hand side */
//PetscErrorCode PetscIJacobianIMEX(TS,PetscReal,Vec,Vec,PetscReal,Mat*,Mat*,MatStructure*,void*);                            
//PetscErrorCode PetscJacobianFunctionIMEX(Mat,Vec,Vec);             

/* Other functions */
PetscErrorCode PetscPreStage(TS,PetscReal);
PetscErrorCode PetscPostStage(TS,PetscReal,PetscInt,Vec*);
PetscErrorCode PetscPostTimeStep(TS);

