#include <DataTypes.h>

PetscErrorCode CreateSVectors(DM da, void *f)
{
  StructField		*field		= (StructField*) 	f;
  PetscErrorCode 	ierr;

    ierr = DMCreateGlobalVector(da,&field->u); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->v); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->w); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->phi); CHKERRQ(ierr);


    ierr = DMCreateGlobalVector(da,&field->uc); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->vc); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->wc); CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(da,&field->Au); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->Av); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->Aw); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->Aphi); CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(da,&field->Du); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->Dv); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->Dw); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->Dphi); CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(da,&field->Su); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->Sv); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->Sw); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->Sphi); CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(da,&field->RHSu); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->RHSv); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->RHSw); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&field->RHSphi); CHKERRQ(ierr);


    ierr = DMCreateGlobalVector(da,&field->RHS_PPE); CHKERRQ(ierr);

   ierr = DMSetMatType(da,MATMPIAIJ);
    //ierr = DMSetMatType(da,MATMPIBAIJ);
   ierr = DMCreateMatrix(da,&field->Amat_PPE); CHKERRQ(ierr);





    return(0);
}
