#include <DataTypes.h>

PetscErrorCode DestroySVectors(DM da, void *f)
{
  StructField		*field		= (StructField*) 	f;
  PetscErrorCode 	ierr;

    ierr = VecDestroy(&field->u); CHKERRQ(ierr);
    ierr = VecDestroy(&field->v); CHKERRQ(ierr);
    ierr = VecDestroy(&field->w); CHKERRQ(ierr);
    ierr = VecDestroy(&field->phi); CHKERRQ(ierr);

    ierr = VecDestroy(&field->uc); CHKERRQ(ierr);
    ierr = VecDestroy(&field->vc); CHKERRQ(ierr);
    ierr = VecDestroy(&field->wc); CHKERRQ(ierr);

    ierr = VecDestroy(&field->Au); CHKERRQ(ierr);
    ierr = VecDestroy(&field->Av); CHKERRQ(ierr);
    ierr = VecDestroy(&field->Aw); CHKERRQ(ierr);
    ierr = VecDestroy(&field->Aphi); CHKERRQ(ierr);

    ierr = VecDestroy(&field->Du); CHKERRQ(ierr);
    ierr = VecDestroy(&field->Dv); CHKERRQ(ierr);
    ierr = VecDestroy(&field->Dw); CHKERRQ(ierr);
    ierr = VecDestroy(&field->Dphi); CHKERRQ(ierr);

    ierr = VecDestroy(&field->Su); CHKERRQ(ierr);
    ierr = VecDestroy(&field->Sv); CHKERRQ(ierr);
    ierr = VecDestroy(&field->Sw); CHKERRQ(ierr);
    ierr = VecDestroy(&field->Sphi); CHKERRQ(ierr);

    ierr = VecDestroy(&field->RHSu); CHKERRQ(ierr);
    ierr = VecDestroy(&field->RHSv); CHKERRQ(ierr);
    ierr = VecDestroy(&field->RHSw); CHKERRQ(ierr);
    ierr = VecDestroy(&field->RHSphi); CHKERRQ(ierr);


    ierr = VecDestroy(&field->RHS_PPE); CHKERRQ(ierr);
    ierr = MatDestroy(&field->Amat_PPE); CHKERRQ(ierr);



    return(0);
}
