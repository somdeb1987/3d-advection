#include <DataTypes.h>
#include <PETSCFunctions.h>
#include <Functions.h>

PetscErrorCode ComputeRHSPPE(DM da, Vec ug, Vec vg, Vec wg, Vec RHS,double dt, void *g)
{
    StructGrid	*grid = (StructGrid*) g;
    PetscErrorCode 	ierr;
    Integer		Nx,Ny,Nz;
    Integer 		i,j,k;
    Integer 		xs,ys,zs,xm,ym,zm;
    double 		hx,hy,hz;
    double 		ux=0, vy=0, wz=0;
    Vec            	ulocal,vlocal,wlocal;
    VarArray 		***rhsarr, ***uarr, ***varr, ***warr;


    hx = grid->dcsi;
    hy = grid->deta;
    hz = grid->dzeta;

    Nx = grid->Nx;
    Ny = grid->Ny;
    Nz = grid->Nz;

    ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);

    DMCreateLocalVector(da,&ulocal);
    VecZeroEntries(ulocal);
    DMGlobalToLocalBegin(da,ug,INSERT_VALUES,ulocal);
    DMGlobalToLocalEnd  (da,ug,INSERT_VALUES,ulocal);

    DMCreateLocalVector(da,&vlocal);
    VecZeroEntries(vlocal);
    DMGlobalToLocalBegin(da,vg,INSERT_VALUES,vlocal);
    DMGlobalToLocalEnd  (da,vg,INSERT_VALUES,vlocal);

    DMCreateLocalVector(da,&wlocal);
    VecZeroEntries(wlocal);
    DMGlobalToLocalBegin(da,wg,INSERT_VALUES,wlocal);
    DMGlobalToLocalEnd  (da,wg,INSERT_VALUES,wlocal);



    ierr = DMDAVecGetArray(da,ulocal,&uarr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,vlocal,&varr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,wlocal,&warr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,RHS,&rhsarr); CHKERRQ(ierr);


    for(i = xs; i < xs+xm; i++){
      for(j = ys; j < ys+ym; j++){
	 for(k = zs; k < zs+zm; k++){
	 // --------------------------------------------------------------------------------------------- //
            if(i == 0){
                ux = 0.5*(-3*uarr[k][j][i] + 4*uarr[k][j][i+1] - uarr[k][j][i+2])/hx;
            }
	    else if(i == Nx-1){
                ux = 0.5*(3*uarr[k][j][i] - 4*uarr[k][j][i-1] + uarr[k][j][i-2])/hx;
            }
            else if(i > 0 && i < Nx-1 ){
                ux = (0.5*(uarr[k][j][i+1] - uarr[k][j][i-1]))/hx;

            }

	 // --------------------------------------------------------------------------------------------- //
            if(j == 0){
                vy = 0.5*(-3*varr[k][j][i] + 4*varr[k][j+1][i] - varr[k][j+2][i])/hy;
            }
            else if (j == Ny-1){
                vy = 0.5*(3*varr[k][j][i] - 4*varr[k][j-1][i] + varr[k][j-2][i])/hy;
            }
            else if(j > 0 && j < Ny-1){
                vy = (0.5*(varr[k][j+1][i] - varr[k][j-1][i]))/hy;

            }
	 // --------------------------------------------------------------------------------------------- //
            if(k == 0){
                wz = 0.5*(-3*warr[k][j][i] + 4*warr[k+1][j][i] - warr[k+2][j][i])/hz;
            }
            else if (k == Nz-1){
                wz = 0.5*( 3*warr[k][j][i] - 4*warr[k-1][j][i] + warr[k-2][j][i])/hz;
            }
            else if(k > 0 && k < Nz-1){
                wz = (0.5*(warr[k+1][j][i] - warr[k-1][j][i]))/hz;

            }
	 // --------------------------------------------------------------------------------------------- //

	
            rhsarr[k][j][i] = (ux + vy + wz)/dt;

         }
       }
    }

    ierr = DMDAVecRestoreArray(da,ulocal,&uarr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,vlocal,&varr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,wlocal,&warr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,RHS,&rhsarr); CHKERRQ(ierr);

    VecDestroy(&ulocal);
    VecDestroy(&vlocal);
    VecDestroy(&wlocal);


    MatNullSpace nullspace;
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace); CHKERRQ(ierr);
    ierr = MatNullSpaceRemove(nullspace,RHS); CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullspace); CHKERRQ(ierr);



    return(0);
}
