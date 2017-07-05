#include <DataTypes.h>
#include <PETSCFunctions.h>
#include <Functions.h>

PetscErrorCode ComputeMatrixPPE(DM da, Mat Amat_PPE, void *g)
{
    StructGrid	*grid = (StructGrid*) g;
    PetscErrorCode 	ierr;
    Integer		Nx,Ny,Nz;
    Integer 		i,j,k;
    Integer 		xs,ys,zs,xm,ym,zm;
    PetscReal 		v[7];
    MatStencil 		row, col[7];
    PetscReal 		hx,hy,hz;
    PetscReal 		hx2,hy2,hz2;
    PetscReal 		a0,a1,a2,a3,a4,a5,a6;


    PetscReal 		***rho;


    hx = grid->dcsi;
    hy = grid->deta;
    hz = grid->dzeta;


    hx2 = hx*hx;
    hy2 = hy*hy;
    hz2 = hz*hz;

    //ierr = MatZeroEntries(Amat_PPE); 

    ierr = DMDAGetInfo(da,0,&Nx,&Ny,&Nz,0,0,0,0,0,0,0,0,0);
    ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
    //ierr = DMDAVecGetArray(da,user->soln.rho,&rho);

    for (k=zs; k<zs+zm; k++) {
       for (j=ys; j<ys+ym; j++) {
          for (i=xs; i<xs+xm; i++) {

           row.i = i;
           row.j = j;        
	   row.k = k;


           if(i == 0){
                a2 = 0;
                a4 = 1.0/hx2;
            }
	   else if(i == Nx-1){
                a2 = 1.0/hx2;
                a4 = 0;
            }
	   else if(i != 0 && i != Nx-1){
                a2 = 1.0/hx2;
                a4 = 1.0/hx2;
            }

           if(j == 0){
                a1 = 0;
                a5 = 1.0/hy2;
            }
           else if(j == Ny-1){
                a1 = 1.0/hy2;
                a5 = 0;
            }
           else if(j != 0 && j != Ny-1){
                a1 = 1.0/hy2;
                a5 = 1.0/hy2;
            }

           if(k == 0){
                a0 = 0;
                a6 = 1.0/hz2;
            }
           else if(k == Nz-1){
                a0 = 1.0/hz2;
                a6 = 0;
            }
            else if(k != 0 && k != Nz-1){
                a0 = 1.0/hz2;
                a6 = 1.0/hz2;
            }

            a3 = -(a0 + a1 + a2 + a4 + a5 + a6);

            /* assemble matrix */
/*-------------------------------------------------------------------------------------------------- */
            /* 8 cases : corner points */
/*-------------------------------------------------------------------------------------------------- */
           if(i == 0 && j == 0 && k == 0){
                v[0] = a3;	col[0].i = i;	col[0].j = j;	col[0].k = k;	
                v[1] = a4;	col[1].i = i+1;	col[1].j = j;	col[1].k = k;
                v[2] = a5;	col[2].i = i;	col[2].j = j+1;	col[2].k = k;
                v[3] = a6;	col[3].i = i;	col[3].j = j;	col[3].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == 0 && j == Ny-1 && k == 0){
                v[0] = a1;	col[0].i = i;	col[0].j = j-1;	col[0].k = k;
                v[1] = a3;	col[1].i = i;	col[1].j = j;	col[1].k = k;
                v[2] = a4;	col[2].i = i+1;	col[2].j = j;	col[2].k = k;
                v[3] = a6;	col[3].i = i;	col[3].j = j;	col[3].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
/*-------------------------------------------------------------------------------------------------- */
           else if(i == Nx-1 && j == 0 && k == 0){
                v[0] = a2;	col[0].i = i-1;	col[0].j = j;	col[0].k = k;
                v[1] = a3;	col[1].i = i;   col[1].j = j;	col[1].k = k;
                v[2] = a5;	col[2].i = i;   col[2].j = j+1;	col[2].k = k;
                v[3] = a6;	col[3].i = i;	col[3].j = j;	col[3].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == Nx-1 && j == Ny-1 && k == 0){
                v[0] = a1;	col[0].i = i;   col[0].j = j-1;	col[0].k = k;
                v[1] = a2;	col[1].i = i-1; col[1].j = j;	col[1].k = k;
                v[2] = a3;	col[2].i = i;   col[2].j = j;	col[2].k = k;
                v[3] = a6;	col[3].i = i;	col[3].j = j;	col[3].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
/*-------------------------------------------------------------------------------------------------- */
           else if(i == 0 && j == 0 && k == Nz-1){
                v[0] = a0;          col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a3;          col[1].i = i;	col[1].j = j;	col[1].k = k;	
                v[2] = a4;          col[2].i = i+1;	col[2].j = j;	col[2].k = k;
                v[3] = a5;          col[3].i = i;	col[3].j = j+1;	col[3].k = k;

                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == 0 && j == Ny-1 && k == Nz-1){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a1;	col[1].i = i;	col[1].j = j-1;	col[1].k = k;
                v[2] = a3;	col[2].i = i;	col[2].j = j;	col[2].k = k;
                v[3] = a4;	col[3].i = i+1;	col[3].j = j;	col[3].k = k;

                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
/*-------------------------------------------------------------------------------------------------- */
           else if(i == Nx-1 && j == 0 && k == Nz-1){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a2;	col[1].i = i-1;	col[1].j = j;	col[1].k = k;
                v[2] = a3;	col[2].i = i;   col[2].j = j;	col[2].k = k;
                v[3] = a5;	col[3].i = i;   col[3].j = j+1;	col[3].k = k;

                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == Nx-1 && j == Ny-1 && k == Nz-1){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a1;	col[1].i = i;   col[1].j = j-1;	col[1].k = k;
                v[2] = a2;	col[2].i = i-1; col[2].j = j;	col[2].k = k;
                v[3] = a3;	col[3].i = i;   col[3].j = j;	col[3].k = k;

                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,4,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
/*-------------------------------------------------------------------------------------------------- */
            /* 12 cases : corner edges */
/*-------------------------------------------------------------------------------------------------- */
           else if(i == 0 && j == 0 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a3;	col[1].i = i;	col[1].j = j;	col[1].k = k;	
                v[2] = a4;	col[2].i = i+1;	col[2].j = j;	col[2].k = k;
                v[3] = a5;	col[3].i = i;	col[3].j = j+1;	col[3].k = k;
                v[4] = a6;	col[4].i = i;	col[4].j = j;	col[4].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == 0 && j == Ny-1){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a1;	col[1].i = i;	col[1].j = j-1;	col[1].k = k;
                v[2] = a3;	col[2].i = i;	col[2].j = j;	col[2].k = k;	
                v[3] = a4;	col[3].i = i+1;	col[3].j = j;	col[3].k = k;
                v[4] = a6;	col[4].i = i;	col[4].j = j;	col[4].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == Nx-1 && j == 0 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a2;	col[1].i = i-1;	col[1].j = j;	col[1].k = k;
                v[2] = a3;	col[2].i = i;	col[2].j = j;	col[2].k = k;	
                v[3] = a4;	col[3].i = i;	col[3].j = j+1;	col[3].k = k;
                v[4] = a6;	col[4].i = i;	col[4].j = j;	col[4].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == Nx-1 && j == Ny-1 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a4;	col[1].i = i;	col[1].j = j-1;	col[1].k = k;
                v[2] = a2;	col[2].i = i-1;	col[2].j = j;	col[2].k = k;
                v[3] = a3;	col[3].i = i;	col[3].j = j;	col[3].k = k;	
                v[4] = a6;	col[4].i = i;	col[4].j = j;	col[4].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
/*-------------------------------------------------------------------------------------------------- */
           else if(i == 0 && k == 0 ){
                v[0] = a2;	col[0].i = i-1;	col[0].j = j;	col[0].k = k;
                v[1] = a3;	col[1].i = i;	col[1].j = j;	col[1].k = k;	
                v[2] = a4;	col[2].i = i+1;	col[2].j = j;	col[2].k = k;
                v[3] = a5;	col[3].i = i;	col[3].j = j+1;	col[3].k = k;
                v[4] = a6;	col[4].i = i;	col[4].j = j;	col[4].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == 0 && k == Nz-1 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a2;	col[1].i = i-1;	col[1].j = j-1;	col[1].k = k;
                v[2] = a3;	col[2].i = i;	col[2].j = j;	col[2].k = k;	
                v[3] = a4;	col[3].i = i+1;	col[3].j = j;	col[3].k = k;
                v[4] = a5;	col[4].i = i;	col[4].j = j+1;	col[4].k = k;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == Nx-1 && k == 0 ){
                v[0] = a1;	col[0].i = i;	col[0].j = j-1;	col[0].k = k;
                v[1] = a2;	col[1].i = i-1;	col[1].j = j;	col[1].k = k;
                v[2] = a3;	col[2].i = i;	col[2].j = j;	col[2].k = k;	
                v[3] = a5;	col[3].i = i;	col[3].j = j+1;	col[3].k = k;
                v[4] = a6;	col[4].i = i;	col[4].j = j;	col[4].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == Nx-1 && k == 0 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a1;	col[1].i = i;	col[1].j = j-1;	col[1].k = k;
                v[2] = a2;	col[2].i = i-1;	col[2].j = j;	col[2].k = k;
                v[3] = a3;	col[3].i = i;	col[3].j = j;	col[3].k = k;	
                v[4] = a5;	col[4].i = i;	col[4].j = j+1;	col[4].k = k;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
/*-------------------------------------------------------------------------------------------------- */
           else if( j == 0 && k == 0 ){
                v[0] = a2;	col[0].i = i-1;	col[0].j = j;	col[0].k = k;
                v[1] = a3;	col[1].i = i;	col[1].j = j;	col[1].k = k;	
                v[2] = a4;	col[2].i = i+1;	col[2].j = j;	col[2].k = k;
                v[3] = a5;	col[3].i = i;	col[3].j = j+1;	col[3].k = k;
                v[4] = a6;	col[4].i = i;	col[4].j = j;	col[4].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if( j == 0 && k == Nz-1 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a2;	col[1].i = i-1;	col[1].j = j;	col[1].k = k;
                v[2] = a3;	col[2].i = i;	col[2].j = j;	col[2].k = k;	
                v[3] = a4;	col[3].i = i+1;	col[3].j = j;	col[3].k = k;
                v[4] = a5;	col[4].i = i;	col[4].j = j+1;	col[4].k = k;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if( j == Ny-1 && k == 0 ){
                v[0] = a1;	col[0].i = i;	col[0].j = j-1;	col[0].k = k;
                v[1] = a2;	col[1].i = i-1;	col[1].j = j;	col[1].k = k;
                v[2] = a3;	col[2].i = i;	col[2].j = j;	col[2].k = k;	
                v[3] = a4;	col[3].i = i+1;	col[3].j = j;	col[3].k = k;
                v[4] = a6;	col[4].i = i;	col[4].j = j;	col[4].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if( j == Ny-1 && k == Nz-1 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a1;	col[1].i = i;	col[1].j = j-1;	col[1].k = k;
                v[2] = a2;	col[2].i = i-1;	col[2].j = j;	col[2].k = k;
                v[3] = a3;	col[3].i = i;	col[3].j = j;	col[3].k = k;	
                v[4] = a4;	col[4].i = i+1;	col[4].j = j;	col[4].k = k;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
/*-------------------------------------------------------------------------------------------------- */
            /* 6 cases : boundary faces */
/*-------------------------------------------------------------------------------------------------- */
           else if(k == 0 ){
                v[0] = a1;	col[0].i = i;	col[0].j = j-1;	col[0].k = k;
                v[1] = a2;	col[1].i = i-1;	col[1].j = j;	col[1].k = k;
                v[2] = a3;	col[2].i = i;	col[2].j = j;	col[2].k = k;	
                v[3] = a4;	col[3].i = i+1;	col[3].j = j;	col[3].k = k;
                v[4] = a5;	col[4].i = i;	col[4].j = j+1;	col[4].k = k;
                v[5] = a6;	col[5].i = i;	col[5].j = j;	col[5].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,6,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(k == Nz-1 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a1;	col[1].i = i;	col[1].j = j-1;	col[1].k = k;
                v[2] = a2;	col[2].i = i-1;	col[2].j = j;	col[2].k = k;	
                v[3] = a3;	col[3].i = i;	col[3].j = j;	col[3].k = k;
                v[4] = a4;	col[4].i = i+1;	col[4].j = j;	col[4].k = k;
                v[5] = a5;	col[5].i = i;	col[5].j = j+1;	col[5].k = k;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,6,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
/*-------------------------------------------------------------------------------------------------- */
           else if(j == 0 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a2;	col[1].i = i-1;	col[1].j = j;	col[1].k = k;
                v[2] = a3;	col[2].i = i;	col[2].j = j;	col[2].k = k;	
                v[3] = a4;	col[3].i = i+1;	col[3].j = j;	col[3].k = k;
                v[4] = a5;	col[4].i = i;	col[4].j = j+1;	col[4].k = k;
                v[5] = a6;	col[5].i = i;	col[5].j = j;	col[5].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,6,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(j == Ny-1  ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a1;	col[1].i = i;	col[1].j = j-1;	col[1].k = k;
                v[2] = a2;	col[2].i = i-1;	col[2].j = j;	col[2].k = k;	
                v[3] = a3;	col[3].i = i;	col[3].j = j;	col[3].k = k;
                v[4] = a4;	col[4].i = i+1;	col[4].j = j;	col[4].k = k;
                v[5] = a6;	col[5].i = i;	col[5].j = j;	col[5].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,6,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
/*-------------------------------------------------------------------------------------------------- */
           else if(i == 0 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k+1;
                v[1] = a1;	col[1].i = i;	col[1].j = j-1;	col[1].k = k;
                v[2] = a3;	col[2].i = i;	col[2].j = j;	col[2].k = k;	
                v[3] = a4;	col[3].i = i+1;	col[3].j = j;	col[3].k = k;
                v[4] = a5;	col[4].i = i;	col[4].j = j+1;	col[4].k = k;
                v[5] = a6;	col[5].i = i;	col[5].j = j;	col[5].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,6,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

           else if(i == Nx-1 ){
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a1;	col[1].i = i;	col[1].j = j-1;	col[1].k = k;
                v[2] = a2;	col[2].i = i-1;	col[2].j = j;	col[2].k = k;	
                v[3] = a3;	col[3].i = i;	col[3].j = j;	col[3].k = k;
                v[4] = a5;	col[4].i = i;	col[4].j = j+1;	col[4].k = k;
                v[5] = a6;	col[5].i = i;	col[5].j = j;	col[5].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,6,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }
/*-------------------------------------------------------------------------------------------------- */
            /*  1 case : Interior points */
/*-------------------------------------------------------------------------------------------------- */
           else {
                v[0] = a0;	col[0].i = i;	col[0].j = j;	col[0].k = k-1;
                v[1] = a1;	col[1].i = i;	col[1].j = j-1;	col[1].k = k;
                v[2] = a2;	col[2].i = i-1;	col[2].j = j;	col[2].k = k;	
                v[3] = a3;	col[3].i = i;	col[3].j = j;	col[3].k = k;
                v[4] = a5;	col[4].i = i+1;	col[4].j = j;	col[4].k = k;
                v[5] = a6;	col[5].i = i;	col[5].j = j+1;	col[5].k = k;
                v[6] = a6;	col[6].i = i;	col[6].j = j;	col[6].k = k+1;
                ierr =  MatSetValuesStencil(Amat_PPE,1,&row,7,col,v,INSERT_VALUES); CHKERRQ(ierr);
            }

	  }
       }
    }

    ierr = MatAssemblyBegin(Amat_PPE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Amat_PPE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


    //ierr = MatSetOption(Amat_PPE,MAT_NO_OFF_PROC_ENTRIES,PETSC_TRUE); CHKERRQ(ierr);
    ierr = MatSetOption(Amat_PPE,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE); CHKERRQ(ierr);


//    ierr = DMDAVecRestoreArray(da,user->soln.rho,&rho); CHKERRQ(ierr);

    MatNullSpace nullspace;
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace); CHKERRQ(ierr);
    ierr = MatSetNullSpace(Amat_PPE,nullspace); CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullspace); CHKERRQ(ierr);


  return(0);
}

