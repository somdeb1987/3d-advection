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
    PetscReal 		a0,a1, a2, a3, a4, a5,a6;
    PetscReal 		***rho;


    hx = grid->dcsi;
    hy = grid->deta;
    hz = grid->dzeta;

    hx2 = hx*hx;
    hy2 = hy*hy;
    hz2 = hz*hz;


    ierr = DMDAGetInfo(da,0,&Nx,&Ny,&Nz,0,0,0,0,0,0,0,0,0);
    ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
    //ierr = DMDAVecGetArray(da,user->soln.rho,&rho);

    for (k=zs; k<zs+zm; k++) {
       for (j=ys; j<ys+ym; j++) {
          for (i=xs; i<xs+xm; i++) {

    	Integer 		num,numi,numj,numk;
    	PetscReal 		HyHzdHx,HxHzdHy,HxHydHz;

    	HyHzdHx = hx*hx;
    	HxHzdHy = hy*hy;
    	HxHydHz = hz*hz;

        //we relate to this point i,j,k
        row.i = i;
        row.j = j;
        row.k = k;


        //do corner cases, at boundary
        if(i==0 || j==0 || k==0 || i==Nx-1 || j==Ny-1 || k==Nz-1)
        {
          //stencil vector position counter
          num = 0;
          numi = 0;         numj = 0;         numk = 0;

          //not on wall?
          if(k!=0)
          {
            v[num] = -HxHydHz;  col[num].i = i;  col[num].j = j;  col[num].k = k-1;

            num++; //increase stencil fill counter
            numk++;
          }

          //not on wall?
          if(j!=0)
          {
            v[num] = -HxHzdHy;  col[num].i = i; col[num].j = j-1; col[num].k = k;

            num++;
            numj++;
          }

          //not on wall?
          if(i!=0)
          {
            v[num] = -HyHzdHx; col[num].i = i-1; col[num].j = j; col[num].k = k;

            num++;
            numi++;
          }

          //not on wall?
          if(i!=Nx-1)
          {
            v[num] = -HyHzdHx; col[num].i = i+1; col[num].j = j; col[num].k = k;

            num++;
            numi++;
          }

          //not on wall?
          if(j!=Ny-1)
          {
            v[num] = -HxHzdHy; col[num].i = i; col[num].j = j+1; col[num].k = k;

            num++;
            numj++;
          }

          //not on wall?
          if(k!=Nz-1)
          {
            v[num] = -HxHydHz; col[num].i = i; col[num].j = j; col[num].k = k+1;

            num++;
            numk++;
          }

          //center point
          v[num] = (PetscReal)(numk)*HxHydHz + (PetscReal)(numj)*HxHzdHy + (PetscReal)(numi)*HyHzdHx;

          col[num].i = i;
          col[num].j = j;
          col[num].k = k;

          num++;

          MatSetValuesStencil(Amat_PPE,1,&row,num,col,v,INSERT_VALUES);
        }
        else //interior points, not lying on any boundary
        {
          v[0] = -HxHydHz;                          col[0].i = i;   col[0].j = j;   col[0].k = k-1;
          v[1] = -HxHzdHy;                          col[1].i = i;   col[1].j = j-1; col[1].k = k;
          v[2] = -HyHzdHx;                          col[2].i = i-1; col[2].j = j;   col[2].k = k;
          v[3] = 2.0*(HyHzdHx + HxHzdHy + HxHydHz); col[3].i = i;   col[3].j = j;   col[3].k = k;
          v[4] = -HyHzdHx;                          col[4].i = i+1; col[4].j = j;   col[4].k = k;
          v[5] = -HxHzdHy;                          col[5].i = i;   col[5].j = j+1; col[5].k = k;
          v[6] = -HxHydHz;                          col[6].i = i;   col[6].j = j;   col[6].k = k+1;

          MatSetValuesStencil(Amat_PPE,1,&row,7,col,v,INSERT_VALUES);
        }


	  }
       }
    }

    ierr = MatAssemblyBegin(Amat_PPE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Amat_PPE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

/*
    ierr = MatSetOption(Amat_PPE,MAT_NO_OFF_PROC_ENTRIES,PETSC_TRUE); CHKERRQ(ierr);
    ierr = MatSetOption(Amat_PPE,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE); CHKERRQ(ierr);


//    ierr = DMDAVecRestoreArray(da,user->soln.rho,&rho); CHKERRQ(ierr);

    MatNullSpace nullspace;
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace); CHKERRQ(ierr);
    ierr = MatSetNullSpace(Amat_PPE,nullspace); CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullspace); CHKERRQ(ierr);
*/

  return(0);
}

