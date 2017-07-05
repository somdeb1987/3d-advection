#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <DataTypes.h>

PetscErrorCode SaveSolution(void *g,int *c, double t, DM da, 
			    Vec ug,
			    Vec vg,
		            Vec wg,
			    Vec phig)
{
   StructGrid	*grid = (StructGrid*) g;

   PetscErrorCode ierr;
   char           filename[32] = "sol";
   PetscMPIInt    rank;
   FILE           *fp;
   Vec            ulocal,vlocal,wlocal, philocal;
   PetscScalar    ***u,***v,***w,***phi;
   DMDALocalInfo	info;
   Integer 		i, is, ie;
   Integer 		j, js, je;
   Integer 		k, ks, ke;
   Integer 		nlocx, nlocy, nlocz;


   ierr = DMDAGetCorners(da, &is, &js, &ks, &nlocx, &nlocy, &nlocz); CHKERRQ(ierr);
   ie=is+nlocx;
   je=js+nlocy;
   ke=ks+nlocz;

   DMGetLocalVector(da,&ulocal);
   VecZeroEntries(ulocal);
   DMGlobalToLocalBegin(da,ug,INSERT_VALUES,ulocal);
   DMGlobalToLocalEnd  (da,ug,INSERT_VALUES,ulocal);

   DMGetLocalVector(da,&vlocal);
   VecZeroEntries(vlocal);
   DMGlobalToLocalBegin(da,vg,INSERT_VALUES,vlocal);
   DMGlobalToLocalEnd  (da,vg,INSERT_VALUES,vlocal);

   DMGetLocalVector(da,&wlocal);
   VecZeroEntries(wlocal);
   DMGlobalToLocalBegin(da,wg,INSERT_VALUES,wlocal);
   DMGlobalToLocalEnd  (da,wg,INSERT_VALUES,wlocal);

   DMGetLocalVector(da,&philocal);
   VecZeroEntries(philocal);
   DMGlobalToLocalBegin(da,phig,INSERT_VALUES,philocal);
   DMGlobalToLocalEnd  (da,phig,INSERT_VALUES,philocal);




   ierr = DMDAVecGetArray(da, ulocal, &u); CHKERRQ(ierr);
   ierr = DMDAVecGetArray(da, vlocal, &v); CHKERRQ(ierr);
   ierr = DMDAVecGetArray(da, wlocal, &w); CHKERRQ(ierr);
   ierr = DMDAVecGetArray(da, philocal, &phi); CHKERRQ(ierr);

   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   sprintf(filename, "sol-%03d-%03d.plt", *c, rank);
   fp = fopen(filename,"w");
   fprintf(fp, "TITLE = \"u_t + u_x + u_y = 0\"\n");
   fprintf(fp, "VARIABLES = x, y,z, u, v, w, phi\n");
   fprintf(fp, "ZONE STRANDID=1, SOLUTIONTIME=%e, I=%d, J=%d, K=%d, DATAPACKING=POINT\n", t, ie-is+1, je-js+1, ke-ks+1);
   for(k=ks; k<ke+1; ++k){
     for(j=js; j<je+1; ++j){
      for(i=is; i<ie+1; ++i){
   		{
      			PetscReal x = grid->csimin + i*grid->dcsi;
      			PetscReal y = grid->etamin + j*grid->deta;
      			PetscReal z = grid->zetamin + k*grid->dzeta;
      			fprintf(fp, "%e %e %e %e %e %e %e\n", 
				x, y, z,
				u[k][j][i], v[k][j][i],w[k][j][i], 
				phi[k][j][i]);
   		}
      }
     }
   }
   fclose(fp);

   ierr = DMDAVecRestoreArray(da, ulocal, &u); CHKERRQ(ierr);
   ierr = DMDAVecRestoreArray(da, vlocal, &v); CHKERRQ(ierr);
   ierr = DMDAVecRestoreArray(da, wlocal, &w); CHKERRQ(ierr);
   ierr = DMDAVecRestoreArray(da, philocal, &phi); CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(da, &ulocal); CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(da, &vlocal); CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(da, &wlocal); CHKERRQ(ierr);
   ierr = DMRestoreLocalVector(da, &philocal); CHKERRQ(ierr);

   ++(*c);
   return(0);
}
