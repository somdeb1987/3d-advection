#include <DataTypes.h>
#include <PETSCFunctions.h>

int TransferVecToTS(DM da,Vec u,Vec v, Vec w,Vec phi,Vec Y,void *ctxt) 
{
  PETScTSContext	*context 	= (PETScTSContext*) ctxt;
  PetscErrorCode  	ierr;
  VarArray		***arru,***arrv,***arrw,***arrphi;
  VarArray    	  	*Yarr;
  DMDALocalInfo		info;
  int             	istart,iend;
  int             	jstart,jend;
  int             	kstart,kend;
  int             	i,j,k,y_index;

  //ierr = VecGetArray(Y,&Yarr); CHKERRQ(ierr);
  ierr = VecGetArrayRead(Y,&Yarr); CHKERRQ(ierr);

  /* ======================================================== */

/*
  DMDAGetLocalInfo(da,&info);
  istart = info.xs;
  iend = info.xs + info.xm;
  jstart = info.ys;
  jend = info.ys + info.ym;
  kstart = info.zs;
  kend = info.zs + info.zm;
*/

  Integer xs,xm,ys,ym,zs,zm;
  DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  istart=xs;	iend=xs+xm;
  jstart=ys;	jend=ys+ym;
  kstart=zs;	kend=zs+zm;

  int offset_u 		= context->offset_u;
  int offset_v 		= context->offset_v;
  int offset_w 		= context->offset_w;
  int offset_phi 	= context->offset_phi;

  /* copy u */
  ierr = DMDAVecGetArray(da,u,&arru); CHKERRQ(ierr);

  y_index = 0;
  for (i=istart; i<iend; i++) {
    for (j=jstart; j<jend; j++) {
    	for (k=kstart; k<kend; k++) {
		Yarr[offset_u+y_index] = arru[k][j][i];
        	y_index++;
	}
    }
  }

  ierr = DMDAVecRestoreArray(da,u,&arru); CHKERRQ(ierr);

  /* copy v */
  ierr = DMDAVecGetArray(da,v,&arrv); CHKERRQ(ierr);

  y_index = 0;
  for (i=istart; i<iend; i++) {
    for (j=jstart; j<jend; j++) {
    	for (k=kstart; k<kend; k++) {
        	Yarr[offset_v+y_index] = arrv[k][j][i];
        	y_index++;
	}
    }
  }

  ierr = DMDAVecRestoreArray(da,v,&arrv); CHKERRQ(ierr);

  /* copy w */
  ierr = DMDAVecGetArray(da,w,&arrw); CHKERRQ(ierr);

  y_index = 0;
  for (i=istart; i<iend; i++) {
    for (j=jstart; j<jend; j++) {
    	for (k=kstart; k<kend; k++) {
        	Yarr[offset_w+y_index] = arrw[k][j][i];
        	y_index++;
	}
    }
  }

  ierr = DMDAVecRestoreArray(da,w,&arrw); CHKERRQ(ierr);

  /* copy phi */
  ierr = DMDAVecGetArray(da,phi,&arrphi); CHKERRQ(ierr);

  y_index = 0;
  for (i=istart; i<iend; i++) {
    for (j=jstart; j<jend; j++) {
    	for (k=kstart; k<kend; k++) {
        	Yarr[offset_phi+y_index] = arrphi[k][j][i];
        	y_index++;
	}
    }
  }
  ierr = DMDAVecRestoreArray(da,phi,&arrphi); CHKERRQ(ierr);

  /* ======================================================== */


  ierr = VecRestoreArrayRead(Y,&Yarr); CHKERRQ(ierr);

  return(0);
}

