#include <stdio.h>
#include <petscsys.h>
#include <petscdmda.h>
#include <petscvec.h>
#include <Definitions.h>
#include <DataTypes.h>
#include <Functions.h>

double absolute(double);
double raiseto(double,double);

/* Calculate R = (1/dx)*[f_x] using fifth order WENO    */
/* R is a global vector, while F is a local vector      */
PetscErrorCode FifthOrderWENO(Vec Fx, Vec Fy, Vec Fz, Vec R, Vec uwX, Vec uwY,Vec uwZ, void *g)
{
  StructGrid    	*grid	= (StructGrid*)	g;
  int     		i, is, ie;
  int     		j, js, je;
  int     		k, ks, ke;
  Vec		  	fxI,fyI,fzI;
  PetscReal  		***fx,***fy,***fz;
  PetscReal  		***fxi,***fyi,***fzi;
  PetscReal  		***upwX,***upwY,***upwZ;
  PetscReal  		***r;
  double		dx,dy,dz;
  double  		dxinv,dyinv,dzinv;
  Integer 		xs,xm,ys,ym,zs,zm;
  DM 			da;
  DMDALocalInfo 	info;
  double 		eps		= 1e-6;
  double 		p_weno		= 2;
  int			ierr;




  dx	=	grid->dcsi;
  dy	=	grid->deta; 
  dz	=	grid->dzeta;

/*
  Integer 		LeftB,RightB;
  Integer 		BottomB,TopB;
  Integer 		FrontB,BackB;

  LeftB		=	grid->LeftB	;
  RightB 	=	grid->RightB	;
  BottomB 	=	grid->BottomB	;
  TopB		=	grid->TopB	;
  FrontB 	=	grid->FrontB	;
  BackB		=	grid->BackB	;

*/

  /* Calculate interface flux */
  VecGetDM(Fx,&da);
/*
  DMDAGetLocalInfo(da,&info);
  is = info.xs;
  ie = info.xs + info.xm;
  js = info.ys;
  je = info.ys + info.ym;
  ks = info.zs;
  ke = info.zs + info.zm;
*/


  DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  is=xs;	ie=xs+xm;
  js=ys;	je=ys+ym;
  ks=zs;	ke=zs+zm;

  DMDAVecGetArray(da,Fx,&fx);
  DMDAVecGetArray(da,uwX,&upwX);
  DMCreateLocalVector(da,&fxI);
  DMDAVecGetArray(da,fxI,&fxi);

  for (i=is; i<ie+1; i++) {
     for (j=js; j<je; j++) {
	  for (k=ks; k<ke; k++) {

    		/* Defining stencil points */
    		double m3, m2, m1, p1, p2;
		if(upwX[k][j][i] > 0){
    			m3 = fx[k][j][i-3];
    			m2 = fx[k][j][i-2];
    			m1 = fx[k][j][i-1];
    			p1 = fx[k][j][i]  ;
    			p2 = fx[k][j][i+1];
		}
		else{

    			m3 = fx[k][j][i+2];
    			m2 = fx[k][j][i+1];
    			m1 = fx[k][j][i]  ;
    			p1 = fx[k][j][i-1];
    			p2 = fx[k][j][i-2];

		}

    		/* Candidate stencils */
    		double f1, f2, f3;
    		f1 = (2*one_sixth)*m1 + (5*one_sixth)*p1 - (one_sixth)*p2;
    		f2 = (-one_sixth)*m2 + (5.0*one_sixth)*m1 + (2*one_sixth)*p1;
    		f3 = (2*one_sixth)*m3 - (7.0*one_sixth)*m2 + (11.0*one_sixth)*m1;


    		/* Candidate stencils and their optimal weights*/
		double w1,w2,w3;
  		WENOParameters	WenoParams;
  		double 		c1,c2,c3;
  		ierr 		= 	SetWENOParams(&WenoParams);if (ierr) return (ierr);
  		eps		=	WenoParams.eps;
  		p_weno		=	WenoParams.p;

  		c1 = 0.3;
  		c2 = 0.6;
  		c3 = 0.1;


  		if (WenoParams.no_limiting) {
			w1 = c1;
   			w2 = c2;
      			w3 = c3;
    		} 
  		else {
      		// Smoothness indicators 
      		double b1, b2, b3;
      		b1 = thirteen_by_twelve*(m1-2*p1+p2)*(m1-2*p1+p2)
           		+ one_fourth*(3*m1-4*p1+p2)*(3*m1-4*p1+p2);
      		b2 = thirteen_by_twelve*(m2-2*m1+p1)*(m2-2*m1+p1)
           		+ one_fourth*(m2-p1)*(m2-p1);
      		b3 = thirteen_by_twelve*(m3-2*m2+m1)*(m3-2*m2+m1) 
           		+ one_fourth*(m3-4*m2+3*m1)*(m3-4*m2+3*m1);
  
      		// This parameter is needed for Borges' or Yamaleev-Carpenter
         	//implementation of non-linear weights or to obtain a fifth-
         	//order interpolation scheme with no limiting. 
      		double tau;
      		if (WenoParams.borges) {
        		tau = absolute(b3 - b1);
      		} 
		else if (WenoParams.yc) {
        		tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);
      		} 
		else {
        		tau = 0;
      		}
  
      		// Definiting the non-linear weights 
      		double a1, a2, a3;
      		if (WenoParams.borges || WenoParams.yc) {
        		a1 = c1 * (1.0 + raiseto(tau/(b1+eps),p_weno));
        		a2 = c2 * (1.0 + raiseto(tau/(b2+eps),p_weno));
        		a3 = c3 * (1.0 + raiseto(tau/(b3+eps),p_weno));
      		} 
		else {
        		a1 = c1 / raiseto(b1+eps,p_weno);
       			a2 = c2 / raiseto(b2+eps,p_weno);
        		a3 = c3 / raiseto(b3+eps,p_weno);
      		}

      		// Convexity 
      		double a_sum_inv;
      		a_sum_inv = 1.0 / (a1 + a2 + a3);
      		w1 = a1 * a_sum_inv;
      		w2 = a2 * a_sum_inv;
     		w3 = a3 * a_sum_inv;
  
      		// Mapping the weights, if needed 
      		if (WenoParams.mapped) {
        		a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1));
        		a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2));
        		a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3));
        		a_sum_inv = 1.0 / (a1 + a2 + a3);
        		w1 = a1 * a_sum_inv;
        		w2 = a2 * a_sum_inv;
        		w3 = a3 * a_sum_inv;

      			}
  
    		}


		fxi[k][j][i] = w1*f1 + w2*f2 + w3*f3;
	   }
	}
  }

  DMDAVecRestoreArray(da,Fx,&fx);
  DMDAVecRestoreArray(da,uwX,&upwX);    


  VecGetDM(Fy,&da);
/*
  DMDAGetLocalInfo(da,&info);
  is = info.xs;
  ie = info.xs + info.xm;
  js = info.ys;
  je = info.ys + info.ym;
  ks = info.zs;
  ke = info.zs + info.zm;
*/

  DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  is=xs;	ie=xs+xm;
  js=ys;	je=ys+ym;
  ks=zs;	ke=zs+zm;

  DMDAVecGetArray(da,Fy,&fy);
  DMDAVecGetArray(da,uwY,&upwY);
  DMCreateLocalVector(da,&fyI);
  DMDAVecGetArray(da,fyI,&fyi);

  for (i=is; i<ie; i++) {
     for (j=js; j<je+1; j++) {
	for (k=ks; k<ke; k++) {
    		/* Defining stencil points */
    		double m3, m2, m1, p1, p2;
		if(upwY[k][j][i] > 0){
    			m3 = fy[k][j-3][i];
    			m2 = fy[k][j-2][i];
    			m1 = fy[k][j-1][i];
    			p1 = fy[k][j][i]  ;
    			p2 = fy[k][j+1][i];
		}
		else{
    			m3 = fy[k][j+2][i];
    			m2 = fy[k][j+1][i];
    			m1 = fy[k][j][i]  ;
    			p1 = fy[k][j-1][i];
    			p2 = fy[k][j-2][i];
		}

    		/* Candidate stencils */
    		double f1, f2, f3;
    		f1 = (2*one_sixth)*m1 + (5*one_sixth)*p1 - (one_sixth)*p2;
    		f2 = (-one_sixth)*m2 + (5.0*one_sixth)*m1 + (2*one_sixth)*p1;
    		f3 = (2*one_sixth)*m3 - (7.0*one_sixth)*m2 + (11.0*one_sixth)*m1;

    		/* Candidate stencils and their optimal weights*/
		double w1,w2,w3;
  		WENOParameters	WenoParams;
  		double 		c1,c2,c3;
  		ierr 		= 	SetWENOParams(&WenoParams);	if (ierr) return (ierr);
  		eps		=	WenoParams.eps;
  		p_weno		=	WenoParams.p; 
  		c1 = 0.3;
  		c2 = 0.6;
  		c3 = 0.1;


  		if (WenoParams.no_limiting) {
			w1 = c1;
   			w2 = c2;
      			w3 = c3;
    		} 
  		else {
      		// Smoothness indicators 
      		double b1, b2, b3;
      		b1 = thirteen_by_twelve*(m1-2*p1+p2)*(m1-2*p1+p2)
           		+ one_fourth*(3*m1-4*p1+p2)*(3*m1-4*p1+p2);
      		b2 = thirteen_by_twelve*(m2-2*m1+p1)*(m2-2*m1+p1)
           		+ one_fourth*(m2-p1)*(m2-p1);
      		b3 = thirteen_by_twelve*(m3-2*m2+m1)*(m3-2*m2+m1) 
           		+ one_fourth*(m3-4*m2+3*m1)*(m3-4*m2+3*m1);
  
      		// This parameter is needed for Borges' or Yamaleev-Carpenter
         	//implementation of non-linear weights or to obtain a fifth-
         	//order interpolation scheme with no limiting. 
      		double tau;
      		if (WenoParams.borges) {
        		tau = absolute(b3 - b1);
      		} 
		else if (WenoParams.yc) {
        		tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);
      		} 
		else {
        		tau = 0;
      		}
  
      		// Definiting the non-linear weights 
      		double a1, a2, a3;
      		if (WenoParams.borges || WenoParams.yc) {
        		a1 = c1 * (1.0 + raiseto(tau/(b1+eps),p_weno));
        		a2 = c2 * (1.0 + raiseto(tau/(b2+eps),p_weno));
        		a3 = c3 * (1.0 + raiseto(tau/(b3+eps),p_weno));
      		} 
		else {
        		a1 = c1 / raiseto(b1+eps,p_weno);
       			a2 = c2 / raiseto(b2+eps,p_weno);
        		a3 = c3 / raiseto(b3+eps,p_weno);
      		}

      		// Convexity 
      		double a_sum_inv;
      		a_sum_inv = 1.0 / (a1 + a2 + a3);
      		w1 = a1 * a_sum_inv;
      		w2 = a2 * a_sum_inv;
     		w3 = a3 * a_sum_inv;
  
      		// Mapping the weights, if needed 
      		if (WenoParams.mapped) {
        		a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1));
        		a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2));
        		a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3));
        		a_sum_inv = 1.0 / (a1 + a2 + a3);
        		w1 = a1 * a_sum_inv;
        		w2 = a2 * a_sum_inv;
        		w3 = a3 * a_sum_inv;

      			}
  
    		}

		fyi[k][j][i] = w1*f1 + w2*f2 + w3*f3;
	    }
	}
  }

  DMDAVecRestoreArray(da,Fy,&fy);
  DMDAVecRestoreArray(da,uwY,&upwY);   

  VecGetDM(Fz,&da);
/*
  DMDAGetLocalInfo(da,&info);
  is = info.xs;
  ie = info.xs + info.xm;
  js = info.ys;
  je = info.ys + info.ym;
  ks = info.zs;
  ke = info.zs + info.zm;
*/

  DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  is=xs;	ie=xs+xm;
  js=ys;	je=ys+ym;
  ks=zs;	ke=zs+zm;

  DMDAVecGetArray(da,Fz,&fz);
  DMDAVecGetArray(da,uwZ,&upwZ);
  DMCreateLocalVector(da,&fzI);
  DMDAVecGetArray(da,fzI,&fzi);

  for (i=is; i<ie; i++) {
     for (j=js; j<je; j++) {
	for (k=ks; k<ke+1; k++) {
    		/* Defining stencil points */
    		double m3, m2, m1, p1, p2;
		if(upwZ[k][j][i] > 0){
    			m3 = fz[k-3][j][i];
    			m2 = fz[k-2][j][i];
    			m1 = fz[k-1][j][i];
    			p1 = fz[k][j][i]  ;
    			p2 = fz[k+1][j][i];
		}
		else{
    			m3 = fz[k+2][j][i];
    			m2 = fz[k+1][j][i];
    			m1 = fz[k][j][i]  ;
    			p1 = fz[k-1][j][i];
    			p2 = fz[k-2][j][i];
		}

    		/* Candidate stencils */
    		double f1, f2, f3;
    		f1 = (2*one_sixth)*m1 + (5*one_sixth)*p1 - (one_sixth)*p2;
    		f2 = (-one_sixth)*m2 + (5.0*one_sixth)*m1 + (2*one_sixth)*p1;
    		f3 = (2*one_sixth)*m3 - (7.0*one_sixth)*m2 + (11.0*one_sixth)*m1;

    		/* Candidate stencils and their optimal weights*/
		double w1,w2,w3;
  		WENOParameters	WenoParams;
  		double 		c1,c2,c3;
  		ierr 		= 	SetWENOParams(&WenoParams);	if (ierr) return (ierr);
  		eps		=	WenoParams.eps;
  		p_weno		=	WenoParams.p; 
  		c1 = 0.3;
  		c2 = 0.6;
  		c3 = 0.1;


  		if (WenoParams.no_limiting) {
			w1 = c1;
   			w2 = c2;
      			w3 = c3;
    		} 
  		else {
      		// Smoothness indicators 
      		double b1, b2, b3;
      		b1 = thirteen_by_twelve*(m1-2*p1+p2)*(m1-2*p1+p2)
           		+ one_fourth*(3*m1-4*p1+p2)*(3*m1-4*p1+p2);
      		b2 = thirteen_by_twelve*(m2-2*m1+p1)*(m2-2*m1+p1)
           		+ one_fourth*(m2-p1)*(m2-p1);
      		b3 = thirteen_by_twelve*(m3-2*m2+m1)*(m3-2*m2+m1) 
           		+ one_fourth*(m3-4*m2+3*m1)*(m3-4*m2+3*m1);
  
      		// This parameter is needed for Borges' or Yamaleev-Carpenter
         	//implementation of non-linear weights or to obtain a fifth-
         	//order interpolation scheme with no limiting. 
      		double tau;
      		if (WenoParams.borges) {
        		tau = absolute(b3 - b1);
      		} 
		else if (WenoParams.yc) {
        		tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);
      		} 
		else {
        		tau = 0;
      		}
  
      		// Definiting the non-linear weights 
      		double a1, a2, a3;
      		if (WenoParams.borges || WenoParams.yc) {
        		a1 = c1 * (1.0 + raiseto(tau/(b1+eps),p_weno));
        		a2 = c2 * (1.0 + raiseto(tau/(b2+eps),p_weno));
        		a3 = c3 * (1.0 + raiseto(tau/(b3+eps),p_weno));
      		} 
		else {
        		a1 = c1 / raiseto(b1+eps,p_weno);
       			a2 = c2 / raiseto(b2+eps,p_weno);
        		a3 = c3 / raiseto(b3+eps,p_weno);
      		}

      		// Convexity 
      		double a_sum_inv;
      		a_sum_inv = 1.0 / (a1 + a2 + a3);
      		w1 = a1 * a_sum_inv;
      		w2 = a2 * a_sum_inv;
     		w3 = a3 * a_sum_inv;
  
      		// Mapping the weights, if needed 
      		if (WenoParams.mapped) {
        		a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1));
        		a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2));
        		a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3));
        		a_sum_inv = 1.0 / (a1 + a2 + a3);
        		w1 = a1 * a_sum_inv;
        		w2 = a2 * a_sum_inv;
        		w3 = a3 * a_sum_inv;

      			}
  
    		}

		fzi[k][j][i] = w1*f1 + w2*f2 + w3*f3;
	    }
	}
  }

  DMDAVecRestoreArray(da,Fz,&fz);
  DMDAVecRestoreArray(da,uwZ,&upwZ);   



  /* Calculating cell-centered flux derivative */

  VecGetDM(R,&da);
/*
  DMDAGetLocalInfo(da,&info);
  is = info.xs;
  ie = info.xs + info.xm;
  js = info.ys;
  je = info.ys + info.ym;
  ks = info.zs;
  ke = info.zs + info.zm;
*/

  DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
  is=xs;	ie=xs+xm;
  js=ys;	je=ys+ym;
  ks=zs;	ke=zs+zm;

  DMDAVecGetArray(da,R,&r);
  dxinv = 1.0/dx;
  dyinv = 1.0/dy;
  dzinv = 1.0/dz;


  for (i=is; i<ie; i++) {
      for (j=js; j<je; j++) {
	  for (k=ks; k<ke; k++) {
    		r[k][j][i] = -dxinv * (fxi[k][j][i+1] - fxi[k][j][i])
			     -dyinv * (fyi[k][j+1][i] - fyi[k][j][i])
			     -dzinv * (fzi[k+1][j][i] - fzi[k][j][i]);

	  }
      }
  }

  DMDAVecRestoreArray(da,R,&r);


  DMDAVecRestoreArray(da,fxI,&fxi);
  VecDestroy(&fxI);

  DMDAVecRestoreArray(da,fyI,&fyi);
  VecDestroy(&fyI);

  DMDAVecRestoreArray(da,fzI,&fzi);
  VecDestroy(&fzI);


  return(0);
}
