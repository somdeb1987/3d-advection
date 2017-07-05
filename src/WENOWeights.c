#include <stdio.h>
#include <Definitions.h>
#include <Functions.h>
#include <DataTypes.h>

double absolute(double);
double raiseto(double,double);

ErrorType WENOWeights(	double w1, double w2, double w3,
			double m3, double m2, double m1, 
			double p1, double p2)
{
  WENOParameters	WenoParams;
  int 			ierr;	
  double 		eps;	//= 1e-6;
  double 		p_weno;	//= 2;
  double 		c1,c2,c3;
	

  ierr 		= 	SetWENOParams(&WenoParams);	if (ierr) return (ierr);


  eps		=	WenoParams.eps;
  p_weno	=	WenoParams.p; 
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



  return(0);
}
