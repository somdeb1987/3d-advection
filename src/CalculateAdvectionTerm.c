#include <stdio.h>
#include <DataTypes.h>
#include <Functions.h>

ErrorType CalculateAdvectionTerm(DM da,Vec Local_u, Vec Local_v, Vec Local_w, 
				Vec Local_p, Vec globalAp, void *g)
{
  StructGrid    	*grid	= (StructGrid*)	g;
  PetscErrorCode  	ierr;
  Vec			Local_uc,Local_vc,Local_wc;
  Vec			upwind_csi,upwind_eta,upwind_zeta;
  Vec			p_Fcsi,p_Feta,p_Fzeta;


  /* Contravarient velocity field */

  DMCreateLocalVector(da,&Local_uc);
  DMCreateLocalVector(da,&Local_vc);
  DMCreateLocalVector(da,&Local_wc);

  /* Cell-centered f value */
  DMCreateLocalVector(da,&p_Fcsi);
  DMCreateLocalVector(da,&p_Feta);
  DMCreateLocalVector(da,&p_Fzeta); 


  /* Upwind Parameter */
  DMCreateLocalVector(da,&upwind_csi);
  DMCreateLocalVector(da,&upwind_eta);
  DMCreateLocalVector(da,&upwind_zeta);


  // Compute advective term //

  // Calculate contravarient velocity //
  ierr = CalculateUContravarient(da,Local_u,Local_uc,grid);
  ierr = CalculateVContravarient(da,Local_v,Local_vc,grid);
  ierr = CalculateWContravarient(da,Local_w,Local_wc,grid);

  // Calculate cell-centered flux f(u) (in u_t + f(u)_x = 0) //
  ierr = CalculateCsiFlux (da,Local_p,Local_uc,p_Fcsi,grid);  if (ierr) return(ierr);
  ierr = CalculateEtaFlux (da,Local_p,Local_vc,p_Feta,grid);  if (ierr) return(ierr);
  ierr = CalculateZetaFlux(da,Local_p,Local_wc,p_Fzeta,grid);  if (ierr) return(ierr);


  // Calculate upwind flag //
  ierr = CalculateCsiUpWind (da,Local_uc,upwind_csi);  if (ierr) return(ierr);
  ierr = CalculateEtaUpWind (da,Local_vc,upwind_eta);  if (ierr) return(ierr);
  ierr = CalculateZetaUpWind(da,Local_wc,upwind_zeta);  if (ierr) return(ierr);




 // Reconstruct the flux derivative //
    // Non-compact scheme: f_x = Bf //
  ierr = FifthOrderWENO(p_Fcsi,p_Feta,p_Fzeta,globalAp,upwind_csi,upwind_eta,upwind_zeta,grid);
  if (ierr) return(ierr);



 // Finished computing advective term //

  VecDestroy(&upwind_csi);
  VecDestroy(&upwind_eta);
  VecDestroy(&upwind_zeta);


  VecDestroy(&Local_uc);
  VecDestroy(&Local_vc);
  VecDestroy(&Local_wc);

  VecDestroy(&p_Fcsi);
  VecDestroy(&p_Feta);
  VecDestroy(&p_Fzeta);


  return(0);
}
