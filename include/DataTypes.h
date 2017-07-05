#include <petscsys.h>
#include <petscdmda.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>
#include <petscts.h>


/* Type Definitions */
typedef PetscReal 	    	Real;
typedef PetscInt 	      	Integer;
typedef PetscBool 	    	Boolean;
typedef Real		        VarArray;
typedef PetscErrorCode		ErrorType;

typedef struct mpi_variable_set {
  int			rank;
  int			nproc;
} MPIVariables;


typedef struct parameters_weno {
  /* Options related to the type of WENO scheme */
  PetscBool 		mapped;		    /* Use mapped weights?                                    */
  PetscBool 		borges;		    /* Use Borges' implementation of weights?                 */
  PetscBool 		yc;		        /* Use Yamaleev-Carpenter implementation of weights?      */
  PetscBool 		no_limiting;  /* Remove limiting -> 5th order polynomial interpolation  */
  double		eps;		      /* epsilon parameter                                      */
  double		p;		/* p parameter                                            */


} WENOParameters;


typedef struct parameters_field {

  Vec 			u,v,w;		    	/* cartesian velocity (solution vectors)	*/
  Vec			uc,vc,wc;		/* contravarient velocity			*/
  Vec			phi;		    	/* scalar solution vector			*/
  Vec 			Au,Av,Aw,Aphi;	    	/* Advection terms				*/
  Vec			Du,Dv,Dw,Dphi;	    	/* Diffusion terms				*/
  Vec 			Su,Sv,Sw,Sphi;	    	/* Source terms					*/
  Vec			RHSu,RHSv,RHSw,RHSphi;	/* contravarient velocity			*/

  Vec			RHS_PPE;
  Mat 			Amat_PPE;

} StructField;

typedef struct parameters_grid {

  Real			csimin,csimax;			/* Domain limits                              */
  Real			etamin,etamax;			/* Domain limits                              */
  Real			zetamin,zetamax;			/* Domain limits                              */
  Integer		Nx,Ny,Nz;			      	/* Number of grid points                      */

  Real			dcsi,deta,dzeta;

  Integer		LeftB,RightB;			/* Left & right boundary indices (not inputs) */
  Integer		BottomB,TopB;			/* Bottom & top boundary indices (not inputs) */
  Integer		FrontB,BackB;			/* Front & back boundary indices (not inputs) */

} StructGrid;

typedef struct parameters_time_integration {

  Integer		T;			      	/* Number of time steps                       */
  Real			tf;			      	/* Final time                                 */
  Real			cfl,dt;				/* Note: not inputs, these will be calculated */

} TSParameters;

typedef struct parameters_input_output_parameters {

  int output_interval,counter;
  /* IO options */
  Boolean		write_initial;		  	/* Write the initial solution to a file?  	*/
  Boolean		write_exact;		    	/* Write the exact solution to a file?    	*/
  Boolean		write_final;		    	/* Write the final solution to a file?    	*/
  char			initial_fname[256];		/* Filename for initial solution          	*/
  char			exact_fname[256];	  	/* Filename for initial solution          	*/
  char			final_fname[256];	  	/* Filename for final solution            	*/

} IOParameters;

typedef struct linear_system_solver {

  KSP	KSP;
  Mat	Amat;
  Vec   rhsX,solnX;

} SolveLinearSystem;


typedef struct parameters_solver {
  /* Some basic options */
  Integer		ic;			      	/* Choice of initial conditions 		*/
  Integer		spatial_scheme;			/* Choice of spatial reconstruction scheme  	*/


  IOParameters		IOParams;
  StructField		Field;
  StructGrid		Grid;
  TSParameters		TSParams;

  SolveLinearSystem	SolveLS;

  MPIVariables   	MPIVars;	    		/* contains MPI related variables         */
  WENOParameters 	WENOParams;        		/* WENO algorithm related options         */

  ErrorType (*ReconstructionRHS)(Vec,Vec, Vec, void*, double,double,Vec,Vec); /* Computes f_x = [B]f (u is the input vector)  */



} StructINS3D;

