#include <stdio.h>
#include <math.h>
#include <DataTypes.h>

ErrorType SetupGrid(DM da,void* g)
{
  StructGrid	*grid = (StructGrid*) g;
  DMDALocalInfo	info;


  Real dcsi,  csimin, csimax,  lengthcsi;
  Real deta,  etamin, etamax,  lengtheta;
  Real dzeta, zetamin,zetamax, lengthzeta;

  DMDAGetLocalInfo(da,&info);


  /* calculating dx */

  csimin 	= grid->csimin;
  csimax 	= grid->csimax;
  lengthcsi	= csimax - csimin;
  dcsi 		= lengthcsi / (grid->Nx -1);
  grid->dcsi 	= dcsi;

  /* calculating dy */

  etamin 	= grid->etamin;
  etamax 	= grid->etamax;
  lengtheta	= etamax - etamin;
  deta 		= lengtheta / (grid->Ny -1);
  grid->deta 	= deta;

  /* calculating dy */

  zetamin 	= grid->zetamin;
  zetamax 	= grid->zetamax;
  lengthzeta	= zetamax - zetamin;
  dzeta		= lengthzeta / (grid->Nz -1);
  grid->dzeta 	= dzeta;

  return(0);
}
