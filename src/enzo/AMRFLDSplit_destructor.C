/*****************************************************************************
 *                                                                           *
 * Copyright 2010 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, Destructor routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:  
/
/  PURPOSE: Frees all memory allocated for the implicit FLD problem.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"



AMRFLDSplit::~AMRFLDSplit()
{

//   if (debug)  printf("Entering AMRFLDSplit::destructor routine\n");

  // delete boundary condition arrays
  int bin, i, j;
  for (bin=0; bin<MAX_FLD_FIELDS; bin++)
    for (i=0; i<3; i++)
      for (j=0; j<2; j++) 
	if (BdryVals[bin][i][j] != NULL)  
	  delete [] BdryVals[bin][i][j];

  // delete amrsolve parameters pointer
  delete amrsolve_params;

}
#endif
