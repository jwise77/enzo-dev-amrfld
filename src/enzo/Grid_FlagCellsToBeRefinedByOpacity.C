/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY SLOPE OF OPACITY)
/
/  written by: Daniel R. Reynolds
/  date:       October, 2011
/
/  PURPOSE: Cell flagging routine for mesh refinement, that computes 
/           the magnitude of the gradient of the opacity.  This 
/           value is compared against a threshold, held in the second
/           entry in MinimumSlopeForRefinement array.  A typical value 
/           for this threshold is 1.0
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


int grid::FlagCellsToBeRefinedByOpacity()
{
  // declarations 
  int i, j, k, index, dim;
 
  // Return if this grid is not on this processor
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  // error check
  if (FlaggingField == NULL)
    ENZO_FAIL("Flagging Field is undefined.");
 
  // Make sure quantities are defined at least to dim 3 (if GridRank < 3)
  for (dim=GridRank; dim<3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }

  // compute total grid size 
  int size=1;
  for (dim=0; dim<GridRank; dim++)   size *= GridDimension[dim];
 
  // initialize return value
  int NumberOfFlaggedCells = 0;
 
  // only take gradient over active dimensions
  int xOffset = (GridDimension[0] > 1) ? 1 : 0;
  int yOffset = (GridDimension[1] > 1) ? GridDimension[0]+1 : 0;
  int zOffset = (GridDimension[2] > 1) ? GridDimension[0]*GridDimension[1]+1 : 0;

  // set opacity constants
  float HIconst   = IntegralRadiationSpectrumHI   / IntegralRadiationSpectrum;
  float HeIconst  = IntegralRadiationSpectrumHeI  / IntegralRadiationSpectrum / 4.0;
  float HeIIconst = IntegralRadiationSpectrumHeII / IntegralRadiationSpectrum / 4.0;

  // get chemistry array pointers
  float *nHI=NULL, *nHeI=NULL, *nHeII=NULL;
  nHI = this->AccessHIDensity();
  if (nHI == NULL)
    ENZO_FAIL("Could not find HIDensity field.");
  nHeI = this->AccessHeIDensity();
  nHeII = this->AccessHeIIDensity();
  int nchem = 1;
  if ((nHeI != NULL) && (nHeII != NULL))  nchem = 3;

  // create and fill temporary opacity array
  float *kappa = new float[GridDimension[0]*GridDimension[1]*GridDimension[2]];
  for (k=0; k<GridDimension[0]*GridDimension[1]*GridDimension[2]; k++)
    kappa[k] = HIconst*nHI[k];
  if (nchem == 3) {
    for (k=0; k<GridDimension[0]*GridDimension[1]*GridDimension[2]; k++)
      kappa[k] += HeIconst*nHeI[k];
    for (k=0; k<GridDimension[0]*GridDimension[1]*GridDimension[2]; k++)
      kappa[k] += HeIIconst*nHeII[k];
  }

  // use second entry in MinimumSlopeForRefinement array to determine refinement 
  // threshold; if left unset, just return
  float MinRefinementSlope = MinimumSlopeForRefinement[1];
  if (MinRefinementSlope == 0.0)  return NumberOfFlaggedCells;

  // iterate over the grid, marking cells as needed
  // float maxgrad=0.0, mingrad=0.0, avggrad=0.0;
  // float maxkappa=0.0, minkappa=0.0, avgkappa=0.0;
  float gradx, grady, gradz, gradnorm;
  float atol = 1.0;
  if (nchem == 1)  atol *= HIconst;
  else             atol *= (HIconst + HeIconst + HeIIconst);
  for (k=GridStartIndex[2]; k<=GridEndIndex[2]; k++)
    for (j=GridStartIndex[1]; j<=GridEndIndex[1]; j++)
      for (i=GridStartIndex[0]; i<=GridEndIndex[0]; i++) {
	index = i + (j + k*GridDimension[1])*GridDimension[0];
	gradx = (kappa[index+xOffset] - kappa[index-xOffset]) /
	  max(fabs(kappa[index]), atol);
	grady = (kappa[index+yOffset] - kappa[index-yOffset]) /
	  max(fabs(kappa[index]), atol);
	gradz = (kappa[index+zOffset] - kappa[index-zOffset]) /
	  max(fabs(kappa[index]), atol);
	gradnorm = sqrt(gradx*gradx + grady*grady + gradz*gradz);
	// maxgrad = (gradnorm > maxgrad) ? gradnorm : maxgrad;
	// mingrad = (gradnorm < mingrad) ? gradnorm : mingrad;
	// avggrad += gradnorm;
	// maxkappa = (kappa[index] > maxkappa) ? kappa[index] : maxkappa;
	// minkappa = (kappa[index] < minkappa) ? kappa[index] : minkappa;
	// avgkappa += kappa[index];
	if (gradnorm > MinRefinementSlope) {
	  FlaggingField[index] += 1;
	  NumberOfFlaggedCells++;
	}
      }  // end loop over cells
  // avggrad /= ((GridEndIndex[2]-GridStartIndex[2]+1) *
  // 	      (GridEndIndex[1]-GridStartIndex[1]+1) * 
  // 	      (GridEndIndex[0]-GridStartIndex[0]+1) );
  // avgkappa /= ((GridEndIndex[2]-GridStartIndex[2]+1) *
  // 	       (GridEndIndex[1]-GridStartIndex[1]+1) * 
  // 	       (GridEndIndex[0]-GridStartIndex[0]+1) );
  // printf("FlagCellsOpacity: grad max = %g, min = %g, avg = %g, NFlagCells = %i\n",
  // 	 maxgrad,mingrad,avggrad,NumberOfFlaggedCells);
  // printf("FlagCellsOpacity: kappa max = %g, min = %g, avg = %g, atol = %g\n",
  // 	 maxkappa,minkappa,avgkappa,atol);

  delete [] kappa;
  return NumberOfFlaggedCells;
}
