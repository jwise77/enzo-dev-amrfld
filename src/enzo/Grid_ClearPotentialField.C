/***********************************************************************
/
/  GRID CLASS (ALLOCATE AND CLEAR THE POTENTIAL FIELD)
/
/  written by: Daniel R. Reynolds
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
 
int grid::ClearPotentialField()
{
 
  // Return is this is not the right processor
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  // Error check
  if (GravitatingMassFieldCellSize == FLOAT_UNDEFINED) 
    ENZO_FAIL("GravitatingMassField uninitialized.\n");
  
  // compute size of the gravitating mass field
  int i, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldDimension[dim];
 
  // allocate and clear the field
  if (PotentialField == NULL) {
    PotentialField = new float[size];
    for (i=0; i<size; i++)  PotentialField[i] = 0.0;
  }
  if (PotentialField == NULL) 
    ENZO_FAIL("malloc error (out of memory?)\n");
 
 
  return SUCCESS;
}
