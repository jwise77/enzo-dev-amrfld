/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY HIERARCHY DEPTH AND LOCATION)
/
/  written by: Daniel R. Reynolds
/  date:       September, 2011
/
/  PURPOSE: Cell flagging routine for mesh refinement, that refines 
/           based on location (creates a hierarchy in the center of 
/           each root grid tile).
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


int grid::FlagCellsToBeRefinedForWeakScaling(int level, HierarchyEntry *MyHierarchyEntry)
{
  // declarations 
  int i, j, k, dim;
 
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

  // determine center of root-grid tile ancestor
  FLOAT RootTileCenter[MAX_DIMENSION];
  FLOAT RootTileWidth[MAX_DIMENSION];
  HierarchyEntry *Ancestor = MyHierarchyEntry;
  while (Ancestor->ParentGrid != NULL) 
    Ancestor = Ancestor->ParentGrid;
  for (i=0; i<MAX_DIMENSION; i++) {
    RootTileCenter[i] = 0.5*(Ancestor->GridData->GetGridRightEdge(i)
			   + Ancestor->GridData->GetGridLeftEdge(i));
    RootTileWidth[i] = Ancestor->GridData->GetGridRightEdge(i)
                     - Ancestor->GridData->GetGridLeftEdge(i);
  }

  // determine refinement region on this level, assuming that each level 
  // refines only the center 50% of the root grid tile in each direction
  FLOAT RefineRegionLeftEdge[MAX_DIMENSION];
  FLOAT RefineRegionRightEdge[MAX_DIMENSION];
  for (i=0; i<MAX_DIMENSION; i++) {
    RefineRegionLeftEdge[i]  = RootTileCenter[i] - RootTileWidth[i]*POW(0.5, level+2);
    RefineRegionRightEdge[i] = RootTileCenter[i] + RootTileWidth[i]*POW(0.5, level+2);
  }
 
  // iterate over the grid, marking cells as needed
  FLOAT cell_center[3];
  for (k=GridStartIndex[2]; k<=GridEndIndex[2]; k++) {
    cell_center[2] = GridLeftEdge[2] + (k-GridStartIndex[2]+0.5)*CellWidth[2][0];
    
    for (j=GridStartIndex[1]; j<=GridEndIndex[1]; j++) {
      cell_center[1] = GridLeftEdge[1] + (j-GridStartIndex[1]+0.5)*CellWidth[1][0];

      for (i=GridStartIndex[0]; i<=GridEndIndex[0]; i++) {
	cell_center[0] = GridLeftEdge[0] + (i-GridStartIndex[1]+0.5)*CellWidth[0][0];

	if ((cell_center[2] >= RefineRegionLeftEdge[2]) &&
	    (cell_center[1] >= RefineRegionLeftEdge[1]) &&
	    (cell_center[0] >= RefineRegionLeftEdge[0]) &&
	    (cell_center[2] <= RefineRegionRightEdge[2]) &&
	    (cell_center[1] <= RefineRegionRightEdge[1]) &&
	    (cell_center[0] <= RefineRegionRightEdge[0])) {

	  FlaggingField[i + (j + k*GridDimension[1])*GridDimension[0]] += 1;
	  NumberOfFlaggedCells++;
	}
      }
    }
  }  // end loop over cells
  
  return NumberOfFlaggedCells;
}
