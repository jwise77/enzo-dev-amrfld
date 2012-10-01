/***********************************************************************
/
/  PREPARE DENSITY FIELD (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       June, 1999
/  modifiedN:  Robert Harkness
/  date:       February, 2008
/
/ ======================================================================= 
/ This routine prepares the density field for all the grids on this level,
/ both particle and baryonic densities.  It also calculates the potential
/ field if this is level 0 (since this involves communication). 
/
/   This is part of a collection of routines called by EvolveLevel.
/   These have been optimized for enhanced message passing
/   performance by performing two passes -- one which generates
/   sends and the second which receives them.
/
/  modified: Robert Harkness, December 2007
/
************************************************************************/

#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "communication.h"
#include "CommunicationUtilities.h"

/* function prototypes */
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
 
#ifdef FAST_SIB
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   SiblingGridList SiblingList[],
				   HierarchyEntry *Grids[], int NumberOfGrids);
#else
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], int NumberOfGrids);
#endif

extern int CopyPotentialFieldAverage;
 
#define GRIDS_PER_LOOP 100000


 
#ifdef FAST_SIB
int CommunicationSolvePotential(LevelHierarchyEntry *LevelArray[],
				SiblingGridList SiblingList[],
				int level, int reallevel, 
				TopGridData *MetaData,
				HierarchyEntry **Grids,
				int NumberOfGrids,
				FLOAT EvaluateTime)
#else   // !FAST_SIB
int CommunicationSolvePotential(LevelHierarchyEntry *LevelArray[],
				int level, int reallevel, 
				TopGridData *MetaData,
				HierarchyEntry **Grids,
				int NumberOfGrids,
				FLOAT EvaluateTime)
#endif  // end FAST_SIB
{

  /* Return if this does not concern us */
  if (!SelfGravity) return SUCCESS;

  int grid1, grid2, StartGrid, EndGrid;
 
  LCAPERF_START("CommunicationSolvePotential");

  /************************************************************************/
  /* Compute the potential for the top grid. */
 
  if (level == 0) {
    TIME_MSG("ComputePotentialFieldLevelZero");
    LCAPERF_START("ComputePotentialFieldLevelZero");
    if (traceMPI) 
      fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): CPFLZero "
	      "(send-receive)\n", MyProcessorNumber);
#ifdef FAST_SIB
    ComputePotentialFieldLevelZero(MetaData, SiblingList,
				   Grids, NumberOfGrids);
#else
    ComputePotentialFieldLevelZero(MetaData, Grids, NumberOfGrids);
#endif
    LCAPERF_STOP("ComputePotentialFieldLevelZero");
  }
       
  /************************************************************************/
  /* Compute a first iteration of the potential and share BV's. */
 
  int iterate;
  if (level > 0) {
    LCAPERF_START("SolveForPotential");
    CopyPotentialFieldAverage = 1;
    for (iterate = 0; iterate < PotentialIterations; iterate++) {
      
      if (iterate > 0)
	CopyPotentialFieldAverage = 2;

 
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
	Grids[grid1]->GridData->SolveForPotential(level, EvaluateTime);
	if (CopyGravPotential)
	  Grids[grid1]->GridData->CopyPotentialToBaryonField();
      }
 
      if (traceMPI) fprintf(tracePtr, "ITPOT post-recv\n");
	
#ifdef FORCE_MSG_PROGRESS 
      CommunicationBarrier();
#endif

      TIME_MSG("CopyPotentialField");
      for (StartGrid = 0; StartGrid < NumberOfGrids; 
	   StartGrid += GRIDS_PER_LOOP) {
	EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);
  
#ifdef BITWISE_IDENTICALITY
	CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
#else
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
#endif
	CommunicationReceiveIndex = 0;
	CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
#ifdef FAST_SIB
	for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
 
	  //fprintf(stderr, "#SIBSend on cpu %"ISYM": %"ISYM"\n", MyProcessorNumber, SiblingList[grid1].NumberOfSiblings);
 
	  // for (grid2 = SiblingList[grid1].NumberOfSiblings-1; grid2 = 0; grid2--)
	  for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	    Grids[grid1]->GridData->
	      CheckForOverlap(SiblingList[grid1].GridList[grid2],
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition,
			      &grid::CopyPotentialField);
	    
	  grid2 = grid1;
	  Grids[grid1]->GridData->
	    CheckForOverlap(Grids[grid2]->GridData,
			    MetaData->LeftFaceBoundaryCondition,
			    MetaData->RightFaceBoundaryCondition,
			    &grid::CopyPotentialField);
	  
	} // ENDFOR grid1
#else
	for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	  for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	    Grids[grid1]->GridData->
	      CheckForOverlap(Grids[grid2]->GridData,
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition,
			      &grid::CopyPotentialField);
#endif

#ifndef BITWISE_IDENTICALITY
#ifdef FORCE_MSG_PROGRESS 
	CommunicationBarrier();
#endif

	if (traceMPI) fprintf(tracePtr, "ITPOT send\n");
 
	CommunicationDirection = COMMUNICATION_SEND;

 
#ifdef FAST_SIB
	for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
 
	  //fprintf(stderr, "#SIBRecv on cpu %"ISYM": %"ISYM"\n", MyProcessorNumber, SiblingList[grid1].NumberOfSiblings);
 
	  // for (grid2 = SiblingList[grid1].NumberOfSiblings-1; grid2 = 0; grid2--)
	  for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	    Grids[grid1]->GridData->
	      CheckForOverlap(SiblingList[grid1].GridList[grid2],
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition,
			      &grid::CopyPotentialField);
 
	  grid2 = grid1;
	  Grids[grid1]->GridData->
	    CheckForOverlap(Grids[grid2]->GridData,
			    MetaData->LeftFaceBoundaryCondition,
			    MetaData->RightFaceBoundaryCondition,
			    &grid::CopyPotentialField);
 
	} // ENDFOR grid1
#else
	for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	  for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	    Grids[grid1]->GridData->
	      CheckForOverlap(Grids[grid2]->GridData,
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition,
			      &grid::CopyPotentialField);
#endif

	CommunicationReceiveHandler();
#endif

      } // ENDFOR grid batches
    } // ENDFOR iterations
    CopyPotentialFieldAverage = 0;
    LCAPERF_STOP("SolveForPotential");
  } // ENDIF level > 0
  
  /* if level > MaximumGravityRefinementLevel, then do final potential
     solve (and acceleration interpolation) here rather than in the main
     EvolveLevel since it involves communications. */
  
  if (reallevel > MaximumGravityRefinementLevel) {
 
    /* compute potential and acceleration on coarser level [LOCAL]
       (but only if there is at least a subgrid -- it should be only
       if there is a subgrrid on reallevel, but this is ok). */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      if (Grids[grid1]->NextGridNextLevel != NULL) {
	Grids[grid1]->GridData->SolveForPotential(MaximumGravityRefinementLevel);
	if (CopyGravPotential)
	  Grids[grid1]->GridData->CopyPotentialToBaryonField();
	else
	  Grids[grid1]->GridData->ComputeAccelerationField
	    ((HydroMethod == Zeus_Hydro) ? DIFFERENCE_TYPE_STAGGERED : 
	     DIFFERENCE_TYPE_NORMAL, MaximumGravityRefinementLevel);
      }
 
    /* Interpolate potential for reallevel grids from coarser grids. */
 
    if (!CopyGravPotential) {
 
      int Dummy, GridCount;
      LevelHierarchyEntry *Temp, *LastTemp;
      HierarchyEntry *Temp3;
      LevelHierarchyEntry *FirstTemp = LevelArray[reallevel];
	
#ifdef FORCE_MSG_PROGRESS 
      CommunicationBarrier();
#endif

      do {

	GridCount = 0;
	CommunicationDirection = COMMUNICATION_POST_RECEIVE;
	CommunicationReceiveIndex = 0;
	CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
	Temp = FirstTemp;
	while (Temp != NULL && GridCount++ < GRIDS_PER_LOOP) {
	  Temp3 = Temp->GridHierarchyEntry;
	  for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	    Temp3 = Temp3->ParentGrid;
	  Temp->GridData->InterpolateAccelerations(Temp3->GridData);
	  Temp = Temp->NextGridThisLevel;
	} // ENDWHILE
	LastTemp = Temp;

	CommunicationDirection = COMMUNICATION_SEND;
	Temp = FirstTemp;
	while (Temp != LastTemp) {
	  Temp3 = Temp->GridHierarchyEntry;
	  for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	    Temp3 = Temp3->ParentGrid;
	  Temp->GridData->InterpolateAccelerations(Temp3->GridData);
	  Temp = Temp->NextGridThisLevel;
	}
	FirstTemp = LastTemp;

	CommunicationReceiveHandler();

      } while (LastTemp != NULL);

    } // end:  if (!CopyGravPotential)
 
  } // end: if (reallevel > MaximumGravityRefinementLevel)

  // --------------------------------------------------

  LCAPERF_STOP("CommunicationSolvePotential");
  return SUCCESS;

}
 
