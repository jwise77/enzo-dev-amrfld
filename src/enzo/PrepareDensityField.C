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
 
int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);

int CommunicationBufferPurge(void);

int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
 
int PrepareGravitatingMassField1(HierarchyEntry *Grid);


#ifdef FAST_SIB
int PrepareGravitatingMassField2a(HierarchyEntry *Grid, int grid1,
				 SiblingGridList SiblingList[],
				 TopGridData *MetaData, int level,
				 FLOAT When);

int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   SiblingGridList SiblingList[],
				   HierarchyEntry *Grids[], int NumberOfGrids);

int CommunicationSolvePotential(LevelHierarchyEntry *LevelArray[], 
				SiblingGridList SiblingList[],
				int level, int reallevel, 
				TopGridData *MetaData, 
				HierarchyEntry **Grids,
				int NumberOfGrids,
				FLOAT EvaluateTime);

int AMRGravitySolve(LevelHierarchyEntry * LevelArray[],	
		    HierarchyEntry **Grids, TopGridData *MetaData,
		    SiblingGridList SiblingList[],
		    int NumberOfGrids, int level_coarse, 
		    int level_fine, FLOAT EvaluateTime);

#else

int PrepareGravitatingMassField2a(HierarchyEntry *Grid, TopGridData *MetaData,
				 LevelHierarchyEntry *LevelArray[], int level,
				 FLOAT When);

int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], int NumberOfGrids);

int CommunicationSolvePotential(LevelHierarchyEntry *LevelArray[], 
				int level, int reallevel, 
				TopGridData *MetaData, 
				HierarchyEntry **Grids,
				int NumberOfGrids,
				FLOAT EvaluateTime);

int AMRGravitySolve(LevelHierarchyEntry * LevelArray[],	
		    HierarchyEntry **Grids, TopGridData *MetaData,
		    int NumberOfGrids, int level_coarse, 
		    int level_fine, FLOAT EvaluateTime);
#endif

int PrepareGravitatingMassField2b(HierarchyEntry *Grid, int level);
 
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);


extern int CopyPotentialFieldAverage;
 
#define GRIDS_PER_LOOP 100000



 
#ifdef FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			SiblingGridList SiblingList[],
			int level, TopGridData *MetaData, FLOAT When)
#else   // !FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			int level, TopGridData *MetaData, FLOAT When)
#endif  // end FAST_SIB
{

  /* Return if this does not concern us */
  if (!SelfGravity) return SUCCESS;
 
  LCAPERF_START("PrepareDensityField");

  int grid1, grid2, StartGrid, EndGrid;
 
  /* Set the time for evaluation of the fields, etc. */
 
  FLOAT EvaluateTime = LevelArray[level]->GridData->ReturnTime() +
    When*LevelArray[level]->GridData->ReturnTimeStep();
 
  /* If level is above MaximumGravityRefinementLevel, then just
     update the gravity at the MaximumGravityRefinementLevel. */
 
  int reallevel = level;
  level = min(level, MaximumGravityRefinementLevel);
 
  /* Create an array (Grids) of all the grids. */
 
  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);

  /************************************************************************/
  /* Grids: Deposit particles in their GravitatingMassFieldParticles.
     (Do a batch of grids at a time; this is a loop over the batches)
  */

  if (traceMPI) 
    fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField (Send)\n");

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  TIME_MSG("Depositing particle mass field");
  LCAPERF_START("DepositParticleMassField");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      DepositParticleMassField(Grids[grid1], EvaluateTime);

#ifdef FORCE_MSG_PROGRESS 
    CommunicationBarrier();
#endif

    if (traceMPI) 
      fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField"
	      " (Receive)\n");
 
    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      DepositParticleMassField(Grids[grid1], EvaluateTime);

    /* Finally, receive the data and process it. */
    
    CommunicationReceiveHandler();

  } // ENDFOR grid batches
  LCAPERF_STOP("DepositParticleMassField");
    

#ifdef FORCE_BUFFER_PURGE
  CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  /******************************************************************/
  /* Grids: compute the GravitatingMassField (baryons & particles). */
  /*   This is now split into two section. */
 
  if (traceMPI) 
    fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF1 (send)\n", 
	    MyProcessorNumber);
 
  TIME_MSG("PrepareGravitatingMassField1");
  LCAPERF_START("PrepareGravitatingMassField1");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* ----- section 1 ---- */
    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
 
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField1(Grids[grid1]);

    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField1(Grids[grid1]);

    /* Finally, receive the data and process it. */
    
    CommunicationReceiveHandler();

  } // ENDFOR grid batches
  LCAPERF_STOP("PrepareGravitatingMassField1");


#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  if (traceMPI) 
    fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF2 (receive)\n", 
	    MyProcessorNumber);
 
  TIME_MSG("PrepareGravitatingMassField2");
  LCAPERF_START("PrepareGravitatingMassField2a");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* ----- section 2 ---- */
    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
#ifdef BITWISE_IDENTICALITY
    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
#else
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
#endif
      
#ifdef FAST_SIB
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2a(Grids[grid1], grid1, SiblingList,
				    MetaData, level, When);
#else
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2a(Grids[grid1], MetaData, LevelArray,
				    level, When);
#endif

#ifndef BITWISE_IDENTICALITY
    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
#ifdef FAST_SIB
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2a(Grids[grid1], grid1, SiblingList,
				   MetaData, level, When);
#else
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2a(Grids[grid1], MetaData, LevelArray,
				   level, When);
#endif

    CommunicationReceiveHandler();
#endif /* BITWISE_IDENTICALITY */

  } // ENDFOR grid batches
  LCAPERF_STOP("PrepareGravitatingMassField2a");

#ifdef FORCE_BUFFER_PURGE
  CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif
 
  /************************************************************************/
  LCAPERF_START("PrepareGravitatingMassField2b");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* ----- section 2 ---- */
    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
      
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2b(Grids[grid1], level);

    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2b(Grids[grid1], level);

    CommunicationReceiveHandler();

  } // ENDFOR grid batches
  LCAPERF_STOP("PrepareGravitatingMassField2b");

  /************************************************************************/
  /* Copy overlapping mass fields to ensure consistency and B.C.'s. */
 
  //  if (level > 0)
 
  if (traceMPI) 
    fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): COMF1 (send)\n", 
	    MyProcessorNumber);
 
  TIME_MSG("CopyOverlappingMassField");
  LCAPERF_START("CopyOverlappingMassField");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
      
#ifdef FAST_SIB
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	Grids[grid1]->GridData->
	  CheckForOverlap(SiblingList[grid1].GridList[grid2],

			  MetaData->LeftFaceBoundaryCondition,
			  MetaData->RightFaceBoundaryCondition,
			  &grid::CopyOverlappingMassField);
#else
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	Grids[grid1]->GridData->
	  CheckForOverlap(Grids[grid2]->GridData,
			  MetaData->LeftFaceBoundaryCondition,
			  MetaData->RightFaceBoundaryCondition,
			  &grid::CopyOverlappingMassField);
#endif

    CommunicationDirection = COMMUNICATION_SEND;
#ifdef FAST_SIB
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	Grids[grid1]->GridData->
	  CheckForOverlap(SiblingList[grid1].GridList[grid2],
			  MetaData->LeftFaceBoundaryCondition,
			  MetaData->RightFaceBoundaryCondition,
			  &grid::CopyOverlappingMassField);
#else
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	Grids[grid1]->GridData->
	  CheckForOverlap(Grids[grid2]->GridData,
			  MetaData->LeftFaceBoundaryCondition,
			  MetaData->RightFaceBoundaryCondition,
			  &grid::CopyOverlappingMassField);
#endif

    CommunicationReceiveHandler();

  } // ENDFOR grid batches
  LCAPERF_STOP("CopyOverlappingMassField");

#ifdef FORCE_BUFFER_PURGE
  CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  /************************************************************************/
  // Compute the potential, either using original Enzo approach, or new
  // approach based on AMRsolve.
  if (SelfGravityConsistent) {
    int level_coarse = 0;
    int level_fine = level;
#ifdef FAST_SIB
    if (AMRGravitySolve(LevelArray, Grids, MetaData, SiblingList, NumberOfGrids, 
			level_coarse, level_fine, EvaluateTime) != SUCCESS) 
      ENZO_FAIL("Error in AMRGravitySolve");
#else
    if (AMRGravitySolve(LevelArray, Grids, MetaData, NumberOfGrids, 
			level_coarse, level_fine, EvaluateTime) != SUCCESS) 
      ENZO_FAIL("Error in AMRGravitySolve");
#endif
  } else {
#ifdef FAST_SIB
    if (CommunicationSolvePotential(LevelArray, SiblingList, level, 
				    reallevel, MetaData, Grids, NumberOfGrids, 
				    EvaluateTime) != SUCCESS)
      ENZO_FAIL("Error in CommunicationSolvePotential");
#else
    if (CommunicationSolvePotential(LevelArray, level, reallevel, MetaData, 
				    Grids, NumberOfGrids, EvaluateTime) != SUCCESS)
      ENZO_FAIL("Error in CommunicationSolvePotential");
#endif 
  }

  // --------------------------------------------------
  // MEMORY LEAK FIX
  //
  // valgrind error: "1,388,304 (67,352 direct, 1,320,952 indirect)
  // bytes in 130 blocks are definitely lost in loss record 22 of 46"
  //
  // Adding missing delete [] () for Grids[] allocated in
  // GenerateGridArray()
  // --------------------------------------------------

  delete [] Grids;

  // --------------------------------------------------

  LCAPERF_STOP("PrepareDensityField");
  return SUCCESS;

}
 
