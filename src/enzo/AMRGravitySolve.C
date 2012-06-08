#ifdef AMR_SOLVE

/// @file     AMRGravitySolve.C
/// @date     2010-12-29
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Daniel Reynolds (reynolds@smu.edu)
/// @brief    Interface between Enzo and the AMRsolve gravity solver
///
/// AMRGravitySolve serves as the single function interface between
/// Enzo and the AMRsolve linear solver package.  It solves the
/// Poisson equation on an Enzo AMR hierarchy between two levels specified
/// by the "level_coarse" and "level_fine" parameters.  
/// The right hand side and solution arrays are held in the Grid class, and
/// are accessed directly by the AMRsolve_Hierarchy::enzo_attach_grav() function.
///
/// Arguments:
///    LevelArray   -- Enzo's current LevelArray
///    level_coarse -- coarsest level involved in hierarchical solve
///    level_fine   -- finest level involved in hierarchical solve

#include "AMRsolve.h"

#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "ExternalBoundary.h"
#include "ProtoSubgrid.h"
#include "GridList.h"
#include "Grid.h"
#include "LevelHierarchy.h"
#include "Hierarchy.h"
#include "communication.h"
#include "CommunicationUtilities.h"


/* function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

#define GRIDS_PER_LOOP 100000




#ifdef FAST_SIB
int AMRGravitySolve(LevelHierarchyEntry * LevelArray[],	
		    HierarchyEntry **Grids, TopGridData *MetaData,
		    SiblingGridList SiblingList[],
		    int NumberOfGrids, int level_coarse, 
		    int level_fine, FLOAT EvaluateTime)
#else
int AMRGravitySolve(LevelHierarchyEntry * LevelArray[],	
		    HierarchyEntry **Grids, TopGridData *MetaData,
		    int NumberOfGrids, int level_coarse, 
		    int level_fine, FLOAT EvaluateTime)
#endif
{
  LCAPERF_START("amr_solve");

  // Ensure that Grid arrays for the PotentialField on this processor 
  // (and for these levels) have been allocated
  LevelHierarchyEntry *Temp;
  int RefinementFactor = RefineBy;
  for (int thislevel=level_coarse; thislevel<=level_fine; thislevel++)
    for (Temp=LevelArray[thislevel]; Temp; Temp=Temp->NextGridThisLevel) 
      Temp->GridHierarchyEntry->GridData->ClearPotentialField();
  




  // Set up AMRsolve MPI object
  pmpi = new AMRsolve_Mpi(MPI_COMM_WORLD);
  AMRsolve_Grid::set_mpi(*pmpi);

  // Set amrsolve parameters
  AMRsolve_Parameters* amrsolve_params = new AMRsolve_Parameters();
  amrsolve_params->set_defaults();
  //  amrsolve_params->set_parameter("solver","fac");
  amrsolve_params->set_parameter("solver","bicgstab");
  //  amrsolve_params->set_parameter("solver","bicgstab-boomer");  // [BROKEN]
  //  amrsolve_params->set_parameter("solver","gmres");
  //  amrsolve_params->set_parameter("solver","pfmg"); // [BROKEN]

  amrsolve_params->set_parameter("dump_a","true");
  amrsolve_params->set_parameter("dump_b","true");
  amrsolve_params->set_parameter("dump_x","true");

  amrsolve_params->set_parameter("solver_itmax","600");
  amrsolve_params->set_parameter("solver_printl", "1");

  // set preconditioning options for BiCGStab and GMRES solvers
  amrsolve_params->set_parameter("prec_itmax",  "20");
  amrsolve_params->set_parameter("prec_restol", "1.0e-6");
  amrsolve_params->set_parameter("prec_rlxtype","1");
  amrsolve_params->set_parameter("prec_npre",   "1");
  amrsolve_params->set_parameter("prec_npost",  "1");
  amrsolve_params->set_parameter("prec_printl", "0");
  amrsolve_params->set_parameter("prec_log",    "1");
 
  // if (debug) {
  //   printf("AMRGravitySolve, customized amrsolve parameters:\n");
  //   amrsolve_params->print();
  // }


  // Insert Enzo grids in this level into a AMRsolve hierarchy
  AMRsolve_Hierarchy* hierarchy = new AMRsolve_Hierarchy;

  LCAPERF_START("amrsolve_attach_grav");
  if (hierarchy->enzo_attach_grav(LevelArray, level_coarse, level_fine) != 0) {
    if (debug)   fprintf(stderr,"AMRGravitySolve: error in enzo_attach_grav\n");
    return FAIL;
  }
  LCAPERF_STOP("amrsolve_attach_grav");

  // Initialize the AMRsolve hierarchy
  AMRsolve_Domain domain(3, DomainLeftEdge, DomainRightEdge);
  //  bool is_periodic[] = {true, true, true};
  bool is_periodic[] = {false, false, false};

  LCAPERF_START("amrsolve_hierarchy");
  hierarchy->initialize(domain,*pmpi,is_periodic);
  LCAPERF_STOP("amrsolvehierarchy");

  // Compute the RHS scaling factor
  FLOAT a=1.0, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(EvaluateTime, &a, &dadt) == FAIL) 
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.\n");
  float f_scale = GravitationalConstant / a;

  // Initialize the AMRsolve linear system
  Eint32 sol_prec   = 0;   // [BROKEN]
  Eint32 zero_guess = 1;
  LCAPERF_START("amrsolve_matrix");
  AMRsolve_Hypre_Grav amrgravsolve(*hierarchy, *amrsolve_params, 
				   sol_prec, zero_guess);
  amrgravsolve.init_hierarchy();
  amrgravsolve.init_stencil();
  amrgravsolve.init_graph();
  amrgravsolve.init_elements(f_scale);
  LCAPERF_STOP("amrsolve_matrix");

  // Solve the linear system
  LCAPERF_START("amrsolve_solve");
  amrgravsolve.solve();
  Eflt64 finalresid = amrgravsolve.residual();
  Eint32 Sits = amrgravsolve.iterations();
  if (debug) printf("   lin resid = %.1e, its = %i\n", finalresid, Sits);
  LCAPERF_STOP("amrsolve_solve");

  // Display solver results
  if (amrgravsolve.evaluate() != 0) {
    fprintf(stderr,"AMRGravitySolve: could not achieve prescribed tolerance!\n");
    fprintf(stderr, "Printing current hierarchy:\n");
    if (debug)   hierarchy->print();
    return FAIL;
  }

  // increment Enzo potential field with amrsolve solution
  amrgravsolve.update_enzo();

  // Clean up
  hierarchy->enzo_detach();

  LCAPERF_STOP("amr_solve");

  delete hierarchy;
  hierarchy = NULL;



  // copy potential values to BaryonField if requested
  int grid1, grid2, StartGrid, EndGrid;
  if (CopyGravPotential)
    for (grid1=0; grid1<NumberOfGrids; grid1++) 
      Grids[grid1]->GridData->CopyPotentialToBaryonField();


  return SUCCESS;


 
  // share boundary values
#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  for (StartGrid=0; StartGrid<NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);
  
#ifdef BITWISE_IDENTICALITY
    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
#else
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
#endif
    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    
#ifdef FAST_SIB
    for (grid1=StartGrid; grid1<EndGrid; grid1++) {
      for (grid2=0; grid2<SiblingList[grid1].NumberOfSiblings; grid2++)
	Grids[grid1]->GridData->CheckForOverlap(SiblingList[grid1].GridList[grid2],
						MetaData->LeftFaceBoundaryCondition,
						MetaData->RightFaceBoundaryCondition,
						&grid::CopyPotentialField);
      grid2 = grid1;
      Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
					      MetaData->LeftFaceBoundaryCondition,
					      MetaData->RightFaceBoundaryCondition,
					      &grid::CopyPotentialField);
    } // ENDFOR grid1
#else
    for (grid1=StartGrid; grid1<EndGrid; grid1++)
      for (grid2=0; grid2<NumberOfGrids; grid2++)
	Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
						MetaData->LeftFaceBoundaryCondition,
						MetaData->RightFaceBoundaryCondition,
						&grid::CopyPotentialField);
#endif

    
#ifndef BITWISE_IDENTICALITY
    
#ifdef FORCE_MSG_PROGRESS 
    CommunicationBarrier();
#endif
    CommunicationDirection = COMMUNICATION_SEND;
 
#ifdef FAST_SIB
    for (grid1=StartGrid; grid1<EndGrid; grid1++) {
      for (grid2=0; grid2<SiblingList[grid1].NumberOfSiblings; grid2++)
	Grids[grid1]->GridData->CheckForOverlap(SiblingList[grid1].GridList[grid2],
						MetaData->LeftFaceBoundaryCondition,
						MetaData->RightFaceBoundaryCondition,
						&grid::CopyPotentialField);
      grid2 = grid1;
      Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
					      MetaData->LeftFaceBoundaryCondition,
					      MetaData->RightFaceBoundaryCondition,
					      &grid::CopyPotentialField);
    } // ENDFOR grid1
#else
    for (grid1=StartGrid; grid1<EndGrid; grid1++)
      for (grid2=0; grid2<NumberOfGrids; grid2++)
	Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
						MetaData->LeftFaceBoundaryCondition,
						MetaData->RightFaceBoundaryCondition,
						&grid::CopyPotentialField);
#endif
    CommunicationReceiveHandler();

#endif   // BITWISE_IDENTICALITY
  } // ENDFOR grid batches


  return SUCCESS;
}

#endif   // AMR_SOLVE
