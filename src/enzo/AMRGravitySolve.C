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

  static float AMRGravSolTime = 0.0;
  static float AMRGravSolTime2 = 0.0;
  // start MPI timer for overall solver
#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif


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
  
  // set solver parameters based on inputs/defaults
  char numstr[80];
  switch (AMRGravitySolve_solver) {
  case 1:  // GMRES
    amrsolve_params->set_parameter("solver","gmres");
    break;
  case 2:  // FAC
    amrsolve_params->set_parameter("solver","fac");
    break;
  default: // BiCGStab
    amrsolve_params->set_parameter("solver","bicgstab");
    break;
  }
  sprintf(numstr, "%i", AMRGravitySolve_maxit);
  amrsolve_params->set_parameter("solver_itmax",numstr);
  amrsolve_params->set_parameter("solver_printl", "0");

  sprintf(numstr, "%e", AMRGravitySolve_restol);
  amrsolve_params->set_parameter("solver_restol",numstr);

  // set preconditioning options for BiCGStab and GMRES solvers
  sprintf(numstr, "%i", AMRGravitySolve_precmaxit);
  amrsolve_params->set_parameter("prec_itmax",numstr);
  sprintf(numstr, "%e", AMRGravitySolve_precrestol);
  amrsolve_params->set_parameter("prec_restol",numstr);
  sprintf(numstr, "%i", AMRGravitySolve_rlxtype);
  amrsolve_params->set_parameter("prec_rlxtype",numstr);
  sprintf(numstr, "%i", AMRGravitySolve_npre);
  amrsolve_params->set_parameter("prec_npre",numstr);
  amrsolve_params->set_parameter("prec_npost",numstr);
  sprintf(numstr, "%i", AMRGravitySolve_Jaciters);
  amrsolve_params->set_parameter("prec_Jaciters",numstr);
  amrsolve_params->set_parameter("prec_printl", "0");
  amrsolve_params->set_parameter("prec_log",    "1");

//   amrsolve_params->set_parameter("dump_x", "true");
 
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
  bool is_periodic[] = {true, true, true};

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
  LCAPERF_START("amrsolve_matrix");
  AMRsolve_Hypre_Grav amrgravsolve(*hierarchy, *amrsolve_params, 
				   AMRGravitySolve_useprec, 
				   AMRGravitySolve_zeroguess);
  amrgravsolve.init_hierarchy();
  amrgravsolve.init_stencil();
  amrgravsolve.init_graph();
  amrgravsolve.init_elements(f_scale);
  LCAPERF_STOP("amrsolve_matrix");


  // start MPI timer for the solve itself
#ifdef USE_MPI
  float stime2 = MPI_Wtime();
#else
  float stime2 = 0.0;
#endif

  // Solve the linear system
  LCAPERF_START("amrsolve_solve");
  amrgravsolve.solve();
  Eflt64 finalresid = amrgravsolve.residual();
  Eint32 Sits = amrgravsolve.iterations();
  if (debug) printf("   lin resid = %.1e, its = %i\n", finalresid, Sits);
  LCAPERF_STOP("amrsolve_solve");

  // stop MPI timer for the solve itself, increment total
#ifdef USE_MPI
  float ftime2 = MPI_Wtime();
#else
  float ftime2 = 0.0;
#endif
  AMRGravSolTime2 += ftime2-stime2;
    

  // Display solver results
  if (amrgravsolve.evaluate() != 0) {
    fprintf(stderr,"AMRGravitySolve: could not achieve prescribed tolerance!\n");
    fprintf(stderr, "Printing current hierarchy:\n");
    if (debug)   hierarchy->print();
    return FAIL;
  }

  // increment Enzo potential field with amrsolve solution
  amrgravsolve.update_enzo();

  LCAPERF_STOP("amr_solve");

  // share boundary values
  int grid1, grid2, StartGrid, EndGrid;
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


//   // write Enzo potential fields (with ghost zones) to disk
//   amrgravsolve.write_potential();

  // copy potential values to BaryonField if requested
  if (CopyGravPotential)
  for (int thislevel=level_coarse; thislevel<=level_fine; thislevel++)
    for (Temp=LevelArray[thislevel]; Temp; Temp=Temp->NextGridThisLevel) 
      Temp->GridHierarchyEntry->GridData->CopyPotentialToBaryonField();
 
  // Clean up
  hierarchy->enzo_detach();
  delete hierarchy;
  delete pmpi;
  delete amrsolve_params;
  hierarchy = NULL;
  amrsolve_params = NULL;

  // stop MPI timer for overall solver itself, increment total
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  flaot ftime = 0.0;
#endif
  AMRGravSolTime += ftime-stime;

  // root node outputs timing information
  if (debug)  printf("Cumulative AMRGravitySolve time = %g (infrastructure = %g)\n\n",
		     AMRGravSolTime2, AMRGravSolTime-AMRGravSolTime2);
    

  return SUCCESS;
}

#endif   // AMR_SOLVE
