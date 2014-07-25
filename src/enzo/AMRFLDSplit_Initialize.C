/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, Initialization routine
/
/  written by: Daniel Reynolds
/  date:       July 2014
/
/  PURPOSE: Allocates all necessary internal memory for problem 
/           definition and associated linear solver.  This begins the
/           interface between Enzo and the FLD solver module, so any
/           and all grid/index transformations must be performed and 
/           stored here.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
#include "BlackbodySED.h"
#include "MonochromaticSED.h"
// #ifdef _OPENMP
// #include <omp.h>
// #endif

// character strings
EXTERN char outfilename[];


// function prototypes
int InitializeRateData(FLOAT Time);
int FreezeRateData(FLOAT Time, HierarchyEntry &TopGrid);
float SED_integral(SED &sed, float a, float b, bool convertHz);
int CosmoIonizationInitialize(FILE *fptr, FILE *Outfptr,
			      HierarchyEntry &TopGrid,
			      TopGridData &MetaData, int local);
int RadHydroConstTestInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);
int RHIonizationTestInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local);
int RHIonizationSteepInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);
int RHIonizationClumpInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);



int AMRFLDSplit::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData) {

#ifdef AMR_SOLVE

  if (debug)  printf("Entering AMRFLDSplit::Initialize routine\n");

  //// Set parallelism information ////

// #ifdef _OPENMP
//   // output number of OpenMP threads that will be used in this run
//   int nthreads = omp_get_max_threads();
//   printf("FLD Initialize: MPI task %"ISYM" has %"ISYM" available OpenMP threads\n",
// 	 MyProcessorNumber,nthreads);
// #endif

#ifndef MPI_INT
  // in case MPI is not included
  int MPI_PROC_NULL = -3;
  int MPI_COMM_WORLD = 0;
#endif

  // initialize the amrsolve Mpi object
  pmpi = new AMRsolve_Mpi(MPI_COMM_WORLD);
  AMRsolve_Grid::set_mpi(*pmpi);

  // set rank of fld problem; error message if not 3 (amrsolve requirement)
  rank = MetaData.TopGridRank;
  if (rank != 3)
    ENZO_FAIL("Error in AMRFLDSplit_Initialize: rank must be 3 (for now)");


  //// input solver parameters ////
  
  // initialize source specification utilities
  int SourceType[MAX_FLD_SOURCES];
  float SourceEnergy[MAX_FLD_SOURCES];
  for (int isrc=0; isrc<MAX_FLD_SOURCES; isrc++) {
    SourceType[isrc] = -1;
    SourceEnergy[isrc] = 0.0;
  }

  // read and check input parameters
  if (this->ReadParameters(MetaData, SourceType, SourceEnergy) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_Initialize: Error in ReadParameters.");


  //// Set local grid information ////

  // find root grid corresponding to this process from the Hierarcy
  HierarchyEntry *RootGrid = &TopGrid;
  bool foundgrid = false;
  for (int i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != RootGrid->GridData->ReturnProcessorNumber()) 
      RootGrid = RootGrid->NextGridThisLevel;
    else {foundgrid=true; break;}
  }
  if (!foundgrid) {
    printf("FLD Initialize ERROR: p%"ISYM" could not locate his grid\n",
	   MyProcessorNumber);
    ENZO_FAIL("Error in AMRFLDSplit_Initialize");
  }

  //   LocDims holds the dimensions of the local domain, 
  //   active cells only (no ghost or boundary cells)
  int dim;
  for (dim=0; dim<rank; dim++)
    LocDims[dim] = RootGrid->GridData->GetGridEndIndex(dim)
                 - RootGrid->GridData->GetGridStartIndex(dim) + 1;

  // initialize internal module units
  if (this->UpdateUnits(MetaData.Time, MetaData.Time) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_Initialize: Error in UpdateUnits.");
  
  // get processor layout from Grid
  int layout[3];     // number of procs in each dim (1-based)
  for (dim=0; dim<rank; dim++) 
    layout[dim] = RootGrid->GridData->GetProcessorLayout(dim);
  
  // get processor location in MPI grid
  int location[3];   // location of this proc in each dim (0-based)
  for (dim=0; dim<rank; dim++) 
    location[dim] = RootGrid->GridData->GetProcessorLocation(dim);

  // set flags denoting if this processor is on the external boundary
  for (dim=0; dim<rank; dim++) {
    if (layout[dim]==0) {
      OnBdry[dim][0] = OnBdry[dim][1] = true;
    }
    else {
      OnBdry[dim][0] = (location[dim] == 0);
      OnBdry[dim][1] = (location[dim] == layout[dim]-1);
    }
  }
  if (debug){
    printf("AMRFLDSplit::Initialize: rank = %"ISYM"\n", rank);
    printf("AMRFLDSplit::Initialize: layout = (%"ISYM",%"ISYM",%"ISYM")\n",
	   layout[0],layout[1],layout[2]);
    printf("AMRFLDSplit::Initialize: BdryType = (%i:%i,%i:%i,%i:%i)\n",
	   BdryType[0][0], BdryType[0][1], BdryType[1][0], 
	   BdryType[1][1], BdryType[2][0], BdryType[2][1]);
  }

  // store number of ghost zones in each direction
  int xghosts = DEFAULT_GHOST_ZONES, yghosts=0, zghosts=0;
  if (rank > 1) {
    yghosts = DEFAULT_GHOST_ZONES;
    if (rank > 2) {
      zghosts = DEFAULT_GHOST_ZONES;
    }
  }


  //// General setup ////

  // dt* gives the time step sizes for each piece of physics
  dtrad = initdt;            // use the input value (scaled units)

  // set a bound on the global initial dt as a factor of the radiation timestep
  dt = initdt*maxsubcycles;

  // set initial time step into TopGrid
  RootGrid->GridData->SetMaxRadiationDt(dt);


  // compute Radiation Energy spectrum integrals to fill 
  // intOpacity*, intIonizing* and intHeating* arrays
  if (this->ComputeRadiationIntegrals() == FAIL) 
    ENZO_FAIL("AMRFLDSplit::Initialize Error in radiation spectrum integrals");


  // fill SourceGroupEnergy for each source and radiation field.
  if (this->SetupSources(RootGrid, SourceType, SourceEnergy) == FAIL) 
    ENZO_FAIL("AMRFLDSplit::Initialize Error in SetupSources");



  //// Setup AMRsolve objects for linear solves ////

#ifdef USE_HYPRE

#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif
  // initialize amrsolve stuff
  //    initialize the diagnostic information
  int ibin;
  for (ibin=0; ibin<NumRadiationFields; ibin++)
    totIters[ibin] = 0;

  //    set amrsolve parameters
  amrsolve_params = new AMRsolve_Parameters();
  amrsolve_params->set_defaults();
  if (sol_type == 0)  amrsolve_params->set_parameter("solver","fac");
  if (sol_type == 1)  amrsolve_params->set_parameter("solver","bicgstab");
  if (sol_type == 2)  amrsolve_params->set_parameter("solver","bicgstab-boomer");
  if (sol_type == 3)  amrsolve_params->set_parameter("solver","gmres");
  if (sol_type == 4)  amrsolve_params->set_parameter("solver","pfmg");
  char numstr[80];
  sprintf(numstr, "%e", sol_tolerance);
  amrsolve_params->set_parameter("solver_restol",numstr);
  sprintf(numstr, "%"ISYM, sol_maxit);
  amrsolve_params->set_parameter("solver_itmax",numstr);
  sprintf(numstr, "%"ISYM, sol_printl);
  amrsolve_params->set_parameter("solver_printl",numstr);
  sprintf(numstr, "%"ISYM, sol_log);
  amrsolve_params->set_parameter("solver_log",numstr);
  sprintf(numstr, "%"ISYM, sol_rlxtype);
  amrsolve_params->set_parameter("solver_rlxtype",numstr);
  sprintf(numstr, "%"ISYM, sol_npre);
  amrsolve_params->set_parameter("solver_npre",numstr);
  sprintf(numstr, "%"ISYM, sol_npost);
  amrsolve_params->set_parameter("solver_npost",numstr);
  sprintf(numstr, "%e", LimiterRmin);
  amrsolve_params->set_parameter("limiter_rmin",numstr);
  sprintf(numstr, "%e", LimiterDmax);
  amrsolve_params->set_parameter("limiter_dmax",numstr);

  // set preconditioning options for BiCGStab and GMRES solvers
  if (sol_type == 1 || sol_type==3) {
    sprintf(numstr, "%"ISYM, sol_precmaxit);
    amrsolve_params->set_parameter("prec_itmax", numstr);
    sprintf(numstr, "%e", sol_precrestol);
    amrsolve_params->set_parameter("prec_restol", numstr);
    sprintf(numstr, "%"ISYM, sol_precnpre);
    amrsolve_params->set_parameter("prec_npre", numstr);
    sprintf(numstr, "%"ISYM, sol_precnpost);
    amrsolve_params->set_parameter("prec_npost", numstr);
    sprintf(numstr, "%"ISYM, sol_precJacit);
    amrsolve_params->set_parameter("prec_Jaciters", numstr);
    sprintf(numstr, "%"ISYM, sol_precrelax);
    amrsolve_params->set_parameter("prec_rlxtype", numstr);
  }
 
  if (debug) {
    printf("AMRFLDSplit::Initialize, customized amrsolve parameters:\n");
    amrsolve_params->print();
  }


#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  AMRSolTime += ftime-stime;

#else  // ifdef USE_HYPRE

  ENZO_FAIL("AMRFLDSplit_Initialize ERROR: module requires USE_HYPRE to be set!");
  
#endif



  //// Setup problem-specific boundary conditions ////

  float ZERO = 0.0;
  float ONE  = 1.0;
  float SMALL = 1.0e-6;
  FILE *fptr = NULL;
  int face;

  // set boundary conditions based on problem type
  // (default to homogeneous Dirichlet)
  switch (ProblemType) {
    
  // ODE test problem, set BCs based on input.  
  // 0 implies periodic, otherwise set to homogeneous Dirichlet
  case 400:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroConstTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroConstTestInitialize.");

    for (ibin=0; ibin<NumRadiationFields; ibin++) {
      for (dim=0; dim<rank; dim++)
	for (face=0; face<2; face++) {
	  if (BdryType[dim][face] != 0) {
	    if (this->SetupBoundary(ibin,dim,face,1,&ZERO) == FAIL) {
	      ENZO_VFAIL("Error setting dim%"ISYM" face%"ISYM" radiation BCs.", dim, face);
	    }
	  }
	}
    }
    break;
    
  // Ionization tests 0 and 1: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 410:
  case 411:
    // first call local problem initializer (to allocate/setup local data)
    if (RHIonizationTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RHIonizationTestInitialize.");
    
    for (ibin=0; ibin<NumRadiationFields; ibin++) {
      for (dim=0; dim<rank; dim++)
	for (face=0; face<2; face++) {
	  if (this->SetupBoundary(ibin,dim,face,1,&ZERO) == FAIL) {
	      ENZO_VFAIL("Error setting dim%"ISYM" face%"ISYM" radiation BCs.", dim, face);
	  }
	}
    }
    break;
    
    
  // Ionization test 2: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all non-periodic faces.
  case 412:
    // first call local problem initializer (to allocate/setup local data)
    if (RHIonizationClumpInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RHIonizationSteepInitialize.");
    
    for (ibin=0; ibin<NumRadiationFields; ibin++) {
      for (dim=0; dim<rank; dim++)
	for (face=0; face<2; face++) {
	  if (BdryType[dim][face] != 0)
	    if (this->SetupBoundary(ibin,dim,face,1,&ZERO) == FAIL) {
	      ENZO_VFAIL("Error setting dim%"ISYM" face%"ISYM" radiation BCs.", dim, face);
	    }
	}
    }
    break;
    
    
  // Ionization test 13: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 413:
    // first call local problem initializer (to allocate/setup local data)
    if (RHIonizationSteepInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RHIonizationSteepInitialize.");
    
    for (ibin=0; ibin<NumRadiationFields; ibin++) {
      for (dim=0; dim<rank; dim++)
	for (face=0; face<2; face++) {
	  if (this->SetupBoundary(ibin,dim,face,1,&ZERO) == FAIL) {
	      ENZO_VFAIL("Error setting dim%"ISYM" face%"ISYM" radiation BCs.", dim, face);
	  }
	}
    }
    break;
    
    
  // Ionization test 14: periodic boundary conditions on all faces (store no data).
  case 414:
    // first call local problem initializer (to allocate/setup local data)
    if (CosmoIonizationInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in CosmoIonizationInitialize.");
    
    break;
    

  // Ionization test 15: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 415:
    // first call local problem initializer (to allocate/setup local data)
    if (CosmoIonizationInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in CosmoIonizationInitialize.");
    
    for (ibin=0; ibin<NumRadiationFields; ibin++) {
      for (dim=0; dim<rank; dim++)
	for (face=0; face<2; face++) {
	  if (this->SetupBoundary(ibin,dim,face,1,&ZERO) == FAIL) {
	      ENZO_VFAIL("Error setting dim%"ISYM" face%"ISYM" radiation BCs.", dim, face);
	  }
	}
    }
    break;
    
    
  // Insert new problem intializers here...



  default:

    // set BCs based on inputs, for non periodic set to 0-valued
    for (ibin=0; ibin<NumRadiationFields; ibin++) {
      for (dim=0; dim<rank; dim++)
	for (face=0; face<2; face++) {
	  if (BdryType[dim][face] != 0)
	    if (this->SetupBoundary(ibin,dim,face,1,&ZERO) == FAIL) {
	      ENZO_VFAIL("Error setting dim%"ISYM" face%"ISYM" radiation BCs.", dim, face);
	    }
	}
    }
    break;
  }
  ////////////////////////////////

  // ensure that CoolData object has been set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) 
      ENZO_FAIL("Error in InitializeRateData.");
  
  // if using an isothermal test, freeze rate data, now that ICs exist
  if (Isothermal) 
    if (FreezeRateData(MetaData.Time, TopGrid) == FAIL) 
      ENZO_FAIL("Error in FreezeRateData.");


  // set 'diags' flag if any processor has 'debug' set
  int glob_debug = 0;
#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  MPI_Allreduce(&debug,&glob_debug,1,DataType,MPI_SUM,MPI_COMM_WORLD);
#else
  glob_debug = debug;
#endif
  diags = (glob_debug != 0);


  //  if (debug) printf("  AMRFLDSplit_Initialize: outputting parameters to log file\n");

  // output solver parameters to output log file 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *outfptr;
    if ((outfptr = fopen(outfilename, "a")) == NULL) {
      fprintf(stderr,"Error opening parameter output file %s!!\n", 
	      outfilename);
      ENZO_FAIL("Error in AMRFLDSplit_Initialize");
    }
    else {
      // write parameters to log file and close
      this->WriteParameters(outfptr);
      fclose(outfptr);
    }
  }

  return SUCCESS;

#else // AMR_SOLVE
  ENZO_FAIL("AMRFLDSplit_Initialize requires AMR_SOLVE in configuration");
#endif // AMR_SOLVE
}
#endif // TRANSFER
