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

  // find root grid corresponding to this process from the Hierarcy
  HierarchyEntry *RootGrid = &TopGrid;
  int i, ibin, isrc, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != RootGrid->GridData->ReturnProcessorNumber()) 
      RootGrid = RootGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    printf("FLD Initialize ERROR: p%"ISYM" could not locate his grid\n",
	   MyProcessorNumber);
    ENZO_FAIL("Error in AMRFLDSplit_Initialize");
  }

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

  // set default module parameters
  int autoscale   = 1;    // enable automatic variable scaling
  int SourceType[MAX_FLD_SOURCES];
  float SourceEnergy[MAX_FLD_SOURCES];
  for (isrc=0; isrc<MAX_FLD_SOURCES; isrc++) {
    SourceType[isrc] = -1;
    SourceEnergy[isrc] = 0.0;
  }

  ////////////////////////////////
  // if input file present, over-write defaults with module inputs
  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int ival, ret;
  float fval, fval2, fval3;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  // check whether input file is non-null
  if (MetaData.RadHydroParameterFname != NULL) {
    if ((fptr = fopen(MetaData.RadHydroParameterFname, "r")) == NULL)
      fprintf(stderr,"Error opening RadHydro parameter file %s, using defaults\n",
	      MetaData.RadHydroParameterFname);
    else {

      // read until out of lines
      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	ret = 0;

	// radiation fields and frequency bands
	ret += sscanf(line, "AMRFLDNumRadiationFields = %"ISYM, &NumRadiationFields);
	if (sscanf(line, "AMRFLDFrequencyBand[%"ISYM"] = %"FSYM" %"FSYM, 
		   &ibin, &fval, &fval2) == 3) {
	   ret++;  
	   if (ibin >= MAX_FLD_FIELDS) {
	     ENZO_VFAIL("AMRFLDFrequencies field %"ISYM" > maximum allowed.\n", ibin);
	   }
	   FrequencyBand[ibin][0] = fval;
	   FrequencyBand[ibin][1] = fval2;
	}

	// radiation scaling factors
	if (sscanf(line, "AMRFLDRadiationScaling[%"ISYM"] = %"FSYM, &ibin, &fval) == 2) {
	  ret++;  
	  if (ibin >= MAX_FLD_FIELDS) {
	    ENZO_VFAIL("AMRFLDRadiationScaling field %"ISYM" > maximum allowed.\n", ibin);
	  }
	  ErScale[ibin] = fval;
	}
	ret += sscanf(line, "AMRFLDAutomaticScaling = %"ISYM, &autoscale);

	// radiation sources
	ret += sscanf(line, "AMRFLDNumSources = %"ISYM, &NumSources);
	if (sscanf(line, "AMRFLDSourceLocation[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM, 
	 	   &isrc, &fval, &fval2, &fval3) == 4) {
	  ret++;  
	  if (isrc >= MAX_FLD_SOURCES) {
	    ENZO_VFAIL("AMRFLDSourceLocation source %"ISYM" > maximum allowed.\n", isrc);
          }
	  SourceLocation[isrc][0] = fval;
	  SourceLocation[isrc][1] = fval2;
	  SourceLocation[isrc][2] = fval3;
	}
	if (sscanf(line, "AMRFLDSourceType[%"ISYM"] = %"ISYM, &isrc, &ival) == 2) {
	  ret++;  
	  if (isrc >= MAX_FLD_SOURCES) {
	    ENZO_VFAIL("AMRFLDSourceType source %"ISYM" > maximum allowed.\n", isrc);
          }
	  SourceType[isrc] = ival;
	}
	if (sscanf(line, "AMRFLDSourceEnergy[%"ISYM"] = %"FSYM, &isrc, &fval) == 2) {
	  ret++;  
	  if (isrc >= MAX_FLD_SOURCES) {
	    ENZO_VFAIL("AMRFLDSourceEnergy source %"ISYM" > maximum allowed.\n", isrc);
          }
	  SourceEnergy[isrc] = fval;
	}
	if (sscanf(line, "AMRFLDSourceGroupEnergy[%"ISYM"][%"ISYM"] = %"FSYM, 
		   &isrc, &ibin, &fval) == 3) {
	  ret++;  
	  if (isrc >= MAX_FLD_SOURCES) {
	    ENZO_VFAIL("AMRFLDSourceGroupEnergy source %"ISYM" > maximum allowed.\n", isrc);
          }
	  if (ibin >= MAX_FLD_FIELDS) {
	    ENZO_VFAIL("AMRFLDSourceGroupEnergy field %"ISYM" > maximum allowed.\n", ibin);
          }
	  SourceGroupEnergy[isrc][ibin] = fval;
	}

	// limiter parameters
	ret += sscanf(line, "AMRFLDLimiterRmin = %"FSYM, &LimiterRmin);
	ret += sscanf(line, "AMRFLDLimiterDmax = %"FSYM, &LimiterDmax);

	// model parameters
	ret += sscanf(line, "AMRFLDIsothermal = %"ISYM, &Isothermal);

	// time-stepping parameters
	ret += sscanf(line, "AMRFLDMaxDt = %"FSYM, &maxdt);
	ret += sscanf(line, "AMRFLDMinDt = %"FSYM, &mindt);
	ret += sscanf(line, "AMRFLDInitDt = %"FSYM, &initdt);
	ret += sscanf(line, "AMRFLDDtControl = %"ISYM, &dt_control);
	ret += sscanf(line, "AMRFLDMaxSubcycles = %"FSYM, &maxsubcycles);
	ret += sscanf(line, "AMRFLDDtNorm = %"FSYM, &dtnorm);
	ret += sscanf(line, "AMRFLDDtGrowth = %"FSYM, &dtgrowth);
	ret += sscanf(line, "AMRFLDTimeAccuracy = %"FSYM, &timeAccuracy);
	ret += sscanf(line, "AMRFLDTheta = %"FSYM, &theta);

	// boundary condition types
	ret += sscanf(line, "AMRFLDRadiationBoundaryX0 = %i %i", 
		      BdryType[0], BdryType[0]+1);
	if (rank > 1) {
	  ret += sscanf(line, "AMRFLDRadiationBoundaryX1 = %i %i",
			BdryType[1], BdryType[1]+1);
	  if (rank > 2) {
	    ret += sscanf(line, "AMRFLDRadiationBoundaryX2 = %i %i",
			  BdryType[2], BdryType[2]+1);
	  }
	}

	// solver parameters
	ret += sscanf(line, "AMRFLDSolType = %"ISYM, &sol_type);
	ret += sscanf(line, "AMRFLDSolTolerance = %"FSYM, &sol_tolerance);
	ret += sscanf(line, "AMRFLDMaxMGIters = %"ISYM, &sol_maxit);
	ret += sscanf(line, "AMRFLDMGRelaxType = %"ISYM, &sol_rlxtype);
	ret += sscanf(line, "AMRFLDMGPreRelax = %"ISYM, &sol_npre);
	ret += sscanf(line, "AMRFLDMGPostRelax = %"ISYM, &sol_npost);

	ret += sscanf(line, "AMRFLDSolPrec = %"ISYM, &sol_prec);
	ret += sscanf(line, "AMRFLDSol_precmaxit = %"ISYM, &sol_precmaxit);
	ret += sscanf(line, "AMRFLDSol_precnpre = %"ISYM, &sol_precnpre);
	ret += sscanf(line, "AMRFLDSol_precnpost = %"ISYM, &sol_precnpost);
	ret += sscanf(line, "AMRFLDSol_precJacit = %"ISYM, &sol_precJacit);
	ret += sscanf(line, "AMRFLDSol_precrelax = %"ISYM, &sol_precrelax);
	ret += sscanf(line, "AMRFLDSol_precrestol = %"FSYM"", &sol_precrestol);

	// flag for setting up weak-scaling runs
	ret += sscanf(line, "AMRFLDWeakScaling = %"ISYM, &WeakScaling);

      }  // end loop over file lines
    }  // end successful file open
  }  // end if file name exists
 
  // clean up
  delete[] dummy;
  rewind(fptr);
  fclose(fptr);


  //// Set local grid information ////

  //   LocDims holds the dimensions of the local domain, 
  //   active cells only (no ghost or boundary cells)
  for (dim=0; dim<rank; dim++)
    LocDims[dim] = RootGrid->GridData->GetGridEndIndex(dim)
                 - RootGrid->GridData->GetGridStartIndex(dim) + 1;


  //// Check input parameters ////
  
  // check for a legal number of radiation fields
  if (NumRadiationFields < 0)
    ENZO_FAIL("AMRFLDSplit_Initialize: negative number of radiation fields requested. Halting run");
  if (NumRadiationFields > MAX_FLD_FIELDS)
    ENZO_FAIL("AMRFLDSplit_Initialize: too many radiation fields requested. Halting run");
  if (debug) 
    printf("AMRFLDSplit::Initialize NumRadiationFields = %"ISYM"\n", NumRadiationFields);

  // warn if no radiation fields are present
  if (NumRadiationFields == 0)
    fprintf(stderr,"WARNING: AMRFLDSplit solver enabled, but NO radiation fields are used!\n");
      
  // check that frequency bands start at postive frequencies
  for (ibin=0; ibin<NumRadiationFields; ibin++)
    if (FrequencyBand[ibin][0] <= 0.0) {
      ENZO_VFAIL("AMRFLDSplit_Initialize: frequency band %"ISYM" must start above 0.0. Halting run.", ibin);
    }

  // set frequency comparison tolerance
  float FreqTol = 1e-4;

  // determine whether fields are monochromatic
  for (ibin=0; ibin<NumRadiationFields; ibin++)
    if (FrequencyBand[ibin][1] - FrequencyBand[ibin][0] < FreqTol) {
      FieldMonochromatic[ibin] = true;
      FrequencyBand[ibin][1] = FrequencyBand[ibin][0]+0.1*FreqTol;  // give it a touch of width
    }

  // ensure that fields are ordered by increasing frequency
  for (ibin=0; ibin<NumRadiationFields-1; ibin++)
    if (FrequencyBand[ibin][0] > FrequencyBand[ibin+1][0]) {
      ENZO_VFAIL("AMRFLDSplit_Initialize: frequencies %"ISYM" and %"ISYM" are out of order. Halting run", 
		 ibin, ibin+1);
    }

  // ensure that groups do not overlap
  for (ibin=0; ibin<NumRadiationFields-1; ibin++)
    if ((!FieldMonochromatic[ibin]) && (FrequencyBand[ibin][1] - FrequencyBand[ibin+1][0] > FreqTol)) {
      ENZO_VFAIL("AMRFLDSplit_Initialize: frequencies %"ISYM" and %"ISYM" overlap. Halting run", 
		 ibin, ibin+1);
    }

  // if we've made it this far, the fields are legal; set neighbor flags
  for (ibin=0; ibin<NumRadiationFields-1; ibin++) {
    if (!FieldMonochromatic[ibin] && !FieldMonochromatic[ibin+1] && 
	(fabs(FrequencyBand[ibin][1] - FrequencyBand[ibin+1][0]) < FreqTol)) {
      FieldNeighbors[ibin][1] = true;
      FieldNeighbors[ibin+1][0] = true;
    }
  }

  // output field band information
  for (ibin=0; ibin<NumRadiationFields; ibin++) {
    if (debug && (FrequencyBand[ibin][1] <= FrequencyBand[ibin][0]))
      printf("   Field %"ISYM" is monochromatic at frequency %g (eV)\n", ibin, FrequencyBand[ibin][0]);
    if (debug && (FrequencyBand[ibin][0] < FrequencyBand[ibin][1]))
      printf("   Field %"ISYM" band = [%g,%g] eV\n", ibin, FrequencyBand[ibin][0], FrequencyBand[ibin][1]);
  }
  for (ibin=0; ibin<NumRadiationFields-1; ibin++)
    if ((debug) && (FieldNeighbors[ibin][1]))
      printf("   Fields %"ISYM" and %"ISYM" are adjacent\n", ibin, ibin+1);

  // check for legal radiation scaling factors
  for (ibin=0; ibin<NumRadiationFields; ibin++) 
    if (ErScale[ibin] <= 0.0) {
      fprintf(stderr,"AMRFLDSplit_Initialize: illegal AMRFLDRadiationScaling[%"ISYM"] = %g\n",
	      ibin, ErScale[ibin]);
      fprintf(stderr,"   re-setting to 1.0\n");
      ErScale[ibin] = 1.0;  // default to no scaling
    }
  for (ibin=0; ibin<NumRadiationFields; ibin++) 
    autoScale[ibin] = (autoscale != 0);  // set bool values based on integer input flag
  if (debug) {
    printf("AMRFLDSplit::Initialize scaling factors:\n");
    for (ibin=0; ibin<NumRadiationFields; ibin++) 
      printf("   ErScale[%"ISYM"] = %g\n", ibin, ErScale[ibin]);
    printf("   autoScale = %"ISYM"\n", autoscale);
  }

  // check for a legal number of sources
  if (NumSources > MAX_FLD_SOURCES) {
    ENZO_VFAIL("AMRFLDSplit_Initialize: too many radiation sources requested (%"ISYM" > %"ISYM"). Halting run\n",
	       NumSources, MAX_FLD_SOURCES);
  }
  if (debug) 
    printf("AMRFLDSplit::Initialize NumSources = %"ISYM"\n", NumSources);

  // check for legal SourceLocation
  for (isrc=0; isrc<NumSources; isrc++) {
    for (dim=0; dim<rank; dim++) 
      if ((SourceLocation[isrc][dim] < DomainLeftEdge[dim]) || 
	  (SourceLocation[isrc][dim] > DomainRightEdge[dim])) {
	ENZO_VFAIL("AMRFLDSplit_Initialize: source %"ISYM" is outside the computational domain. Halting run\n",
		   dim);
      }
    if (debug) {
      printf("AMRFLDSplit::Initialize sources %"ISYM" at location ", isrc);
      for (dim=0; dim<rank; dim++)  printf(" %g", SourceLocation[isrc][dim]);
      printf("\n");
    }
  }

  // check for legal SourceType/SourceEnergy or SourceGroupEnergy
  for (isrc=0; isrc<NumSources; isrc++) {
    bool sourceOK = false;

    // check for SourceGroupEnergy specification approach
    bool sourceBinsOK = true;
    for (ibin=0; ibin<NumRadiationFields; ibin++) 
      if (SourceGroupEnergy[isrc][ibin] < 0.0)  sourceBinsOK = false;
    if (sourceBinsOK) {
      sourceOK = true;
      if (debug) {
	printf("AMRFLDSplit::Initialize source %"ISYM" has group energies", isrc);
	for (ibin=0; ibin<NumRadiationFields; ibin++)  printf(" %g", SourceGroupEnergy[isrc][ibin]);
	printf("\n");
      }
    }

    // check for SourceType/SourceEnergy specification approach
    if ((SourceType[isrc] >= 0) && (SourceType[isrc] < NUM_FLD_SOURCE_TYPES) && (SourceEnergy[isrc] >= 0.0)) {

      // warn if the source set up using both approaches
      if (sourceOK) {
	fprintf(stderr,"AMRFLDSplit_Initialize Warning: source %"ISYM" set up using both approaches; choosing SourceGroupEnergy values\n",isrc);
	SourceType[isrc] = -1;
      } else {
	sourceOK = true;
	if (debug) 
	  printf("AMRFLDSplit::Initialize source %"ISYM" has type %"ISYM" and energy %g\n", 
		 isrc, SourceType[isrc], SourceEnergy[isrc]);
      }
    }

    // error out if the source is not set up
    if (!sourceOK) {
      ENZO_VFAIL("AMRFLDSplit_Initialize: source %"ISYM" is not properly configured. Halting run\n",isrc);
    }
  }

  // check for legal limiter parameters
  if (LimiterRmin < 0.0) {
    ENZO_VFAIL("AMRFLDSplit_Initialize: illegal LimiterRmin = %g.\n", LimiterRmin);
  }
  if (debug) 
    printf("AMRFLDSplit::Initialize LimiterRmin = %g\n", LimiterRmin);
  if (LimiterDmax <= 0.0) {
    ENZO_VFAIL("AMRFLDSplit_Initialize: illegal LimiterDmax = %g.\n", LimiterDmax);
  }
  if (debug) 
    printf("AMRFLDSplit::Initialize LimiterDmax = %g\n", LimiterDmax);

  // ensure that Enzo was called with RadiativeCooling enabled 
  // (since AMRFLDSplit doesn't handle chemistry/cooling)
  if (!RadiativeCooling) 
    ENZO_FAIL("AMRFLDSplit_Initialize: RadiativeCooling must be on!  Halting run");
  
  // maxdt gives the maximum radiation time step size
  if (maxdt <= 0.0) {
    fprintf(stderr,"AMRFLDSplit::Initialize: illegal MaxDt = %g\n",maxdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    maxdt = huge_number;  // default is no limit
  }

  // mindt gives the minimum radiation time step size
  if (mindt < 0.0) {
    fprintf(stderr,"AMRFLDSplit::Initialize: illegal MinDt = %g\n",mindt);
    fprintf(stderr,"   re-setting to %g\n",0.0);
    mindt = 0.0;  // default is 0.0
  }

  // initdt gives the initial time step size
  if (initdt <= 0.0) {
    fprintf(stderr,"AMRFLDSplit::Initialize: illegal InitDt = %g\n",initdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    initdt = huge_number;  // default is no limit
  }

  // dt_control gives the time step controller algorithm
  if ((dt_control < -1) || (dt_control > NUM_FLD_DT_CONTROLLERS-2)) {
    fprintf(stderr,"AMRFLDSplit::Initialize: illegal DtControl = %"ISYM"\n",
	    dt_control);
    fprintf(stderr,"   re-setting to -1 (original controller)\n");
    dt_control = -1;
  }
  
  // maxsubcycles gives the maximum desired ratio between hydro time step 
  // size and radiation time step size (dt_rad <= dt_hydro)
  // ***warn if subcycling radiation***
  if (maxsubcycles < 1.0) {
    if (debug) {
      fprintf(stderr,"AMRFLDSplit::Initialize: illegal MaxSubcycles = %g\n",maxsubcycles);
      fprintf(stderr,"   re-setting to 1.0\n");
    }
    maxsubcycles = 1.0;    // default is to synchronize steps
  }
  if (maxsubcycles > 1.0) {
    if (debug) {
      fprintf(stderr,"\n*********************************************************\n");
      fprintf(stderr," WARNING: radiation subcycling (MaxSubcycles = %g > 1.0)\n",
	      maxsubcycles);
      fprintf(stderr,"          may not work properly with Enzo chemistry module!\n");
      fprintf(stderr,"*********************************************************\n\n");
    }
  }

  // dtnorm gives the norm for calculating allowed relative change per step
  if (dtnorm < 0.0) {
    fprintf(stderr,"AMRFLDSplit::Initialize: illegal DtNorm = %g\n",dtnorm);
    fprintf(stderr,"   re-setting to 2.0 (2-norm)\n");
    dtnorm = 2.0;  // default is 2-norm
  }

  // dtgrowth gives the maximum growth factor in dt per step
  if (dtgrowth < 1.0 || dtgrowth > 10.0) {
    fprintf(stderr,"AMRFLDSplit::Initialize: illegal DtGrowth = %g\n",dtgrowth);
    fprintf(stderr,"   re-setting to 1.1\n");
    dtgrowth = 1.1;
  }

  // timeAccuracy gives the desired percent change in values per step
  if (timeAccuracy <= 0.0) {
    fprintf(stderr,"AMRFLDSplit::Initialize: illegal TimeAccuracy = %g\n",timeAccuracy);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    timeAccuracy = huge_number;  // default is no limit
  }

  // theta gives the implicit time-stepping method (1->BE, 0.5->CN, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"AMRFLDSplit::Initialize: illegal theta = %g\n",theta);
    fprintf(stderr,"   re-setting theta to 1.0 (Backwards Euler)\n");
    theta = 1.0;  // default is backwards Euler
  }

  // check for appropriate BdryType values, otherwise set dim to periodic
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++)
      if ((BdryType[dim][face] < 0) || (BdryType[dim][face] >= NUM_FLD_BDRY_TYPES)) {
	fprintf(stderr,"AMRFLDSplit_Initialize Warning: re-setting BC to periodic, dim %"ISYM", face %"ISYM"\n",
		dim,face);
	BdryType[dim][face] = 0;
      }

  // check that periodic faces match
  for (dim=0; dim<rank; dim++) 
    if ((BdryType[dim][0]*BdryType[dim][1] == 0) && 
	(BdryType[dim][0]+BdryType[dim][1] != 0)) {
      fprintf(stderr,"AMRFLDSplit_Initialize Warning: non-matching periodic BCs, dim %"ISYM"\n",dim);
      BdryType[dim][0] = 0;
      BdryType[dim][1] = 0;
    }

  //   check linear solver parameters
  if ((sol_type < 0) || (sol_type >= NUM_FLD_SOL_TYPES)) {
    fprintf(stderr,"Illegal SolType = %"ISYM",  Setting to 1 (BiCGStab)\n", sol_type);
    sol_type = 1;
  }
  if ((sol_tolerance < 1.0e-10) || (sol_tolerance > 1.0)) {
    fprintf(stderr,"Illegal SolTolerance = %g. Setting to 1e-5\n",
	    sol_tolerance);
    sol_tolerance = 1e-5;
  }
  if (sol_maxit <= 0) {
    fprintf(stderr,"Illegal MaxMGIters = %"ISYM". Setting to 200\n",
	    sol_maxit);
    sol_maxit = 200;
  }
  if ((sol_rlxtype < 0) || (sol_rlxtype >= NUM_HYPRE_RLX_TYPES)) {
    fprintf(stderr,"Illegal MGRelaxType = %"ISYM". Setting to 1\n",
	    sol_rlxtype);
    sol_rlxtype = 1;
  }
  if (sol_npre < 1) {
    fprintf(stderr,"Illegal MGPreRelax = %"ISYM", Setting to 1\n",
	    sol_npre);
    sol_npre = 1;
  }
  if (sol_npost < 1) {
    fprintf(stderr,"Illegal MGPostRelax = %"ISYM". Setting to 1\n",
	    sol_npost);
    sol_npost = 1;
  }
  if ((sol_prec < 0) || (sol_prec > 1)) {
    fprintf(stderr,"Illegal SolPrec = %"ISYM".  Setting to 1 (enabled)\n", 
	    sol_prec);
    sol_prec = 1;
  }
  if (sol_precmaxit <= 0) {
    fprintf(stderr,"Illegal Sol_precmaxit = %"ISYM". Setting to 1\n",
	    sol_precmaxit);
    sol_precmaxit = 1;
  }
  if (sol_precnpre < 1) {
    fprintf(stderr,"Illegal Sol_precnpre = %"ISYM". Setting to 1\n",
	    sol_precnpre);
    sol_precnpre = 1;
  }
  if (sol_precnpost < 1) {
    fprintf(stderr,"Illegal Sol_precnpost = %"ISYM". Setting to 1\n",
	    sol_precnpost);
    sol_precnpost = 1;
  }
  if (sol_precJacit < 1) {
    fprintf(stderr,"Illegal Sol_precJacit = %"ISYM". Setting to 2\n",
	    sol_precJacit);
    sol_precJacit = 2;
  }
  if ((sol_precrelax < 0) || (sol_precrelax >= NUM_HYPRE_RLX_TYPES)) {
    fprintf(stderr,"Illegal Sol_precrelax = %"ISYM". Setting to 1\n",
	    sol_precrelax);
    sol_precrelax = 1;
  }
  if ((sol_precrestol < 0.0) || (sol_precrestol > 1.0)) {
    fprintf(stderr,"Illegal Sol_precrestol = %g. Setting to 0.0\n",
	    sol_precrestol);
    sol_precrestol = 0.0;
  }

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

  // dt* gives the time step sizes for each piece of physics
  dtrad = initdt;            // use the input value (scaled units)

  // set a bound on the global initial dt as a factor of the radiation timestep
  dt = initdt*maxsubcycles;

  // set initial time step into TopGrid
  RootGrid->GridData->SetMaxRadiationDt(dt);
  
  // store number of ghost zones in each direction
  int xghosts = DEFAULT_GHOST_ZONES, yghosts=0, zghosts=0;
  if (rank > 1) {
    yghosts = DEFAULT_GHOST_ZONES;
    if (rank > 2) {
      zghosts = DEFAULT_GHOST_ZONES;
    }
  }

  // if this is a weak scaling test, overwrite SourceLocation for each 
  // source to replicate source setup on each root grid tile
  if (WeakScaling)
    for (isrc=0; isrc<NumSources; isrc++) 
      for (dim=0; dim<rank; dim++) {
	float frac = (SourceLocation[isrc][dim] - DomainLeftEdge[dim]) 
	           / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	OriginalSourceLocation[isrc][dim] = SourceLocation[isrc][dim];
	SourceLocation[isrc][dim] = RootGrid->GridData->GetGridLeftEdge(dim) +
                                    frac*(RootGrid->GridData->GetGridRightEdge(dim) -
					  RootGrid->GridData->GetGridLeftEdge(dim));
      }
  
  // compute Radiation Energy spectrum integrals to fill 
  // intOpacity*, intIonizing* and intHeating* arrays
  if (this->ComputeRadiationIntegrals() == FAIL) 
    ENZO_FAIL("AMRFLDSplit::Initialize Error in radiation spectrum integrals");


  // fill SourceGroupEnergy for sources that are specified by type/energy
  for (isrc=0; isrc<NumSources; isrc++) {
    if (SourceType[isrc] == 0) {   // monochromatic source at hnu = 13.6 eV
      MonochromaticSED tmp_src(13.6);
      for (ibin=0; ibin<NumRadiationFields; ibin++) {
	if (FieldMonochromatic[ibin]) {   // skip monochromatic fields
	  SourceGroupEnergy[isrc][ibin] = 0.0;
	  continue;
	}
	SourceGroupEnergy[isrc][ibin] = SourceEnergy[isrc]
	  * SED_integral(tmp_src, FrequencyBand[ibin][0], FrequencyBand[ibin][1], true);
      }
      if (debug) {
	printf("AMRFLDSplit::Initialize monochromatic source %"ISYM" has group energies", isrc);
	for (ibin=0; ibin<NumRadiationFields; ibin++)  printf(" %g", SourceGroupEnergy[isrc][ibin]);
	printf("\n");
      }
    }
    if (SourceType[isrc] == 1) {  // blackbody spectrum at T = 1e5 K
      BlackbodySED tmp_src(1.0e5);
      float total_integral = SED_integral(tmp_src, 13.6, -1.0, true);
      for (ibin=0; ibin<NumRadiationFields; ibin++) {
	if (FieldMonochromatic[ibin]) {   // skip monochromatic fields
	  SourceGroupEnergy[isrc][ibin] = 0.0;
	  continue;
	}
	SourceGroupEnergy[isrc][ibin] = SourceEnergy[isrc] / total_integral
	  * SED_integral(tmp_src, FrequencyBand[ibin][0], FrequencyBand[ibin][1], true);
      }
      if (debug) {
	printf("AMRFLDSplit::Initialize blackbody source %"ISYM" has group energies", isrc);
	for (ibin=0; ibin<NumRadiationFields; ibin++)  printf(" %g", SourceGroupEnergy[isrc][ibin]);
	printf("\n");
      }
    }
  }


#ifdef USE_HYPRE

#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif
  // initialize amrsolve stuff
  //    initialize the diagnostic information
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

  ////////////////////////////////
  // set up the boundary conditions on the radiation field, 
  // depending on the ProblemType
  float ZERO = 0.0;
  float ONE  = 1.0;
  float SMALL = 1.0e-6;
  fptr = NULL;

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
