/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, ReadParameters routine
/
/  written by: Daniel Reynolds
/  date:       July 2014
/
/  PURPOSE: Reads the AMRFLD-specific solver parameters from the
/           input file.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"



int AMRFLDSplit::ReadParameters(TopGridData &MetaData, 
				int SourceType[MAX_FLD_SOURCES], 
				float SourceEnergy[MAX_FLD_SOURCES]) {

  if (debug)  printf("Entering AMRFLDSplit::ReadParameters routine\n");

  // set default module parameters
  int autoscale = 1;    // enable automatic variable scaling

  // if input file present, over-write defaults with module inputs
  int ibin, isrc, dim, face;
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

	// analytical opacity parameters
	ret += sscanf(line, "AMRFLDAnalyticOpacity = %"ISYM, &AnalyticOpacity);
	if (sscanf(line, "AMRFLDOpacityC0[%"ISYM"] = %"FSYM, &ibin, &fval) == 2) {
	  ret++;  
	  if (ibin >= MAX_FLD_FIELDS) {
	    ENZO_VFAIL("AMRFLDOpacityC0 field %"ISYM" > maximum allowed.\n", ibin);
	  }
	  OpacityC0[ibin] = fval;
	}
	if (sscanf(line, "AMRFLDOpacityC1[%"ISYM"] = %"FSYM, &ibin, &fval) == 2) {
	  ret++;  
	  if (ibin >= MAX_FLD_FIELDS) {
	    ENZO_VFAIL("AMRFLDOpacityC1 field %"ISYM" > maximum allowed.\n", ibin);
	  }
	  OpacityC1[ibin] = fval;
	}
	if (sscanf(line, "AMRFLDOpacityC2[%"ISYM"] = %"FSYM, &ibin, &fval) == 2) {
	  ret++;  
	  if (ibin >= MAX_FLD_FIELDS) {
	    ENZO_VFAIL("AMRFLDOpacityC2 field %"ISYM" > maximum allowed.\n", ibin);
	  }
	  OpacityC2[ibin] = fval;
	}
	if (sscanf(line, "AMRFLDOpacityC3[%"ISYM"] = %"FSYM, &ibin, &fval) == 2) {
	  ret++;  
	  if (ibin >= MAX_FLD_FIELDS) {
	    ENZO_VFAIL("AMRFLDOpacityC3 field %"ISYM" > maximum allowed.\n", ibin);
	  }
	  OpacityC3[ibin] = fval;
	}
	if (sscanf(line, "AMRFLDOpacityC4[%"ISYM"] = %"FSYM, &ibin, &fval) == 2) {
	  ret++;  
	  if (ibin >= MAX_FLD_FIELDS) {
	    ENZO_VFAIL("AMRFLDOpacityC4 field %"ISYM" > maximum allowed.\n", ibin);
	  }
	  OpacityC4[ibin] = fval;
	}

	// flag for setting up weak-scaling runs
	ret += sscanf(line, "AMRFLDWeakScaling = %"ISYM, &WeakScaling);

      }  // end loop over file lines
    }  // end successful file open
  }  // end if file name exists
 
  // clean up
  delete[] dummy;
  rewind(fptr);
  fclose(fptr);


  //// Check input parameters ////
  
  // check for a legal number of radiation fields
  if (NumRadiationFields < 0)
    ENZO_FAIL("AMRFLDSplit_ReadParameters: negative number of radiation fields requested. Halting run");
  if (NumRadiationFields > MAX_FLD_FIELDS)
    ENZO_FAIL("AMRFLDSplit_ReadParameters: too many radiation fields requested. Halting run");
  if (debug) 
    printf("AMRFLDSplit::ReadParameters NumRadiationFields = %"ISYM"\n", NumRadiationFields);

  // warn if no radiation fields are present
  if (NumRadiationFields == 0)
    fprintf(stderr,"WARNING: AMRFLDSplit solver enabled, but NO radiation fields are used!\n");
      
  // check that frequency bands start at postive frequencies
  for (ibin=0; ibin<NumRadiationFields; ibin++)
    if (FrequencyBand[ibin][0] <= 0.0) {
      ENZO_VFAIL("AMRFLDSplit_ReadParameters: frequency band %"ISYM" must start above 0.0. Halting run.", ibin);
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
      ENZO_VFAIL("AMRFLDSplit_ReadParameters: frequencies %"ISYM" and %"ISYM" are out of order. Halting run", 
		 ibin, ibin+1);
    }

  // ensure that groups do not overlap
  for (ibin=0; ibin<NumRadiationFields-1; ibin++)
    if ((!FieldMonochromatic[ibin]) && (FrequencyBand[ibin][1] - FrequencyBand[ibin+1][0] > FreqTol)) {
      ENZO_VFAIL("AMRFLDSplit_ReadParameters: frequencies %"ISYM" and %"ISYM" overlap. Halting run", 
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
      fprintf(stderr,"AMRFLDSplit_ReadParameters: illegal AMRFLDRadiationScaling[%"ISYM"] = %g\n",
	      ibin, ErScale[ibin]);
      fprintf(stderr,"   re-setting to 1.0\n");
      ErScale[ibin] = 1.0;  // default to no scaling
    }
  for (ibin=0; ibin<NumRadiationFields; ibin++) 
    autoScale[ibin] = (autoscale != 0);  // set bool values based on integer input flag
  if (debug) {
    printf("AMRFLDSplit::ReadParameters scaling factors:\n");
    for (ibin=0; ibin<NumRadiationFields; ibin++) 
      printf("   ErScale[%"ISYM"] = %g\n", ibin, ErScale[ibin]);
    printf("   autoScale = %"ISYM"\n", autoscale);
  }

  // check for a legal number of sources
  if (NumSources > MAX_FLD_SOURCES) {
    ENZO_VFAIL("AMRFLDSplit_ReadParameters: too many radiation sources requested (%"ISYM" > %"ISYM"). Halting run\n",
	       NumSources, MAX_FLD_SOURCES);
  }
  if (debug) 
    printf("AMRFLDSplit::ReadParameters NumSources = %"ISYM"\n", NumSources);

  // check for legal SourceLocation
  for (isrc=0; isrc<NumSources; isrc++) {
    for (dim=0; dim<rank; dim++) 
      if ((SourceLocation[isrc][dim] < DomainLeftEdge[dim]) || 
	  (SourceLocation[isrc][dim] > DomainRightEdge[dim])) {
	ENZO_VFAIL("AMRFLDSplit_ReadParameters: source %"ISYM" is outside the computational domain. Halting run\n", isrc);
      }
    if (debug) {
      printf("AMRFLDSplit::ReadParameters source %"ISYM" at location ", isrc);
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
	printf("AMRFLDSplit::ReadParameters source %"ISYM" has group energies", isrc);
	for (ibin=0; ibin<NumRadiationFields; ibin++)  printf(" %g", SourceGroupEnergy[isrc][ibin]);
	printf("\n");
      }
    }

    // check for SourceType/SourceEnergy specification approach
    if ((SourceType[isrc] >= 0) && (SourceType[isrc] < NUM_FLD_SOURCE_TYPES) && (SourceEnergy[isrc] >= 0.0)) {

      // warn if the source set up using both approaches
      if (sourceOK) {
	fprintf(stderr,"AMRFLDSplit_ReadParameters Warning: source %"ISYM" set up using both approaches; choosing SourceGroupEnergy values\n",isrc);
	SourceType[isrc] = -1;
      } else {
	sourceOK = true;
	if (debug) 
	  printf("AMRFLDSplit::ReadParameters source %"ISYM": type %"ISYM", %g photons/sec\n", 
		 isrc, SourceType[isrc], SourceEnergy[isrc]);
      }
    }

    // error out if the source is not set up
    if (!sourceOK) {
      ENZO_VFAIL("AMRFLDSplit_ReadParameters: source %"ISYM" is not properly configured. Halting run\n",isrc);
    }
  }

  // check for legal limiter parameters
  if (LimiterRmin < 0.0) {
    ENZO_VFAIL("AMRFLDSplit_ReadParameters: illegal LimiterRmin = %g.\n", LimiterRmin);
  }
  if (debug) 
    printf("AMRFLDSplit::ReadParameters LimiterRmin = %g\n", LimiterRmin);
  if (LimiterDmax <= 0.0) {
    ENZO_VFAIL("AMRFLDSplit_ReadParameters: illegal LimiterDmax = %g.\n", LimiterDmax);
  }
  if (debug) 
    printf("AMRFLDSplit::ReadParameters LimiterDmax = %g\n", LimiterDmax);

  // ensure that Enzo was called with RadiativeCooling enabled 
  // (since AMRFLDSplit doesn't handle chemistry/cooling)
  if (!RadiativeCooling) 
    ENZO_FAIL("AMRFLDSplit_ReadParameters: RadiativeCooling must be on!  Halting run");
  
  // maxdt gives the maximum radiation time step size
  if (maxdt <= 0.0) {
    fprintf(stderr,"AMRFLDSplit::ReadParameters: illegal MaxDt = %g\n",maxdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    maxdt = huge_number;  // default is no limit
  }

  // mindt gives the minimum radiation time step size
  if (mindt < 0.0) {
    fprintf(stderr,"AMRFLDSplit::ReadParameters: illegal MinDt = %g\n",mindt);
    fprintf(stderr,"   re-setting to %g\n",0.0);
    mindt = 0.0;  // default is 0.0
  }

  // initdt gives the initial time step size
  if (initdt <= 0.0) {
    fprintf(stderr,"AMRFLDSplit::ReadParameters: illegal InitDt = %g\n",initdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    initdt = huge_number;  // default is no limit
  }

  // dt_control gives the time step controller algorithm
  if ((dt_control < -1) || (dt_control > NUM_FLD_DT_CONTROLLERS-2)) {
    fprintf(stderr,"AMRFLDSplit::ReadParameters: illegal DtControl = %"ISYM"\n",
	    dt_control);
    fprintf(stderr,"   re-setting to -1 (original controller)\n");
    dt_control = -1;
  }
  
  // maxsubcycles gives the maximum desired ratio between hydro time step 
  // size and radiation time step size (dt_rad <= dt_hydro)
  // ***warn if subcycling radiation***
  if (maxsubcycles < 1.0) {
    if (debug) {
      fprintf(stderr,"AMRFLDSplit::ReadParameters: illegal MaxSubcycles = %g\n",maxsubcycles);
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
    fprintf(stderr,"AMRFLDSplit::ReadParameters: illegal DtNorm = %g\n",dtnorm);
    fprintf(stderr,"   re-setting to 2.0 (2-norm)\n");
    dtnorm = 2.0;  // default is 2-norm
  }

  // dtgrowth gives the maximum growth factor in dt per step
  if (dtgrowth < 1.0 || dtgrowth > 10.0) {
    fprintf(stderr,"AMRFLDSplit::ReadParameters: illegal DtGrowth = %g\n",dtgrowth);
    fprintf(stderr,"   re-setting to 1.1\n");
    dtgrowth = 1.1;
  }

  // timeAccuracy gives the desired percent change in values per step
  if (timeAccuracy <= 0.0) {
    fprintf(stderr,"AMRFLDSplit::ReadParameters: illegal TimeAccuracy = %g\n",timeAccuracy);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    timeAccuracy = huge_number;  // default is no limit
  }

  // theta gives the implicit time-stepping method (1->BE, 0.5->CN, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"AMRFLDSplit::ReadParameters: illegal theta = %g\n",theta);
    fprintf(stderr,"   re-setting theta to 1.0 (Backwards Euler)\n");
    theta = 1.0;  // default is backwards Euler
  }

  // check for appropriate BdryType values, otherwise set dim to periodic
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++)
      if ((BdryType[dim][face] < 0) || (BdryType[dim][face] >= NUM_FLD_BDRY_TYPES)) {
	fprintf(stderr,"AMRFLDSplit_ReadParameters Warning: re-setting BC to periodic, dim %"ISYM", face %"ISYM"\n",
		dim,face);
	BdryType[dim][face] = 0;
      }

  // check that periodic faces match
  for (dim=0; dim<rank; dim++) 
    if ((BdryType[dim][0]*BdryType[dim][1] == 0) && 
	(BdryType[dim][0]+BdryType[dim][1] != 0)) {
      fprintf(stderr,"AMRFLDSplit_ReadParameters Warning: non-matching periodic BCs, dim %"ISYM"\n",dim);
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


  return SUCCESS;

}
#endif // TRANSFER
