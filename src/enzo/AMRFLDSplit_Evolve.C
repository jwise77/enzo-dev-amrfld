/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, Evolve Routine
/
/  written by: Daniel Reynolds
/  date:       July 2014
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
#ifdef TRANSFER
 
#include "AMRFLDSplit.h"
#include "phys_constants.h"

//#define FAIL_ON_NAN
#define NO_FAIL_ON_NAN

// function prototypes
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
#ifdef FAST_SIB
			  SiblingGridList SiblingList[],
#endif
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
float TotalFieldValue(LevelHierarchyEntry *LevelArray[], 
		      int levelStart, int levelFinish, int field);



// This routine evolves the radiation field in an operator-split fashion, 
// subcycling the physics in the following manner: 
//     dt_rad <= dt_hydro
// Prior to completion, the routine also updates the maximum time step the 
// overall Grid module can take to meet a maximum subcycling ratio of 
// radiation to hydrodynamics.
int AMRFLDSplit::Evolve(LevelHierarchyEntry *LevelArray[], int level, 
			HierarchyEntry *Grids[], int NumberOfGrids,
			TopGridData *MetaData, ExternalBoundary *Exterior, 
#ifdef FAST_SIB
			SiblingGridList SiblingList[],
#endif
			float dthydro)
{

#ifdef AMR_SOLVE

  // start MPI timer for overall solver
#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif

  ////////////////////////////////////
  // Problem Setup Phase

  if (debug)  printf("\n AMRFLDSplit Evolve:\n");

  // Find a grid on this level to get current time
  HierarchyEntry* ThisGrid=LevelArray[level]->GridHierarchyEntry;
  tnew = ThisGrid->GridData->ReturnTime();
  told = tnew;

  // update internal units for current times
  if (this->UpdateUnits(told, tnew) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_Evolve: Error in UpdateUnits.");

  // scale radiation fields on all relevant grids to solver units
  if (ScaleFields(LevelArray, level, false) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_Evolve: Error in ScaleFields.");
    
  // collect and output statistics on current radiation field
  float Etyp[MAX_FLD_FIELDS];
  float Emax[MAX_FLD_FIELDS];
  if (PrintStats(LevelArray, level, Etyp, Emax) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_Evolve: PrintStats error");

  // if autoScale enabled, determine scaling factor updates here
  float ScaleCorrTol = 1.e-2;
  float ErScaleCorr[NumRadiationFields];
  for (int ibin=0; ibin<NumRadiationFields; ibin++) {
    ErScaleCorr[ibin] = 1.0;
    if (StartAutoScale[ibin] && autoScale[ibin])
      if ((Emax[ibin] - Etyp[ibin]) > ScaleCorrTol*Emax[ibin])
	ErScaleCorr[ibin] = Emax[ibin];
  }

  // insert Enzo grids into an AMRsolve hierarchy
  AMRsolve_Hierarchy* hierarchy = new AMRsolve_Hierarchy;
  hierarchy->enzo_attach_fld(LevelArray, NumRadiationFields, level, MAX_DEPTH_OF_HIERARCHY);

  // Initialize the AMRsolve domain & hierarchy
  AMRsolve_Domain domain(3, DomainLeftEdge, DomainRightEdge);
  bool is_periodic[] = {BdryType[0][0]==0, BdryType[1][0]==0, BdryType[2][0]==0};
  hierarchy->initialize(domain, *pmpi, is_periodic);

  // initialize variables that we'll use throughout the time subcycling
  float  stime2, ftime2;     // radiation, chemistry timers
  Eflt64 Eerror, errtmp;     // temporary used for computing time step
  int    radstep, radstop;   // subcycle iterators


  ////////////////////////////////////
  // Problem Solve Phase

  // add redshifting sources to emissivity fields for each bin
  if (this->Redshifting(LevelArray, level) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_Evolve: Redshifting error");

  // add parameter-specified sources to emissivity fields for each bin
  if (this->RadiationSource(LevelArray, level, tnew) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_Evolve: RadiationSource error");

  // internal time-stepping loop to catch up with Hydro time
  float end_time = tnew + dthydro;
  radstop = 0;
  for (radstep=0; radstep<=maxsubcycles*100; radstep++) {
      
    // start MPI timer for radiation solver
#ifdef USE_MPI
    stime2 = MPI_Wtime();
#else
    stime2 = 0.0;
#endif

    // update time-step information
    told = tnew;

    // keep trying time steps until radiation solver succeeds. 
    // Note: if we reach the minimum time step size, RadStep will call ENZO_FAIL
    int recompute_step = 1;
    while (recompute_step) {

      // update time-step information.  Note: dtrad was set on previous 
      // iteration of solver, or by user input for first iteration
      tnew = told + dtrad;
      if ((tnew - end_time)/end_time > -1.0e-14) {   // don't exceed synchronization time
	tnew = end_time;
	radstop = 1;
      }
      dt = tnew - told;
      if (debug) 
	printf("\n subcycled rad %"ISYM": dt=%7.1e, t=%7.1e (hydro dt=%7.1e, t=%7.1e)\n",
	       radstep, dt, tnew, dthydro, end_time);
      
      if (debug)  printf(" ----------------------------------------------------------------------\n");

      // loop over radiation bins, taking a step of each
      Eerror = 0.0;
      recompute_step = 0;
      for (int ibin=0; ibin<NumRadiationFields; ibin++) {
	recompute_step += this->RadStep(ibin, LevelArray, level, hierarchy, 
					Etyp[ibin], Emax[ibin], &errtmp);
	if (recompute_step) break;
	Eerror = max(Eerror, errtmp);  // accumulate the maximum error over all bins
      }

      // if the radiation step was unsuccessful, back-track to previous 
      // step and pull back on dtrad
      if (recompute_step) {
	dtrad = max(dtrad*0.5, mindt);
	tnew = told;
	radstop = 0;
      }

    }
    if (debug)  printf("\n");

	
    // stop MPI timer for radiation solver, increment total
#ifdef USE_MPI
    ftime2 = MPI_Wtime();
#else
    ftime2 = 0.0;
#endif
    AMRSolTime += ftime2-stime2;
    
    // update the radiation time step size for next time step
    //   (limit growth at each cycle)
    float dt_est = this->ComputeTimeStep(Eerror);
    dtrad = min(dt_est, dtgrowth*dtrad);


    //////////////////////////////////////////////////////////////////////
    // Initial attempt at just filling ghost zones and overlap regions 
    // throughout the hierarchy using the SetBoundaryConditions routine.
    // We only need to communicate the radiation fields, but Enzo has no 
    // current routines for exchanging a subset of fields at a time.
    // Replace this with the new 'FieldObjects' infrastructure when it 
    // is ready.

    if (SetBoundaryConditions(Grids, NumberOfGrids, 
#ifdef FAST_SIB
			      SiblingList, 
#endif
			      level, MetaData,
			      Exterior, NULL) == FAIL)
      ENZO_FAIL("SetBoundaryConditions() failed!\n");


    //////////////////////////////////////////////////////////////////////


    // break out of time-stepping loop if we've reached the end
    if (radstop)  break;
	
  } // end outer radiation time-stepping loop
  

  ////////////////////////////////////
  // Problem Cleanup and Preparation for Next Call Phase

  // clean up amrsolve structures
  hierarchy->enzo_detach();
  delete hierarchy;
  hierarchy = NULL;

  // fill the chemistry and cooling rates
  if (this->FillRates(LevelArray, level) == FAIL)
    ENZO_FAIL("AMRFLDSplit Evolve: FillRates error");

  // if diagnostics enabled, output RMS of photo-ionization and photo-heating rates
  if (diags) {
    float phHI_rms       = TotalFieldValue(LevelArray, level, MAX_DEPTH_OF_HIERARCHY, kphHI);
    float phHeI_rms      = TotalFieldValue(LevelArray, level, MAX_DEPTH_OF_HIERARCHY, kphHeI);
    float phHeII_rms     = TotalFieldValue(LevelArray, level, MAX_DEPTH_OF_HIERARCHY, kphHeII);
    float photogamma_rms = TotalFieldValue(LevelArray, level, MAX_DEPTH_OF_HIERARCHY, PhotoGamma);
    if (debug) {
      printf("   photo-ionization and photo-heating rates:\n");
      printf("            phHI = %13.7e\n", phHI_rms);
      if (!RadiativeTransferHydrogenOnly) {
	printf("           phHeI = %13.7e\n", phHeI_rms);
	printf("          phHeII = %13.7e\n", phHeII_rms);
      }
      printf("      photogamma = %13.7e\n", photogamma_rms);
    }    
  }

  // update the radiation time step size at this level for next time step 
  if (dtrad != huge_number) {
    for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
      for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	   Temp=Temp->NextGridThisLevel)
	Temp->GridHierarchyEntry->GridData->SetMaxRadiationDt(dtrad*maxsubcycles);
  }

  // scale radiation fields on all relevant grids back to enzo units
  if (ScaleFields(LevelArray, level, true) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_Evolve: Error in ScaleFields.");

  // zero out emissivity fields for all grids on this processor
  if (ResetEmissivityFields(LevelArray, level) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_Evolve: Error in ResetEmissivityFields.");

  // update scaling factors to account for new values
  for (int ibin=0; ibin<NumRadiationFields; ibin++)
    if (StartAutoScale[ibin] && autoScale[ibin]) 
      ErScale[ibin] *= ErScaleCorr[ibin];

  // stop MPI timer, add to cumulative clock, output to stdout
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  RTtime += ftime-stime;
  if (debug)  printf("AMRFLDSplit cumulative time = %g (AMRSolve = %g)\n\n",
		     RTtime, AMRSolTime);

  return SUCCESS;

#else
  ENZO_FAIL("AMRFLDSplit_Evolve ERROR: module requires AMR_SOLVE to be set!");
  return FAIL;
#endif   // AMR_SOLVE
 
}



#ifdef AMR_SOLVE

// This routine evolves the radiation subsystem within the AMRFLDSplit module.  
// This is performed in a robust manner; if the current radiation subsystem 
// cannot be solved to the desired tolerance (meaning the time step is 
// likely much too large), it will return a flag to the calling routine to 
// have the step recomputed with a smaller time step size.
int AMRFLDSplit::RadStep(int ibin, LevelHierarchyEntry *LevelArray[], 
			 int level, AMRsolve_Hierarchy *hierarchy, 
			 float Etyp, float Emax, Eflt64 *Eerror)
{

  // update internal units for current times
  if (this->UpdateUnits(told, tnew) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in UpdateUnits.");

  // rescale dt, told, tnew to physical values for use within solver
  dt   *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;

  //   compute current opacity for this bin
  if (this->Opacity(ibin, LevelArray, level, tnew) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_RadStep: Opacity failure!!");
    
  //   enforce boundary conditions on new radiation field
  if (this->EnforceBoundary(ibin, LevelArray) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: EnforceBoundary failure!!");
    
  //   copy current radiation field into temporary field (KPhHI)
  //   on all grids owned by this processor, all levels here and down;
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	  // set grid dimension information
	  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	  int ghXl = DEFAULT_GHOST_ZONES;
	  int n3[] = {1, 1, 1};
	  for (int dim=0; dim<rank; dim++)
	    n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	            - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	  int x0len = n3[0] + 2*ghXl;
	  int x1len = n3[1] + 2*ghYl;
	  int x2len = n3[2] + 2*ghZl;

	  // access old/new radiation fields (old stored in KPhHI)
	  float *Eold = Temp->GridHierarchyEntry->GridData->AccessKPhHI();
	  float *Enew = AccessRadiationField(ibin, Temp->GridHierarchyEntry);
	  if (Eold == NULL)
	    ENZO_FAIL("AMRFLDSplit_RadStep: cannot access KPhHI field!!\n");
	  if (Enew == NULL) {
	    ENZO_VFAIL("AMRFLDSplit_RadStep: cannot access radiation field %"ISYM"!!\n", ibin);
	  }

	  // copy new radiation field into into old field
	  for (int k=0; k<x0len*x1len*x2len; k++)  Eold[k] = Enew[k];

      }  // end loop over grids on this proc
  
  // calculate emissivity source norm
  float glob_eta = TotalFieldValue(LevelArray, level, MAX_DEPTH_OF_HIERARCHY, Emissivity0+ibin);
  if (debug)   printf("   total emissivity %"ISYM" = %g\n", ibin, glob_eta);

  // determine if the radiation system needs to be solved in this step:
  //   if not, set error to 0, reset units, return with a successful step; 
  float srctol=1.0e-30;
  float radtol=1.0e-8;
  if ((glob_eta*dt < srctol) && (fabs(Etyp-Emax) < radtol)) {
    *Eerror = 0.1;
    dt   /= TimeUnits;
    told /= TimeUnits;
    tnew /= TimeUnits;
    return 0;
  }

  // if we've made it here, turn on automatic scaling for next time step
  StartAutoScale[ibin] = true;

  // set up and solve this radiation field's equation via amrsolve

  // Initialize the amrsolve FLD solver
  AMRsolve_Hypre_FLD amrfldsolve(*hierarchy, *amrsolve_params, ibin, sol_prec);
  amrfldsolve.init_hierarchy();
  // if (debug)   hierarchy->print();
  amrfldsolve.init_stencil();
  amrfldsolve.init_graph();

  //    initialize amrfldsolve system
  amrfldsolve.init_elements(dt, theta, NiUnits, NiUnits0, LenUnits, LenUnits0, 
			    ErUnits[ibin], ErUnits0[ibin], BdryType);


  // /////////  AMRsolve_Hypre_FLD testing routine ////////
  // amrfldsolve.tester();
  // MPI_Barrier(MPI_COMM_WORLD);
  // ENZO_FAIL("Stopping run (testing)");
  // //////////////////////////////////////////////////////


  //    solve amrfldsolve system
  amrfldsolve.solve();
  Eflt64 finalresid = amrfldsolve.residual();
  Eint32 Sits = amrfldsolve.iterations();
  if (debug) printf("   lin resid = %.1e (tol = %.1e), its = %i\n",
		    finalresid, sol_tolerance, Sits);


  // check solution, and if the solver fails:
  //  * if we have room to decrease the time step size, do so and allow remainder 
  //    of function to complete (to reset units, etc.), but have calling routine 
  //    update dt and compute step again.
  //  * if we have don't have room to decrease the time step size, output 
  //    parameters and error message.
  int recompute_step = 0;
  if (amrfldsolve.evaluate() != 0) {
    if (dt > mindt*TimeUnits) {
      recompute_step = 1;
    }
    else {
      fprintf(stderr,"AMRFLDSplit_RadStep: could not achieve prescribed tolerance!\n");
      fprintf(stderr, "Printing current hierarchy:\n");
      if (debug)   hierarchy->print();
      
      // dump amrsolve matrices, module parameters to disk
      amrfldsolve.abort_dump();
      if (debug) {
	fprintf(stderr,"Dumping AMRFLDSplit module parameters to file RTdump.params\n");
	FILE *fptr = fopen("RTdump.params", "w");
	this->WriteParameters(fptr);
	fclose(fptr);
      }
      
      ENZO_FAIL("Error in AMRFLDSplit_RadStep");
    }
  }
  if (debug)  printf(" ======================================================================\n");

  // increment Enzo radiation field with amrsolve solution
  amrfldsolve.update_enzo();

  // enforce a solution floor on radiation
  float epsilon=1.0;      // radiation floor, set to 100*roundoff
  while ((1.0 + epsilon) > 1.0)  epsilon*=0.5;
  epsilon *= 100.0;
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	  // set grid dimension information
	  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	  int ghXl = DEFAULT_GHOST_ZONES;
	  int n3[] = {1, 1, 1};
	  for (int dim=0; dim<rank; dim++)
	    n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	            - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	  int x0len = n3[0] + 2*ghXl;
	  int x1len = n3[1] + 2*ghYl;
	  int x2len = n3[2] + 2*ghZl;

	  // access new radiation fields
	  float *Enew = AccessRadiationField(ibin, Temp->GridHierarchyEntry);
	  if (Enew == NULL) {
	    ENZO_VFAIL("AMRFLDSplit_RadStep: cannot access radiation field %"ISYM"!!\n",
		       ibin);
	  }
	  
	  // enforce floor on new field
	  for (int k=0; k<x0len*x1len*x2len; k++)  
	    Enew[k] = max(Enew[k], epsilon);

      }  // end loop over grids on this proc


  // If this was a successful solve, compute relative change in radiation field 
  // solution over time step
  if (recompute_step == 0) 
    *Eerror = amrfldsolve.rdiff_norm(dtnorm, 0.1);

  // rescale dt, told, tnew back to normalized values
  dt   /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;

  return recompute_step;

}

#endif   // AMR_SOLVE


// This routine scales the current radiation fields on this process, 
// from level down, either to/from solver units based on the "Rescale" 
// argument (false -> solver units, true -> Enzo units)
int AMRFLDSplit::ScaleFields(LevelHierarchyEntry *LevelArray[], int level, bool Rescale) {

  // scale radiation fields on all relevant grids to/from solver units
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry *Temp=LevelArray[thislevel]; Temp; Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	// set grid dimension information
	int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	int ghXl = DEFAULT_GHOST_ZONES;
	int n3[] = {1, 1, 1};
	for (int dim=0; dim<rank; dim++)
	  n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	          - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	int x0len = n3[0] + 2*ghXl;
	int x1len = n3[1] + 2*ghYl;
	int x2len = n3[2] + 2*ghZl;
	  
	// scale radiation fields
	for (int ibin=0; ibin<NumRadiationFields; ibin++) {
	  float *Enew = AccessRadiationField(ibin, Temp->GridHierarchyEntry);
	  if (Enew == NULL) {
	    ENZO_VFAIL("AMRFLDSplit_ScaleFields error: cannot access radiation bin %"ISYM"\n",ibin);
	  }
	  if (Rescale)
	    for (int k=0; k<x0len*x1len*x2len; k++)  Enew[k] *= ErScale[ibin];
	  else
	    for (int k=0; k<x0len*x1len*x2len; k++)  Enew[k] /= ErScale[ibin];
	}

      }  // end loop over grids on this proc
  
  return SUCCESS;
}



// This routine zeros out the emissivity fields from this level down
// (used after the step is over and energy has been deposited)
int AMRFLDSplit::ResetEmissivityFields(LevelHierarchyEntry *LevelArray[], int level) {

  // zero out emissivity fields for all grids on this processor
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	// set grid dimension information
	int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	int ghXl = DEFAULT_GHOST_ZONES;
	int n3[] = {1, 1, 1};
	for (int dim=0; dim<rank; dim++)
	  n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	          - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	int x0len = n3[0] + 2*ghXl;
	int x1len = n3[1] + 2*ghYl;
	int x2len = n3[2] + 2*ghZl;
	  
	// reset emissivity fields
	for (int ibin=0; ibin<NumRadiationFields; ibin++) {
	  float *eta = AccessEmissivityField(ibin, Temp->GridHierarchyEntry);
	  if (eta == NULL) {
	    ENZO_VFAIL("AMRFLDSplit_Evolve error: cannot access emissivity bin %"ISYM"\n", ibin);
	  }
	  for (int k=0; k<x0len*x1len*x2len; k++)  eta[k]=0.0;
	}

      }  // end loop over grids on this proc
    
  return SUCCESS;
}



// This routine compiles and outputs statistics on the current 
// radiation field, returning the RMS and maximum values of the
// radiation fields, and printing statistics to the screen
int AMRFLDSplit::PrintStats(LevelHierarchyEntry *LevelArray[], 
			    int level, float *Etyp, float *Emax) {

  // generate local statistics on this process
  int ibin;
  for (ibin=0; ibin<MAX_FLD_FIELDS; ibin++) {
    Etyp[ibin] = 0.0;
    Emax[ibin] = 0.0;
  }
  float Enum=0.0;
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry *Temp=LevelArray[thislevel]; Temp; Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	  // set grid dimension information
	  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	  int ghXl = DEFAULT_GHOST_ZONES;
	  int n3[] = {1, 1, 1};
	  for (int dim=0; dim<rank; dim++)
	    n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	            - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	  int x0len = n3[0] + 2*ghXl;
	  int x1len = n3[1] + 2*ghYl;
	  int x2len = n3[2] + 2*ghZl;
	  Enum += n3[0]*n3[1]*n3[2];
	  
	  // scale radiation fields, and accumulate field statistics
	  for (ibin=0; ibin<NumRadiationFields; ibin++) {
	    float *Enew = AccessRadiationField(ibin, Temp->GridHierarchyEntry);
	    if (Enew == NULL) {
	      ENZO_VFAIL("AMRFLDSplit_PrintStats error: cannot access radiation bin %"ISYM"\n",
			 ibin);
	    }

	    for (int k=ghZl; k<n3[2]+ghZl; k++) 
	      for (int j=ghYl; j<n3[1]+ghYl; j++) 
		for (int i=ghXl; i<n3[0]+ghXl; i++) {
		  Etyp[ibin] += Enew[(k*x1len+j)*x0len+i]*Enew[(k*x1len+j)*x0len+i];
		  Emax[ibin] = max(Emax[ibin], fabs(Enew[(k*x1len+j)*x0len+i]));
		}
	  }

      }  // end loop over grids on this proc

  // communicate among all procs to get statistics of current radiation field
#ifdef USE_MPI
  MPI_Datatype FDataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg nbins_p1 = NumRadiationFields+1;
  MPI_Arg nbins = NumRadiationFields;
  float dtmpIn[NumRadiationFields+1];
  for (ibin=0; ibin<NumRadiationFields; ibin++)  dtmpIn[ibin] = Etyp[ibin];
  dtmpIn[NumRadiationFields] = Enum;
  float dtmpOut[NumRadiationFields+1];
  if (MPI_Allreduce(Emax, dtmpOut, nbins, FDataType, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS)
    ENZO_FAIL("AMRFLDSplit_PrintStats error in MPI_Allreduce\n");
  for (ibin=0; ibin<NumRadiationFields; ibin++)  
    Emax[ibin] = dtmpOut[ibin];
  MPI_Allreduce(dtmpIn, dtmpOut, nbins_p1, FDataType, MPI_SUM, MPI_COMM_WORLD);
  for (ibin=0; ibin<NumRadiationFields; ibin++)  
    Etyp[ibin] = sqrt(dtmpOut[ibin]/dtmpOut[NumRadiationFields]);
#else
  for (ibin=0; ibin<NumRadiationFields; ibin++)  
    Etyp[ibin] = sqrt(Etyp[ibin]/Enum);
#endif
  if (debug) {
    printf("   current internal quantities:\n");
    for (ibin=0; ibin<NumRadiationFields; ibin++) 
      printf("      E%"ISYM": rms = %13.7e, max = %13.7e\n", 
	     ibin, Etyp[ibin], Emax[ibin]);
  }    

  return SUCCESS;
}


#endif   // TRANSFER
