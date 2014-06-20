/*****************************************************************************
 *                                                                           *
 * Copyright 2010 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Single-Group, Multi-species, AMR, Gray Flux-Limited Diffusion 
/  Split Implicit Problem Class, Evolve Routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
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
#include "CosmologyParameters.h"


//#define FAIL_ON_NAN
#define NO_FAIL_ON_NAN


/* function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
#ifdef FAST_SIB
			  SiblingGridList SiblingList[],
#endif
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);


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


  // scale radiation fields on all relevant grids to solver units; output statistics
  LevelHierarchyEntry *Temp;
  float Etyp[MAX_RADIATION_BINS];
  float Emax[MAX_RADIATION_BINS];
  int ibin;
  for (ibin=0; ibin<MAX_RADIATION_BINS; ibin++) {
    Etyp[ibin] = 0.0;
    Emax[ibin] = 0.0;
  }
  float Enum=0.0;
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (Temp=LevelArray[thislevel]; Temp; Temp=Temp->NextGridThisLevel)
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
	  Enum += x0len*x1len*x2len;
	  
	  // scale radiation fields
	  for (ibin=0; ibin<NumBins; ibin++) {
	    float *Enew = AccessRadiationBin(ibin, Temp->GridHierarchyEntry);
	    if (Enew == NULL) {
	      ENZO_VFAIL("AMRFLDSplit_Evolve error: cannot access radiation bin %"ISYM"\n",
			 ibin);
	    }
	    for (int k=0; k<x0len*x1len*x2len; k++)  Enew[k]/=ErScale[ibin];

	    for (int k=0; k<x0len*x1len*x2len; k++)  {
	      Etyp[ibin] += Enew[k]*Enew[k];
	      Emax[ibin] = max(Emax[ibin], fabs(Enew[k]));
	    }
	  }

      }  // end loop over grids on this proc

  // communicate among all procs to get information on current radiation field
#ifdef USE_MPI
  MPI_Datatype FDataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg nbins_p1 = NumBins+1;
  MPI_Arg nbins = NumBins;
  float dtmpIn[NumBins+1];
  for (ibin=0; ibin<NumBins; ibin++)  
    dtmpIn[ibin] = Etyp[ibin];
  dtmpIn[NumBins] = Enum;
  float dtmpOut[NumBins+1];
  MPI_Allreduce(&Emax, &dtmpOut, nbins, FDataType, MPI_MAX, MPI_COMM_WORLD);
  for (ibin=0; ibin<NumBins; ibin++)  
    Emax[ibin] = dtmpOut[ibin];
  MPI_Allreduce(&dtmpIn, &dtmpOut, nbins_p1, FDataType, MPI_SUM, MPI_COMM_WORLD);
  for (ibin=0; ibin<NumBins; ibin++)  
    Etyp[ibin] = sqrt(dtmpOut[ibin]/dtmpOut[NumBins]);
#else
  for (ibin=0; ibin<NumBins; ibin++)  
    Etyp[ibin] = sqrt(Etyp[ibin]/Enum);
#endif
  if (debug) {
    printf("   current internal quantities:\n");
    for (ibin=0; ibin<NumBins; ibin++) 
      printf("      E%"ISYM": rms = %13.7e, max = %13.7e\n", 
	     ibin, Etyp[ibin], Emax[ibin]);
  }    

  // if autoScale enabled, determine scaling factor updates here
  float ScaleCorrTol = 1.e-2;
  float ErScaleCorr[NumBins];
  for (ibin=0; ibin<NumBins; ibin++)  ErScaleCorr[ibin] = 1.0;
  if (StartAutoScale && autoScale)
    for (ibin=0; ibin<NumBins; ibin++)
      if ((Emax[ibin] - Etyp[ibin]) > ScaleCorrTol*Emax[ibin])
	ErScaleCorr[ibin] = Emax[ibin];

  // insert Enzo grids into an AMRsolve hierarchy
  AMRsolve_Hierarchy* hierarchy = new AMRsolve_Hierarchy;
  hierarchy->enzo_attach_fld(LevelArray, NumBins, level, MAX_DEPTH_OF_HIERARCHY);

  // Initialize the AMRsolve domain & hierarchy
  AMRsolve_Domain domain(3, DomainLeftEdge, DomainRightEdge);
  bool is_periodic[] = {BdryType[0][0]==0, BdryType[1][0]==0, BdryType[2][0]==0};
  hierarchy->initialize(domain, *pmpi, is_periodic);

  // initialize variables that we'll use throughout the time subcycling
  float stime2, ftime2;   // radiation, chemistry timers
  Eflt64 Eerror, errtmp;  // temporary used for computing time step
  int radstep, radstop;   // subcycle iterators


  ////////////////////////////////////
  // Problem Solve Phase

  // Find a grid on this level, owned by this processor, to get current time
  //// NOTE: IS IT SAFE TO ASSUME THAT ALL PROCS WILL OWN AT LEAST ONE GRID ON EACH LEVEL?  IS THIS EVEN NECESSARY?
  HierarchyEntry* ThisGrid=NULL;
  for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
    ThisGrid = Temp->GridHierarchyEntry;
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber())
      break;
  }

  // safety check: if ThisGrid was not set, just pick the first one on this level
  if (ThisGrid == NULL) 
    ThisGrid = LevelArray[level]->GridHierarchyEntry;

  // Get general time-related information
  tnew = ThisGrid->GridData->ReturnTime();

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
      for (ibin=0; ibin<NumBins; ibin++) {
	recompute_step += this->RadStep(ibin, LevelArray, level, hierarchy, 
					Etyp[ibin], Emax[ibin], &errtmp);
	if (recompute_step) break;
	Eerror = max(Eerror, errtmp);  // accumulate the maximum change over all bins
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

  // update the radiation time step size at this level for next time step 
  if (dtrad != huge_number) {
    for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
      for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	   Temp=Temp->NextGridThisLevel)
	Temp->GridHierarchyEntry->GridData->SetMaxRadiationDt(dtrad*maxsubcycles);
  }

  // scale radiation fields on all relevant grids back to enzo units
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
	  
	  // scale radiation fields
	  for (ibin=0; ibin<NumBins; ibin++) {
	    float *Enew = AccessRadiationBin(ibin, Temp->GridHierarchyEntry);
	    if (Enew == NULL) {
	      ENZO_VFAIL("AMRFLDSplit_Evolve error: cannot access radiation bin %"ISYM"\n",
			 ibin);
	    }
	    for (int k=0; k<x0len*x1len*x2len; k++)  Enew[k]*=ErScale[ibin];
	  }

      }  // end loop over grids on this proc
    
  // update scaling factors to account for new values
  if (StartAutoScale && autoScale) 
    for (ibin=0; ibin<NumBins; ibin++)
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
int AMRFLDSplit::RadStep(int Bin, LevelHierarchyEntry *LevelArray[], 
			 int level, AMRsolve_Hierarchy *hierarchy, 
			 float Etyp, float Emax, Eflt64 *Eerror)
{

  // update internal Enzo units for current times
  float TempUnits, RadUnits;
  double MassUnits;
  DenUnits0 = LenUnits0 = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits0, &LenUnits0, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, told) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, told) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in RadiationGetUnits.");
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  ErUnits0[Bin] = RadUnits*ErScale[Bin];
  NiUnits0 = (Nchem == 0) ? 1.0 : DenUnits0/mp;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(told, &a0, &adot0) != SUCCESS) 
      ENZO_FAIL("AMRFLDSplit_RadStep: Error in CosmologyComputeExpansionFactor.");
  adot0 /= TimeUnits;  // rescale to physical units

  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tnew) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tnew) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in RadiationGetUnits.");
  ErUnits[Bin] = RadUnits*ErScale[Bin];
  NiUnits = (Nchem == 0) ? 1.0 : DenUnits/mp;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) != SUCCESS) 
      ENZO_FAIL("AMRFLDSplit_RadStep: Error in CosmologyComputeExpansionFactor.");
  adot /= TimeUnits;  // rescale to physical units
    
  // rescale dt, told, tnew to physical values for use within solver
  dt   *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;

  //   compute emissivity at this internal time step (if provided internally)
  int eta_set = 0;
  float srcNorm = 0.0;
#ifdef EMISSIVITY 
  if (StarMakerEmissivityField > 0)  eta_set = 1;
#endif
  if (eta_set == 0) 
    srcNorm = this->RadiationSource(Bin, LevelArray, level, tnew);
    
  //   compute opacity at this internal time step
  if (this->Opacity(Bin, LevelArray, level, tnew) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_RadStep: Opacity failure!!");
    
  //   enforce boundary conditions on old/new radiation fields
  if (this->EnforceBoundary(Bin, LevelArray) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: EnforceBoundary failure!!");
    
  //   copy current radiation field into temporary field (KPhHI)
  //   on all grids owned by this processor, all levels here and down;
  //   calculate emissivty source norm (if not already calculated)
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
	  float *Enew = AccessRadiationBin(Bin, Temp->GridHierarchyEntry);
	  if (Eold == NULL)
	    ENZO_FAIL("AMRFLDSplit_RadStep: cannot access KPhHI field!!\n");
	  if (Enew == NULL) {
	    ENZO_VFAIL("AMRFLDSplit_RadStep: cannot access radiation field %"ISYM"!!\n",
		       Bin);
	  }

	  // update emissivity norm (if needed)
	  if (eta_set) {
	    float dVscale = 1;
	    for (int dim=0; dim<rank; dim++)
	      dVscale *= (Temp->GridData->GetGridRightEdge(dim) 
			- Temp->GridData->GetGridLeftEdge(dim)) 
		      / n3[dim] / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	    float *eta = AccessEmissivityBin(Bin, Temp->GridHierarchyEntry);
	    if (eta == NULL) {
	      ENZO_VFAIL("AMRFLDSplit_RadStep: cannot access emissivity field %"ISYM"!!\n",
			 Bin);
	    }
	    for (int k=ghZl; k<n3[2]+ghZl; k++) 
	      for (int j=ghYl; j<n3[1]+ghYl; j++) 
		for (int i=ghXl; i<n3[0]+ghXl; i++) 
		  srcNorm += eta[(k*x1len+j)*x0len+i]*eta[(k*x1len+j)*x0len+i]*dVscale;
	  }
	  
	  // copy new radiation field into into old field
	  for (int k=0; k<x0len*x1len*x2len; k++)  Eold[k] = Enew[k];

      }  // end loop over grids on this proc

  // finish off emissivity norm (if needed)
  if (eta_set) {
    float glob_eta;
#ifdef USE_MPI
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg ONE = 1;
    MPI_Allreduce(&srcNorm,&glob_eta,ONE,DataType,MPI_SUM,MPI_COMM_WORLD);
#else
    glob_eta = srcNorm;
#endif
    srcNorm = sqrt(glob_eta);
  }
  if (debug)   printf("   emissivity %"ISYM" norm = %g\n", Bin, srcNorm);

  // determine if the radiation system needs to be solved in this step:
  //   if not, set error = 0.0, reset units, return with a successful step; 
  //   if so, continue
  float srctol=1.0e-30;
  float radtol=1.0e-8;
  if ((srcNorm*dt < srctol) && (fabs(Etyp-Emax) < radtol)) {
    *Eerror = 0.0;
    dt   /= TimeUnits;
    told /= TimeUnits;
    tnew /= TimeUnits;
    return 0;
  }

  // if we've made it here, turn on automatic scaling for next step
  StartAutoScale = true;


  // set up and solve radiation equation via amrsolve

  // Initialize the amrsolve FLD solver
  AMRsolve_Hypre_FLD amrfldsolve(*hierarchy, *amrsolve_params, Bin, sol_prec);
  amrfldsolve.init_hierarchy();
  // if (debug)   hierarchy->print();
  amrfldsolve.init_stencil();
  amrfldsolve.init_graph();

  //    initialize amrfldsolve system
  Eflt64 a_ = a;
  Eflt64 a0_ = a0;
  Eflt64 adot_  = (ESpectrum[Bin] < 0) ? 0.0 : adot;
  Eflt64 adot0_ = (ESpectrum[Bin] < 0) ? 0.0 : adot0;
  amrfldsolve.init_elements(dt, Nchem, theta, a_, a0_, adot_, 
			    adot0_, NiUnits, NiUnits0, LenUnits, 
			    LenUnits0, ErUnits[Bin], ErUnits0[Bin], BdryType);


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
  float epsilon=1.0;      // radiation floor
  //  epsilon = 1.0e-8;
  while ((1.0 + epsilon) > 1.0)  epsilon*=0.5;
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
	  float *Enew = AccessRadiationBin(Bin, Temp->GridHierarchyEntry);
	  if (Enew == NULL) {
	    ENZO_VFAIL("AMRFLDSplit_RadStep: cannot access radiation field %"ISYM"!!\n",
		       Bin);
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


#endif   // TRANSFER
