/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, Parameter output routine
/
/  written by: Daniel Reynolds
/  date:       July 2014
/
/  PURPOSE: Writes all necessary internal parameters for problem restart.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"

int AMRFLDSplit::WriteParameters(FILE *fptr)
{

  // all non-root processors return immediately
  if (MyProcessorNumber != ROOT_PROCESSOR)  return SUCCESS;

//   if (debug)  printf("Entering AMRFLDSplit::WriteParameters routine\n");
  int ibin, isrc, dim;

  // radiation fields and frequency bands
  fprintf(fptr, "AMRFLDNumRadiationFields = %"ISYM"\n", NumRadiationFields);
  for (ibin=0; ibin<NumRadiationFields; ibin++)
    fprintf(fptr, "AMRFLDFrequencyBand[%"ISYM"] = %22.16e %22.16e\n", 
   	    ibin, FrequencyBand[ibin][0], FrequencyBand[ibin][1]);

  // radiation scaling factors
  for (ibin=0; ibin<NumRadiationFields; ibin++)
    fprintf(fptr, "AMRFLDRadiationScaling["ISYM"] = %22.16e\n", 
	    ibin, ErScale[ibin]);
  if (autoScale[0]) {
    fprintf(fptr, "AMRFLDAutomaticScaling = 1\n");
  } else {
    fprintf(fptr, "AMRFLDAutomaticScaling = 0\n");
  }

  // radiation sources
  fprintf(fptr, "AMRFLDNumSources = %"ISYM"\n", NumSources);
  if (WeakScaling) {
    for (isrc=0; isrc<NumSources; isrc++) 
      fprintf(fptr, "AMRFLDSourceLocation[%"ISYM"] = %22.16e %22.16e %22.16e\n", isrc, 
	      OriginalSourceLocation[isrc][0], OriginalSourceLocation[isrc][1], 
	      OriginalSourceLocation[isrc][2]);
  } else {
    for (isrc=0; isrc<NumSources; isrc++)
      fprintf(fptr, "AMRFLDSourceLocation[%"ISYM"] = %22.16e %22.16e %22.16e\n", isrc, 
	      SourceLocation[isrc][0], SourceLocation[isrc][1], SourceLocation[isrc][2]);
  }
  for (isrc=0; isrc<NumSources; isrc++)
    for (ibin=0; ibin<NumRadiationFields; ibin++)
      fprintf(fptr, "AMRFLDSourceGroupEnergy[%"ISYM"][%"ISYM"] = %22.16e\n", 
	      isrc, ibin, SourceGroupEnergy[isrc][ibin]);
      
  // limiter parameters
  fprintf(fptr, "AMRFLDLimiterRmin = %22.16e\n", LimiterRmin);
  fprintf(fptr, "AMRFLDLimiterDmax = %22.16e\n", LimiterDmax);

  // model parameters
  fprintf(fptr, "AMRFLDIsothermal = %"ISYM"\n", Isothermal);

  // time-stepping parameters
  fprintf(fptr, "AMRFLDMaxDt = %22.16e\n", maxdt);
  fprintf(fptr, "AMRFLDMinDt = %22.16e\n", mindt);
  if (dt == 0.0) {
    fprintf(fptr, "AMRFLDInitDt = %22.16e\n", initdt);
  } else {      // set restart initial time step to current time step
    fprintf(fptr, "AMRFLDInitDt = %22.16e\n", dt);
  }
  fprintf(fptr, "AMRFLDDtControl = %"ISYM"\n", dt_control);
  fprintf(fptr, "AMRFLDMaxSubcycles = %22.16e\n", maxsubcycles);
  fprintf(fptr, "AMRFLDDtNorm = %22.16e\n", dtnorm);
  fprintf(fptr, "AMRFLDDtGrowth = %22.16e\n", dtgrowth);
  fprintf(fptr, "AMRFLDTimeAccuracy = %22.16e\n", timeAccuracy);
  fprintf(fptr, "AMRFLDTheta = %22.16e\n", theta);

  // boundary condition types
  fprintf(fptr, "AMRFLDRadiationBoundaryX0 = %i %i\n", 
	  BdryType[0][0], BdryType[0][1]);
  if (rank > 1) 
    fprintf(fptr, "AMRFLDRadiationBoundaryX1 = %i %i\n", 
	    BdryType[1][0], BdryType[1][1]);
  if (rank > 2) 
    fprintf(fptr, "AMRFLDRadiationBoundaryX2 = %i %i\n", 
	    BdryType[2][0], BdryType[2][1]);

  // solver parameters
  fprintf(fptr, "AMRFLDSolType = %"ISYM"\n", sol_type);
  fprintf(fptr, "AMRFLDSolTolerance = %22.16e\n", sol_tolerance);
  fprintf(fptr, "AMRFLDMaxMGIters = %"ISYM"\n", sol_maxit);    
  fprintf(fptr, "AMRFLDMGRelaxType = %"ISYM"\n", sol_rlxtype);    
  fprintf(fptr, "AMRFLDMGPreRelax = %"ISYM"\n", sol_npre);    
  fprintf(fptr, "AMRFLDMGPostRelax = %"ISYM"\n", sol_npost);    

  fprintf(fptr, "AMRFLDSolPrec = %"ISYM"\n", sol_prec);
  fprintf(fptr, "AMRFLDSol_precmaxit = %"ISYM"\n", sol_precmaxit);
  fprintf(fptr, "AMRFLDSol_precnpre = %"ISYM"\n", sol_precnpre);
  fprintf(fptr, "AMRFLDSol_precnpost = %"ISYM"\n", sol_precnpost);
  fprintf(fptr, "AMRFLDSol_precJacit = %"ISYM"\n", sol_precJacit);
  fprintf(fptr, "AMRFLDSol_precrelax = %"ISYM"\n", sol_precrelax);
  fprintf(fptr, "AMRFLDSol_precrestol = %g\n", sol_precrestol);

  // flag for setting up weak-scaling runs
  fprintf(fptr, "AMRFLDWeakScaling = %"ISYM"\n", WeakScaling);

  // output relevant units: although these aren't required for restart, 
  // cosmology runs never output the units (why?), making data analysis tough
  fprintf(fptr, "DensityUnits = %22.16e\n", DenUnits);
  fprintf(fptr, "LengthUnits = %22.16e\n",  LenUnits);
  fprintf(fptr, "TimeUnits = %22.16e\n",    TimeUnits);

  return SUCCESS;
}
#endif
