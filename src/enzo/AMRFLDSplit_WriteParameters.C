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
/  Split Implicit Problem Class, Parameter output routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/
/  PURPOSE: Writes all necessary internal parameters for problem 
/           restart.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"

int AMRFLDSplit::WriteParameters(FILE *fptr)
{

//   if (debug)  printf("Entering AMRFLDSplit::WriteParameters routine\n");
  
  fprintf(fptr, "RadHydroESpectrum = %"ISYM"\n", ESpectrum);
  fprintf(fptr, "RadHydroChemistry = %"ISYM"\n", Nchem);
  fprintf(fptr, "RadHydroModel = %"ISYM"\n", Model);

  fprintf(fptr, "RadHydroMaxDt = %22.16e\n", maxdt);
  fprintf(fptr, "RadHydroMinDt = %22.16e\n", mindt);

  // set restart initial time step to current time step
  if (dt == 0.0) {
    fprintf(fptr, "RadHydroInitDt = %22.16e\n", initdt);
  } else {
    fprintf(fptr, "RadHydroInitDt = %22.16e\n", dt);
  }
  fprintf(fptr, "RadHydroDtControl = %"ISYM"\n", dt_control);
  fprintf(fptr, "RadHydroMaxSubcycles = %22.16e\n", maxsubcycles);
  fprintf(fptr, "RadHydroDtNorm = %22.16e\n", dtnorm);
  fprintf(fptr, "RadHydroDtRadFac = %22.16e\n", dtfac);
  fprintf(fptr, "RadiationScaling = %22.16e\n", ErScale);
  fprintf(fptr, "RadHydroTheta = %22.16e\n", theta);
  fprintf(fptr, "RadiationBoundaryX0Faces = %i %i\n", 
	  BdryType[0][0], BdryType[0][1]);
  if (rank > 1) 
    fprintf(fptr, "RadiationBoundaryX1Faces = %i %i\n", 
	    BdryType[1][0], BdryType[1][1]);
  if (rank > 2) 
    fprintf(fptr, "RadiationBoundaryX2Faces = %i %i\n", 
	    BdryType[2][0], BdryType[2][1]);
  fprintf(fptr, "RadHydroSolType = %i\n", sol_type);
  fprintf(fptr, "RadHydroSolTolerance = %22.16e\n", sol_tolerance);
  fprintf(fptr, "RadHydroMaxMGIters = %i\n", sol_maxit);    
  fprintf(fptr, "RadHydroMGRelaxType = %i\n", sol_rlxtype);    
  fprintf(fptr, "RadHydroMGPreRelax = %i\n", sol_npre);    
  fprintf(fptr, "RadHydroMGPostRelax = %i\n", sol_npost);    

  fprintf(fptr, "RadHydroSolPrec = %i\n", sol_prec);
  fprintf(fptr, "RadHydroSol_precmaxit = %i\n", sol_precmaxit);
  fprintf(fptr, "RadHydroSol_precnpre = %i\n", sol_precnpre);
  fprintf(fptr, "RadHydroSol_precnpost = %i\n", sol_precnpost);
  fprintf(fptr, "RadHydroSol_precJacit = %i\n", sol_precJacit);
  fprintf(fptr, "RadHydroSol_precrelax = %i\n", sol_precrelax);
  fprintf(fptr, "RadHydroSol_precrestol = %g\n", sol_precrestol);

  fprintf(fptr, "WeakScaling = %i\n", &WeakScaling);    

  // if doing an ionization problem (ProblemTypes 410-415),
  // output additional parameters 
  if ((ProblemType >= 410) && (ProblemType <= 415)) {
    fprintf(fptr, "NGammaDot = %22.16e\n", NGammaDot);
    fprintf(fptr, "EtaRadius = %22.16e\n", EtaRadius);
    fprintf(fptr, "EtaCenter = %22.16e %22.16e %22.16e\n",
	    EtaCenter[0], EtaCenter[1], EtaCenter[2]);
  }

  // output relevant units: although these aren't required for restart, 
  // cosmology runs never output the units (why?), making data analyisis tough
  fprintf(fptr, "DensityUnits = %22.16e\n", DenUnits);
  fprintf(fptr, "LengthUnits = %22.16e\n",  LenUnits);
  fprintf(fptr, "TimeUnits = %22.16e\n",    TimeUnits);

  return SUCCESS;
}
#endif
