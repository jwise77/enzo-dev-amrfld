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
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, Constructor routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:  
/
/  PURPOSE: Initializes all values to illegal numbers, and sets all 
/           arrays to NULL;  Requires call to Initialize to actually 
/           set up these values.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"


AMRFLDSplit::AMRFLDSplit()
{

  int bin, dim, src, face;
#ifndef MPI_INT
  int MPI_PROC_NULL = -3;
#endif

  // initialize solver parameters
  sol_tolerance  = 1e-5;  // HYPRE solver tolerance
  sol_maxit      = 200;   // HYPRE max linear iters
  sol_type       = 1;     // BiCGSTab solver
  sol_printl     = 1;     // HYPRE print level
  sol_log        = 1;     // enable HYPRE logging
  sol_prec       = 1;     // Enable HG preconditioner by default
  sol_precmaxit  = 1;     // one preconditioning sweep per Krylov iteration
  sol_precnpre   = 1;     // one pre-relaxation sweep
  sol_precnpost  = 1;     // one post-relaxation sweep
  sol_precJacit  = 2;     // two Jacobi iterations per HG call
  sol_precrelax  = 1;     // weighted Jacobi
  sol_precrestol = 0.0;   // use an iteration-based stop criteria
  sol_rlxtype    = 1;     // weighted Jacobi
  sol_npre       = 1;     // one pre-relaxation sweep
  sol_npost      = 1;     // one post-relaxation sweep

  // initialize limiter parameters
  LimiterRmin = 1.e-20;
  LimiterDmax = 1.e-2;

  // initialize solver diagnostics
  RTtime     = 0.0;
  AMRSolTime = 0.0;
  for (bin=0; bin<MAX_FLD_FIELDS; bin++)
    totIters[bin] = 0;

  // initialize problem grid information (update to defaults)
  rank = -1;
  for (dim=0; dim<3; dim++) {
    LocDims[dim] = 1;
    for (face=0; face<2; face++) {
      OnBdry[dim][face] = false;   // in domain interior
      BdryType[dim][face] = 0;     // periodic
      for (bin=0; bin<MAX_FLD_FIELDS; bin++) 
	BdryVals[bin][dim][face] = NULL;
    }
  }
  
  // initialize time-stepping data
  initdt       = 1.0e20;   // use Hydro step
  maxdt        = 1.0e20;   // no maximum
  mindt        = 0.0;      // no minimum
  maxsubcycles = 1.0;      // no subcycling
  timeAccuracy = 1.0e20;   // no accuracy tolerance
  dt_control   = 2;        // PID controller
  Err_cur      = 1.0;      // perfect error in current step
  Err_new      = 1.0;      // perfect error in previous step
  dtnorm       = 2.0;      // 2-norm for time error estimate
  dtgrowth     = 1.1;      // 10% allowed dt growth per step
  tnew         = -1.0;
  told         = -1.0;
  dt           = -1.0;
  dtrad        = -1.0;
  theta        = 1.0;      // backward Euler


  // initialize radiation and chemistry problem-defining data
  NumRadiationFields = 0;  // no radiation by default
  for (bin=0; bin<MAX_FLD_FIELDS; bin++) {
    FrequencyBand[bin][0] = -1.0;
    FrequencyBand[bin][1] = -1.0;
    FieldMonochromatic[bin] = false;
    FieldNeighbors[bin][0] = false;
    FieldNeighbors[bin][1] = false;
  }
  Isothermal = 0;          // temperature-dependent run
  
  // initialize ionization source parameters
  NumSources = 0;
  for (src=0; bin<MAX_FLD_SOURCES; src++) {
    for (dim=0; dim<3; dim++) {
      SourceLocation[src][dim] = 0.0;
      OriginalSourceLocation[src][dim] = 0.0;
    }
    for (bin=0; bin<MAX_FLD_FIELDS; bin++)
      SourceGroupEnergy[src][bin] = 0.0;
  }
  WeakScaling = 0;         // standard run, do not replicate input sources
  
  // initialize cosmology and scaling constants
  a              = 1.0;
  a0             = 1.0;
  adot           = 0.0;
  adot0          = 0.0;
  aUnits         = 1.0;
  autoScale      = true;
  StartAutoScale = false;
  for (bin=0; bin<MAX_FLD_FIELDS; bin++) 
    ErScale[bin] = 1.0;
  for (bin=0; bin<MAX_FLD_FIELDS; bin++) 
    ErUnits[bin] = 1.0;
  for (bin=0; bin<MAX_FLD_FIELDS; bin++) 
    ErUnits0[bin] = 1.0;
  NiUnits   = 1.0;
  NiUnits0  = 1.0;
  DenUnits  = 1.0;
  DenUnits0 = 1.0;
  LenUnits  = 1.0;
  LenUnits0 = 1.0;
  TimeUnits = 1.0;
  VelUnits  = 1.0;

  // initialize integrals over frequency space
  for (bin=0; bin<MAX_FLD_FIELDS; bin++) {
    intOpacity_HI[bin]    = 0.0;
    intOpacity_HeI[bin]   = 0.0;
    intOpacity_HeII[bin]  = 0.0;
    intIonizing_HI[bin]   = 0.0;
    intIonizing_HeI[bin]  = 0.0;
    intIonizing_HeII[bin] = 0.0;
    intHeating_HI[bin]    = 0.0;
    intHeating_HeI[bin]   = 0.0;
    intHeating_HeII[bin]  = 0.0;
  }

}
#endif
