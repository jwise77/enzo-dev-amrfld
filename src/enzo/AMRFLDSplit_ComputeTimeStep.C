/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, Time Step Computation Routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:  
/
/  PURPOSE: Computes the radiation time step size.  We note that this 
/           value affects the global hydrodynamics time step: 
/                 dt = min(dt_hydro,dt_radiation).
/           This routine is called with scaled arguments.
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
 
 
float AMRFLDSplit::ComputeTimeStep(Eflt64 Eerror)
{

  // If timeAccuracy is set, compute maximum time step as estimate 
  // allowing timeAccuracy relative error (error estimated elsewhere)
  float dt_est = huge_number;    // max time step (normalized units)
  if (timeAccuracy != huge_number) {

    // local variables
    float safety = 1.0;
    float dt_growth = 1.1;
    float den = ((theta < 0.52) && (theta > 0.48)) ? 2.0 : 1.0;

    // compute volume factor to account for fact that error estimate 
    // integrates over domain volume
    float Vol = 1.0;
    for (int i=0; i<rank; i++)
      Vol *= (DomainRightEdge[i]-DomainLeftEdge[i]);
    float VolFac = (dtnorm > 0.0) ? pow(Vol,1.0/dtnorm) : 1.0;
      
    // update old error estimates
    float Err_old = Err_cur;
    Err_cur = Err_new;
    
    // compute new error estimate ratio (include relative tolerance, floor)
    Err_new = max(Eerror/timeAccuracy/VolFac, 1.0e-8);
      
    // Set time step depending on how it has been set up by the user.
    //    dt_control determines the time adaptivity algorithm:
    //             0 -> I controller
    //             1 -> PI controller
    //             2 -> PID controller
    //          else -> original time controller
    if (dt_control == 0) {

      float k1 = -1.0/den;
      dt_est = safety * dt * pow(Err_new,k1);
      if (debug)
	printf(" ComputeTimeStep (I): Err = %g, dt_est = %g (%.1fx change)\n", 
	       Err_new, dt_est, dt_est/dt);

    } else if (dt_control == 1) {

      float k1 = -0.7/den;
      float k2 =  0.4/den;
      dt_est = safety * dt * pow(Err_new,k1) * pow(Err_cur,k2);
      if (debug)
	printf(" ComputeTimeStep (PI): Errs = %g %g, dt_est = %g (%.1fx change)\n", 
	       Err_new, Err_cur, dt_est, dt_est/dt);

    } else if (dt_control == 2) {

      float k1 = -0.49/den;
      float k2 =  0.34/den;
      float k3 = -0.1/den;
      dt_est = safety * dt * pow(Err_new,k1) * pow(Err_cur,k2) * pow(Err_old,k3);
      if (debug)
	printf(" ComputeTimeStep (PID): Errs = %g %g %g, dt_est = %g (%.1fx change)\n", 
	       Err_new, Err_cur, Err_old, dt_est, dt_est/dt);

    } else {

      dt_est = safety * dt / Err_new;
      if (debug)
	printf(" ComputeTimeStep: Err = %g, dt_est = %g (%.1fx change)\n", 
	       Err_new, dt_est, dt_est/dt);

    }

    // enforce growth factor.  Note: dtrad is the last estimate; since dt 
    // could be limited by output times, we use dtrad when limiting growth
    dt_est = min(dt_growth*dtrad, dt_est);

  }

  // account for min/max time step size (according to user)
  dt_est = max(dt_est, mindt);
  dt_est = min(dt_est, maxdt);

  return dt_est;
}

#endif
