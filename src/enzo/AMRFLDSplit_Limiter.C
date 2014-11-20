/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class -- Flux Limiter function
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Computes the flux limiter for a given cell face.
/
************************************************************************/
#ifdef TRANSFER
#include "phys_constants.h"
#include "AMRFLDSplit.h"


float AMRFLDSplit::Limiter(float E1, float E2, float k1, float k2, 
			   float nUn, float lUn, float dxi)
{

  // limiter bounds
  float Rmin = min(LimiterRmin/lUn, 1e-20);            // 1st is astro/cosmo, 2nd is lab frame
  float Dmax = max(LimiterDmax * clight * lUn, 1e20);  // 1st is astro/cosmo, 2nd is lab frame
  float Eavg = max((E1 + E2)*0.5, 1e-30);              // arithmetic mean, bound from zero
  k1 = max(k1, 1e-20);                                 // bound away from zero since we'll divide
  k2 = max(k2, 1e-20);
  float kap = 2.0*k1*k2/(k1+k2)*nUn;                   // harmonic mean
  float R = max(dxi*fabs(E1-E2)/Eavg, Rmin);

  // compute limiter
  return (min(clight/sqrt(9.0*kap*kap + R*R), Dmax));

}

#endif
