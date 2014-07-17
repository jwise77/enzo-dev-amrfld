/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Adaptive integration of a given SED function over a prescribed interval
/
/  written by: Daniel Reynolds
/  date:       July 2014
/  modified1:  
/
/  PURPOSE: Computes the integral over frequency space
/                  int_{a}^{b} sed(hnu) d hnu
/           where [a,b] is a prescribed interval in frequency 
/           space, and sed(hnu) is a given frequency-dependent function.
/
/           These are computed using a high-order composite Gauss-Lobatto
/           quadrature rule [O(h^{16}) accurate].
/
/           The SED may be monochromatic, in which if the monochromatic 
/           frequency is within the interval we evaluate the SED; 
/           otherwise we return 0.
/
/           Since frequencies must be strictly greater than zero, 
/           if (a <= 0.0) we set a to zero, and if (b <= 0.0) we set 
/           b to infinity.  For problems that require integration to 
/           infinity, we note that it is typically more efficient 
/           to set b to infinity than to arbitrarily select a large 
/           value.  If the resulting interval does not satisfy the
/           requirement that (a < b), an error flag will be returned, 
/           and the integral value will be set to 0.0.
/
/  NOTE:    Both the SED and integration limits must operate on 
/           frequencies given in eV.  However, if the integral over 
/           true frequency space is desired (in Hz), we merely scale the 
/           resulting eV integral value by an appropriate unit factor.
/
/  ARGUMENTS: 
/           sed -- object of type SED
/           [a,b] -- integration interval, in units of eV
/           convertHz -- flag denoting whether the integral (performed 
/                        over hnu in eV) should be converted to Hz.
/           R -- the resulting value of the integral
/
/  RETURN:
/           returns SUCCESS on successful integration
/           returns FAIL on an illegal interval
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "phys_constants.h"
#include "SED.h"


// utility function
float SED_composite_integral(SED &sed, float a, float b, int n);
float SED_composite_integral2(SED &sed, float a, int n);


// adaptive numerical integration of a specified integrand over 
// a specified interval to a specified accuracy
int SED_integral(SED &sed, float a, float b, bool convertHz, double &R) {

  // initialize result
  R = 0.0;

  // return error on illegal interval
  if ((a > b) && (a > 0) && (b > 0.0))  return FAIL;

  // if the sed is monochromatic:
  // return 1.0 if the sed frequency lies inside the interval, otherwise return 0.0
  if (sed.monochromatic()) {
    if (sed.lower_bound() <= 0.0) {
      R = 0.0;
    } else if ((a > 0) && (sed.lower_bound() < a)) {
      R = 0.0;
    } else if ((b > 0) && (sed.lower_bound() > b)) {
      R = 0.0;
    } else {
      R = 1.0;
    }
    return SUCCESS;
  }

  // set local variables for integration bounds
  float nu_L=a, nu_R=b;
  nu_L = (nu_L > 0.0) ? nu_L : 0.0;
  nu_L = (sed.lower_bound() > nu_L) ? sed.lower_bound() : nu_L;

  bool b_infinite = ((b <= 0.0) && (sed.upper_bound() <= 0.0));
  if (!b_infinite && (b > 0.0) && (sed.upper_bound() > 0.0)) 
    nu_R = (sed.upper_bound() < b) ? sed.upper_bound() : b;
  if (!b_infinite && (b > 0.0) && (sed.upper_bound() <= 0.0)) 
    nu_R = b;
  if (!b_infinite && (b <= 0.0) && (sed.upper_bound() > 0.0)) 
    nu_R = sed.upper_bound();
	    
  // if the integral has now disappeared, return 0.0
  if (!b_infinite && (nu_L >= nu_R))  return SUCCESS;

  // set the number of subintervals and compute approximation
  int N = 5000;
  if (b_infinite)   // improper integral (infinite upper bound)
    R = SED_composite_integral2(sed, nu_L, N);
  else              // definite integral
    R = SED_composite_integral(sed, nu_L, nu_R, N);

  // R is computed based on integration in eV; scale result if Hz was desired
  if (convertHz)  R *= ev2erg/hplanck;

  return SUCCESS;
}

 

// composite numerical integration of a specified integrand over a specified interval
float SED_composite_integral(SED &sed, float a, float b, int n) {

  // nodes/weights defining quadrature method
  float hwid = (b-a)/n;      // subinterval width
  int nodes = 8;             // num quadrature nodes
  float x[] = { -0.18343464249564980493,
		 0.18343464249564980493,
		-0.52553240991632898581,
		 0.52553240991632898581,
		-0.79666647741362673959,
		 0.79666647741362673959,
		-0.96028985649753623168,
		 0.96028985649753623168 };
  float w[] = {  0.36268378337836198296,
		 0.36268378337836198296,
		 0.31370664587788728733,
		 0.31370664587788728733,
		 0.22238103445337447054,
		 0.22238103445337447054,
		 0.10122853629037625915,
		 0.10122853629037625915 };

  float F = 0.0;                        // initialize result
  for (int i=0; i<n; i++) {             // loop over subintervals
    float xmid = a + (i+0.5)*hwid;     // subinterval midpoint
    for (int j=0; j<nodes; j++) {       // loop over quadrature nodes
      float z = xmid + 0.5*hwid*x[j];
      F += w[j]*sed.value(z);           // add contribution to result
    }
  } // end subinterval loop
  return (0.5*hwid*F);    // return final result, scaled by hwid/2

}

 

// remapped composite numerical integration of a specified 
// integrand over the interval [a,infinity]
float SED_composite_integral2(SED &sed, float a, int n) {

  // remapped interval
  float a2 = 0.0;
  float b2 = 1.0;

  // nodes/weights defining quadrature method
  float hwid = (b2-a2)/n;    // subinterval width
  int nodes = 8;             // num quadrature nodes
  float x[] = { -0.18343464249564980493,
		 0.18343464249564980493,
		-0.52553240991632898581,
		 0.52553240991632898581,
		-0.79666647741362673959,
		 0.79666647741362673959,
		-0.96028985649753623168,
		 0.96028985649753623168 };
  float w[] = {  0.36268378337836198296,
		 0.36268378337836198296,
		 0.31370664587788728733,
		 0.31370664587788728733,
		 0.22238103445337447054,
		 0.22238103445337447054,
		 0.10122853629037625915,
		 0.10122853629037625915 };

  float F = 0.0;                        // initialize result
  for (int i=0; i<n; i++) {             // loop over subintervals
    float xmid = a2 + (i+0.5)*hwid;     // subinterval midpoint
    for (int j=0; j<nodes; j++) {       // loop over quadrature nodes
      float z = xmid + 0.5*hwid*x[j];
      F += w[j]*sed.value(a/z)/z/z;     // add contribution to result
    }
  } // end subinterval loop
  return (a*0.5*hwid*F);    // return final result, scaled by a*hwid/2

}

