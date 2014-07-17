/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Species cross-section SED instantiation of SED base class.
/
/  written by: Daniel Reynolds
/  date:       July, 2014
/  modified1:  
/
/  PURPOSE: Implementation of species cross-sections (and related 
/  functions) used in integrals that couple radiation fields with 
/  heating, photo-ionization and chemistry-dependent opacities.
/
************************************************************************/

#ifndef CROSS_SECTION_SED_DEFINED__
#define CROSS_SECTION_SED_DEFINED__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "phys_constants.h"
#include "SED.h"

class CrossSectionSED : public virtual SED {

 private:

  // SED-defining data
  float threshold;        // ionization frequency threshold
  int Species;            /* choice of chemical species:
                                0 -> HI
                                1 -> HeI
                                2 -> HeII */

  int SED_type;           /* choice of integrand:
                                0 -> sigma(nu)
                                1 -> sigma(nu)/nu
                                2 -> sigma(nu)*(1-nu0/nu)
                             where sigma(nu) is the cross-section, and nu0 
                             is the ionization threshold of the species */

 public:

  // constructor
  CrossSectionSED(int species, int sed_type) {
    this->Species = species;
    this->SED_type = sed_type;
    if (species == 0)              // HI lower frequency cutoff
      this->threshold = 13.6;
    if (species == 1)              // HeI lower frequency cutoff
      this->threshold = 24.6;
    if (species == 2)              // HeII lower frequency cutoff
      this->threshold = 54.4;
  };

  // monochromatic return function
  bool monochromatic() { return false; };    // not monochromatic

  // lower bound function
  float lower_bound() { return threshold; }; // threshold of species

  // upper bound function
  float upper_bound() { return -1.0; };      // no upper bound

  // SED function
  float value(float hnu) {
    if (hnu < threshold)  return 0.0;

    // set some constants
    float nu0_HI   = 13.6*ev2erg/hplanck;   // ionization threshold of HI (hz)
    float nu0_HeI  = 24.6*ev2erg/hplanck;   // ionization threshold of HeI (hz)
    float nu0_HeII = 54.4*ev2erg/hplanck;   // ionization threshold of HeII (hz)
    float nuscaled;                         // normalized frequency
    float eps;                              // constant in cross section definition
    float sigma=0.0;

    // compute sigma based on species choice
    switch (Species) {

    case 0:  // HI
      if (hnu == 13.6)  
	sigma = 6.30e-18;
      else {
	nuscaled = hnu/13.6;
	eps = sqrt(nuscaled - 1.0);
	sigma = (6.30e-18)*POW(nuscaled,-4.0)
	  *exp(4.0-4.0*atan(eps)/eps)/(1.0-exp(-2.0*pi/eps));
      }
      break;

    case 1:  // HeI
      if (hnu == 24.6)
	sigma = 7.42e-18;
      else {
	nuscaled = hnu/24.6;
	sigma = (7.42e-18)*(1.66*POW(nuscaled,-2.05) - 0.66*POW(nuscaled,-3.05));
 	/* // updated cross-section from Verland et al (1996) -- see calc_rates.src
	   float x = hnu*ev2erg/hplanck/13.61 - 0.4434 - 1.0;
	   float y = sqrt(x*x + 2.136*2.136);
	   sigma = 9.492e-16 * (x*x + 2.039*2.039) * POW(y, -3.906) 
                 * POW(1.0+sqrt(y/1.469), -3.188); */
      }
      break;
      
    case 2:  // HeII

      if (hnu == 54.4) 
	sigma = 1.575e-18;
      else {
	nuscaled = hnu/54.4;
	eps = sqrt(nuscaled - 1.0);
	sigma = (1.575e-18)*POW(nuscaled,-4.0)
	  *exp(4.0-4.0*atan(eps)/eps)/(1.0-exp(-2.0*pi/eps));
      }
      break;

    }

    // compute output based on SED choice
    switch (SED_type) {

    case 0:  // sigma(nu)
      return sigma;
      break;

    case 1:  // sigma(nu)/nu
      return (sigma*hplanck/hnu/ev2erg);
      break;

    case 2:  // sigma(nu)*(1-nu0/nu)
      return (sigma*(1.0 - threshold/hnu));
      break;

    }

  };

};
  
#endif
