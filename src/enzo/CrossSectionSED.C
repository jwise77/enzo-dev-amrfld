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

#include "CrossSectionSED.h"

// constructor
CrossSectionSED::CrossSectionSED(int species) {
  this->Species = species;
  if (species == 0)              // HI lower frequency cutoff
    this->threshold = 13.6;
  if (species == 1)              // HeI lower frequency cutoff
    this->threshold = 24.6;
  if (species == 2)              // HeII lower frequency cutoff
    this->threshold = 54.4;
}

// monochromatic return function
bool CrossSectionSED::monochromatic() { 
  return false;       // not monochromatic
}

// lower bound function
float CrossSectionSED::lower_bound() { 
  return threshold;   // threshold of species
}

// upper bound function
float CrossSectionSED::upper_bound() { 
  return -1.0;        // no upper bound
}

// SED function
float CrossSectionSED::value(float hnu) {

  // all cross-sections are zero below their ionization threshold
  if (hnu < this->threshold)  return 0.0;

  // set some constants
  float nuscaled;      // normalized frequency
  float eps;           // constant in cross section definition
  float sigma=0.0;
  
  // compute sigma based on species choice
  switch (Species) {
    
  case 0:  // HI
    if (hnu == this->threshold)  
      sigma = 6.30e-18;
    else {
      nuscaled = hnu/this->threshold;
      eps = sqrt(nuscaled - 1.0);
      sigma = (6.30e-18)*POW(nuscaled,-4.0)
	*exp(4.0-4.0*atan(eps)/eps)/(1.0-exp(-2.0*pi/eps));
    }
    break;

  case 1:  // HeI
    if (hnu == this->threshold)
      sigma = 7.42e-18;
    else {
      nuscaled = hnu/this->threshold;
      sigma = (7.42e-18)*(1.66*POW(nuscaled,-2.05) - 0.66*POW(nuscaled,-3.05));
      /* // updated cross-section from Verland et al (1996) -- see calc_rates.src
	 float x = hnu*ev2erg/hplanck/13.61 - 0.4434 - 1.0;
	 float y = sqrt(x*x + 2.136*2.136);
	 sigma = 9.492e-16 * (x*x + 2.039*2.039) * POW(y, -3.906) 
	 * POW(1.0+sqrt(y/1.469), -3.188); */
    }
    break;
      
  case 2:  // HeII

    if (hnu == this->threshold) 
      sigma = 1.575e-18;
    else {
      nuscaled = hnu/this->threshold;
      eps = sqrt(nuscaled - 1.0);
      sigma = (1.575e-18)*POW(nuscaled,-4.0)
	*exp(4.0-4.0*atan(eps)/eps)/(1.0-exp(-2.0*pi/eps));
    }
    break;

  }

  return sigma;
}
