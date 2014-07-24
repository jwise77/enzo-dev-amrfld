/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Blackbody radiation SED instantiation of the SED base class.
/
/  written by: Daniel Reynolds
/  date:       July, 2014
/  modified1:  
/
/  PURPOSE: This is an example SED class implementation, that is used 
/  for test problems.
/
************************************************************************/

#include "BlackbodySED.h"

// constructor
BlackbodySED::BlackbodySED(float Temp) { this->Temperature = Temp; };

// monochromatic return function
bool BlackbodySED::monochromatic() { 
  return false;  // not monochromatic
}

// lower bound function
float BlackbodySED::lower_bound() { 
  return 0.0;    // lower bound of hnu=0
}

// upper bound function
float BlackbodySED::upper_bound() { 
  return -1.0;   // no upper bound
}

// main SED function
float BlackbodySED::value(float hnu) {
  float nu = hnu*ev2erg/hplanck;    // convert frequency to Hz
  return (8.0*pi*hplanck*POW(nu/clight,3)/(exp(hplanck*nu/kboltz/Temperature)-1.0));
}
