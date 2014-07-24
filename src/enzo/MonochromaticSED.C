/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Monochromatic radiation SED instantiation of the SED base class.
/
/  written by: Daniel Reynolds
/  date:       July, 2014
/  modified1:  
/
/  PURPOSE: This is an example SED class implementation, that is used 
/  for test problems.
/
************************************************************************/

#include "MonochromaticSED.h"

// constructor
MonochromaticSED::MonochromaticSED(float Frequency) { 
  this->frequency = Frequency; 
};

// monochromatic return function
bool MonochromaticSED::monochromatic() { 
  return true; 
};

// lower bound function
float MonochromaticSED::lower_bound() { 
  return frequency; 
};

// upper bound function
float MonochromaticSED::upper_bound() { 
  return frequency; 
};

// SED function
float MonochromaticSED::value(float hnu) { 
  return 1.0; 
};
