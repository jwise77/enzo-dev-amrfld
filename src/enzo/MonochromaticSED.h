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

#ifndef MONOCHROMATIC_SED_DEFINED__
#define MONOCHROMATIC_SED_DEFINED__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "SED.h"

class MonochromaticSED : public virtual SED {

 private:

  // SED-defining data
  float frequency;

 public:

  // constructor
  MonochromaticSED(float Frequency) { this->frequency = Frequency; };

  // monochromatic return function
  bool monochromatic() { return true; };

  // lower bound function
  float lower_bound() { return frequency; };

  // upper bound function
  float upper_bound() { return frequency; };

  // SED function
  float value(float hnu) { return 1.0; };

};
  
#endif
