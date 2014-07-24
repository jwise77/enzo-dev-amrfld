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

#ifndef BLACKBODY_SED_DEFINED__
#define BLACKBODY_SED_DEFINED__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "phys_constants.h"
#include "SED.h"

class BlackbodySED : public virtual SED {

 private:

  // SED-defining data
  float Temperature;      // Temperature (in Kelvin) of the blackbody source

 public:

  // constructor
  BlackbodySED(float Temp);

  // required functions
  bool monochromatic();
  float lower_bound();
  float upper_bound();
  float value(float hnu);

};
  
#endif
