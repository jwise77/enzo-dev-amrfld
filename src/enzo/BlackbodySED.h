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
//#include "preincludes.h"
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
  BlackbodySED(float Temp) { this->Temperature = Temp; };

  // monochromatic return function
  bool monochromatic() { return false; };  // not monochromatic

  // lower bound function
  float lower_bound() { return 0.0; };     // lower bound of hnu=0

  // upper bound function
  float upper_bound() { return -1.0; };    // no upper bound

  // SED function
  float value(float hnu) {
    float nu = hnu*ev2erg/hplanck;    // convert frequency to Hz
    return (8.0*pi*hplanck*POW(nu/clight,3)/(exp(hplanck*nu/kboltz/Temperature)-1.0));
  };

};
  
#endif
