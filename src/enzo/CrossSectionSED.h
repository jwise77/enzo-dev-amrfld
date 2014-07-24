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

 public:

  // constructor
  CrossSectionSED(int species);
  
  // required functions
  bool monochromatic();
  float lower_bound();
  float upper_bound();
  float value(float hnu);

};
  
#endif
