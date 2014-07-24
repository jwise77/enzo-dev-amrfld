/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Spectral Energy Distribution Abstract Base Class, used in conjunction
/  with AMRFLDSplit radiation solver
/
/  written by: Daniel Reynolds
/  date:       July, 2014
/  modified1:  
/
/  PURPOSE: This class defines SED-specific functions required for 
/  interfacing between desired radiation emission sources (e.g. star 
/  particle types), species cross-section functions, and the 
/  multi-frequency/group radiation solver, AMRFLDSplit.  For a given 
/  emission spectrum, the user must create a derived class that fills 
/  the relevant data structures and instantiates their relevant 
/  SED-specific functions.
/
/  Two example SED intantiations are provided:
/      BlackbodySED.h -- black-body SED at a given temperature
/      MonochromaticSED.h -- monochromatic SED at a given frequency
/
************************************************************************/

#ifndef SED_BASE_CLASS_DEFINED__
#define SED_BASE_CLASS_DEFINED__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

class SED
{

 public:

  /* Should return boolean depending on whether SED exists at only a specific 
     frequency (true), or if it exists across a frequency range (false). */
  virtual bool monochromatic() = 0;

  /* If the SED cuts off below a given frequency, return that cutoff frequency here.  
     If there is no lower bound, return a negative value. 
     If the SED is monochromatic, return that frequency here.
     The value should be given in eV, e.g. if the SED exists from the 
     ionization threshold of Hydrogen, this should be 13.6. */
  virtual float lower_bound() = 0;    

  /* If the SED cuts off above a given frequency, return that cutoff frequency here.  
     If there is no upper bound, return a negative value. 
     If the SED is monochromatic, this is unused.
     The value should be given in eV, e.g. if the SED exists from the 
     ionization threshold of Hydrogen, this should be 13.6. */
  virtual float upper_bound() = 0;    

  /* The main SED function: implements the non-monochromatic SED, and should 
     return the SED value at a given frequency (specified in units of eV).  
     If the SED is monochromatic, this function must still be implemented, but 
     will only be evaluated at the specified frequency. */
  virtual float value(float hnu) = 0;

};
  
#endif
