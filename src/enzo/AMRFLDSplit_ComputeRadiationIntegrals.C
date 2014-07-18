/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, Radiation numerical integral routine 
/
/  written by: Daniel Reynolds
/  date:       July 2014
/  modified1:  
/
/  PURPOSE: Computes the frequency-space integrals
/                  1/|binwidth| * int_{bin} sigmaHI(nu) d nu
/                  1/|binwidth| * int_{bin} sigmaHeI(nu) d nu
/                  1/|binwidth| * int_{bin} sigmaHeII(nu) d nu
/                  1/|binwidth| * int_{bin} sigmaHI(nu)/nu d nu
/                  1/|binwidth| * int_{bin} sigmaHeI(nu)/nu d nu
/                  1/|binwidth| * int_{bin} sigmaHeII(nu)/nu d nu
/                  1/|binwidth| * int_{bin} sigmaHI(nu)*(1-nuHI/nu) d nu
/                  1/|binwidth| * int_{bin} sigmaHeI(nu)*(1-nuHeI/nu) d nu
/                  1/|binwidth| * int_{bin} sigmaHeII(nu)*(1-nuHeII/nu) d nu
/           for each radiation group (and associated radiation bin); or
/                  sigmaHI(nu_freq)
/                  sigmaHeI(nu_freq)
/                  sigmaHeII(nu_freq)
/                  sigmaHI(nu_freq) / nu_freq^2
/                  sigmaHeI(nu_freq) / nu_freq^2
/                  sigmaHeII(nu_freq) / nu_freq^2
/                  sigmaHI(nu_freq) * (1/nu_freq - nuHI / nu_freq^2)
/                  sigmaHeI(nu_freq) * (1/nu_freq - nuHeI / nu_freq^2)
/                  sigmaHeII(nu_freq) * (1/nu_freq - nuHeII / nu_freq^2)
/           for each monochromatic radiation field.  Here nuHI, nuHeI and 
/           nuHeII are the ionization thresholds of each species, sigmaHI(nu),
/           sigmaHeI(nu) and sigmaHeII(nu) are the species cross-sections.  
/           For the monochromatic fields, nu_freq is the frequency.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
#include "phys_constants.h"
#include "CrossSectionSED.h"


// function prototypes
int SED_integral(SED &sed, float a, float b, bool convertHz, double &R);



int AMRFLDSplit::ComputeRadiationIntegrals() {

  //  if (debug)  printf("Entering AMRFLDSplit::ComputeRadiationIntegrals\n");

  // create SED objects
  CrossSectionSED sHI(0, 0);
  CrossSectionSED sHeI(1, 0);
  CrossSectionSED sHeII(2, 0);
  CrossSectionSED sHI_ion(0, 1);
  CrossSectionSED sHeI_ion(1, 1);
  CrossSectionSED sHeII_ion(2, 1);
  CrossSectionSED sHI_heat(0, 2);
  CrossSectionSED sHeI_heat(1, 2);
  CrossSectionSED sHeII_heat(2, 2);

  // loop over radiation fields
  for (int ibin=0; ibin<NumRadiationFields; ibin++) {

    // radiation frequency
    if (FieldMonochromatic[ibin]) {

      float hnu = FrequencyBand[ibin][0];   // eV
      float nu  = hnu * ev2erg / hplanck;   // Hz

      // sigma(nu_freq)
      intOpacity_HI[ibin]   = sHI.value(hnu);
      intOpacity_HeI[ibin]  = sHeI.value(hnu);
      intOpacity_HeII[ibin] = sHeII.value(hnu);

      // sigma(nu_freq) / nu_freq^2
      intIonizing_HI[ibin]   = sHI_ion.value(hnu)/nu;
      intIonizing_HeI[ibin]  = sHeI_ion.value(hnu)/nu;
      intIonizing_HeII[ibin] = sHeII_ion.value(hnu)/nu;

      // sigma(nu_freq) * (1/nu_freq - nu_species / nu_freq^2)
      intHeating_HI[ibin]   = sHI_heat.value(hnu)/nu;
      intHeating_HeI[ibin]  = sHeI_heat.value(hnu)/nu;
      intHeating_HeII[ibin] = sHeII_heat.value(hnu)/nu;
      
    // radiation group
    } else {

      float hnu_L = FrequencyBand[ibin][0];  // eV
      float hnu_R = FrequencyBand[ibin][1];  // eV
      float binwidth = (hnu_R-hnu_L) * ev2erg / hplanck;   // Hz
      float IntVal;

      // 1/|binwidth| * int_{bin} sigmaHI(nu) d nu
      if (SED_integral(sHI, hnu_L, hnu_R, true, IntVal) != SUCCESS)
	ENZO_FAIL("Error in computing HI opacity integral")
      intOpacity_HI[ibin] = IntVal / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeI(nu) d nu
      if (SED_integral(sHeI, hnu_L, hnu_R, true, IntVal) != SUCCESS)
	ENZO_FAIL("Error in computing HeI opacity integral")
      intOpacity_HeI[ibin] = IntVal / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeII(nu) d nu
      if (SED_integral(sHeII, hnu_L, hnu_R, true, IntVal) != SUCCESS)
	ENZO_FAIL("Error in computing HeII opacity integral")
      intOpacity_HeII[ibin] = IntVal / binwidth;
      
      // 1/|binwidth| * int_{bin} sigmaHI(nu)/nu d nu
      if (SED_integral(sHI_ion, hnu_L, hnu_R, true, IntVal) != SUCCESS)
	ENZO_FAIL("Error in computing HI ionization integral")
      intIonizing_HI[ibin] = IntVal / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeI(nu)/nu d nu
      if (SED_integral(sHeI_ion, hnu_L, hnu_R, true, IntVal) != SUCCESS)
	ENZO_FAIL("Error in computing HeI ionization integral")
      intIonizing_HeI[ibin] = IntVal / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeII(nu)/nu d nu
      if (SED_integral(sHeII_ion, hnu_L, hnu_R, true, IntVal) != SUCCESS)
	ENZO_FAIL("Error in computing HeII ionization integral")
      intIonizing_HeII[ibin] = IntVal / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHI(nu)*(1-nuHI/nu) d nu
      if (SED_integral(sHI_heat, hnu_L, hnu_R, true, IntVal) != SUCCESS)
	ENZO_FAIL("Error in computing HI photoheating integral")
      intHeating_HI[ibin] = IntVal / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeI(nu)*(1-nuHeI/nu) d nu
      if (SED_integral(sHeI_heat, hnu_L, hnu_R, true, IntVal) != SUCCESS)
	ENZO_FAIL("Error in computing HeI photoheating integral")
      intHeating_HeI[ibin] = IntVal / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeII(nu)*(1-nuHeII/nu) d nu
      if (SED_integral(sHeII_heat, hnu_L, hnu_R, true, IntVal) != SUCCESS)
	ENZO_FAIL("Error in computing HeII photoheating integral")
      intHeating_HeII[ibin] = IntVal / binwidth;

    }

  }  // end loop over radiation fields

  return SUCCESS;
}
 
#endif
