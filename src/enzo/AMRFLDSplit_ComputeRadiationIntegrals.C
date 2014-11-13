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
float SED_integral(SED &sed, float a, float b, bool convertHz);

// SED types for computing cross-section integrals 
class CrossSection_Ionizing : public virtual SED {
 private:
  CrossSectionSED *baseSED;     // base SED
 public:
  CrossSection_Ionizing(CrossSectionSED &base) { this->baseSED = &base; };
  bool monochromatic() { return this->baseSED->monochromatic(); };
  float lower_bound() { return this->baseSED->lower_bound(); };
  float upper_bound() { return this->baseSED->upper_bound(); };
  float value(float hnu) { return this->baseSED->value(hnu)*hplanck/hnu/ev2erg; };
};

class CrossSection_Heating : public virtual SED {
 private:
  CrossSectionSED *baseSED;     // base SED
 public:
  CrossSection_Heating(CrossSectionSED &base) { this->baseSED = &base; };
  bool monochromatic() { return this->baseSED->monochromatic(); };
  float lower_bound() { return this->baseSED->lower_bound(); };
  float upper_bound() { return this->baseSED->upper_bound(); };
  float value(float hnu) { 
    return this->baseSED->value(hnu)*(1.0 - this->baseSED->lower_bound()/hnu); 
  };
};



int AMRFLDSplit::ComputeRadiationIntegrals() {

  //  if (debug)  printf("Entering AMRFLDSplit::ComputeRadiationIntegrals\n");

  // create SED objects for each chemical species and integrand
  CrossSectionSED sHI(0);
  CrossSectionSED sHeI(1);
  CrossSectionSED sHeII(2);
  CrossSection_Ionizing sHI_i(sHI);
  CrossSection_Ionizing sHeI_i(sHeI);
  CrossSection_Ionizing sHeII_i(sHeII);
  CrossSection_Heating sHI_h(sHI);
  CrossSection_Heating sHeI_h(sHeI);
  CrossSection_Heating sHeII_h(sHeII);

  // loop over radiation fields
  for (int ibin=0; ibin<NumRadiationFields; ibin++) {

    // radiation frequency
    if (FieldMonochromatic[ibin]) {

      // sigma(nu_freq)
      intOpacity_HI[ibin]   = sHI.value(FrequencyBand[ibin][0]);
      intOpacity_HeI[ibin]  = sHeI.value(FrequencyBand[ibin][0]);
      intOpacity_HeII[ibin] = sHeII.value(FrequencyBand[ibin][0]);

      // sigma(nu_freq) / nu_freq
      intIonizing_HI[ibin]   = sHI_i.value(FrequencyBand[ibin][0]);
      intIonizing_HeI[ibin]  = sHeI_i.value(FrequencyBand[ibin][0]);
      intIonizing_HeII[ibin] = sHeII_i.value(FrequencyBand[ibin][0]);

      // sigma(nu_freq) * (1 - nu_species / nu_freq)
      intHeating_HI[ibin]   = sHI_h.value(FrequencyBand[ibin][0]);
      intHeating_HeI[ibin]  = sHeI_h.value(FrequencyBand[ibin][0]);
      intHeating_HeII[ibin] = sHeII_h.value(FrequencyBand[ibin][0]);
      
    // radiation group
    } else {

      float hnu_L = FrequencyBand[ibin][0];  // eV
      float hnu_R = FrequencyBand[ibin][1];  // eV
      float binwidth = (hnu_R-hnu_L) * ev2erg / hplanck;   // Hz

      // 1/|binwidth| * int_{bin} sigmaHI(nu) d nu
      intOpacity_HI[ibin] = SED_integral(sHI, hnu_L, hnu_R, true) / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeI(nu) d nu
      intOpacity_HeI[ibin] = SED_integral(sHeI, hnu_L, hnu_R, true) / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeII(nu) d nu
      intOpacity_HeII[ibin] = SED_integral(sHeII, hnu_L, hnu_R, true) / binwidth;
      
      // 1/|binwidth| * int_{bin} sigmaHI(nu)/nu d nu
      intIonizing_HI[ibin] = SED_integral(sHI_i, hnu_L, hnu_R, true) / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeI(nu)/nu d nu
      intIonizing_HeI[ibin] = SED_integral(sHeI_i, hnu_L, hnu_R, true) / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeII(nu)/nu d nu
      intIonizing_HeII[ibin] = SED_integral(sHeII_i, hnu_L, hnu_R, true) / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHI(nu)*(1-nuHI/nu) d nu
      intHeating_HI[ibin] = SED_integral(sHI_h, hnu_L, hnu_R, true) / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeI(nu)*(1-nuHeI/nu) d nu
      intHeating_HeI[ibin] = SED_integral(sHeI_h, hnu_L, hnu_R, true) / binwidth;

      // 1/|binwidth| * int_{bin} sigmaHeII(nu)*(1-nuHeII/nu) d nu
      intHeating_HeII[ibin] = SED_integral(sHeII_h, hnu_L, hnu_R, true) / binwidth;

    }

    // output bin integrals
    if (debug) {
      printf("Computed Integrals for Radiation Field %"ISYM":\n", ibin);
      printf("  intOpacity_HI    = %22.16e\n", intOpacity_HI[ibin]);
      printf("  intOpacity_HeI   = %22.16e\n", intOpacity_HeI[ibin]);
      printf("  intOpacity_HeII  = %22.16e\n", intOpacity_HeII[ibin]);
      printf("  intIonizing_HI   = %22.16e\n", intIonizing_HI[ibin]);
      printf("  intIonizing_HeI  = %22.16e\n", intIonizing_HeI[ibin]);
      printf("  intIonizing_HeII = %22.16e\n", intIonizing_HeII[ibin]);
      printf("  intHeating_HI    = %22.16e\n", intHeating_HI[ibin]);
      printf("  intHeating_HeI   = %22.16e\n", intHeating_HeI[ibin]);
      printf("  intHeating_HeII  = %22.16e\n", intHeating_HeII[ibin]);
    }

  }  // end loop over radiation fields

  return SUCCESS;
}

#endif
