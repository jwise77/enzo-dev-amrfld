// inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "BlackbodySED.h"
#include "MonochromaticSED.h"
#include "CrossSectionSED.h"
#include "SED.h"
#include "phys_constants.h"

// namespace
using namespace std;

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

class Blackbody_Moment : public virtual SED {
 private:
  BlackbodySED *baseSED;     // base SED
 public:
  Blackbody_Moment(BlackbodySED &base) { this->baseSED = &base; };
  bool monochromatic() { return this->baseSED->monochromatic(); };
  float lower_bound() { return this->baseSED->lower_bound(); };
  float upper_bound() { return this->baseSED->upper_bound(); };
  float value(float hnu) { return this->baseSED->value(hnu)*hnu*ev2erg/hplanck; };
};






// main routine
Eint32 main(Eint32 argc, char* argv[]) {
  
  // reset default output precision
  cout << setprecision(16);

  // First try some monochromatic spectra
  MonochromaticSED m1(13.6);
  cout << "Monochromatic spectrum at 13.6 eV:\n";
  cout << "  monochromatic = " << m1.monochromatic() << endl;
  cout << "  lower_bound = " << m1.lower_bound() << endl;
  cout << "  upper_bound = " << m1.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(m1, 1.0, -1.0, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(m1, 20.0, 100.0, true) << endl << endl;

  MonochromaticSED m2(24.6);
  cout << "Monochromatic spectrum at 24.6 eV:\n";
  cout << "  monochromatic = " << m2.monochromatic() << endl;
  cout << "  lower_bound = " << m2.lower_bound() << endl;
  cout << "  upper_bound = " << m2.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(m2, 1.0, -1.0, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(m2, 20.0, 100.0, true) << endl << endl;



  // Second try some blackbody spectra
  BlackbodySED b1(1.0e5);
  cout << "Blackbody spectrum at T=1e5 K:\n";
  cout << "  monochromatic = " << b1.monochromatic() << endl;
  cout << "  lower_bound = " << b1.lower_bound() << endl;
  cout << "  upper_bound = " << b1.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(b1, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(b1, 1.0, 1.e5, true) << endl;
  cout << "  integral [13.6,1000] = " << SED_integral(b1, 13.6, 1000.0, true) << endl;
  cout << "  integral [13.6,infty] = " << SED_integral(b1, 13.6, -1.0, true) << endl << endl;

  BlackbodySED b2(2.0e5);
  cout << "Blackbody spectrum at T=2e5 K:\n";
  cout << "  monochromatic = " << b2.monochromatic() << endl;
  cout << "  lower_bound = " << b2.lower_bound() << endl;
  cout << "  upper_bound = " << b2.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(b2, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(b2, 1.0, 1.e5, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(b2, 20.0, 100.0, true) << endl << endl;



  // Third try all cross-section spectra
  CrossSectionSED c0(0);
  cout << "Cross-section spectrum, HI:\n";
  cout << "  monochromatic = " << c0.monochromatic() << endl;
  cout << "  lower_bound = " << c0.lower_bound() << endl;
  cout << "  upper_bound = " << c0.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(c0, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(c0, 1.0, 1.e5, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(c0, 20.0, 100.0, true) << endl << endl;

  CrossSection_Ionizing c0I(c0);
  cout << "Cross-section spectrum, HI, ionizing:\n";
  cout << "  monochromatic = " << c0I.monochromatic() << endl;
  cout << "  lower_bound = " << c0I.lower_bound() << endl;
  cout << "  upper_bound = " << c0I.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(c0I, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(c0I, 1.0, 1.e5, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(c0I, 20.0, 100.0, true) << endl << endl;

  CrossSection_Heating c0H(c0);
  cout << "Cross-section spectrum, HI, heating:\n";
  cout << "  monochromatic = " << c0H.monochromatic() << endl;
  cout << "  lower_bound = " << c0H.lower_bound() << endl;
  cout << "  upper_bound = " << c0H.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(c0H, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(c0H, 1.0, 1.e5, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(c0H, 20.0, 100.0, true) << endl << endl;

  CrossSectionSED c1(1);
  cout << "Cross-section spectrum, HeI:\n";
  cout << "  monochromatic = " << c1.monochromatic() << endl;
  cout << "  lower_bound = " << c1.lower_bound() << endl;
  cout << "  upper_bound = " << c1.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(c1, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(c1, 1.0, 1.e5, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(c1, 20.0, 100.0, true) << endl << endl;

  CrossSection_Ionizing c1I(c1);
  cout << "Cross-section spectrum, HeI, ionizing:\n";
  cout << "  monochromatic = " << c1I.monochromatic() << endl;
  cout << "  lower_bound = " << c1I.lower_bound() << endl;
  cout << "  upper_bound = " << c1I.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(c1I, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(c1I, 1.0, 1.e5, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(c1I, 20.0, 100.0, true) << endl << endl;

  CrossSection_Heating c1H(c1);
  cout << "Cross-section spectrum, HeI, heating:\n";
  cout << "  monochromatic = " << c1H.monochromatic() << endl;
  cout << "  lower_bound = " << c1H.lower_bound() << endl;
  cout << "  upper_bound = " << c1H.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(c1H, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(c1H, 1.0, 1.e5, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(c1H, 20.0, 100.0, true) << endl << endl;

  CrossSectionSED c2(2);
  cout << "Cross-section spectrum, HeII:\n";
  cout << "  monochromatic = " << c2.monochromatic() << endl;
  cout << "  lower_bound = " << c2.lower_bound() << endl;
  cout << "  upper_bound = " << c2.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(c2, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(c2, 1.0, 1.e5, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(c2, 20.0, 100.0, true) << endl << endl;

  CrossSection_Ionizing c2I(c2);
  cout << "Cross-section spectrum, HeII, ionizing:\n";
  cout << "  monochromatic = " << c2I.monochromatic() << endl;
  cout << "  lower_bound = " << c2I.lower_bound() << endl;
  cout << "  upper_bound = " << c2I.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(c2I, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(c2I, 1.0, 1.e5, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(c2I, 20.0, 100.0, true) << endl << endl;

  CrossSection_Heating c2H(c2);
  cout << "Cross-section spectrum, HeII, heating:\n";
  cout << "  monochromatic = " << c2H.monochromatic() << endl;
  cout << "  lower_bound = " << c2H.lower_bound() << endl;
  cout << "  upper_bound = " << c2H.upper_bound() << endl;
  cout << "  integral [1,infty] = " << SED_integral(c2H, 1.0, -1.0, true) << endl;
  cout << "  integral [1,1e5] = " << SED_integral(c2H, 1.0, 1.e5, true) << endl;
  cout << "  integral [20,100] = " << SED_integral(c2H, 20.0, 100.0, true) << endl << endl;


  Blackbody_Moment b1M(b1);
  float nu0 = 13.6 * ev2erg / hplanck;
  cout << "Blackbody scaling factors?:\n";
  cout << "  test 1 = " 
       << SED_integral(b1M, 1.0, -1.0, true) / SED_integral(b1, 13.6, -1.0, true) / nu0
       << endl;


  
  return 0;
} // end main

