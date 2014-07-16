// inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "BlackbodySED.h"
#include "MonochromaticSED.h"
#include "CrossSectionSED.h"
#include "BlackbodySED.h"

// namespace
using namespace std;

// function prototypes
int SED_integral(SED &sed, float a, float b, bool convertHz, double &R);

// main routine
Eint32 main(Eint32 argc, char* argv[]) {
  
  // temporary variables
  float R;

  // reset default output precision
  cout << setprecision(16);

  // First try some monochromatic spectra
  MonochromaticSED m1(13.6);
  cout << "Monochromatic spectrum at 13.6 eV:\n";
  cout << "  monochromatic = " << m1.monochromatic() << endl;
  cout << "  lower_bound = " << m1.lower_bound() << endl;
  cout << "  upper_bound = " << m1.upper_bound() << endl;
  SED_integral(m1, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(m1, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;

  MonochromaticSED m2(24.6);
  cout << "Monochromatic spectrum at 24.6 eV:\n";
  cout << "  monochromatic = " << m2.monochromatic() << endl;
  cout << "  lower_bound = " << m2.lower_bound() << endl;
  cout << "  upper_bound = " << m2.upper_bound() << endl;
  SED_integral(m2, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(m2, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;



  // Second try some blackbody spectra
  BlackbodySED b1(1.0e5);
  cout << "Blackbody spectrum at T=1e5 K:\n";
  cout << "  monochromatic = " << b1.monochromatic() << endl;
  cout << "  lower_bound = " << b1.lower_bound() << endl;
  cout << "  upper_bound = " << b1.upper_bound() << endl;
  SED_integral(b1, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(b1, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(b1, 13.6, 1000.0, true, R);
  cout << "  integral [13.6,1000] = " << R << endl;
  SED_integral(b1, 13.6, -1.0, true, R);
  cout << "  integral [13.6,infty] = " << R << endl << endl;

  BlackbodySED b2(2.0e5);
  cout << "Blackbody spectrum at T=2e5 K:\n";
  cout << "  monochromatic = " << b2.monochromatic() << endl;
  cout << "  lower_bound = " << b2.lower_bound() << endl;
  cout << "  upper_bound = " << b2.upper_bound() << endl;
  SED_integral(b2, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(b2, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(b2, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;



  // Third try all cross-section spectra
  CrossSectionSED c00(0, 0);
  cout << "Cross-section spectrum, HI, type 0:\n";
  cout << "  monochromatic = " << c00.monochromatic() << endl;
  cout << "  lower_bound = " << c00.lower_bound() << endl;
  cout << "  upper_bound = " << c00.upper_bound() << endl;
  SED_integral(c00, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(c00, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(c00, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;

  CrossSectionSED c01(0, 1);
  cout << "Cross-section spectrum, HI, type 1:\n";
  cout << "  monochromatic = " << c01.monochromatic() << endl;
  cout << "  lower_bound = " << c01.lower_bound() << endl;
  cout << "  upper_bound = " << c01.upper_bound() << endl;
  SED_integral(c01, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(c01, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(c01, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;

  CrossSectionSED c02(0, 2);
  cout << "Cross-section spectrum, HI, type 2:\n";
  cout << "  monochromatic = " << c02.monochromatic() << endl;
  cout << "  lower_bound = " << c02.lower_bound() << endl;
  cout << "  upper_bound = " << c02.upper_bound() << endl;
  SED_integral(c02, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(c02, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(c02, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;

  CrossSectionSED c10(1, 0);
  cout << "Cross-section spectrum, HeI, type 0:\n";
  cout << "  monochromatic = " << c10.monochromatic() << endl;
  cout << "  lower_bound = " << c10.lower_bound() << endl;
  cout << "  upper_bound = " << c10.upper_bound() << endl;
  SED_integral(c10, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(c10, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(c10, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;

  CrossSectionSED c11(1, 1);
  cout << "Cross-section spectrum, HeI, type 1:\n";
  cout << "  monochromatic = " << c11.monochromatic() << endl;
  cout << "  lower_bound = " << c11.lower_bound() << endl;
  cout << "  upper_bound = " << c11.upper_bound() << endl;
  SED_integral(c11, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(c11, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(c11, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;

  CrossSectionSED c12(1, 2);
  cout << "Cross-section spectrum, HeI, type 2:\n";
  cout << "  monochromatic = " << c12.monochromatic() << endl;
  cout << "  lower_bound = " << c12.lower_bound() << endl;
  cout << "  upper_bound = " << c12.upper_bound() << endl;
  SED_integral(c12, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(c12, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(c12, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;

  CrossSectionSED c20(2, 0);
  cout << "Cross-section spectrum, HeII, type 0:\n";
  cout << "  monochromatic = " << c20.monochromatic() << endl;
  cout << "  lower_bound = " << c20.lower_bound() << endl;
  cout << "  upper_bound = " << c20.upper_bound() << endl;
  SED_integral(c20, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(c20, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(c20, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;

  CrossSectionSED c21(2, 1);
  cout << "Cross-section spectrum, HeII, type 1:\n";
  cout << "  monochromatic = " << c21.monochromatic() << endl;
  cout << "  lower_bound = " << c21.lower_bound() << endl;
  cout << "  upper_bound = " << c21.upper_bound() << endl;
  SED_integral(c21, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(c21, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(c21, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;

  CrossSectionSED c22(2, 2);
  cout << "Cross-section spectrum, HeII, type 2:\n";
  cout << "  monochromatic = " << c22.monochromatic() << endl;
  cout << "  lower_bound = " << c22.lower_bound() << endl;
  cout << "  upper_bound = " << c22.upper_bound() << endl;
  SED_integral(c22, 1.0, -1.0, true, R);
  cout << "  integral [1,infty] = " << R << endl;
  SED_integral(c22, 1.0, 1.e5, true, R);
  cout << "  integral [1,1e5] = " << R << endl;
  SED_integral(c22, 20.0, 100.0, true, R);
  cout << "  integral [20,100] = " << R << endl << endl;

  
  return 0;
} // end main

