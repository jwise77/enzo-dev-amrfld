/*****************************************************************************
 *                                                                           *
 * Copyright 2010 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Single-Group, Multi-species, AMR, Gray Flux-Limited Diffusion 
/  Split Implicit Problem Class, Radiation spectrum evaluation routine 
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:  
/
/  PURPOSE: Takes in a frequency and returns the assumed radiation 
/  energy spectrum at that frequency.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"

 
 

float AMRFLDSplit::RadiationSpectrum(int Bin, float nu)
{

  // set necessary constants
  float h = 6.6260693e-27;          // Planck's constant [ergs*s]
  float kb = 1.3806504e-16;         // Boltzmann's constant [ergs/K]
  float pi = 4.0*atan(1.0);         // pi
  float c = 2.99792458e10;          // speed of light [cm/s]
  float ev2erg = 1.60217653e-12;    // conversion constant from eV to ergs
  float nu0 = hnu0_HI*ev2erg/h;     // ionization threshold of Hydrogen (hz)
  float nu1 = 2.5*nu0;              // ionization of Wolf Reyet stars + HeliumI
  float nu2 = 4.0*nu0;              // ionization threshold of HeliumII
  float nu3 = 100.0*nu0;            // parameter used to characterize PopII SED cutoff frequency
  float GBbrk = 2.5;                // parameters in Ricotti 2002 fig.4 fit
  float GXray = 2.0e-3;             // parameters in Ricotti 2002 fig.4 fit
  float sigma;

  // check that frequency is within the allowed range
  if (nu < nu0) {
    fprintf(stderr,"AMRFLDSplit::RadiationSpectrum Error: (nu = %g) < (nu0 = %g)\n",
	    nu,nu0);
    return -1.0;
  }

  // evaluate the radiation spectrum based on the internal ESpectrum parameter
  switch (ESpectrum[Bin]) {

  // T = 1e5 K blackbody spectrum
  case 1:
    sigma = 8.0*pi*h*POW(nu/c,3)/(exp(h*nu/kb/1e5)-1.0);
    break;

  // SED with photons above 4 Ryd truncated
  case 2:
    if (nu < nu1)
      sigma = 1.0/nu0/POW(nu/nu0, 1.8);
    else if (nu < nu2)
      sigma = 1.0/nu0/GBbrk/POW(nu/nu0, 1.8);
    else if (nu >= nu2)
      sigma = 0.0;
    break;

  // // T = 1e5 K blackbody spectrum, bin 0: nu0 <= nu <= 5e15
  // case 3:
  //   if ((nu >= nu0) && (nu <= 5.e15))
  //     sigma = 8.0*pi*h*POW(nu/c,3)/(exp(h*nu/kb/1e5)-1.0);
  //   break;

  // T = 1e5 K blackbody spectrum, bin 0: 0 < nu <= 5e15
  case 3:
    if (nu <= 5.e15)
      sigma = 8.0*pi*h*POW(nu/c,3)/(exp(h*nu/kb/1e5)-1.0);
    break;

  // T = 1e5 K blackbody spectrum, bin 2: 5e15 < nu <= 1e16
  case 4:
    if ((nu > 5.e15) && (nu <= 1.e16))
      sigma = 8.0*pi*h*POW(nu/c,3)/(exp(h*nu/kb/1e5)-1.0);
    break;

  // T = 1e5 K blackbody spectrum, bin 3: 1e16 < nu
  case 5:
    if (nu > 1.e16)
      sigma = 8.0*pi*h*POW(nu/c,3)/(exp(h*nu/kb/1e5)-1.0);
    break;

  // Add new spectrum choices here
  // case 6:
  //   ...
  //   break;

  default:
    // simple power law spectrum with power -1.5
    sigma = POW((nu/nu0),-1.5);
  }

  return sigma;
}

#endif
