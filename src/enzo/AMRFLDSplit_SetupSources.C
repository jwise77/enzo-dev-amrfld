/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, SetupSources routine
/
/  written by: Daniel Reynolds
/  date:       July 2014
/
/  PURPOSE: Handles user inputs to set up relevant desired emissivity
/           sources for each radiation frequency/group.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
#include "BlackbodySED.h"
#include "MonochromaticSED.h"


// function prototypes
float SED_integral(SED &sed, float a, float b, bool convertHz);


// SED type for computing radiation energy contribution integrals
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



int AMRFLDSplit::SetupSources(HierarchyEntry *RootGrid,
			      int SourceType[MAX_FLD_SOURCES], 
			      float SourceEnergy[MAX_FLD_SOURCES]) {

  if (debug)  printf("Entering AMRFLDSplit::SetupSources routine\n");

  // if this is a weak scaling test, overwrite SourceLocation for each 
  // source to replicate source setup on each root grid tile
  int isrc, ibin, dim;
  if (WeakScaling)
    for (isrc=0; isrc<NumSources; isrc++) 
      for (dim=0; dim<rank; dim++) {
	float frac = (SourceLocation[isrc][dim] - DomainLeftEdge[dim]) 
	           / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	OriginalSourceLocation[isrc][dim] = SourceLocation[isrc][dim];
	SourceLocation[isrc][dim] = RootGrid->GridData->GetGridLeftEdge(dim) +
                                    frac*(RootGrid->GridData->GetGridRightEdge(dim) -
					  RootGrid->GridData->GetGridLeftEdge(dim));
      }


  // set up some constants
  float hnu0 = 13.6;                      // HI ionization threshold [eV]
  float nu0 = hnu0 * ev2erg / hplanck;    // [Hz]
  

  // fill SourceGroupEnergy for sources that are specified by type/energy
  for (isrc=0; isrc<NumSources; isrc++) {

    // monochromatic source at hnu = 13.6 eV
    if (SourceType[isrc] == 0) {
      MonochromaticSED tmp_src(13.6);
      for (ibin=0; ibin<NumRadiationFields; ibin++) {
	if (FieldMonochromatic[ibin]) {   // skip monochromatic fields
	  SourceGroupEnergy[isrc][ibin] = 0.0;
	  continue;
	}
	SourceGroupEnergy[isrc][ibin] = SourceEnergy[isrc] * hplanck * nu0
	  * SED_integral(tmp_src, FrequencyBand[ibin][0], FrequencyBand[ibin][1], true);
      }
      if (debug) {
	printf("AMRFLDSplit::Initialize monochromatic source %"ISYM" has group energies", isrc);
	for (ibin=0; ibin<NumRadiationFields; ibin++)  printf(" %g", SourceGroupEnergy[isrc][ibin]);
	printf("\n");
      }
    }

    // blackbody spectrum at T = 1e5 K
    if (SourceType[isrc] == 1) {
      BlackbodySED bSED(1.0e5);
      Blackbody_Moment bSED_M(bSED);
      float total_integral = SED_integral(bSED, 13.6, -1.0, true);
      for (ibin=0; ibin<NumRadiationFields; ibin++) {
	if (FieldMonochromatic[ibin]) {   // skip monochromatic fields
	  SourceGroupEnergy[isrc][ibin] = 0.0;
	  continue;
	}
	SourceGroupEnergy[isrc][ibin] = SourceEnergy[isrc] * hplanck / total_integral
	  * SED_integral(bSED_M, FrequencyBand[ibin][0], FrequencyBand[ibin][1], true);
      }
      if (debug) {
	printf("AMRFLDSplit::Initialize blackbody source %"ISYM" has group energies", isrc);
	for (ibin=0; ibin<NumRadiationFields; ibin++)  printf(" %g", SourceGroupEnergy[isrc][ibin]);
	printf("\n");
      }
    }

  }   // end isrc loop

  return SUCCESS;
}
#endif // TRANSFER
