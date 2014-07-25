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
#include "SED.h"
#include "BlackbodySED.h"
#include "MonochromaticSED.h"


// function prototypes
float SED_integral(SED &sed, float a, float b, bool convertHz);


// SED "moment" class for computing radiation field contributions
class SED_Moment : public virtual SED {
 private:
  SED *baseSED;     // base SED
 public:
  SED_Moment(SED &base) { this->baseSED = &base; };
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

    // skip sources with SourceType of -1 (indicates alternate source strategy)
    if (SourceType[isrc] == -1)  continue;

    // set up source SED
    SED *tmp_src;
    switch (SourceType[isrc]) {
    case 0:      // monochromatic source at hnu = 13.6 eV
      tmp_src = new MonochromaticSED(13.6);
      break;
    case 1:      // blackbody spectrum at T = 1e5 K
      tmp_src = new BlackbodySED(1.0e5);
      break;
    default:
      ENZO_VFAIL("AMRFLDSplit::SetupSources error, source type %"ISYM" does not exist\n", 
		 SourceType[isrc]);
    }

    // set up source SED moment
    SED_Moment tmp_src_m(*tmp_src);

    // compute integral of SED over all ionizing frequencies
    float total_integral = SED_integral(*tmp_src, hnu0, -1.0, true);

    // iterate over radiation fields
    for (ibin=0; ibin<NumRadiationFields; ibin++) {

      // skip monochromatic fields
      if (FieldMonochromatic[ibin]) {
	SourceGroupEnergy[isrc][ibin] = 0.0;
	continue;
      }

      // computing source contributions to this radiation field
      SourceGroupEnergy[isrc][ibin] = SourceEnergy[isrc] * hplanck / total_integral
	* SED_integral(tmp_src_m, FrequencyBand[ibin][0], FrequencyBand[ibin][1], true);
    }

    // report results to stdout
    if (debug) {
      printf("AMRFLDSplit::Initialize source %"ISYM" has group energies", isrc);
      for (ibin=0; ibin<NumRadiationFields; ibin++)  
	printf(" %g", SourceGroupEnergy[isrc][ibin]);
      printf("\n");
    }

    // clean up
    delete tmp_src;

  }   // end isrc loop

  return SUCCESS;
}
#endif // TRANSFER
