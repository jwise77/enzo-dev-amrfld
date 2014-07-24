/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, Unit Update Routine
/
/  written by: Daniel Reynolds
/  date:       July 2014
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
#ifdef TRANSFER

#include "AMRFLDSplit.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"



// function prototypes
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);



// This routine updates the internal AMRFLDSplit unit values for a given time.
int AMRFLDSplit::UpdateUnits(FLOAT Told, FLOAT Tnew) {

  // local variables
  int ibin;
  float TempUnits, RadUnits;
  double MassUnits;

  // update internal units for old time
  DenUnits0 = LenUnits0 = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits0, &LenUnits0, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, Told) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_UpdateUnits: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, Told) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_UpdateUnits: Error in RadiationGetUnits.");
  for (ibin=0; ibin<NumRadiationFields; ibin++)
    ErUnits0[ibin] = RadUnits*ErScale[ibin];
  NiUnits0 = DenUnits0/mh;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(Told, &a0, &adot0) != SUCCESS) 
      ENZO_FAIL("AMRFLDSplit_UpdateUnits: Error in CosmologyComputeExpansionFactor.");
  adot0 /= TimeUnits;    // rescale to physical units

  // update internal units for new time
  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, Tnew) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_UpdateUnits: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, Tnew) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_UpdateUnits: Error in RadiationGetUnits.");
  for (ibin=0; ibin<NumRadiationFields; ibin++)
    ErUnits[ibin] = RadUnits*ErScale[ibin];
  NiUnits = DenUnits/mh;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(Tnew, &a, &adot) != SUCCESS) 
      ENZO_FAIL("AMRFLDSplit_UpdateUnits: Error in CosmologyComputeExpansionFactor.");
  adot /= TimeUnits;     // rescale to physical units

  return SUCCESS;
}

#endif   // TRANSFER
