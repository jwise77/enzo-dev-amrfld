/*****************************************************************************
 *                                                                           *
 * Copyright 2014 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Single-Group, Multi-species, AMR, Gray Flux-Limited Diffusion 
/  Split Implicit Problem Class, EnforceBoundary routine
/
/  written by: Daniel Reynolds
/  date:       May 2014
/  modified1:  
/
/  PURPOSE: Accesses a specific radiation bin out of a Grid object.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"



float* AMRFLDSplit::AccessRadiationBin(int Bin, HierarchyEntry *ThisGrid)
{
//   if (debug)
//     printf("Entering AMRFLDSplit::AccessRadiationBin routine\n");

  // access old/new radiation fields (old stored in KPhHI)
  float *Efield = NULL;
  switch (Bin) {
  case 0:
    Efield = ThisGrid->GridData->AccessRadiationFrequency0();
    break;
  case 1:
    Efield = ThisGrid->GridData->AccessRadiationFrequency1();
    break;
  case 2:
    Efield = ThisGrid->GridData->AccessRadiationFrequency2();
    break;
  case 3:
    Efield = ThisGrid->GridData->AccessRadiationFrequency3();
    break;
  case 4:
    Efield = ThisGrid->GridData->AccessRadiationFrequency4();
    break;
  case 5:
    Efield = ThisGrid->GridData->AccessRadiationFrequency5();
    break;
  case 6:
    Efield = ThisGrid->GridData->AccessRadiationFrequency6();
    break;
  case 7:
    Efield = ThisGrid->GridData->AccessRadiationFrequency7();
    break;
  case 8:
    Efield = ThisGrid->GridData->AccessRadiationFrequency8();
    break;
  case 9:
    Efield = ThisGrid->GridData->AccessRadiationFrequency9();
    break;
  default:
    ENZO_FAIL("AMRFLDSplit::AccessRadiationBin error: illegal Bin value");
  }

  // return success
  return Efield;

}
#endif
