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
/  Split Implicit Problem Class, AccessEmissivityBin routine
/
/  written by: Daniel Reynolds
/  date:       May 2014
/  modified1:  
/
/  PURPOSE: Accesses a specific emissivity bin out of a Grid object.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"



float* AMRFLDSplit::AccessEmissivityBin(int Bin, HierarchyEntry *ThisGrid)
{
//   if (debug)
//     printf("Entering AMRFLDSplit::AccessRadiationBin routine\n");

  // access emissivity fields
  float *Efield = NULL;
  switch (Bin) {
  case 0:
    Efield = ThisGrid->GridData->AccessEmissivity0();
    break;
  case 1:
    Efield = ThisGrid->GridData->AccessEmissivity1();
    break;
  case 2:
    Efield = ThisGrid->GridData->AccessEmissivity2();
    break;
  case 3:
    Efield = ThisGrid->GridData->AccessEmissivity3();
    break;
  case 4:
    Efield = ThisGrid->GridData->AccessEmissivity4();
    break;
  case 5:
    Efield = ThisGrid->GridData->AccessEmissivity5();
    break;
  case 6:
    Efield = ThisGrid->GridData->AccessEmissivity6();
    break;
  case 7:
    Efield = ThisGrid->GridData->AccessEmissivity7();
    break;
  case 8:
    Efield = ThisGrid->GridData->AccessEmissivity8();
    break;
  case 9:
    Efield = ThisGrid->GridData->AccessEmissivity9();
    break;
  default:
    ENZO_FAIL("AMRFLDSplit::AccessEmissivityBin error: illegal Bin value");
  }

  // return success
  return Efield;

}
#endif
