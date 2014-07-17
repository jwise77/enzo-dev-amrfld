/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, radiation field access routine
/
/  written by: Daniel Reynolds
/  date:       July 2014
/  modified1:  
/
/  PURPOSE: Accesses a specific radiation field out of a Grid object.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"



float* AMRFLDSplit::AccessRadiationField(int Field, HierarchyEntry *ThisGrid)
{
//   if (debug)
//     printf("Entering AMRFLDSplit::AccessRadiationField routine\n");

  // access desired radiation field
  float *Efield = NULL;
  switch (Field) {
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
    ENZO_FAIL("AMRFLDSplit::AccessRadiationField error: illegal Field value");
  }

  // return success
  return Efield;

}
#endif
