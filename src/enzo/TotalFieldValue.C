/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Accessibility routine for debugging: accumulates the total integrated 
/  amount of a given field.
/
/  written by: Daniel Reynolds
/  date:       July 2014
/  modified1:
/
************************************************************************/
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"




// This routine computes the total integrated value of a given field over
// the Enzo hierarchy from levelStart down through levelFinish.  Here, 
// "field" may be any allowable fieldtype (see typedefs.h)
float TotalFieldValue(LevelHierarchyEntry *LevelArray[], 
		      int levelStart, int levelFinish, int field)
{

  // initialize output
  float result = 0.0;

  // ensure that levelFinish >= levelStart, otherwise return
  if (levelStart > levelFinish)  return result;

  // ensure that levelStart >= 0, otherwise return
  if (levelStart < 0) return result;

  // iterate over all grids on this processor from levelStart through levelFinish
  for (int thislevel=levelStart; thislevel<levelFinish; thislevel++)
    for (LevelHierarchyEntry *Temp=LevelArray[thislevel]; Temp; Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {
	
	// access field data
	float *data = NULL;
	switch (field) {
	case Density:     //  Density 
	  data = Temp->GridHierarchyEntry->GridData->AccessDensity();
	  break;
	case TotalEnergy:     //  Total Energy 
	  data = Temp->GridHierarchyEntry->GridData->AccessTotalEnergy();
	  break;
	case InternalEnergy:     //  Gas Energy 
	  data = Temp->GridHierarchyEntry->GridData->AccessGasEnergy();
	  break;
	case Velocity1:     //  Velocity1 
	  data = Temp->GridHierarchyEntry->GridData->AccessVelocity1();
	  break;
	case Velocity2:     //  Velocity2
	  data = Temp->GridHierarchyEntry->GridData->AccessVelocity2();
	  break;
	case Velocity3:     //  Velocity3
	  data = Temp->GridHierarchyEntry->GridData->AccessVelocity3();
	  break;
	case ElectronDensity:     //  Electron Density
	  data = Temp->GridHierarchyEntry->GridData->AccessElectronDensity();
	  break;
	case HIDensity:     //  Hydrogen-I Density
	  data = Temp->GridHierarchyEntry->GridData->AccessHIDensity();
	  break;
	case HIIDensity:     //  Hydrogen-II Density 
	  data = Temp->GridHierarchyEntry->GridData->AccessHIIDensity();
	  break;
	case HeIDensity:     //  Helium-I Density 
	  data = Temp->GridHierarchyEntry->GridData->AccessHeIDensity();
	  break;
	case HeIIDensity:     //  Helium-II Density 
	  data = Temp->GridHierarchyEntry->GridData->AccessHeIIDensity();
	  break;
	case HeIIIDensity:     //  Helium-III Density
	  data = Temp->GridHierarchyEntry->GridData->AccessHeIIIDensity();
	  break;
	case HMDensity:     //  Molecular Hydrogen Density
	  data = Temp->GridHierarchyEntry->GridData->AccessHMDensity();
	  break;
	case H2IDensity:     //  H2I Density 
	  data = Temp->GridHierarchyEntry->GridData->AccessH2IDensity();
	  break;
	case H2IIDensity:     //  H2II Density 
	  data = Temp->GridHierarchyEntry->GridData->AccessH2IIDensity();
	  break;
	case DIDensity:     //  DI Density 
	  data = Temp->GridHierarchyEntry->GridData->AccessDIDensity();
	  break;
	case DIIDensity:     //  DII Density 
	  data = Temp->GridHierarchyEntry->GridData->AccessDIIDensity();
	  break;
	case HDIDensity:     //  HDI Density 
	  data = Temp->GridHierarchyEntry->GridData->AccessHDIDensity();
	  break;
	case SNColour:     //  SN Colour 
	  data = Temp->GridHierarchyEntry->GridData->AccessSNColour();
	  break;
	case Metallicity:     //  Metallicity 
	  data = Temp->GridHierarchyEntry->GridData->AccessMetallicity();
	  break;
	case ExtraType0:     //  ExtraType0 
	  data = Temp->GridHierarchyEntry->GridData->AccessExtraType0();
	  break;
	case ExtraType1:     //  ExtraType1 
	  data = Temp->GridHierarchyEntry->GridData->AccessExtraType1();
	  break;
	case kphHI:     //  kphHI 
	  data = Temp->GridHierarchyEntry->GridData->AccessKPhHI();
	  break;
	case PhotoGamma:     //  PhotoGamma 
	  data = Temp->GridHierarchyEntry->GridData->AccessPhotoGamma();
	  break;
	case kphHeI:     //  kphHeI 
	  data = Temp->GridHierarchyEntry->GridData->AccessKPhHeI();
	  break;
	case gammaHeI:     //  Gamma HeI 
	  data = Temp->GridHierarchyEntry->GridData->AccessGammaHeI();
	  break;
	case kphHeII:     //  kphHeII 
	  data = Temp->GridHierarchyEntry->GridData->AccessKPhHeII();
	  break;
	case gammaHeII:     //  Gamma HeII 
	  data = Temp->GridHierarchyEntry->GridData->AccessGammaHeII();
	  break;
	case kdissH2I:     //  kdissH2I 
	  data = Temp->GridHierarchyEntry->GridData->AccessKDissH2I();
	  break;
	case GravPotential:     //  Gravitational Potential 
	  data = Temp->GridHierarchyEntry->GridData->AccessGravPotential();
	  break;
	case Acceleration0:     //  Acceleration0 
	  data = Temp->GridHierarchyEntry->GridData->AccessAcceleration0();
	  break;
	case Acceleration1:     //  Acceleration1 
	  data = Temp->GridHierarchyEntry->GridData->AccessAcceleration1();
	  break;
	case Acceleration2:     //  Acceleration2 Density 
	  data = Temp->GridHierarchyEntry->GridData->AccessAcceleration2();
	  break;
	case RadPressure0:     //  Radiation Pressure0 
	  data = Temp->GridHierarchyEntry->GridData->AccessRadPressure0();
	  break;
	case RadPressure1:     //  Radiation Pressure1 
	  data = Temp->GridHierarchyEntry->GridData->AccessRadPressure1();
	  break;
	case RadPressure2:     //  Radiation Pressure2 
	  data = Temp->GridHierarchyEntry->GridData->AccessRadPressure2();
	  break;
	case RadiationFreq0:     //  RadiationFreq0
	  data = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency0();
	  break;
	case RadiationFreq1:     //  RadiationFreq1
	  data = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency1();
	  break;
	case RadiationFreq2:     //  RadiationFreq2
	  data = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency2();
	  break;
	case RadiationFreq3:     //  RadiationFreq3
	  data = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency3();
	  break;
	case RadiationFreq4:     //  RadiationFreq4
	  data = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency4();
	  break;
	case RadiationFreq5:     //  RadiationFreq5
	  data = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency5();
	  break;
	case RadiationFreq6:     //  RadiationFreq6
	  data = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency6();
	  break;
	case RadiationFreq7:     //  RadiationFreq7
	  data = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency7();
	  break;
	case RadiationFreq8:     //  RadiationFreq8
	  data = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency8();
	  break;
	case RadiationFreq9:     //  RadiationFreq9
	  data = Temp->GridHierarchyEntry->GridData->AccessRadiationFrequency9();
	  break;
	case Emissivity0:     //  Emissivity0 
	  data = Temp->GridHierarchyEntry->GridData->AccessEmissivity0();
	  break;
	case Emissivity1:     //  Emissivity1
	  data = Temp->GridHierarchyEntry->GridData->AccessEmissivity1();
	  break;
	case Emissivity2:     //  Emissivity2 
	  data = Temp->GridHierarchyEntry->GridData->AccessEmissivity2();
	  break;
	case Emissivity3:     //  Emissivity3 
	  data = Temp->GridHierarchyEntry->GridData->AccessEmissivity3();
	  break;
	case Emissivity4:     //  Emissivity4 
	  data = Temp->GridHierarchyEntry->GridData->AccessEmissivity4();
	  break;
	case Emissivity5:     //  Emissivity5 
	  data = Temp->GridHierarchyEntry->GridData->AccessEmissivity5();
	  break;
	case Emissivity6:     //  Emissivity6 
	  data = Temp->GridHierarchyEntry->GridData->AccessEmissivity6();
	  break;
	case Emissivity7:     //  Emissivity7 
	  data = Temp->GridHierarchyEntry->GridData->AccessEmissivity7();
	  break;
	case Emissivity8:     //  Emissivity8 
	  data = Temp->GridHierarchyEntry->GridData->AccessEmissivity8();
	  break;
	case Emissivity9:     //  Emissivity9 
	  data = Temp->GridHierarchyEntry->GridData->AccessEmissivity9();
	  break;
	default:
	  data = NULL;
	}

	// skip this grid if data not present
	if (data == NULL) continue;

	// set grid dimension information
	int rank = Temp->GridHierarchyEntry->GridData->GetGridRank();
	int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	int ghXl = DEFAULT_GHOST_ZONES;
	int n3[] = {1, 1, 1};
	for (int dim=0; dim<rank; dim++)
	  n3[dim] = Temp->GridHierarchyEntry->GridData->GetGridEndIndex(dim)
	          - Temp->GridHierarchyEntry->GridData->GetGridStartIndex(dim) + 1;
	int x0len = n3[0] + 2*ghXl;
	int x1len = n3[1] + 2*ghYl;
	int x2len = n3[2] + 2*ghZl;

	// compute volume scaling factor for this grid/level
	float dVscale = 1;
	for (int dim=0; dim<rank; dim++)
	  dVscale *= (Temp->GridData->GetGridRightEdge(dim) 
		    - Temp->GridData->GetGridLeftEdge(dim)) / n3[dim]
		   / (DomainRightEdge[dim] - DomainLeftEdge[dim]);

	// accumulate result
	for (int k=ghZl; k<n3[2]+ghZl; k++) 
	  for (int j=ghYl; j<n3[1]+ghYl; j++) 
	    for (int i=ghXl; i<n3[0]+ghXl; i++) 
	      result += data[(k*x1len+j)*x0len+i]*dVscale;
	
      }  // end loop over grids on this proc
  
  // communicate among all procs to get statistics of current field
  float glob_result;
#ifdef USE_MPI
  MPI_Datatype FDataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Allreduce(&result, &glob_result, 1, FDataType, MPI_SUM, MPI_COMM_WORLD);
#else
  glob_result = result;
#endif
  
  // return with result
  return glob_result;
}

