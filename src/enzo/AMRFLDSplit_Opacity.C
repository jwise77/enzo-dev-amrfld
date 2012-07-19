/*****************************************************************************
 *                                                                           *
 * Copyright 2012 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Single-Group, Multi-species, AMR, Gray Flux-Limited Diffusion 
/  Split Implicit Problem Class, opacity field calculation routine 
/
/  written by: Daniel Reynolds
/  date:       July 2012
/  modified1:  
/
/  PURPOSE: Computes the temperature, density and chemistry-dependent 
/  opacity field throughout the domain, storing it in the photo-heating
/  baryon field (since that is not used by the code during this phase 
/  of the calculation).
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"

 
 

int AMRFLDSplit::Opacity(LevelHierarchyEntry *LevelArray[], int level, float time)
{

  // initialize local variables to be reused
  int i, j, k;
  float HIconst   = intSigESigHI   / intSigE;
  float HeIconst  = intSigESigHeI  / intSigE / 4.0;
  float HeIIconst = intSigESigHeII / intSigE / 4.0;
  
  // iterate over grids owned by this processor (this level down)
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridData->ReturnProcessorNumber()) {
	
	// set dimension information
	int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
	int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
	int ghXl = DEFAULT_GHOST_ZONES;
	int n3[] = {1, 1, 1};
	for (int dim=0; dim<rank; dim++)
	  n3[dim] = Temp->GridData->GetGridEndIndex(dim)
  	          - Temp->GridData->GetGridStartIndex(dim) + 1;
	int x0len = n3[0] + 2*ghXl;
	int x1len = n3[1] + 2*ghYl;
	int x2len = n3[2] + 2*ghZl;

	// access chemistry fields
	float *HI=NULL, *HeI=NULL, *HeII=NULL;
	HI   = Temp->GridData->AccessHIDensity();
	HeI  = Temp->GridData->AccessHeIDensity();
	HeII = Temp->GridData->AccessHeIIDensity();

	// access photo-heating field to store opacity
	float *kap=NULL;
	kap = Temp->GridData->AccessPhotoGamma();

	// check that required field data exists
	if (HI==NULL)
	  ENZO_FAIL("AMRFLDSplit_Opacity ERROR: no HI array!");
	if (Nchem==3 && HeI==NULL)
	  ENZO_FAIL("AMRFLDSplit_Opacity ERROR: no HeI array!");
	if (Nchem==3 && HeII==NULL)
	  ENZO_FAIL("AMRFLDSplit_Opacity ERROR: no HeII array!");
	if (kap==NULL)
	  ENZO_FAIL("AMRFLDSplit_Opacity ERROR: no PhotoGamma array!");


	/////////////
	// compute opacity field based on ProblemType
	switch (ProblemType) {


	// Insert user-defined opacity fields here (don'te forget to "break")...


	// Standard chemistry-dependent opacity field
	default:

	  // Hydrogen-only calculation
	  if (Nchem == 1)
	    for (i=0; i<x0len*x1len*x2len; i++) 
	      kap[i] = HI[i]*HIconst;

	  // Hydrogen + Helium calculation
	  if (Nchem == 3)
	    for (i=0; i<x0len*x1len*x2len; i++) 
	      kap[i] = HI[i]*HIconst + HeI[i]*HeIconst + HeII[i]*HeIIconst;
	  
	  break;

	}
	
      }  // end iteration over grids on this processor


  return SUCCESS;
}

#endif
