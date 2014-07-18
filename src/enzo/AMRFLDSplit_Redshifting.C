/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, Redshifting Routine
/
/  written by: Daniel Reynolds
/  date:       July 2014
/  modified1:  
/
/  PURPOSE: Computes the redshifting fluxes between neighboring groups 
/  within a cosmology simulation, and adds these fluxes to the current 
/  emissivity fields for each group.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
#include "phys_constants.h"
 

int AMRFLDSplit::Redshifting(LevelHierarchyEntry *LevelArray[], int level)
{

  // if not running a cosmology simulation, just return
  if (!ComovingCoordinates)  return SUCCESS;

  // iterate over grids owned by this processor (this level down)
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridData->ReturnProcessorNumber()) {
	
	// set dimension/grid information
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

	// iterate over all radiation fields
	for (int ibin=0; ibin<NumRadiationFields; ibin++) {

	  // skip monochromatic fields
	  if (FieldMonochromatic[ibin])  continue;

	  // access emissivity field to place redshifting contributions
	  float *eta = AccessEmissivityField(ibin, Temp->GridHierarchyEntry);
	  if (eta == NULL) {
	    ENZO_VFAIL("AMRFLDSplit_Redshifting error: emissivity field %"ISYM" missing\n", ibin);
	  }

	  // access relevant radiation fields
	  float *E=NULL, *E_l=NULL, *E_r=NULL;
	  E = AccessRadiationField(ibin, Temp->GridHierarchyEntry);
	  if (E == NULL) {
	    ENZO_VFAIL("AMRFLDSplit_Redshifting error: radiation field %"ISYM" missing\n", ibin);
	  }
	  if (FieldNeighbors[ibin][0]) {
	    E_l = AccessRadiationField(ibin-1, Temp->GridHierarchyEntry);
	    if (E_l == NULL) {
	      ENZO_VFAIL("AMRFLDSplit_Redshifting error: radiation field %"ISYM" missing\n", ibin-1);
	    }
	  }
	  if (FieldNeighbors[ibin][1]) {
	    E_r = AccessRadiationField(ibin+1, Temp->GridHierarchyEntry);
	    if (E_r == NULL) {
	      ENZO_VFAIL("AMRFLDSplit_Redshifting error: radiation field %"ISYM" missing\n", ibin+1);
	    }
	  }

	  // set bin widths
	  float bwid, bwid_l, bwid_r;
	  bwid = FrequencyBand[ibin][1] - FrequencyBand[ibin][0];
	  if (FieldNeighbors[ibin][0]) 
	    bwid_l = FrequencyBand[ibin-1][1] - FrequencyBand[ibin-1][0];
	  if (FieldNeighbors[ibin][1]) 
	    bwid_r = FrequencyBand[ibin+1][1] - FrequencyBand[ibin+1][0];

	  // time-centered unit/expansion factors
	  float adot_c  = 0.5*(adot+adot0);
	  float a_c     = 0.5*(a+a0);
	  float units_c = 0.5*(ErUnits[ibin]/TimeUnits + ErUnits0[ibin]/TimeUnits);

	  // set reusable constants
	  float rconst=0.0, rconst_l=0.0, rconst_r=0.0;
	  if (FieldNeighbors[ibin][0]) {
	    rconst   -= 0.5*adot_c*FrequencyBand[ibin][0]/a_c/bwid*units_c;
	    rconst_l -= 0.5*adot_c*FrequencyBand[ibin][0]/a_c/bwid_l*units_c;
	  } else {
	    rconst -= adot_c*FrequencyBand[ibin][0]/a_c/bwid*units_c;
	  }
	  if (FieldNeighbors[ibin][1]) {
	    rconst   += 0.5*adot_c*FrequencyBand[ibin][1]/a_c/bwid*units_c;
	    rconst_r += 0.5*adot_c*FrequencyBand[ibin][1]/a_c/bwid_r*units_c;
	  } else {
	    rconst += adot_c*FrequencyBand[ibin][1]/a_c/bwid*units_c;
	  }
	  rconst -= adot_c/a_c*units_c;
	  
	  // add redshifting source terms for this equation to current 
	  // emissivity, based on neighbor types
	  //    field has neighbors on both sides
	  if (FieldNeighbors[ibin][0] && FieldNeighbors[ibin][1]) {
	    for (int k=ghZl; k<n3[2]+ghZl; k++) 
	      for (int j=ghYl; j<n3[1]+ghYl; j++) 
		for (int i=ghXl; i<n3[0]+ghXl; i++) 
		  eta[(k*x1len + j)*x0len + i] += rconst   * E[(k*x1len + j)*x0len + i]
		                                + rconst_l * E_l[(k*x1len + j)*x0len + i]
                                                + rconst_r * E_r[(k*x1len + j)*x0len + i];

	  //    field has left neighbor only
	  } else if (FieldNeighbors[ibin][0] && !FieldNeighbors[ibin][1]) {
	    for (int k=ghZl; k<n3[2]+ghZl; k++) 
	      for (int j=ghYl; j<n3[1]+ghYl; j++) 
		for (int i=ghXl; i<n3[0]+ghXl; i++) 
		  eta[(k*x1len + j)*x0len + i] += rconst   * E[(k*x1len + j)*x0len + i]
		                                + rconst_l * E_l[(k*x1len + j)*x0len + i];

	  //    field has right neighbor only
	  } else if (!FieldNeighbors[ibin][0] && FieldNeighbors[ibin][1]) {
	    for (int k=ghZl; k<n3[2]+ghZl; k++) 
	      for (int j=ghYl; j<n3[1]+ghYl; j++) 
		for (int i=ghXl; i<n3[0]+ghXl; i++) 
		  eta[(k*x1len + j)*x0len + i] += rconst   * E[(k*x1len + j)*x0len + i]
                                                + rconst_r * E_r[(k*x1len + j)*x0len + i];

	  //    field has no neighbors
	  } else {
	    for (int k=ghZl; k<n3[2]+ghZl; k++) 
	      for (int j=ghYl; j<n3[1]+ghYl; j++) 
		for (int i=ghXl; i<n3[0]+ghXl; i++) 
		  eta[(k*x1len + j)*x0len + i] += rconst * E[(k*x1len + j)*x0len + i];

	  }

	} // end loop over radiation fields
	
      }  // end iteration over grids on this processor

  return SUCCESS;  
}

#endif
