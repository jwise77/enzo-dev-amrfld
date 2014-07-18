/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, EnforceBoundary routine
/
/  written by: Daniel Reynolds
/  date:       December 2010
/  modified1:  
/
/  PURPOSE: Enforces boundary conditions on the radiation field (top 
/           grid only).  Assumes that the stored radiation units are 
/           updated appropriately (see 2nd note below).
/
/           Note: Neumann values are enforced on the first 
/                 layer of ghost zones using a first-order central 
/                 difference approximation to the first (outward-normal) 
/                 derivative.
/           Note: Since the internal radiation variables are comoving 
/                 and normalized, we renormalize the boundary conditions 
/                 as they are enforced to match the internal units.
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"



int AMRFLDSplit::EnforceBoundary(int ibin, LevelHierarchyEntry *LevelArray[])
{
//   if (debug)
//     printf("Entering AMRFLDSplit::EnforceBoundary routine\n");

  // find this processor's TopGrid
  LevelHierarchyEntry* Temp;
  HierarchyEntry* ThisGrid;
  for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel) {
    ThisGrid = Temp->GridHierarchyEntry;
    if (MyProcessorNumber == ThisGrid->GridData->ReturnProcessorNumber())
      break;
  }

  // set grid information
  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
  int ghXl = DEFAULT_GHOST_ZONES;
  int n3[] = {1, 1, 1};
  for (int dim=0; dim<rank; dim++)
    n3[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
            - ThisGrid->GridData->GetGridStartIndex(0) + 1;
  int x0len = n3[0] + 2*ghXl;
  int x1len = n3[1] + 2*ghYl;
  int x2len = n3[2] + 2*ghZl;
  float dx[3];
  for (int dim=0; dim<rank; dim++)
    dx[dim] = (ThisGrid->GridData->GetGridRightEdge(dim) 
	     - ThisGrid->GridData->GetGridLeftEdge(dim)) 
             / n3[dim];
      
  // access current radiation field
  float *Enew = AccessRadiationField(ibin, ThisGrid);
  
  // set some shortcuts
  float dxa, dya, dza, unit;
  int i, i2, j, j2, k, k2, idx, idx2, idxbc;
  dxa = dx[0] * LenUnits;
  
  // x0 left boundary
  //   Dirichlet
  if (OnBdry[0][0] && (BdryType[0][0]==1)) {
    for (k=0; k<n3[2]; k++)
      for (j=0; j<n3[1]; j++)
	for (i=0; i<ghXl; i++) {
	  idxbc = k*n3[1] + j;
	  idx = ((k+ghZl)*x1len + j+ghYl)*x0len + i;
	  Enew[idx] = BdryVals[ibin][0][0][idxbc]/ErUnits[ibin];
	}
  }
  //   Neumann
  if (OnBdry[0][0] && (BdryType[0][0]==2)) {
    i = -1;  i2 = i+1;
    for (k=0; k<n3[2]; k++)
      for (j=0; j<n3[1]; j++) {
	idx = ((k+ghZl)*x1len + j+ghYl)*x0len + i+ghXl;
	idx2 = ((k+ghZl)*x1len + j+ghYl)*x0len + i2+ghXl;
	idxbc = k*n3[1] + j;
	Enew[idx] = Enew[idx2] + dxa*BdryVals[ibin][0][0][idxbc]/ErUnits[ibin];
      }
  }
    
  // x0 right boundary
  //   Dirichlet
  if (OnBdry[0][1] && (BdryType[0][1]==1)) {
    for (k=0; k<n3[2]; k++)
      for (j=0; j<n3[1]; j++)
	for (i=x0len-ghXl; i<x0len; i++) {
	  idxbc = k*n3[1] + j;
	  idx = ((k+ghZl)*x1len + j+ghYl)*x0len + i;
	  Enew[idx] = BdryVals[ibin][0][1][idxbc]/ErUnits[ibin];
	}
  }
  //   Neumann
  if (OnBdry[0][1] && (BdryType[0][1]==2)) {
    i = n3[0];  i2 = i-1;
    for (k=0; k<n3[2]; k++)
      for (j=0; j<n3[1]; j++) {
	idx = ((k+ghZl)*x1len + j+ghYl)*x0len + i+ghXl;
	idx2 = ((k+ghZl)*x1len + j+ghYl)*x0len + i2+ghXl;
	idxbc = k*n3[1] + j;
	Enew[idx] = Enew[idx2] + dxa*BdryVals[ibin][0][1][idxbc]/ErUnits[ibin];
      }
  }
    
  if (rank > 1) {
    dya = dx[1]*LenUnits;
    // x1 left boundary
    //   Dirichlet
    if (OnBdry[1][0] && (BdryType[1][0]==1)) {
      for (k=0; k<n3[2]; k++)
	for (j=0; j<ghYl; j++)
	  for (i=0; i<n3[0]; i++) {
	    idx = ((k+ghZl)*x1len + j)*x0len + i+ghXl;
	    idxbc = i*n3[2] + k;
	    Enew[idx] = BdryVals[ibin][1][0][idxbc]/ErUnits[ibin];
	  }
    }
    //   Neumann
    if (OnBdry[1][0] && (BdryType[1][0]==2)) {
      j = -1;  j2 = j+1;
      for (k=0; k<n3[2]; k++)
	for (i=0; i<n3[0]; i++) {
	  idx = ((k+ghZl)*x1len + j+ghYl)*x0len + i+ghXl;
	  idx2 = ((k+ghZl)*x1len + j2+ghYl)*x0len + i+ghXl;
	  idxbc = i*n3[2] + k;
	  Enew[idx] = Enew[idx2] + dya*BdryVals[ibin][1][0][idxbc]/ErUnits[ibin];
	}
    }
      
    // x1 right boundary
    //   Dirichlet
    if (OnBdry[1][1] && (BdryType[1][1]==1)) {
      for (k=0; k<n3[2]; k++)
	for (j=x1len-ghYl; j<x1len; j++)
	  for (i=0; i<n3[0]; i++) {
	    idx = ((k+ghZl)*x1len + j)*x0len + i+ghXl;
	    idxbc = i*n3[2] + k;
	    Enew[idx] = BdryVals[ibin][1][1][idxbc]/ErUnits[ibin];
	  }
    }
    //   Neumann
    if (OnBdry[1][1] && (BdryType[1][1]==2)) {
      j = n3[1];  j2 = j-1;
      for (k=0; k<n3[2]; k++)
	for (i=0; i<n3[0]; i++) {
	  idx = ((k+ghZl)*x1len + j+ghYl)*x0len + i+ghXl;
	  idx2 = ((k+ghZl)*x1len + j2+ghYl)*x0len + i+ghXl;
	  idxbc = i*n3[2] + k;
	  Enew[idx] = Enew[idx2] + dya*BdryVals[ibin][1][1][idxbc]/ErUnits[ibin];
	}
    }
  }  // end if rank > 1
    
  if (rank > 2) {
    dza = dx[2]*LenUnits;
    // x2 left boundary
    //   Dirichlet
    if (OnBdry[2][0] && (BdryType[2][0]==1)) {
      for (k=0; k<ghZl; k++)
	for (j=0; j<n3[1]; j++)
	  for (i=0; i<n3[0]; i++) {
	    idx = (k*x1len + j+ghYl)*x0len + i+ghXl;
	    idxbc = j*n3[0] + i;
	    Enew[idx] = BdryVals[ibin][2][0][idxbc]/ErUnits[ibin];
	  }
    }
    //   Neumann
    if (OnBdry[2][0] && (BdryType[2][0]==2)) {
      k = -1;  k2 = k+1;
      for (j=0; j<n3[1]; j++)
	for (i=0; i<n3[0]; i++) {
	  idx = ((k+ghZl)*x1len + j+ghYl)*x0len + i+ghXl;
	  idx2 = ((k2+ghZl)*x1len + j+ghYl)*x0len + i+ghXl;
	  idxbc = j*n3[0] + i;
	  Enew[idx] = Enew[idx2] + dza*BdryVals[ibin][2][0][idxbc]/ErUnits[ibin];
	}
    }
      
    // x2 right boundary
    //   Dirichlet
    if (OnBdry[2][1] && (BdryType[2][1]==1)) {
      for (k=x2len-ghZl; k<x2len; k++)
	for (j=0; j<n3[1]; j++)
	  for (i=0; i<n3[0]; i++) {
	    idx = (k*x1len + j+ghYl)*x0len + i+ghXl;
	    idxbc = j*n3[0] + i;
	    Enew[idx] = BdryVals[ibin][2][1][idxbc]/ErUnits[ibin];
	  }
    }
    //   Neumann
    if (OnBdry[2][1] && (BdryType[2][1]==2)) {
      k = n3[2];  k2 = k-1;
      for (j=0; j<n3[1]; j++)
	for (i=0; i<n3[0]; i++) {
	  idx = ((k+ghZl)*x1len + j+ghYl)*x0len + i+ghXl;
	  idx2 = ((k2+ghZl)*x1len + j+ghYl)*x0len + i+ghXl;
	  idxbc = j*n3[0] + i;
	  Enew[idx] = Enew[idx2] + dza*BdryVals[ibin][2][1][idxbc]/ErUnits[ibin];
	}
    }
  }  // end if rank > 2

//   if (debug)
//     printf("Exiting AMRFLDSplit::EnforceBoundary routine\n");

  // return success
  return SUCCESS;

}
#endif
