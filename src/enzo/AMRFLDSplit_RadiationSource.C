/*****************************************************************************
 *                                                                           *
 * Copyright 2010 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, Emissivity Field Computation Routine
/
/  written by: Daniel Reynolds
/  date:       January 2011
/  modified1:  
/
/  PURPOSE: Computes the emissivity field to enforce on the radiation 
/           energy equation.  This is only called if StarMaker sources 
/           are not handled elsewhere (i.e. for test problems).
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
#include "phys_constants.h"
#include "BlackbodySED.h"
 

// function prototypes
int SED_integral(SED &sed, float a, float b, bool convertHz, double &R);
 


float AMRFLDSplit::RadiationSource(int ibin, LevelHierarchyEntry *LevelArray[], 
				   int level, float time)
{

  // initialize local variables to be reused
  int i, j, k;
  float total_eta = 0.0;
  float cellZl, cellZr, cellYl, cellYr, cellXl, cellXr, cellXc, cellYc, cellZc;

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
	float x0L = Temp->GridData->GetGridLeftEdge(0);
	float x1L = Temp->GridData->GetGridLeftEdge(1);
	float x2L = Temp->GridData->GetGridLeftEdge(2);
	float dx[3] = {1, 1, 1};
	for (int dim=0; dim<rank; dim++)
	  dx[dim] = (Temp->GridData->GetGridRightEdge(dim) 
		     - Temp->GridData->GetGridLeftEdge(dim)) 
	          / n3[dim];
	float lUn = (LenUnits + LenUnits0)*0.5;
	float dV = dx[0]*dx[1]*dx[2]*lUn*lUn*lUn;

	// set a cell "normalized volume"
	float dVscale = 1;
	for (int dim=0; dim<rank; dim++)
	  dVscale *= (Temp->GridData->GetGridRightEdge(dim) 
 		    - Temp->GridData->GetGridLeftEdge(dim)) 
	            / n3[dim] / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	
	// access emissivity field
	float *eta = AccessEmissivityField(ibin, Temp->GridHierarchyEntry);
	if (eta == NULL)
	  ENZO_VFAIL("AMRFLDSplit_RadiationSource error: emissivity field %"ISYM" missing\n",
		     ibin);

	// iterate over sources, adding emissivity to this bin at specified location
	for (int isrc=0; isrc<NumSources; isrc++) {

	  // all sources have radius of one cell; count number of cells to receive source
	  int num_cells = 0;
	  for (k=ghZl; k<n3[2]+ghZl; k++) {
	    cellZc = x2L + (k-ghZl+0.5)*dx[2];	      // z-center (comoving) for this cell
	    for (j=ghYl; j<n3[1]+ghYl; j++) {
	      cellYc = x1L + (j-ghYl+0.5)*dx[1];      // y-center (comoving) for this cell
	      for (i=ghXl; i<n3[0]+ghXl; i++) {
		cellXc = x0L + (i-ghXl+0.5)*dx[0];    // x-center (comoving) for this cell
		if ( (fabs(cellXc-SourceLocation[ibin][0]) < dx[0]) &&
		     (fabs(cellYc-SourceLocation[ibin][1]) < dx[1]) &&
		     (fabs(cellZc-SourceLocation[ibin][2]) < dx[2]) )
		  num_cells++;                        // cell is within source region
	      } // x-loop
	    } // y-loop
	  } // z-loop

	  // equi-partition energy among affected cells
	  float cell_energy = SourceGroupEnergy[isrc][ibin] / num_cells / dV;
	  for (k=ghZl; k<n3[2]+ghZl; k++) {
	    cellZc = x2L + (k-ghZl+0.5)*dx[2];	      // z-center (comoving) for this cell
	    for (j=ghYl; j<n3[1]+ghYl; j++) {
	      cellYc = x1L + (j-ghYl+0.5)*dx[1];      // y-center (comoving) for this cell
	      for (i=ghXl; i<n3[0]+ghXl; i++) {
		cellXc = x0L + (i-ghXl+0.5)*dx[0];    // x-center (comoving) for this cell
		if ( (fabs(cellXc-SourceLocation[ibin][0]) < dx[0]) &&
		     (fabs(cellYc-SourceLocation[ibin][1]) < dx[1]) &&
		     (fabs(cellZc-SourceLocation[ibin][2]) < dx[2]) )
		  eta[(k*x1len + j)*x0len + i] += cell_energy;
	      } // x-loop
	    } // y-loop
	  } // z-loop

	} // end loop over sources from input parameters

	
	/////////////
	// add emissivity based on ProblemType, if desired
	
	switch (ProblemType) {

	case 412:    // emissivity flux along x=0 wall (NGammaDot photons/s/cm^2)
	  
	  // only relevant for radiation groups (not individual frequencies)
	  if (!FieldMonochromatic[ibin]) {

	    // place ionization sources along left wall (if on this subdomain)
	    if (x0L == 0.0) {

	      // determine this group's portion of total blackbody emissivity
	      BlackbodySED tmp_src(1.0e5);
	      float total_integral;
	      if (SED_integral(tmp_src, 13.6, -1.0, true, total_integral) != SUCCESS) {
		ENZO_FAIL("ERROR in integrating the overall blackbody SED\n");
	      }
	      float integral;
	      if (SED_integral(tmp_src, FrequencyBand[ibin][0], FrequencyBand[ibin][1], 
			       true, integral) != SUCCESS) {
		ENZO_VFAIL("ERROR in integrating the blackbody SED for group %"ISYM"\n", ibin);
	      }
	      float wall_energy = 1e6 * 13.6 * ev2erg * integral / total_integral / dx[0] / lUn;

	      // place along wall (i=ghXl)
	      for (k=ghZl; k<n3[2]+ghZl; k++)
		for (j=ghYl; j<n3[1]+ghYl; j++)
		  eta[(k*x1len + j)*x0len + ghXl] = wall_energy;
	    }
	  } 
	  break;
		
	  
	}
	
	// accumulate L2 norm of total active emissivity on this grid
	for (k=ghZl; k<n3[2]+ghZl; k++) 
	  for (j=ghYl; j<n3[1]+ghYl; j++) 
	    for (i=ghXl; i<n3[0]+ghXl; i++) 
 	      total_eta += eta[(k*x1len+j)*x0len+i]*eta[(k*x1len+j)*x0len+i]*dVscale;
	
      }  // end iteration over grids on this processor
 
  // communicate to obtain overall emissivity
  float glob_eta = 0.0;
#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg ONE = 1;
  MPI_Allreduce(&total_eta,&glob_eta,ONE,DataType,MPI_SUM,MPI_COMM_WORLD);
#else
  glob_eta = total_eta;
#endif
  
  return sqrt(glob_eta);

}

#endif
