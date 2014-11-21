/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class -- SetupSystem routine
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Fills HYPRE matrix and right-hand side objects corresponding
/           to a unigrid FLD radiation linear system of equations.
/              d_t E - Div(D(E)*Grad(E)) = -adot/a*E - c*kappa*E + eta + src
/           where D(E) is a nonlinear flux-limiter depending on E0
/           (time lagged).  We pose the linear system in predictor-corrector 
/           form, so that a solution of zero corresponds to a correct 
/           predictor.  This helps with both (a) determination of whether 
/           a solve is actually necessary, and (b) enforcement of boundary 
/           conditions, since they may be directly placed onto the current 
/           state, and the correction need only refrain from interfering.
/
/           In addition to returning the SUCCESS/FAIL flag, we also fill 
/           the rhsnorm value, corresponding to the RMS norm of the 
/           right-hand side vector.
/
/  NOTE: This routine is only called when the AMRFLDSplit solver is run
/        on a unigrid problem, in order to bypass the full AMRsolve 
/        infrastructure (that does not scale well at present).
/
************************************************************************/
#ifdef TRANSFER
#include "AMRFLDSplit.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif


int AMRFLDSplit::SetupSystem(HierarchyEntry *ThisGrid, 
			     int ibin, 
			     int SolvIndices[3][2],
			     HYPRE_StructMatrix &P, 
			     HYPRE_StructVector &rhsvec, 
			     HYPRE_StructVector &solvec, 
			     float &rhsnorm) {

  // check for valid bin index
  if (ibin<0 || ibin>=NumRadiationFields) {
    ENZO_VFAIL("AMRFLDSplit::SetupSystem: Error in ibin = %"ISYM"\n", ibin)
  }

  // set reusable constants
  float rUn  = ErUnits[ibin];
  float rUn0 = ErUnits0[ibin];
  rhsnorm = 0.0;
  float dtfac  = dt*theta;       // time step conversion factor
  float dtfac0 = dt*(1.0-theta); // time step conversion factor
  float afac  = 0.0;
  float afac0 = 0.0;
  if (!FieldMonochromatic[ibin]) {
    afac  = adot/a;      // expansion factor (new time)
    afac0 = adot0/a0;    // expansion factor (old time)
  }

  // set stencil indices based on dimensionality
  int stsize, s_xl, s_xr, s_yl, s_yr, s_zl, s_zr, s_c;
  stsize = 0; s_xl = s_xr = s_yl = s_yr = s_zl = s_zr = s_c = -1;
  if (rank > 2)  s_zl = stsize++;
  if (rank > 1)  s_yl = stsize++;
  s_xl = stsize++;
  s_c = stsize++;
  s_xr = stsize++;
  if (rank > 1)  s_yr = stsize++;
  if (rank > 2)  s_zr = stsize++;
  Eint32 entries[] = {0, 1, 2, 3, 4, 5, 6, 7};


  // set grid dimension information
  int dim, i;
  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
  int ghXl = DEFAULT_GHOST_ZONES;
  int n3[] = {1, 1, 1};
  for (dim=0; dim<rank; dim++)
    n3[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
            - ThisGrid->GridData->GetGridStartIndex(dim) + 1;
  int x0len = n3[0] + 2*ghXl;
  int x1len = n3[1] + 2*ghYl;
  int x2len = n3[2] + 2*ghZl;
  int size = x0len*x1len*x2len;

  // set dimension information
  int dxinv[] = {1.0, 1.0, 1.0};
  for (int dim=0; dim<rank; dim++) 
    dxinv[dim] = n3[dim] / ( ThisGrid->GridData->GetGridRightEdge(dim) -
			     ThisGrid->GridData->GetGridLeftEdge(dim) );

  // access relevant radiation, opacity and emissivity fields
  float *E       = AccessRadiationField(ibin, ThisGrid);
  float *Eta     = AccessEmissivityField(ibin, ThisGrid);
  float *Opacity = ThisGrid->GridData->AccessPhotoGamma();  // temporarily stored here
  if (E == NULL)
    ENZO_FAIL("AMRFLDSplit::SetupSystem ERROR: no radiation array!");
  if (Eta == NULL)
    ENZO_FAIL("AMRFLDSplit::SetupSystem ERROR: no emissivity array!");
  if (Opacity == NULL)
    ENZO_FAIL("AMRFLDSplit::SetupSystem ERROR: no opacity array!");
  
  
  ////////////////////////////
  // compute matrix entries, place within HYPRE objects
  
  // initialize constants and reusable local variables
  float dxi  = dxinv[0]/LenUnits;
  float dyi  = dxinv[1]/LenUnits;
  float dzi  = dxinv[2]/LenUnits;
  float dxi0 = dxinv[0]/LenUnits0;
  float dyi0 = dxinv[1]/LenUnits0;
  float dzi0 = dxinv[2]/LenUnits0;
  float dxfac = dtfac*dxi*dxi;
  float dyfac = dtfac*dyi*dyi;
  float dzfac = dtfac*dzi*dzi;
  float dxfac0 = dtfac0*dxi0*dxi0;
  float dyfac0 = dtfac0*dyi0*dyi0;
  float dzfac0 = dtfac0*dzi0*dzi0;

  // iterate over the domain, computing the RHS array
/*#pragma omp for reduction(+:rhsnorm) schedule(static) default(shared)*/
  for (int i2=0; i2<n3[2]; i2++) {
    int k2 = ghZl + i2;
    for (int i1=0; i1<n3[1]; i1++) { 
      int k1 = ghYl + i1;
      for (int i0=0; i0<n3[0]; i0++) {
	int k0 = ghXl + i0;
	
	// compute indices of neighboring Enzo cells
	int k_l00 = k0-1 + x0len*(k1   + x1len*k2);
	int k_0l0 = k0   + x0len*(k1-1 + x1len*k2);
	int k_00l = k0   + x0len*(k1   + x1len*(k2-1));
	int k_000 = k0   + x0len*(k1   + x1len*k2);
	int k_00r = k0   + x0len*(k1   + x1len*(k2+1));
	int k_0r0 = k0   + x0len*(k1+1 + x1len*k2);
	int k_r00 = k0+1 + x0len*(k1   + x1len*k2);
	    
	// declare loop-local variables
	float Ediff[7], D[7], D0[7], vals[7];
	    
	// z-directional limiter, lower face
	if (rank > 2) {
	  Ediff[s_zl] = E[k_000] - E[k_00l];
	  D0[s_zl] = Limiter(E[k_000], E[k_00l], Opacity[k_000], 
			     Opacity[k_00l], NiUnits0, LenUnits0, dzi0);
	  D[s_zl] = Limiter(E[k_000], E[k_00l], Opacity[k_000], 
			    Opacity[k_00l], NiUnits, LenUnits, dzi);
	}
	    
	// y-directional limiter, lower face
	if (rank > 1) {
	  Ediff[s_yl] = E[k_000] - E[k_0l0];
	  D0[s_yl] = Limiter(E[k_000], E[k_0l0], Opacity[k_000], 
			     Opacity[k_0l0], NiUnits0, LenUnits0, dyi0);
	  D[s_yl] = Limiter(E[k_000], E[k_0l0], Opacity[k_000], 
			    Opacity[k_0l0], NiUnits, LenUnits, dyi);
	}
	    
	// x-directional limiter, lower face
	Ediff[s_xl] = E[k_000] - E[k_l00];
	D0[s_xl] = Limiter(E[k_000], E[k_l00], Opacity[k_000], 
			   Opacity[k_l00], NiUnits0, LenUnits0, dxi0);
	D[s_xl] = Limiter(E[k_000], E[k_l00], Opacity[k_000], 
			  Opacity[k_l00], NiUnits, LenUnits, dxi);
	    
	// x-directional limiter, upper face
	Ediff[s_xr] = E[k_r00] - E[k_000];
	D0[s_xr] = Limiter(E[k_000], E[k_r00], Opacity[k_000], 
			   Opacity[k_r00], NiUnits0, LenUnits0, dxi0);
	D[s_xr] = Limiter(E[k_000], E[k_r00], Opacity[k_000], 
			  Opacity[k_r00], NiUnits, LenUnits, dxi);
	    
	// y-directional limiter, upper face
	if (rank > 1) {
	  Ediff[s_yr] = E[k_0r0] - E[k_000];
	  D0[s_yr] = Limiter(E[k_000], E[k_0r0], Opacity[k_000], 
			     Opacity[k_0r0], NiUnits0, LenUnits0, dyi0);
	  D[s_yr] = Limiter(E[k_000], E[k_0r0], Opacity[k_000], 
			    Opacity[k_0r0], NiUnits, LenUnits, dyi);
	}
	    
	// z-directional limiter, upper face
	if (rank > 2) {
	  Ediff[s_zr] = E[k_00r] - E[k_000];
	  D0[s_zr] = Limiter(E[k_000], E[k_00r], Opacity[k_000], 
			     Opacity[k_00r], NiUnits0, LenUnits0, dzi0);
	  D[s_zr] = Limiter(E[k_000], E[k_00r], Opacity[k_000], 
			    Opacity[k_00r], NiUnits, LenUnits, dzi);
	}
	    
	// opacity values in this cell
	float kap  = Opacity[k_000]*NiUnits;
	float kap0 = Opacity[k_000]*NiUnits0;

	// set the matrix entries.  Note: the diffusive component 
	// need not be rescaled, since scaling and chain rule cancel 
	if (rank > 2)  vals[s_zl] = -dzfac*D[s_zl];      // z-left
	if (rank > 1)  vals[s_yl] = -dyfac*D[s_yl];      // y-left
	vals[s_xl]                = -dxfac*D[s_xl];      // x-left
	vals[s_c] = 1.0 + dtfac*(afac + clight*kap)      // self
	  + dxfac*(D[s_xl]+D[s_xr])
	  + dyfac*(D[s_yl]+D[s_yr])
	  + dzfac*(D[s_zl]+D[s_zr]);
	vals[s_xr]                = -dxfac*D[s_xr];      // x-right
	if (rank > 1)  vals[s_yr] = -dyfac*D[s_yr];      // y-right
	if (rank > 2)  vals[s_zr] = -dzfac*D[s_zr];      // z-right

	// set the rhs entries
	float rhs = (dtfac/rUn + dtfac0/rUn0)*Eta[k_000]
	  - (dtfac*(afac+clight*kap) + dtfac0*(afac0+clight*kap0))*E[k_000]
	  + dxfac0*(D0[s_xr]*Ediff[s_xr] - D0[s_xl]*Ediff[s_xl])
	  + dxfac*(D[s_xr]*Ediff[s_xr] - D[s_xl]*Ediff[s_xl]);
	if (rank > 1)
	  rhs += (dyfac0*(D0[s_yr]*Ediff[s_yr] - D0[s_yl]*Ediff[s_yl]) +
		  dyfac*(D[s_yr]*Ediff[s_yr] - D[s_yl]*Ediff[s_yl]));
	if (rank > 2)
	  rhs += (dzfac0*(D0[5]*Ediff[s_zr] - D0[0]*Ediff[s_zl]) +
		  dzfac*(D[5]*Ediff[s_zr] - D[0]*Ediff[s_zl]));
	rhsnorm += (rhs*rhs);

	// check for legal matrix/rhs values, otherwise error with useful message
	for (int idx=0; idx<6; idx++)
	  if (isinf(vals[idx]) || isnan(vals[idx])) {
	    fprintf(stderr,"AMRFLDSplit::SetupSystem ERROR: illegal matrix value (%g)\n   "
		    "proc = %i, i* = %i %i %i %i, LocDims = %i %i %i, eta = %g, "
		    "E = %g %g %g %g %g %g %g, dtfac = %g, dtfac0 = %g, "
		    "kap = %g %g %g %g %g %g %g\n   D* = %g, %g, %g, %g, %g, %g\n   "
		    "Ed* = %g %g %g %g %g %g\n\n",
		    vals[idx], MyProcessorNumber, idx, i0, i1, i2, n3[0], n3[1], n3[2], 
		    Eta[k_000], E[k_000], E[k_00l], E[k_0l0], E[k_l00], E[k_r00], E[k_0r0], 
		    E[k_00r], dtfac, dtfac0, Opacity[k_000], Opacity[k_00l], 
		    Opacity[k_0l0], Opacity[k_l00], Opacity[k_r00], Opacity[k_0r0], 
		    Opacity[k_00r], D[0], D[1], D[2], D[3], D[4], D[5], Ediff[0], 
		    Ediff[1], Ediff[2], Ediff[3], Ediff[4], Ediff[5]);
	    ENZO_FAIL("NaN error in AMRFLDSplit::SetupSystem\n");
	  }
	if (isinf(rhs) || isnan(rhs)) {
	  fprintf(stderr,"AMRFLDSplit::SetupSystem ERROR: illegal vector value (%g)\n   "
		  "proc = %i, i* = %i %i %i, LocDims = %i %i %i, eta = %g, "
		  "E = %g %g %g %g %g %g %g, dtfac = %g, dtfac0 = %g, "
		  "kap = %g %g %g %g %g %g %g\n   D* = %g, %g, %g, %g, %g, %g\n   "
		  "Ed* = %g %g %g %g %g %g\n\n",
		  rhs, MyProcessorNumber, i0, i1, i2, n3[0], n3[1], n3[2], 
		  Eta[k_000], E[k_000], E[k_00l], E[k_0l0], E[k_l00], E[k_r00], E[k_0r0], 
		  E[k_00r], dtfac, dtfac0, Opacity[k_000], Opacity[k_00l], 
		  Opacity[k_0l0], Opacity[k_l00], Opacity[k_r00], Opacity[k_0r0], 
		  Opacity[k_00r], D[0], D[1], D[2], D[3], D[4], D[5], Ediff[0], 
		  Ediff[1], Ediff[2], Ediff[3], Ediff[4], Ediff[5]);
	  ENZO_FAIL("NaN error in AMRFLDSplit::SetupSystem\n");
	}
	    
	    
	// update the matrix entries based on boundary conditions
	//    z-left face
	if ((rank > 2) && (OnBdry[2][0]) && (i2==0)) {
	  if (BdryType[2][0] == 1)         // Dirichlet
	    vals[s_zl] = 0.0;
	  else if (BdryType[2][0] == 2) {  // Neumann
	    vals[s_c] += vals[s_zl];
	    vals[s_zl] = 0.0;
	  }
	}
	    
	//    y-left face
	if ((rank > 1) && (OnBdry[1][0]) && (i1==0)) {
	  if (BdryType[1][0] == 1)         // Dirichlet
	    vals[s_yl] = 0.0;
	  else if (BdryType[1][0] == 2) {  // Neumann
	    vals[s_c] += vals[s_yl];
	    vals[s_yl] = 0.0;
	  }
	}
	    
	//    x-left face
	if ((OnBdry[0][0]) && (i0==0)) {
	  if (BdryType[0][0] == 1)         // Dirichlet
	    vals[s_xl] = 0.0;
	  else if (BdryType[0][0] == 2) {  // Neumann
	    vals[s_c] += vals[s_xl];
	    vals[s_xl] = 0.0;
	  }
	}
	    
	//    x-right face
	if ((OnBdry[0][1]) && (i0==n3[0]-1)) {
	  if (BdryType[0][1] == 1)         // Dirichlet
	    vals[s_xr] = 0.0;
	  else if (BdryType[0][1] == 2) {  // Neumann
	    vals[s_c] += vals[s_xr];
	    vals[s_xr] = 0.0;
	  }
	}
	    
	//    y-right face
	if ((rank > 1) && (OnBdry[1][1]) && (i1==n3[1]-1)) {
	  if (BdryType[1][1] == 1)         // Dirichlet
	    vals[s_yr] = 0.0;
	  else if (BdryType[1][1] == 2) {  // Neumann
	    vals[s_c] += vals[s_yr];
	    vals[s_yr] = 0.0;
	  }
	}
	    
	//    z-right face
	if ((rank > 2) && (OnBdry[2][1]) && (i2==n3[2]-1)) {
	  if (BdryType[2][1] == 1)         // Dirichlet
	    vals[s_zr] = 0.0;
	  else if (BdryType[2][1] == 2) {  // Neumann
	    vals[s_c] += vals[s_zr];
	    vals[s_zr] = 0.0;
	  }
	}
	    
	    
	// insert matrix, rhs and initial solution entries into HYPRE objects
	Eint32 iloc[] = {SolvIndices[0][0]+i0, SolvIndices[1][0]+i1, SolvIndices[2][0]+i2};
	float zero=0.0;
	HYPRE_StructVectorSetBoxValues(rhsvec, iloc, iloc, &rhs);
	HYPRE_StructVectorSetBoxValues(solvec, iloc, iloc, &zero);
	HYPRE_StructMatrixSetBoxValues(P, iloc, iloc, stsize, entries, vals); 
	    
      }  // for i
    }  // for j
  }  // for k
      
  // assemble HYPRE matrix and vectors
  HYPRE_StructVectorAssemble(solvec);
  HYPRE_StructVectorAssemble(rhsvec);
  HYPRE_StructMatrixAssemble(P);
  
  // accumulate full rhsnorm
  // combine the processor-local rhsnorm values together before returning 
  float rhssum=0.0;
#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Allreduce(&rhsnorm, &rhssum, 1, DataType, MPI_SUM, MPI_COMM_WORLD);
#else
  rhssum = rhsnorm;
#endif
  rhsnorm = sqrt(rhssum);


  return SUCCESS;
}

#endif  /* TRANSFER */
