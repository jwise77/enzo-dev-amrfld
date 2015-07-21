/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, RadStepUnigrid routine
/
/  written by: Daniel Reynolds
/  date:       November 2014
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
#include "phys_constants.h"

//#define FAIL_ON_NAN
#define NO_FAIL_ON_NAN

// function prototypes
float TotalFieldValue(LevelHierarchyEntry *LevelArray[], 
		      int levelStart, int levelFinish, int field);



// This routine evolves the radiation subsystem within the AMRFLDSplit module.
// This is performed in a robust manner; if the current radiation subsystem 
// cannot be solved to the desired tolerance (meaning the time step is 
// likely much too large), it will return a flag to the calling routine to 
// have the step recomputed with a smaller time step size.
//
// This is the unigrid version of the routine; it does not rely on anything 
// from the AMRsolve infrastructure 
int AMRFLDSplit::RadStepUnigrid(int ibin, LevelHierarchyEntry *LevelArray[], Eflt64 &Echange)
{

#ifdef USE_HYPRE

  // initialize return value
  int recompute_step = 0;

  // find this MPI task's root-grid tile
  HierarchyEntry* ThisGrid=NULL;
  for (LevelHierarchyEntry* Temp=LevelArray[0]; Temp; Temp=Temp->NextGridThisLevel)
    if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {
      ThisGrid = Temp->GridHierarchyEntry;
      break;
    }
  if (ThisGrid == NULL) {
    ENZO_VFAIL("AMRFLDSplit::RadStepUnigrid: MPI task %i failed to find root-grid tile",
	       MyProcessorNumber);
  }

  // set grid dimension information
  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
  int ghXl = DEFAULT_GHOST_ZONES;
  int n3[] = {1, 1, 1};
  for (int dim=0; dim<rank; dim++)
    n3[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
            - ThisGrid->GridData->GetGridStartIndex(dim) + 1;
  int x0len = n3[0] + 2*ghXl;
  int x1len = n3[1] + 2*ghYl;
  int x2len = n3[2] + 2*ghZl;
 
  // update internal units for current times
  if (this->UpdateUnits(told, tnew) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit::RadStepUnigrid: Error in UpdateUnits.");

  // rescale dt, told, tnew to physical values for use within solver
  dt   *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;

  // compute current opacity for this bin
  if (this->Opacity(ibin, LevelArray, 0, tnew) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit::RadStepUnigrid: Opacity failure!!");
    
  // enforce boundary conditions on new radiation field
  if (this->EnforceBoundary(ibin, LevelArray) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit::RadStepUnigrid: EnforceBoundary failure!!");
    
  // copy current radiation field into temporary field (KPhHI)
  float *Eold = ThisGrid->GridData->AccessKPhHI();
  float *Enew = AccessRadiationField(ibin, ThisGrid);
  if (Eold == NULL)
    ENZO_FAIL("AMRFLDSplit::RadStepUnigrid: cannot access KPhHI field!!\n");
  if (Enew == NULL) {
    ENZO_VFAIL("AMRFLDSplit::RadStepUnigrid: cannot access radiation field %"ISYM"!!\n", ibin);
  }
  for (int k=0; k<x0len*x1len*x2len; k++)  Eold[k] = Enew[k];

  // calculate emissivity source norm
  float glob_eta = TotalFieldValue(LevelArray, 0, MAX_DEPTH_OF_HIERARCHY, Emissivity0+ibin);
  if (debug)   printf("   total emissivity %"ISYM" = %g\n", ibin, glob_eta);


  ////////////////// set up the radiation system //////////////////

  // allocate HYPRE structures
  HYPRE_StructGrid grid;         // HYPRE grid object for setup
  HYPRE_StructStencil stencil;   // stencil object
  HYPRE_StructMatrix P;          // holds radiation matrix
  HYPRE_StructVector rhsvec;     // holds radiation rhs vector
  HYPRE_StructVector solvec;     // holds radiation solution vector
  Eflt64 *HYPREbuff;             // holds contiguous sections of rhs/sol

  // set up the grid
  //    create the grid object
  HYPRE_StructGridCreate(MPI_COMM_WORLD, rank, &grid);

  // compute global index information for this subdomain
  int SolvIndices[3][2] = {{0, 0}, {0, 0}, {0, 0}};    // global indices of this root grid
  int GlobDims[3] = {1, 1, 1};                         // global domain grid size
  float fCellsLeft;
  float dV = 1.0;                                      // volume of each cell
  for (int dim=0; dim<rank; dim++) {

    // grid mesh size in this direction
    float dx = (ThisGrid->GridData->GetGridRightEdge(dim) -
		ThisGrid->GridData->GetGridLeftEdge(dim)) / n3[dim];
    dV *= dx;

    // the global indexing is easy if we're at the left edge
    if (ThisGrid->GridData->GetProcessorLocation(dim)==0)  
      SolvIndices[dim][0]=0;

    // otherwise we compute the number of intervening cells to left edge
    else {

      // get floating point value for number of cells
      float fCellsLeft = (ThisGrid->GridData->GetGridLeftEdge(dim) - 
			  DomainLeftEdge[dim]) / dx;

      // round floating point value to closest integer
      SolvIndices[dim][0] =  (long) (fCellsLeft >= 0.0) ?
	(trunc(fCellsLeft+0.5)) : (trunc(fCellsLeft-0.5));
    }

    // add on local subdomain size to obtain right edge indices
    SolvIndices[dim][1] = SolvIndices[dim][0] + n3[dim] - 1;

    // compute global domain size in this direction
    float fCells = (DomainRightEdge[dim] - DomainLeftEdge[dim]) / dx;
    GlobDims[dim] = (long) (trunc(fCells+0.5));
  }

  //    set my grid extents as if we have one part with multiple boxes.
  //    Have each processor describe it's own global extents
  Eint32 ilower[3] = {SolvIndices[0][0], SolvIndices[1][0], SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1], SolvIndices[1][1], SolvIndices[2][1]};
  HYPRE_StructGridSetExtents(grid, ilower, iupper);

  //    set grid periodicity
  Eint32 periodicity[3] = {0, 0, 0};
  if (BdryType[0][0] == 0)  periodicity[0] = GlobDims[0];
  if (BdryType[1][0] == 0)  periodicity[1] = GlobDims[1];
  if (BdryType[2][0] == 0)  periodicity[2] = GlobDims[2];
  HYPRE_StructGridSetPeriodic(grid, periodicity);
  
  //    assemble the grid
  HYPRE_StructGridAssemble(grid);

  // set up the stencil
  int stSize;
  if (rank == 1) 
    stSize = 3;
  else if (rank == 2)
    stSize = 5;
  else 
    stSize = 7;
  HYPRE_StructStencilCreate(rank, stSize, &stencil);

  //   set stencil entries
  Eint32 offset[3];
  Eint32 stentry=0;
  //      dependency to x2 left
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //      dependency to x1 left
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //      dependency to x0 left
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //      dependency to self
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //      dependency to x0 right
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //      dependency to x1 right
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //      dependency to x2 right
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }

  // allocate temporary HYPRE matrix and vectors
  HYPREbuff = new Eflt64[n3[0]];
  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &P);
  HYPRE_StructMatrixInitialize(P);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &rhsvec);
  HYPRE_StructVectorInitialize(rhsvec);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &solvec);
  HYPRE_StructVectorInitialize(solvec);


  // fill in the radiation system (assembles matrix and vectors)
  float rhsnorm;
  if (this->SetupSystem(ThisGrid, ibin, SolvIndices, P, rhsvec, solvec, rhsnorm) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit::RadStepUnigrid: SetupSystem Error!!");


  // skip solve if ||rhs|| is small enough  (i.e. old solution is fine)
  if (rhsnorm < 0.01*sol_tolerance) {

    if (debug)
      printf("   no solve required: |rhs| = %.1e  <  tol = %.1e\n", rhsnorm, sol_tolerance);

    // set small error value, rescale dt, told, tnew, adot back to normalized values
    Echange = 1e-16;
    dt   /= TimeUnits;
    told /= TimeUnits;
    tnew /= TimeUnits;

    // destroy HYPRE structures
    HYPRE_StructVectorDestroy(rhsvec);
    HYPRE_StructVectorDestroy(solvec);
    HYPRE_StructMatrixDestroy(P);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructGridDestroy(grid);

    // return with a successful step
    return 0;
  }

  // if we've made it here, turn on automatic scaling for next time step
  StartAutoScale[ibin] = true;


  ////////////////// solve this radiation field's equation //////////////////

  // set linear solver tolerance (rescale to relative residual and not actual)
  bool use_abs_resid = (rhsnorm < 1.e-8);
  Eflt64 delta = (use_abs_resid) ? sol_tolerance : sol_tolerance/rhsnorm;
  delta = min(delta, 1.0e-2);

  // set up the solver and preconditioner [PFMG]
  //    create the solver & preconditioner
  HYPRE_StructSolver solver;            // HYPRE solver structure
  HYPRE_StructSolver preconditioner;    // HYPRE preconditioner structure
  switch (sol_type) {
  case 4:   // PCG
    HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
    break;
  case 3:   // GMRES
    HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &solver);
    break;
  default:  // BiCGStab
    HYPRE_StructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
    break;
  }
  HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &preconditioner);
  
  // Multigrid solver: for periodic dims, only coarsen until grid no longer divisible by 2
  Eint32 max_levels, level=-1;
  int Ndir;
  if (BdryType[0][0] == 0) {
    level = 0;
    Ndir = GlobDims[0];
    while ( Ndir%2 == 0 ) {
      level++;
      Ndir /= 2;
    }
  }
  max_levels = level;
  if (rank > 1) {
    if (BdryType[1][0] == 0) {
      level = 0;
      Ndir = GlobDims[1];
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = min(level,max_levels);
  }
  if (rank > 2) {
    if (BdryType[2][0] == 0) {
      level = 0;
      Ndir = GlobDims[2];
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = min(level,max_levels);
  }

  //    set preconditioner options
  if (max_levels > -1) 
    HYPRE_StructPFMGSetMaxLevels(preconditioner, max_levels);
  HYPRE_StructPFMGSetMaxIter(preconditioner, sol_precmaxit);
  HYPRE_StructPFMGSetRelaxType(preconditioner, sol_rlxtype);
  HYPRE_StructPFMGSetNumPreRelax(preconditioner, sol_npre);
  HYPRE_StructPFMGSetNumPostRelax(preconditioner, sol_npost);
  
  //    set solver options
  switch (sol_type) {
  case 4:   // PCG
    HYPRE_StructPCGSetPrintLevel(solver, sol_printl);
    HYPRE_StructPCGSetLogging(solver, sol_log);
    HYPRE_StructPCGSetRelChange(solver, 1);
    if (rank > 1) {
      HYPRE_StructPCGSetMaxIter(solver, sol_maxit);
      HYPRE_StructPCGSetPrecond(solver, 
				(HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
				(HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
				preconditioner);
    } else {    // ignore preconditioner for 1D tests (bug); increase CG its
      HYPRE_StructPCGSetMaxIter(solver, sol_maxit*500);
    }
    if (use_abs_resid) {
      HYPRE_StructPCGSetTol(solver, 0.0);
      HYPRE_StructPCGSetAbsoluteTol(solver, delta);
    } else {
      HYPRE_StructPCGSetTol(solver, delta);
    }
    HYPRE_StructPCGSetup(solver, P, rhsvec, solvec);
    break;
  case 3:   // GMRES
    //  HYPRE_StructGMRESSetPrintLevel(solver, sol_printl);
    HYPRE_StructGMRESSetLogging(solver, sol_log);
    //  HYPRE_StructGMRESSetRelChange(solver, 1);
    if (rank > 1) {
      HYPRE_StructGMRESSetMaxIter(solver, sol_maxit);
      HYPRE_StructGMRESSetKDim(solver, sol_maxit);
      HYPRE_StructGMRESSetPrecond(solver, 
				  (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
				  (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
				  preconditioner);
    } else {    // ignore preconditioner for 1D tests (bug); increase its
      HYPRE_StructGMRESSetMaxIter(solver, sol_maxit*50);
      HYPRE_StructGMRESSetKDim(solver, sol_maxit*50);
    }
    if (use_abs_resid) {
      HYPRE_StructGMRESSetTol(solver, 0.0);
      HYPRE_StructGMRESSetAbsoluteTol(solver, delta);
    } else {
      HYPRE_StructGMRESSetTol(solver, delta);
    }
    HYPRE_StructGMRESSetup(solver, P, rhsvec, solvec);
    break;
  default:  // BiCGStab
    //  HYPRE_StructBiCGSTABSetPrintLevel(solver, sol_printl);
    HYPRE_StructBiCGSTABSetLogging(solver, sol_log);
    if (rank > 1) {
      HYPRE_StructBiCGSTABSetMaxIter(solver, sol_maxit);
      HYPRE_StructBiCGSTABSetPrecond(solver, 
				     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
				     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
				     preconditioner);
    } else {    // ignore preconditioner for 1D tests (bug); increase its
      HYPRE_StructBiCGSTABSetMaxIter(solver, sol_maxit*500);
    }
    if (use_abs_resid) {
      HYPRE_StructBiCGSTABSetTol(solver, 0.0);
      HYPRE_StructBiCGSTABSetAbsoluteTol(solver, delta);
    } else {
      HYPRE_StructBiCGSTABSetTol(solver, delta);
    }
    HYPRE_StructBiCGSTABSetup(solver, P, rhsvec, solvec);
    break;
  }
  
  // solve the linear system
  switch (sol_type) {
  case 4:   // PCG
    HYPRE_StructPCGSolve(solver, P, rhsvec, solvec);
    break;
  case 3:   // GMRES
    HYPRE_StructGMRESSolve(solver, P, rhsvec, solvec);
    break;
  default:  // BiCGStab
    HYPRE_StructBiCGSTABSolve(solver, P, rhsvec, solvec);
    break;
  }
  
  // extract solver & preconditioner statistics
  Eflt64 finalresid=1.0;  // HYPRE solver statistics
  Eint32 Sits=0, Pits=0;  // HYPRE solver statistics
  switch (sol_type) {
  case 4:   // PCG
    if (use_abs_resid)
      finalresid = 0.0;
    else
      HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
    HYPRE_StructPCGGetNumIterations(solver, &Sits);
    break;
  case 3:   // GMRES
    if (use_abs_resid)
      finalresid = 0.0;
    else
      HYPRE_StructGMRESGetFinalRelativeResidualNorm(solver, &finalresid);
    HYPRE_StructGMRESGetNumIterations(solver, &Sits);
    break;
  default:  // BiCGStab
    if (use_abs_resid)
      finalresid = 0.0;
    else
      HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(solver, &finalresid);
    HYPRE_StructBiCGSTABGetNumIterations(solver, &Sits);
    break;
  }
  HYPRE_StructPFMGGetNumIterations(preconditioner, &Pits);
  *totIters += Sits;
  if (debug) printf("   lin resid = %.1e (tol = %.1e, |rhs| = %.1e), its = (%i,%i)\n",
		    finalresid*rhsnorm, delta, rhsnorm, Sits, Pits);

  // check if residual is NaN
  if (finalresid != finalresid) {

#ifndef FAIL_ON_NAN
    if (dt > mindt*TimeUnits*1.00000001) {
      // allow remainder of function to complete (to reset units, etc.), 
      // but have calling routine update dt and compute step again.
      recompute_step = 1;
    } else {
#endif
      fprintf(stderr,"AMRFLDSplit::RadStepUnigrid: could not solve problem at minimum step size!\n");
      
      // output linear system to disk
      if (debug)  printf("Writing out matrix to file P.mat\n");
      HYPRE_StructMatrixPrint("P.mat",P,0);
      if (debug)  printf("Writing out rhs to file b.vec\n");
      HYPRE_StructVectorPrint("b.vec",rhsvec,0);
      if (debug)  printf("Writing out current solution to file x.vec\n");
      HYPRE_StructVectorPrint("x.vec",solvec,0);
      
      ENZO_FAIL("Error in AMRFLDSplit::RadStepUnigrid");
#ifndef FAIL_ON_NAN
    }
#endif
  }

  ////////////////// postprocess results and clean up //////////////////

  // if solve was successful: extract values and add to current solution; 
  // accumulate overall relative change in radiation field
  float atol = 0.1;
  Eflt64 loc_diff=0.0;
  Eflt64 weight, diff;
  dV = (dtnorm > 0.0) ? dV : 1.0;
  if (!recompute_step) {
    Eint32 ilower[3] = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
    Eint32 iupper[3] = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
    int xBuff = ghXl;
    int yBuff = ghYl-SolvIndices[1][0];
    int zBuff = ghZl-SolvIndices[2][0];
    for (int iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
      int Zbl = (iz+zBuff)*x0len*x1len;
      ilower[2] = iz;  iupper[2] = iz;
      for (int iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
	int Ybl = (iy+yBuff)*x0len;
	ilower[1] = iy;  iupper [1] = iy;
	HYPRE_StructVectorGetBoxValues(solvec, ilower, iupper, HYPREbuff);
	for (int ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) {
	  int idx = Zbl+Ybl+xBuff+ix;
	  weight = sqrt(fabs(Enew[idx]*Eold[idx])) + atol;
	  diff = fabs(HYPREbuff[ix])/weight;
	  if (dtnorm > 0.0) {
	    loc_diff += pow(diff,dtnorm)*dV;
	  } else {
	    loc_diff = max(loc_diff, diff);
	  }
	  Enew[idx] += HYPREbuff[ix];
	}
      }
    }
    // communicate to obtain overall relative change
    if (dtnorm > 0.0) {
      if (MPI_Allreduce(&loc_diff, &Echange, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != 0)
	ENZO_FAIL("AMRFLDSplit::RadStepUnigrid: Error in MPI_Allreduce");
      Echange = pow(Echange, 1.0/dtnorm);
    } else {
      if (MPI_Allreduce(&loc_diff, &Echange, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) != 0)
	ENZO_FAIL("AMRFLDSplit::RadStepUnigrid: Error in MPI_Allreduce");
    }
  }

  // destroy HYPRE structures
  switch (sol_type) {
  case 4:   // PCG
    HYPRE_StructPCGDestroy(solver);
    break;
  case 3:   // GMRES
    HYPRE_StructGMRESDestroy(solver);
    break;
  default:  // BiCGStab
    HYPRE_StructBiCGSTABDestroy(solver);
    break;
  }
  HYPRE_StructPFMGDestroy(preconditioner);
  HYPRE_StructVectorDestroy(rhsvec);
  HYPRE_StructVectorDestroy(solvec);
  HYPRE_StructMatrixDestroy(P);
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructGridDestroy(grid);

  // enforce a solution floor on radiation
  float epsilon=1.0;      // radiation floor
  while (epsilon*0.25 > 0.0)  epsilon*=0.5;
  epsilon *= 100.0;
  for (int i=0; i<x0len*x1len*x2len; i++)  
    Enew[i] = max(Enew[i], epsilon);

  // rescale dt, told, tnew back to normalized values
  dt   /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;

  return recompute_step;

#else   // USE_HYPRE

  ENZO_FAIL("AMRFLDSplit::RadStepUnigrid requires HYPRE to be enabled; halting run");
  return 1;

#endif
}

#endif   // TRANSFER
