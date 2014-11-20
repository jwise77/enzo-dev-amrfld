/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class, RadStep Routine
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



#ifdef AMR_SOLVE

// This routine evolves the radiation subsystem within the AMRFLDSplit module.  
// This is performed in a robust manner; if the current radiation subsystem 
// cannot be solved to the desired tolerance (meaning the time step is 
// likely much too large), it will return a flag to the calling routine to 
// have the step recomputed with a smaller time step size.
int AMRFLDSplit::RadStep(int ibin, LevelHierarchyEntry *LevelArray[], 
			 int level, AMRsolve_Hierarchy *hierarchy, 
			 float Etyp, float Emax, Eflt64 &Eerror)
{

  // update internal units for current times
  if (this->UpdateUnits(told, tnew) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_RadStep: Error in UpdateUnits.");

  // rescale dt, told, tnew to physical values for use within solver
  dt   *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;

  //   compute current opacity for this bin
  if (this->Opacity(ibin, LevelArray, level, tnew) != SUCCESS)
    ENZO_FAIL("AMRFLDSplit_RadStep: Opacity failure!!");
    
  //   enforce boundary conditions on new radiation field
  if (this->EnforceBoundary(ibin, LevelArray) != SUCCESS) 
    ENZO_FAIL("AMRFLDSplit_RadStep: EnforceBoundary failure!!");
    
  //   copy current radiation field into temporary field (KPhHI)
  //   on all grids owned by this processor, all levels here and down;
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	  // set grid dimension information
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

	  // access old/new radiation fields (old stored in KPhHI)
	  float *Eold = Temp->GridHierarchyEntry->GridData->AccessKPhHI();
	  float *Enew = AccessRadiationField(ibin, Temp->GridHierarchyEntry);
	  if (Eold == NULL)
	    ENZO_FAIL("AMRFLDSplit_RadStep: cannot access KPhHI field!!\n");
	  if (Enew == NULL) {
	    ENZO_VFAIL("AMRFLDSplit_RadStep: cannot access radiation field %"ISYM"!!\n", ibin);
	  }

	  // copy new radiation field into into old field
	  for (int k=0; k<x0len*x1len*x2len; k++)  Eold[k] = Enew[k];

      }  // end loop over grids on this proc
  
  // calculate emissivity source norm
  float glob_eta = TotalFieldValue(LevelArray, level, MAX_DEPTH_OF_HIERARCHY, Emissivity0+ibin);
  if (debug)   printf("   total emissivity %"ISYM" = %g\n", ibin, glob_eta);

  // determine if the radiation system needs to be solved in this step:
  //   if not, set error to 0, reset units, return with a successful step; 
  float srctol=1.0e-30;
  float radtol=1.0e-8;
  if ((glob_eta*dt < srctol) && (fabs(Etyp-Emax) < radtol)) {
    Eerror = 1e-16;
    dt   /= TimeUnits;
    told /= TimeUnits;
    tnew /= TimeUnits;
    return 0;
  }

  // if we've made it here, turn on automatic scaling for next time step
  StartAutoScale[ibin] = true;

  // set up and solve this radiation field's equation via amrsolve

  // Initialize the amrsolve FLD solver
  AMRsolve_Hypre_FLD amrfldsolve(*hierarchy, *amrsolve_params, ibin, sol_prec);
  amrfldsolve.init_hierarchy();
  // if (debug)   hierarchy->print();
  amrfldsolve.init_stencil();
  amrfldsolve.init_graph();

  //    initialize amrfldsolve system
  amrfldsolve.init_elements(dt, theta, NiUnits, NiUnits0, LenUnits, LenUnits0, 
			    ErUnits[ibin], ErUnits0[ibin], BdryType);


  // /////////  AMRsolve_Hypre_FLD testing routine ////////
  // amrfldsolve.tester();
  // MPI_Barrier(MPI_COMM_WORLD);
  // ENZO_FAIL("Stopping run (testing)");
  // //////////////////////////////////////////////////////


  //    solve amrfldsolve system
  amrfldsolve.solve();
  Eflt64 finalresid = amrfldsolve.residual();
  Eint32 Sits = amrfldsolve.iterations();
  if (debug) printf("   lin resid = %.1e (tol = %.1e), its = %i\n",
		    finalresid, sol_tolerance, Sits);


  // check solution, and if the solver fails:
  //  * if we have room to decrease the time step size, do so and allow remainder 
  //    of function to complete (to reset units, etc.), but have calling routine 
  //    update dt and compute step again.
  //  * if we have don't have room to decrease the time step size, output 
  //    parameters and error message.
  int recompute_step = 0;
  if (amrfldsolve.evaluate() != 0) {
    if (dt > mindt*TimeUnits) {
      recompute_step = 1;
    }
    else {
      fprintf(stderr,"AMRFLDSplit_RadStep: could not achieve prescribed tolerance!\n");
      fprintf(stderr, "Printing current hierarchy:\n");
      if (debug)   hierarchy->print();
      
      // dump amrsolve matrices, module parameters to disk
      amrfldsolve.abort_dump();
      if (debug) {
	fprintf(stderr,"Dumping AMRFLDSplit module parameters to file RTdump.params\n");
	FILE *fptr = fopen("RTdump.params", "w");
	this->WriteParameters(fptr);
	fclose(fptr);
      }
      
      ENZO_FAIL("Error in AMRFLDSplit_RadStep");
    }
  }

  // increment Enzo radiation field with amrsolve solution
  amrfldsolve.update_enzo();

  // enforce a solution floor on radiation
  float epsilon=1.0;      // radiation floor, set to 100*roundoff
  while ((1.0 + epsilon) > 1.0)  epsilon*=0.5;
  epsilon *= 100.0;
  for (int thislevel=level; thislevel<MAX_DEPTH_OF_HIERARCHY; thislevel++)
    for (LevelHierarchyEntry* Temp=LevelArray[thislevel]; Temp; 
	 Temp=Temp->NextGridThisLevel)
      if (MyProcessorNumber == Temp->GridHierarchyEntry->GridData->ReturnProcessorNumber()) {

	  // set grid dimension information
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

	  // access new radiation fields
	  float *Enew = AccessRadiationField(ibin, Temp->GridHierarchyEntry);
	  if (Enew == NULL) {
	    ENZO_VFAIL("AMRFLDSplit_RadStep: cannot access radiation field %"ISYM"!!\n",
		       ibin);
	  }
	  
	  // enforce floor on new field
	  for (int k=0; k<x0len*x1len*x2len; k++)  
	    Enew[k] = max(Enew[k], epsilon);

      }  // end loop over grids on this proc


  // If this was a successful solve, compute relative change in radiation field 
  // solution over time step
  if (recompute_step == 0) 
    Eerror = amrfldsolve.rdiff_norm(dtnorm, 0.1);

  // rescale dt, told, tnew back to normalized values
  dt   /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;

  return recompute_step;

}

#endif   // AMR_SOLVE

#endif   // TRANSFER
