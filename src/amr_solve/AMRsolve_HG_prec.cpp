/// @file      AMRsolve_HG_prec.cpp
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Implementation of the AMRsolve_HG_prec class

#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>

#include "HYPRE_sstruct_ls.h"
#include "HYPRE_krylov.h"

#include "hdf5.h"

#include "AMRsolve_defs.h"

//----------------------------------------------------------------------

#include "AMRsolve_mpi.h"
#include "AMRsolve_scalar.h"
#include "AMRsolve_error.h"
#include "AMRsolve_constants.h"
#include "AMRsolve_faces.h"
#include "AMRsolve_mpi.h"
#include "AMRsolve_domain.h"
#include "AMRsolve_grid.h"
#include "AMRsolve_level.h"
#include "AMRsolve_hierarchy.h"
#include "AMRsolve_parameters.h"
#include "AMRsolve_HG_prec.h"
#include "AMRsolve_error.h"

//======================================================================
// PUBLIC MEMBER FUNCTIONS
//======================================================================

/// AMRsolve_HG_prec constructor
AMRsolve_HG_prec::AMRsolve_HG_prec(AMRsolve_Hierarchy& hierarchy, 
				   int BdryType[3][2]) :
  hierarchy_(&hierarchy), csolver_(0), citer_(0), resid_(-1.0)
{
  // set initialized flag to false
  initialized = false;

  // initialize remaining pointers to NULL
  Ac_ = NULL;
  Bc_ = NULL;
  Xc_ = NULL;
  Y_  = NULL; 

  // insert array-valued argument into object
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      BdryType_[i][j] = BdryType[i][j];
}

//----------------------------------------------------------------------

/// AMRsolve_Hypre_FLD destructor
AMRsolve_HG_prec::~AMRsolve_HG_prec()
{
  // destroy HYPRE coarse grid solver
  if (HYPRE_StructPFMGDestroy(csolver_))  ERROR("could not destroy csolver_\n");
}

//----------------------------------------------------------------------

/// Initialize routine for 2-level hierarchical-grid preconditioner.
/// Need to create and set parameters into coarse grid PFMG solver.
int AMRsolve_HG_prec::Initialize_(AMRsolve_Parameters *parameters,
				  HYPRE_StructMatrix  *A,
				  HYPRE_StructVector  *X,
				  HYPRE_StructVector  *B,
				  HYPRE_SStructVector *Y)
{
  // set output flag to success
  int ierr = 0;

  // insert arguments into object
  Ac_ = A;
  Xc_ = X;
  Bc_ = B;
  Y_  = Y;

  // Create and initialize the solver
  ierr = HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &csolver_);
  if (ierr != 0) ERROR("could not create PFMG solver_\n");

  // extract preconditioner parameters
  std::string sitmax   = parameters->value("prec_itmax");
  std::string srestol  = parameters->value("prec_restol");
  std::string srlxtype = parameters->value("prec_rlxtype");
  std::string snpre    = parameters->value("prec_npre");
  std::string snpost   = parameters->value("prec_npost");
  std::string sprintl  = parameters->value("prec_printl");
  std::string slog     = parameters->value("prec_log");
  std::string Jitmax   = parameters->value("prec_Jaciters");

  //   if not defined, then define them
  if (sitmax   == "")  parameters->add_parameter("prec_itmax","20"); 
  if (srestol  == "")  parameters->add_parameter("prec_restol","1.0e-6");
  if (srlxtype == "")  parameters->add_parameter("prec_rlxtype","1");
  if (snpre    == "")  parameters->add_parameter("prec_npre","1"); 
  if (snpost   == "")  parameters->add_parameter("prec_npost","1");
  if (sprintl  == "")  parameters->add_parameter("prec_printl","1");
  if (slog     == "")  parameters->add_parameter("prec_log","1");
  if (Jitmax   == "")  parameters->add_parameter("prec_Jaciters","1");

  // re-extract solver parameters now that everything is set
  sitmax   = parameters->value("prec_itmax");
  srestol  = parameters->value("prec_restol");
  srlxtype = parameters->value("prec_rlxtype");
  snpre    = parameters->value("prec_npre");
  snpost   = parameters->value("prec_npost");
  sprintl  = parameters->value("prec_printl");
  slog     = parameters->value("prec_log");
  Jitmax   = parameters->value("prec_Jaciters");

  //   set local and saved variables
  double restol = atof(srestol.c_str());
  int itmax     = atoi(sitmax.c_str());
  int rlxtype   = atoi(srlxtype.c_str());
  int npre      = atoi(snpre.c_str());
  int npost     = atoi(snpost.c_str());
  int printl    = atoi(sprintl.c_str());
  int log       = atoi(slog.c_str());
  Jacobi_iters  = atoi(Jitmax.c_str());

  // determine the maximum number of MG levels (due to periodicity)
  int max_levels, Ndir, level=-1;
  max_levels = 100000;
  for (int idim=0; idim<3; idim++) {
    if (BdryType_[idim][0] == 0) {
      level = 0;
      Ndir = hierarchy_->period_index(idim,0);
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = MIN(max_levels,level);
  }
  
  // set solver options
  if (max_levels > -1) 
    HYPRE_StructPFMGSetMaxLevels(csolver_, max_levels);
  if (itmax != 0 ) {
    ierr = HYPRE_StructPFMGSetMaxIter(csolver_,itmax);
    if (ierr != 0) ERROR("could not set itmax\n");
  }
  if (restol != 0.0) {
    ierr = HYPRE_StructPFMGSetTol(csolver_,restol);
    if (ierr != 0) ERROR("could not set restol\n");
  }
  ierr = HYPRE_StructPFMGSetRelaxType(csolver_,rlxtype);
  if (ierr != 0) ERROR("could not set rlxtype\n");
  ierr = HYPRE_StructPFMGSetNumPreRelax(csolver_,npre);
  if (ierr != 0) ERROR("could not set npre\n");
  ierr = HYPRE_StructPFMGSetNumPostRelax(csolver_,npost);
  if (ierr != 0) ERROR("could not set npost\n");
  ierr = HYPRE_StructPFMGSetPrintLevel(csolver_,printl);
  if (ierr != 0) ERROR("could not set printl\n");
  ierr = HYPRE_StructPFMGSetLogging(csolver_,log);
  if (ierr != 0) ERROR("could not set log\n");

  // set initialized flag to true
  initialized = true;

  return ierr;
}  // AMRsolve_HG_prec::Initialize_

//------------------------------------------------------------------------

/// Setup routine for 2-level hierarchical-grid preconditioner.
int AMRsolve_HG_prec::Setup_(HYPRE_SStructMatrix A, 
			     HYPRE_SStructVector b,
			     HYPRE_SStructVector x)
{
  // check whether coarse solver has been initialized
  if (!initialized) {
    printf("AMRsolve_HG_prec::Setup_ error, module is not initialized!\n");
    return 1;
  }

  // Only needs to set up the coarse grid solver
  return (HYPRE_StructPFMGSetup(csolver_, *Ac_, *Bc_, *Xc_));
}  // AMRsolve_HG_prec::Setup_


// wrapper routine for passing HG_prec_setup_ to HYPRE solvers
HYPRE_Int HG_prec_setup(HYPRE_SStructSolver solver, 
			HYPRE_SStructMatrix A, 
			HYPRE_SStructVector b, 
			HYPRE_SStructVector x) {

  // unpack "solver" and call module setup routine
  AMRsolve_HG_prec *prec = (AMRsolve_HG_prec *) solver;
  return (prec->Setup_(A,b,x));
}

//------------------------------------------------------------------------

/// Solve routine for 2-level hierarchical-grid preconditioner.
int AMRsolve_HG_prec::Solve_(HYPRE_SStructMatrix A, 
			     HYPRE_SStructVector b,
			     HYPRE_SStructVector x)
{
  // check whether coarse solver has been initialized
  if (!initialized) {
    printf("AMRsolve_HG_prec::Setup_ error, module is not initialized!\n");
    return 1;
  }

  // set output flag to success
  int ierr = 0;

  // copy b data into AMRsolve grid u vectors
  ierr = HYPRE_to_AMRsolve_(&b, 1);
  if (ierr != 0)  ERROR("could not convert b into AMRsolve grid u vectors\n");
  
  // restrict u onto coarse grid
  ierr = restrict(hierarchy_->num_levels()-1, 0);
  if (ierr != 0)  ERROR("could not restrict u to coarse grid\n");
  
  // copy coarse grid u into Bc_ vector
  ierr = AMRsolve_to_HYPRE_coarse_(Bc_, 1);

  // clear out solution vector of old data
  ierr = HYPRE_StructVectorSetConstantValues(*Xc_, 0.0);
  if (ierr != 0)  ERROR("could not initialize Xc_ to 0.0\n");

  // Solve the coarse-grid linear system, update statistics
  int iters;
  ierr = HYPRE_StructPFMGSolve(csolver_,*Ac_,*Bc_,*Xc_);
  if (ierr != 0)  ERROR("could not solve coarse problem with PFMG\n");
  ierr = HYPRE_StructPFMGGetFinalRelativeResidualNorm(csolver_,&resid_);
  if (ierr != 0) ERROR("could not get resid_\n");
  ierr = HYPRE_StructPFMGGetNumIterations(csolver_,&iters);
  if (ierr != 0) ERROR("could not get iters\n");
  citer_ += iters;

  // copy Xc_ vector back to coarse grid u
  ierr = HYPRE_to_AMRsolve_coarse_(Xc_, 1);

  // prolong u onto full hierarchy
  ierr = prolong(0, hierarchy_->num_levels()-1, 0);
  if (ierr != 0)  ERROR("could not prolong u to hierarchy\n");

  // copy u into output vector x
  ierr = AMRsolve_to_HYPRE_(&x, 1);
  if (ierr != 0)  ERROR("could not convert AMRsolve grid u vector into x\n");

  // perform Jacobi smoothing steps on full hierarchy, using temporary vector Y_
  for (int Jit=0; Jit<Jacobi_iters; Jit++) {
    ierr = Jacobi_smooth_(A, x, b, *Y_);
    if (ierr != 0)  ERROR("could not perform Jacobi smoother\n");
  }

  return ierr;
}  // AMRsolve_HG_prec::Solve_

// wrapper routine for passing HG_prec_solve_ to HYPRE solvers
HYPRE_Int HG_prec_solve(HYPRE_SStructSolver solver, 
			HYPRE_SStructMatrix A, 
			HYPRE_SStructVector b, 
			HYPRE_SStructVector x) {

  // unpack "solver" object and call corresponding solve routine
  AMRsolve_HG_prec *prec = (AMRsolve_HG_prec *) solver;
  return (prec->Solve_(A,b,x));
}

//----------------------------------------------------------------------

/// Utility routine to perform one Jacobi iteration,  
///                v = v + D^{-1}*(b - A*v).
/// Arguments: v holds the current guess, b holds the RHS, and y is a 
///            temporary.  The vector b is not modified.
int AMRsolve_HG_prec::Jacobi_smooth_(HYPRE_SStructMatrix A,
				     HYPRE_SStructVector v, 
				     HYPRE_SStructVector b,
				     HYPRE_SStructVector y)
{
  // set output flag to success
  int ierr = 0;

  // copy b into y  [y = b]
  ierr = HYPRE_SStructVectorCopy(b, y);
  if (ierr != 0) {
    fprintf(stderr,"Jacobi_smooth_ error in HYPRE_SStructVectorCopy = %i\n",ierr);
    return ierr;
  }

  // perform matrix axpy  [y = b - A*v]
  ierr = HYPRE_SStructMatrixMatvec(-1.0, A, v, 1.0, y);
  if (ierr != 0) {
    fprintf(stderr,"Jacobi_smooth_ error in HYPRE_SStructMatrixMatvec = %i\n",ierr);
    return ierr;
  }

  // gather y
  ierr = HYPRE_SStructVectorGather(y);
  if (ierr != 0) {
    fprintf(stderr,"Jacobi_smooth_ error in HYPRE_SStructVectorGather = %i\n",ierr);
    return ierr;
  }

  // perform Jacobi solve on all grids owned by this processor  [y = D^{-1}*(b - A*v)]
  int n0, n1, n2;
  Scalar *u, *f;
  int lower[3], upper[3];
  int diagonal[1] = {0};
  ItHierarchyGridsLocal itgl(*hierarchy_);
  while (AMRsolve_Grid* grid = itgl++) {

    // access y and D, storing them in local arrays u and f, respectively
    grid->get_limits(lower,upper);
    u = grid->get_u(&n0,&n1,&n2);
    f = grid->get_f(&n0,&n1,&n2);
    ierr = HYPRE_SStructVectorGetBoxValues(y, grid->level(), lower, upper, 0, u);
    if (ierr != 0) {
      fprintf(stderr,"Jacobi_smooth_ error in HYPRE_SStructVectorGetBoxValues = %i\n",ierr);
      return ierr;
    }
    ierr = HYPRE_SStructMatrixGetBoxValues(A, grid->level(), lower, upper, 0, 1, diagonal, f);
    if (ierr != 0) {
      fprintf(stderr,"Jacobi_smooth_ error in HYPRE_SStructMatrixGetBoxValues = %i\n",ierr);
      return ierr;
    }

    // perform Jacobi update
    for (int i=0; i<n0*n1*n2; i++)  u[i] /= f[i];
    
    // place result back into y
    ierr = HYPRE_SStructVectorSetBoxValues(y, grid->level(), lower, upper, 0, u);
    if (ierr != 0) {
      fprintf(stderr,"Jacobi_smooth_ error in HYPRE_SStructVectorSetBoxValues = %i\n",ierr);
      return ierr;
    }

  } // while itgl

  // assemble y
  ierr = HYPRE_SStructVectorAssemble(y);
  if (ierr != 0) {
    fprintf(stderr,"Jacobi_smooth_ error in HYPRE_SStructVectorAssemble = %i\n",ierr);
    return ierr;
  }

  // apply Jacobi update to v  [v = v + D^{-1}*(b - A*v)]
  ierr = HYPRE_SStructAxpy(1.0, y, v);
  if (ierr != 0) {
    fprintf(stderr,"Jacobi_smooth_ error in HYPRE_SStructAxpy = %i\n",ierr);
    return ierr;
  }

  return ierr;
}  // AMRsolve_HG_prec::Jacobi_smooth_

//------------------------------------------------------------------------

/// Utility routine to convert from u to HYPRE_SStructVector over full hierarchy. 
int AMRsolve_HG_prec::AMRsolve_to_HYPRE_(HYPRE_SStructVector *V, int u_vs_f)
{
  // set output flag to success
  int ierr = 0;

  // set boolean for choice of u or f
  int use_u = (u_vs_f != 0);

  // iterate over all local grids, copying from HYPRE vector to u or f
  int n0, n1, n2;
  Scalar *values;
  int lower[3], upper[3];
  ItHierarchyGridsLocal itg(*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // get level, grid information
    int level = grid->level();
    grid->get_limits(lower,upper);

    // copy AMRsolve vector into HYPRE vector
    if (use_u)
      values = grid->get_u(&n0,&n1,&n2);
    else
      values = grid->get_f(&n0,&n1,&n2);
    ierr = HYPRE_SStructVectorSetBoxValues(*V, level, lower, upper, 0, values);
    if (ierr != 0)  return ierr;
  }

  // Assemble HYPRE SStruct vector
  ierr = HYPRE_SStructVectorAssemble(*V);
  return ierr;

}  // AMRsolve_HG_prec::AMRsolve_to_HYPRE_

//------------------------------------------------------------------------

/// Utility routine to convert from HYPRE_SStructVector to u over full hierarchy. 
int AMRsolve_HG_prec::HYPRE_to_AMRsolve_(HYPRE_SStructVector *V, int u_vs_f)
{
  // set output flag to success
  int ierr = 0;

  // set boolean for choice of u or f
  int use_u = (u_vs_f != 0);

  // gather HYPRE SStruct vector
  ierr = HYPRE_SStructVectorGather(*V);
  if (ierr != 0)  return ierr;

  // iterate over all local grids, copying from HYPRE vector to u or f
  int n0, n1, n2;
  Scalar *values;
  int lower[3], upper[3];
  ItHierarchyGridsLocal itg(*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // get level, grid information
    int level = grid->level();
    grid->get_limits(lower,upper);

    // extract HYPRE vector into AMRsolve vector
    if (use_u)
      values = grid->get_u(&n0,&n1,&n2);
    else
      values = grid->get_f(&n0,&n1,&n2);
    ierr = HYPRE_SStructVectorGetBoxValues(*V, level, lower, upper, 0, values);
    if (ierr != 0)  return ierr;
  }
  return ierr;

}  // AMRsolve_HG_prec::HYPRE_to_AMRsolve_

//------------------------------------------------------------------------

/// Utility routine to convert from u to HYPRE_StructVector over coarse grid only. 
int AMRsolve_HG_prec::AMRsolve_to_HYPRE_coarse_(HYPRE_StructVector *V, int u_vs_f)
{
  // set output flag to success
  int ierr = 0;

  // set boolean for choice of u or f
  int use_u = (u_vs_f != 0);

  // iterate over local coarse grids, copying from HYPRE vector to u or f
  int n0, n1, n2;
  Scalar *values;
  int lower[3], upper[3];
  AMRsolve_Level *coarse = &(hierarchy_->level(0));
  ItLevelGridsLocal itlg (*coarse);
  while (AMRsolve_Grid* grid = itlg++) {

    // get grid information
    grid->get_limits(lower,upper);

    // copy AMRsolve vector into HYPRE vector
    if (use_u)
      values = grid->get_u(&n0,&n1,&n2);
    else
      values = grid->get_f(&n0,&n1,&n2);
    ierr = HYPRE_StructVectorSetBoxValues(*V, lower, upper, values);
    if (ierr != 0)  return ierr;
  }

  // Assemble HYPRE Struct vector
  ierr = HYPRE_StructVectorAssemble(*V);
  return ierr;

}  // AMRsolve_HG_prec::AMRsolve_to_HYPRE_coarse_

//------------------------------------------------------------------------

/// Utility routine to convert from HYPRE SStructVector to u over coarse grid only.
int AMRsolve_HG_prec::HYPRE_to_AMRsolve_coarse_(HYPRE_StructVector *V, int u_vs_f)
{
  // set output flag to success
  int ierr = 0;

  // set boolean for choice of u or f
  int use_u = (u_vs_f != 0);

  // iterate over local coarse grids, copying from HYPRE vector to u or f
  int n0, n1, n2;
  Scalar *values;
  int lower[3], upper[3];
  AMRsolve_Level *coarse = &(hierarchy_->level(0));
  ItLevelGridsLocal itlg (*coarse);
  while (AMRsolve_Grid* grid = itlg++) {

    // get grid information
    grid->get_limits(lower,upper);

    // extract HYPRE vector into AMRsolve vector
    if (use_u)
      values = grid->get_u(&n0,&n1,&n2);
    else
      values = grid->get_f(&n0,&n1,&n2);
    ierr = HYPRE_StructVectorGetBoxValues(*V, lower, upper, values);
    if (ierr != 0)  return ierr;
  }
  return ierr;

}  // AMRsolve_HG_prec::HYPRE_to_AMRsolve_coarse_

//------------------------------------------------------------------------

/// Restriction operator between levels
/** Restricts all grid's u_ data arrays from level_fine through 
  * level_coarse, using the requested restriction operator; returns
  *  integer denoting pass (0) or fail (1) */
int AMRsolve_HG_prec::restrict(int level_fine, int level_coarse) throw()
{
  // check for legal input arguments
  if (level_coarse < 0) {
    printf("AMRsolve_HG_prec::restrict error -- illegal level_coarse (%i < 0)\n",
	   level_coarse);
    return 1;
  }
  if (level_fine > hierarchy_->num_levels()-1) {
    printf("AMRsolve_HG_prec::restrict error -- illegal level_fine (%i > %i)\n",
	   level_fine, hierarchy_->num_levels()-1);
    return 1;
  }

  // count total number of receive buffers needed on this processor
  int total_buffers=0;
  ItHierarchyGridsLocal itg (*hierarchy_);
  while (AMRsolve_Grid * g = itg++) {
    if (g->level() < level_fine) {
      ItGridChildren itgc (*g);
      while (AMRsolve_Grid * c = itgc++) {
	total_buffers++;
      }
    }
  }

  // allocate/open buffers for all local grids to receive childrens' overlap data
  Scalar **rbuffs = new Scalar*[total_buffers];
  MPI_Request *reqs = new MPI_Request[total_buffers];
  int *tags = new int[total_buffers];
  int g1_ilo[3], g1_ihi[3], g2_ilo[3], g2_ihi[3];
  ItHierarchyGridsLocal itg1 (*hierarchy_);
  int ibuff = 0;
  MPI_Datatype Stype = (sizeof(Scalar) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  while (AMRsolve_Grid * g = itg1++) {
    if (g->level() < level_fine) {
      ItGridChildren itgc (*g);
      while (AMRsolve_Grid * c = itgc++) {

	// check for overlap, and get region indices
	if (g->overlap_indices(c, g1_ilo, g1_ihi, g2_ilo, g2_ihi)) {

	  // allocate receive buffer
	  int size=1;
	  for (int i=0; i<3; i++)  size *= (g2_ihi[i]-g2_ilo[i]+1);
	  rbuffs[ibuff] = new Scalar[size];

	  // set unique exchange tag
	  tags[ibuff] = g->id() + c->id() * 1000000;

	  // open receive channel, and increment buffer index
	  MPI_Comm com = g->get_mpi()->comm();
	  if (MPI_Irecv(&(rbuffs[ibuff][0]), size, Stype, c->ip(), tags[ibuff], 
			com, &(reqs[ibuff])) != 0) 
	    ERROR("Error in MPI_Irecv call");
	  ibuff++;
	}
      }
    }
  }

  
  /* starting at the top-most level:
     (a) wait for childrens' overlap region data
     (b) restrict fine grid data in overlap region onto u_ grid array
     (c) send data to parent
  */
  int tag;
  MPI_Status status;
  for (int ilevel=level_fine; ilevel >= level_coarse; ilevel--) {

    // iterate over grids local to this level
    ItLevelGridsLocal itlg (hierarchy_->level(ilevel));
    while (AMRsolve_Grid *g = itlg++) {

      // if I'm not on level_fine, get data and perform restrictions
      if (ilevel < level_fine) {

	// iterate over children
	ItGridChildren itgc (*g);
	while (AMRsolve_Grid * c = itgc++) {
	  
	  // check for overlap, and get region indices
	  if (g->overlap_indices(c, g1_ilo, g1_ihi, g2_ilo, g2_ihi)) {

	    // determine unique exchange tag
	    tag = g->id() + c->id() * 1000000;
	    
	    // find buffer matching this tag
	    ibuff = -1;
	    for (int j=0; j<total_buffers; j++) {
	      if (tag == tags[j])  ibuff = j;
	    }
	    if (ibuff == -1) 
	      ERROR("Cannot find matching buffer");

	    // wait for receive to finish
	    if (MPI_Wait(&(reqs[ibuff]), &status) != 0)
	      ERROR("Error in MPI_Wait call");

	    // access my u array
	    int nx[3];
	    Scalar *mydata = g->get_u(&nx[0],&nx[1],&nx[2]);

	    // perform restriction
	    int ovsize[3] = { g2_ihi[0] - g2_ilo[0] + 1, 
			      g2_ihi[1] - g2_ilo[1] + 1, 
			      g2_ihi[2] - g2_ilo[2] + 1};
	    do_restrict(mydata, rbuffs[ibuff], nx, g1_ilo, g1_ihi, ovsize);

	    // delete this receive buffer
	    delete[] rbuffs[ibuff];
	  }  // if overlap_indices
	} // while itgc
      } // if ilevel

      // if I'm not on level_coarse, send data to parent
      if (ilevel > level_coarse) {

	// set pointer to parent grid, through searching over grids on coarser level
	ItLevelGridsAll itpg (hierarchy_->level(ilevel-1));
	AMRsolve_Grid *parent;
	while (AMRsolve_Grid * p = itpg++) {
	  if (g->id_parent() == p->id()) {
	    parent = p;
	    break;
	  }
	}

	// determine overlap region
	if (g->overlap_indices(parent, g1_ilo, g1_ihi, g2_ilo, g2_ihi)) {

	  // allocate and fill send buffer
	  int size=1;
	  for (int i=0; i<3; i++)  size *= (g1_ihi[i]-g1_ilo[i]+1);
	  Scalar *sbuff = new Scalar[size];
	  int nx[3];
	  Scalar *mydata = g->get_u(&nx[0],&nx[1],&nx[2]);
	  int loc=0;
	  for (int k=g1_ilo[2]; k<=g1_ihi[2]; k++)
	    for (int j=g1_ilo[1]; j<=g1_ihi[1]; j++)
	      for (int i=g1_ilo[0]; i<=g1_ihi[0]; i++) 
		sbuff[loc++] = mydata[(k*nx[1] + j)*nx[0] + i];

	  // send buffer to parent
	  tag = parent->id() + g->id() * 1000000;
	  MPI_Request req;
	  MPI_Comm com = g->get_mpi()->comm();
	  if (MPI_Isend(sbuff, size, Stype, parent->ip(), tag, 
			com, &req) != 0) 
	    ERROR("Error in MPI_Isend call");

	  // wait for send to finish, and then delete buffer
	  if (MPI_Wait(&req, &status) != 0)
	    ERROR("Error in MPI_Wait call");
	  delete[] sbuff;
	  
	} // if overlap_indices
      } // if ilevel
    } // while itlg
  } // for ilevel


  // clean up
  delete[] rbuffs;
  delete[] reqs;
  delete[] tags;

  return 0;
}

//----------------------------------------------------------------------

/// Routine to actually perform the restriction of the overlap data
/** Performs the restriction of the data in overlap region (that has 
  * size ovsize[3]) into the local data array mydata (of size mysize[3]), 
  * where the overlap data should be placed into the region defined by 
  * g1_ilo[3] nad g1_ihi[3]. */
void AMRsolve_HG_prec::do_restrict(Scalar *mydata, Scalar *overlap, 
				     int *mysize, int *g1_ilo, 
				     int *g1_ihi, int *ovsize) throw()
{

  // ensure that overlap region fits within mydata
  for (int i=0; i<3; i++) 
    if (g1_ihi[i]-g1_ilo[i]+1 > mysize[i])
      ERROR("Error, overlap region larger than local data array");

  // determine refinement factor (should be uniform in all dims)
  int rfac = ovsize[0]/(g1_ihi[0]-g1_ilo[0]+1);
  Scalar rsc = 1.0/rfac/rfac/rfac;

  // Perform restriction -- iterate over local overlap region
  for (int k=g1_ilo[2]; k<=g1_ihi[2]; k++)
    for (int j=g1_ilo[1]; j<=g1_ihi[1]; j++)
      for (int i=g1_ilo[0]; i<=g1_ihi[0]; i++) {

	// iterate over refined cells to accumulate average
	mydata[(k*mysize[1] + j)*mysize[0] + i] = 0.0;
	for (int k2=rfac*(k-g1_ilo[2]); k2<rfac*(k-g1_ilo[2]+1); k2++)
	  for (int j2=rfac*(j-g1_ilo[1]); j2<rfac*(j-g1_ilo[1]+1); j2++)
	    for (int i2=rfac*(i-g1_ilo[0]); i2<rfac*(i-g1_ilo[0]+1); i2++)
	      mydata[(k*mysize[1] + j)*mysize[0] + i] += 
		overlap[(k2*ovsize[1] + j2)*ovsize[0] + i2]*rsc;
      }
  
}

//----------------------------------------------------------------------

/// Prolongation operator between levels
/** Prolongs all grid's u_ data arrays from level_coarse up to 
  * level_fine, using the requested prolongation operator; returns 
  * integer denoting pass (0) or fail (1) */
int AMRsolve_HG_prec::prolong(int level_coarse, int level_fine, 
				int method) throw()
{

  // check for legal input arguments
  if (level_coarse < 0) {
    printf("AMRsolve_HG_prec::prolong error -- illegal level_coarse (%i < 0)\n",
	   level_coarse);
    return 1;
  }
  if (level_fine > hierarchy_->num_levels()-1) {
    printf("AMRsolve_HG_prec::prolong error -- illegal level_fine (%i > %i)\n",
	   level_fine,hierarchy_->num_levels()-1);
    return 1;
  }

  // count total number of receive buffers needed on this processor
  int total_buffers=0;
  ItHierarchyGridsLocal itg (*hierarchy_);
  while (AMRsolve_Grid * g = itg++) {
    if (g->level() > level_coarse)  total_buffers++;
  }

  // allocate/open buffers for all local grids to receive parent's overlap data
  Scalar **rbuffs = new Scalar*[total_buffers];
  MPI_Request *reqs = new MPI_Request[total_buffers];
  int *tags = new int[total_buffers];
  int g1_ilo[3], g1_ihi[3], g2_ilo[3], g2_ihi[3];
  ItHierarchyGridsLocal itg1 (*hierarchy_);
  int ibuff = 0;
  MPI_Datatype Stype = (sizeof(Scalar) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  while (AMRsolve_Grid * g = itg1++) {
    if (g->level() > level_coarse) {

      // set pointer to parent grid, through searching over grids on coarser level
      ItLevelGridsAll itpg (hierarchy_->level(g->level()-1));
      AMRsolve_Grid *parent;
      while (AMRsolve_Grid * p = itpg++) {
	if (g->id_parent() == p->id()) {
	  parent = p;
	  break;
	}
      }
      
      // determine overlap region
      if (g->overlap_indices(parent, g1_ilo, g1_ihi, g2_ilo, g2_ihi)) {

	// allocate receive buffer
	int size=1;
	for (int i=0; i<3; i++)  size *= (g1_ihi[i]-g1_ilo[i]+1);
	rbuffs[ibuff] = new Scalar[size];
	
	// set unique exchange tag
	tags[ibuff] = g->id() + parent->id() * 1000000;

	// open receive channel, and increment buffer index
	MPI_Comm com = g->get_mpi()->comm();
	if (MPI_Irecv(rbuffs[ibuff], size, Stype, parent->ip(), 
		      tags[ibuff], com, &(reqs[ibuff])) != 0) 
	  ERROR("Error in MPI_Irecv call");
	ibuff++;
      } // if overlap_indices
    } // if g->level()
  } // while itg1


  /* starting at the bottom-most level:
     (a) wait for parent's overlap region data
     (b) prolong coarse grid data in overlap region onto u_ grid array
     (c) send data to children
  */
  int tag;
  MPI_Status status;
  for (int ilevel=level_coarse; ilevel <= level_fine; ilevel++) {

    // iterate over grids local to this level
    ItLevelGridsLocal itlg (hierarchy_->level(ilevel));
    while (AMRsolve_Grid *g = itlg++) {

      // if I'm not on level_coarse, get data and perform prolongation
      if (ilevel > level_coarse) {

	// set pointer to parent grid, through searching over grids on coarser level
	ItLevelGridsAll itpg (hierarchy_->level(ilevel-1));
	AMRsolve_Grid *parent;
	while (AMRsolve_Grid * p = itpg++) {
	  if (g->id_parent() == p->id()) {
	    parent = p;
	    break;
	  }
	}
	
	// determine overlap region
	if (g->overlap_indices(parent, g1_ilo, g1_ihi, g2_ilo, g2_ihi)) {

	  // determine unique exchange tag
	  tag = g->id() + parent->id() * 1000000;
	  
	  // find buffer matching this tag
	  ibuff = -1;
	  for (int j=0; j<total_buffers; j++) {
	    if (tag == tags[j])  ibuff = j;
	  }
	  if (ibuff == -1) 
	    ERROR("Cannot find matching buffer");


	  // wait for receive to finish
	  if (MPI_Wait(&(reqs[ibuff]), &status) != 0)
	    ERROR("Error in MPI_Wait call");
	  
	  // access my u array
	  int nx[3];
	  Scalar *mydata = g->get_u(&nx[0],&nx[1],&nx[2]);
	  
	  // perform prolongation
	  int ovsize[3] = { g2_ihi[0] - g2_ilo[0] + 1, 
			    g2_ihi[1] - g2_ilo[1] + 1, 
			    g2_ihi[2] - g2_ilo[2] + 1};
	  do_prolong(mydata, rbuffs[ibuff], nx, g1_ilo, g1_ihi, ovsize, method);
	  
	  // delete this receive buffer
	  delete[] rbuffs[ibuff];
	  
	}  // if overlap_indices
      } // if ilevel


      // if I'm not on level_fine, send data to children
      if (ilevel < level_fine) {

	// iterate over children
	ItGridChildren itgc (*g);
	while (AMRsolve_Grid * c = itgc++) {
	  
	  // check for overlap, and get region indices
	  if (g->overlap_indices(c, g1_ilo, g1_ihi, g2_ilo, g2_ihi)) {

	    // allocate and fill send buffer
	    int size=1;
	    for (int i=0; i<3; i++)  size *= (g1_ihi[i]-g1_ilo[i]+1);
	    Scalar *sbuff = new Scalar[size];
	    int nx[3];
	    Scalar *mydata = g->get_u(&nx[0],&nx[1],&nx[2]);
	    int loc=0;
	    for (int k=g1_ilo[2]; k<=g1_ihi[2]; k++)
	      for (int j=g1_ilo[1]; j<=g1_ihi[1]; j++)
		for (int i=g1_ilo[0]; i<=g1_ihi[0]; i++)
		  sbuff[loc++] = mydata[(k*nx[1] + j)*nx[0] + i];

	    // send buffer to child
	    tag = c->id() + g->id() * 1000000;
	    MPI_Request req;
	    MPI_Comm com = g->get_mpi()->comm();
	    if (MPI_Isend(sbuff, size, Stype, c->ip(), tag, com, &req) != 0) 
	      ERROR("Error in MPI_Isend call");

	    // wait for send to finish, and then delete buffer
	    if (MPI_Wait(&req, &status) != 0)
	      ERROR("Error in MPI_Wait call");
	    delete[] sbuff;
	  
	  } // if overlap_indices
	} // while itgc
      } // if ilevel
    } // while itlg
  } // for ilevel


  // clean up
  delete[] rbuffs;
  delete[] reqs;
  delete[] tags;

  return 0;
}

//----------------------------------------------------------------------

/// Routine to actually perform the prolongation of the overlap data
/** Performs the prolongation of the data in overlap region (that has 
  * size ovsize[3]) into the local data array mydata (of size mysize[3]), 
  * where the overlap data should be placed into the region defined by 
  * g1_ilo[3] and g1_ihi[3]. */
void AMRsolve_HG_prec::do_prolong(Scalar *mydata, Scalar *overlap, 
				    int *mysize, int *g1_ilo, 
				    int *g1_ihi, int *ovsize, 
				    int method) throw()
{

  // ensure that overlap region fits within mydata
  for (int i=0; i<3; i++) 
    if (g1_ihi[i]-g1_ilo[i]+1 > mysize[i])
      ERROR("Error, overlap region larger than local data array");

  // determine refinement factor (should be uniform in all dims)
  int rfac = (g1_ihi[0]-g1_ilo[0]+1)/ovsize[0];
  
  // Perform restriction based on 'method'
  switch (method) {

  case 0:    // piecewise constant

    // Perform prolongation -- iterate over local overlap region (fine cells)
    for (int k=g1_ilo[2]; k<=g1_ihi[2]; k++)
      for (int j=g1_ilo[1]; j<=g1_ihi[1]; j++)
	for (int i=g1_ilo[0]; i<=g1_ihi[0]; i++) {
	  
	  // indices of coarse grid cell covering this cell
	  int k2 = (k-g1_ilo[2])/rfac;   // integer division rounds down, as desired
	  int j2 = (j-g1_ilo[1])/rfac;
	  int i2 = (i-g1_ilo[0])/rfac;
	  
	  // insert coarse grid value into fine cell
	  mydata[(k*mysize[1] + j)*mysize[0] + i] = 
	    overlap[(k2*ovsize[1] + j2)*ovsize[0] + i2];
	}
    
    break;
    
  default: 
      ERROR("Error, prolongation method undefined");
  }
    
  
}

//======================================================================

