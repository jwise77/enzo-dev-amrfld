/// @file      AMRsolve_fld.cpp
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Implementation of the AMRsolve_FLD class

#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>

#include "HYPRE_struct_ls.h"

#include "hdf5.h"

#include "AMRsolve_defs.h"

//----------------------------------------------------------------------

// defines

#define WHERE printf ("%s:%d ",__FILE__,__LINE__);

//----------------------------------------------------------------------

// Constants

const bool debug = false;
const bool trace = false;

//----------------------------------------------------------------------

// Typedefs

typedef int int3[3];

//----------------------------------------------------------------------

#include "AMRsolve_mpi.h"
#include "AMRsolve_scalar.h"
#include "AMRsolve_error.h"
#include "AMRsolve_constants.h"
#include "AMRsolve_faces.h"
#include "AMRsolve_domain.h"
#include "AMRsolve_grid.h"
#include "AMRsolve_level.h"
#include "AMRsolve_hierarchy.h"
#include "AMRsolve_parameters.h"
#include "AMRsolve_HG_prec.h"
#include "AMRsolve_fld.h"
#include "AMRsolve_error.h"


//======================================================================
// PUBLIC MEMBER FUNCTIONS
//======================================================================

/// AMRsolve_FLD constructor
AMRsolve_FLD::AMRsolve_FLD(AMRsolve_Hierarchy& hierarchy, 
				       AMRsolve_Parameters& parameters,
				       int precflag)
  : Ac_(0), Bc_(0), Xc_(0), cgrid_(0), cstencil_(0), 
    parameters_(&parameters), hierarchy_(&hierarchy), resid_(-1.0), 
    iter_(-1), citer_(-1), r_factor_(const_r_factor), Nchem_(-1), theta_(-1.0), 
    dt_(-1.0), aval_(-1.0), aval0_(-1.0), adot_(-1.0), adot0_(-1.0), 
    nUn_(-1.0), nUn0_(-1.0), lUn_(-1.0), lUn0_(-1.0), rUn_(-1.0), rUn0_(-1.0)
{
  // set preconditioner flag
  use_prec = (precflag != 0);

  // set array-valued items
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      BdryType_[i][j] = -1;
}

//----------------------------------------------------------------------

/// AMRsolve_FLD destructor
AMRsolve_FLD::~AMRsolve_FLD()
{
  // destroy HYPRE objects that we created along the way
  int ierr;
  ierr = HYPRE_SStructVectorDestroy(B_);
  if (ierr != 0)  ERROR("could not destroy B_\n");

  ierr = HYPRE_SStructVectorDestroy(X_);
  if (ierr != 0)  ERROR("could not destroy X_\n");

  ierr = HYPRE_SStructMatrixDestroy(A_);
  if (ierr != 0)  ERROR("could not destroy A_\n");

  ierr = HYPRE_SStructGraphDestroy(graph_);
  if (ierr != 0)  ERROR("could not destroy graph_\n");

  ierr = HYPRE_SStructStencilDestroy(stencil_);
  if (ierr != 0)  ERROR("could not destroy stencil_\n");

  ierr = HYPRE_SStructGridDestroy(grid_);
  if (ierr != 0)  ERROR("could not destroy grid_\n");

  if (use_prec) {
    ierr = HYPRE_SStructVectorDestroy(Y_);
    if (ierr != 0)  ERROR("could not destroy Y_\n");

    ierr = HYPRE_StructVectorDestroy(Bc_);
    if (ierr != 0)  ERROR("could not destroy Bc_\n");

    ierr = HYPRE_StructVectorDestroy(Xc_);
    if (ierr != 0)  ERROR("could not destroy Xc_\n");

    ierr = HYPRE_StructMatrixDestroy(Ac_);
    if (ierr != 0)  ERROR("could not destroy Ac_\n");

    ierr = HYPRE_StructStencilDestroy(cstencil_);
    if (ierr != 0)  ERROR("could not destroy cstencil_\n");

    ierr = HYPRE_StructGridDestroy(cgrid_);
    if (ierr != 0)  ERROR("could not destroy cgrid_\n");
  }
}

//----------------------------------------------------------------------

/// Initialize the Grid Hierarchy
/** Creates a hypre grid, with one part per level and one box per 
    Grid patch object, for an AMR problem.  Sets grid box extents, 
    grid part variables, and periodicity. */
void AMRsolve_FLD::init_hierarchy()
{
  if (!use_prec)  return;
  int ierr;
  int dim       = hierarchy_->dimension();

  // Create the hypre grids
  _TRACE_;
  ierr = HYPRE_StructGridCreate(MPI_COMM_WORLD, dim, &cgrid_);
  if (ierr != 0) ERROR("could not create cgrid_\n");

  _TRACE_;
  ItHierarchyLevels itl (*hierarchy_);
  while (AMRsolve_Level* level = itl++) {

    _TRACE_;
    if (level->index() == 0) { // only set up the operator on coarse grids

      // Set extents for boxes that comprise the hypre grid
      ItLevelGridsLocal itgl (*level);
      while (AMRsolve_Grid* grid = itgl++) {
	int lower[3] = {grid->index_lower(0),
			grid->index_lower(1),
			grid->index_lower(2)};
	int upper[3] = {grid->index_upper(0),
			grid->index_upper(1),
			grid->index_upper(2)};
	// set up cgrid for coarse level
	ierr = HYPRE_StructGridSetExtents(cgrid_, lower, upper);
	if (ierr != 0) ERROR("could not set cgrid_ extents\n");
      } // while grid = itgl++

      _TRACE_;
      // Set periodicity of the grid
      int period[3] = { hierarchy_->period_index(0,0),
			hierarchy_->period_index(1,0),
			hierarchy_->period_index(2,0) };
      
      _TRACE_;
      ierr = HYPRE_StructGridSetPeriodic(cgrid_, period);
      if (ierr != 0) ERROR("could not set cgrid_ periodicity\n");
    } // if (level->index() == 0)
  } // while level = itl++

  // When finished, assemble the hypre grids
  _TRACE_;
  ierr = HYPRE_StructGridAssemble(cgrid_);
  if (ierr != 0) ERROR("could not assemble cgrid_\n");
  _TRACE_;
  
} // AMRsolve_FLD::init_hierarchy()

//----------------------------------------------------------------------

/// Initialize the discretization stencils.  
/** Creates and initializes a hypre stencil object. */
void AMRsolve_FLD::init_stencil()
{
  int ierr;
  if (!use_prec) return;

  _TRACE_;
  int dim = hierarchy_->dimension();

  _TRACE_;
  ierr = HYPRE_StructStencilCreate(dim, dim*2+1, &cstencil_);
  if (ierr != 0) ERROR("could not initialize stencil_\n");

  int entries[][3] = { {  0, 0, 0 },     // center
		       {  1, 0, 0 },     // X+
		       { -1, 0, 0 },     // X-
		       {  0, 1, 0 },     // Y+
		       {  0,-1, 0 },     // Y-
		       {  0, 0, 1 },     // Z+
		       {  0, 0,-1 } };   // Z-

  _TRACE_;
  if (dim >= 1) {
    ierr = HYPRE_StructStencilSetElement(cstencil_, 0, entries[0]);
    if (ierr != 0) ERROR("could not set cstencil_ entry\n");
  }
  if (dim >= 1) {
    ierr = HYPRE_StructStencilSetElement(cstencil_, 1, entries[1]);
    if (ierr != 0) ERROR("could not set cstencil_ entry\n");
  }
  if (dim >= 1) {
    ierr = HYPRE_StructStencilSetElement(cstencil_, 2, entries[2]);
    if (ierr != 0) ERROR("could not set cstencil_ entry\n");
  }
  if (dim >= 2) {
    ierr = HYPRE_StructStencilSetElement(cstencil_, 3, entries[3]);
    if (ierr != 0) ERROR("could not set cstencil_ entry\n");
  }
  if (dim >= 2) {
    ierr = HYPRE_StructStencilSetElement(cstencil_, 4, entries[4]);
    if (ierr != 0) ERROR("could not set cstencil_ entry\n");
  }
  if (dim >= 3) {
    ierr = HYPRE_StructStencilSetElement(cstencil_, 5, entries[5]);
    if (ierr != 0) ERROR("could not set cstencil_ entry\n");
  }
  if (dim >= 3) {
    ierr = HYPRE_StructStencilSetElement(cstencil_, 6, entries[6]);
    if (ierr != 0) ERROR("could not set cstencil_ entry\n");
  }
  _TRACE_;

} // AMRsolve_FLD::init_stencil()

//----------------------------------------------------------------------

/// Initialize the matrix A and right-hand-side vector b
/* Creates a matrix with a given nonzero structure, and sets nonzero
   values. */
void AMRsolve_FLD::init_elements(double dt, int Nchem, double theta, 
				 double aval, double aval0, 
				 double adot, double adot0, 
				 double nUn, double nUn0, 
				 double lUn, double lUn0, double rUn, 
				 double rUn0, int BdryType[3][2])
{
  // set input arguments into AMRsolve_FLD object
  dt_        = dt;
  Nchem_     = Nchem;
  theta_     = theta;
  aval_      = aval;
  aval0_     = aval0;
  adot_      = adot;
  adot0_     = adot0;
  nUn_       = nUn;
  nUn0_      = nUn0;
  lUn_       = lUn;
  lUn0_      = lUn0;
  rUn_       = rUn;
  rUn0_      = rUn0;
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      BdryType_[i][j] = BdryType[i][j];
 
  // Set up the coarse grid matrix if requested
  if (!use_prec)  return;

  // Create the hypre matrices and vectors
  int ierr;
  ierr = HYPRE_StructMatrixCreate(MPI_COMM_WORLD, cgrid_, cstencil_,  &Ac_);
  if (ierr != 0) ERROR("could not create Ac_\n");
  ierr = HYPRE_StructVectorCreate(MPI_COMM_WORLD, cgrid_, &Xc_);
  if (ierr != 0) ERROR("could not create Xc_\n");
  ierr = HYPRE_StructVectorCreate(MPI_COMM_WORLD, cgrid_, &Bc_);
  if (ierr != 0) ERROR("could not create Bc_\n");

  // Initialize the hypre matrices and vector objects
  ierr = HYPRE_StructMatrixInitialize(Ac_);
  if (ierr != 0) ERROR("could not initialize Ac_\n");
  ierr = HYPRE_StructVectorInitialize(Xc_);
  if (ierr != 0) ERROR("could not initialize Xc_\n");
  ierr = HYPRE_StructVectorInitialize(Bc_);
  if (ierr != 0) ERROR("could not initialize Bc_\n");

  //--------------------------------------------------
  // Initialize the Ac_ matrix elements
  //--------------------------------------------------
  init_elements_matrix_();

  // Assemble the matrix and vectors
  ierr = HYPRE_StructMatrixAssemble(Ac_);
  if (ierr != 0) ERROR("could not assemble Ac_\n");
  ierr = HYPRE_StructVectorAssemble(Bc_);
  if (ierr != 0) ERROR("could not assemble Ac_\n");
  ierr = HYPRE_StructVectorAssemble(Xc_);
  if (ierr != 0) ERROR("could not assemble Ac_\n");

  // Optionally write the matrices and RHS vector to files for debugging
  if (parameters_->value("dump_ac") == "true") {
    ierr = HYPRE_StructMatrixPrint("Ac-hypre",Ac_,0);
    if (ierr != 0) ERROR("could not print Ac_\n");
  }

} // AMRsolve_FLD::init_elements()

//----------------------------------------------------------------------

/// Initialize and solve the linear system
void AMRsolve_FLD::solve()
{
  int    ierr;
  int    itmax  = 0;
  double restol = 0.0;

  // Check solver parameters
  std::string sitmax  = parameters_->value("solver_itmax");
  std::string srestol = parameters_->value("solver_restol");

  // If not defined, then define them
  if (sitmax == "")  parameters_->add_parameter("solver_itmax","200");
  if (srestol == "") parameters_->add_parameter("solver_restol","1e-6");

  // recheck solver parameters
  sitmax  = parameters_->value("solver_itmax");
  srestol = parameters_->value("solver_restol");

  // Set local variables
  itmax  = atoi(sitmax.c_str());
  restol = atof(srestol.c_str());

  // call solver
  solve_bicgstab_(itmax,restol);

} // AMRsolve_FLD::solve()

//----------------------------------------------------------------------

/// Evaluate the success of the solve, return values (0=success, 1=failure)
int AMRsolve_FLD::evaluate()
{

  // output solution/rhs output if requested
  if (parameters_->value("dump_x") == "true" || 
      parameters_->value("dump_b") == "true") {

    // iterate over processor-local grids
    ItHierarchyGridsLocal itg(*hierarchy_);
    while (AMRsolve_Grid* grid = itg++) {

      char filename[80];

      // save solution
      if (parameters_->value("dump_x") == "true") {
	sprintf(filename,"X.%d",grid->id());
	grid->write("header",filename);
	grid->write("u",filename);
      }
      
      // save rhs
      if (parameters_->value("dump_b") == "true") {
	sprintf(filename,"B.%d",grid->id());
	grid->write("header",filename);
	grid->write("f",filename);
      }
      
    } // grid = itg++
  } // if dump_x or dump_b


  // check for error flags; output info to stdout, if necessary
  int err_flag = 0;

  // Residual too high
  double restol = atof(parameters_->value("solver_restol").c_str());
  if (resid_ > restol) {
    if (pmpi->is_root()) 
      printf("Diverged: %g > %g\n", resid_,restol);
    err_flag = 1;
  }

  // Iterations reached limit
  int itmax = atoi(parameters_->value("solver_itmax").c_str());
  if (iter_ >= itmax) {
    if (pmpi->is_root()) 
      printf("Stalled: %d >= %d\n", iter_,itmax);
    err_flag = 1;
  }

  return err_flag;

} // AMRsolve_FLD::evaluate()

//----------------------------------------------------------------------

/// Extracts HYPRE solution and updates Enzo radiation field
void AMRsolve_FLD::update_enzo()
{
  int ierr;

  // iterate over grids on this processor
  ItHierarchyGridsLocal itg(*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // access solution
    int n3[3];
    Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);

    // access Enzo radiation field
    Scalar* E = grid->get_E();

    // get buffering information on relating amrsolve grid to Enzo data
    int ghosts[3][2]; 
    grid->get_Ghosts(ghosts);
    int en0 = n3[0] + ghosts[0][0] + ghosts[0][1];
    int en1 = n3[1] + ghosts[1][0] + ghosts[1][1];
    int en2 = n3[2] + ghosts[2][0] + ghosts[2][1];

    // update Enzo data with amrsolve solution
    int k0, k1, k2, k, i0, i1, i2, i;
    for (i2=0; i2<n3[2]; i2++) {
      k2 = ghosts[2][0] + i2;

      for (i1=0; i1<n3[1]; i1++) {
	k1 = ghosts[1][0] + i1;

	for (i0=0; i0<n3[0]; i0++) {
	  k0 = ghosts[0][0] + i0;

	  // compute indices of amrsolve, enzo cells
	  i = i0 + n3[0]*(i1 + n3[1]*i2);
	  k = k0 + en0*(k1 + en1*k2);

	  // update Enzo solution with amrsolve solution (corrector)
	  E[k] += u[i];
	}
      }
    }

  } // grid = itg++

} // AMRsolve_FLD::update_enzo()


//----------------------------------------------------------------------

/// dumps HYPRE matrix and RHS (called when aborting solve)
void AMRsolve_FLD::abort_dump()
{
  int ierr;

  // have HYPRE dump out coarse matrix to disk
  if (use_prec) {
    ierr = HYPRE_StructMatrixPrint("Ac.mat",Ac_,0);
    if (ierr != 0) ERROR("could not print Ac_\n");
  }

} // AMRsolve_FLD::abort_dump()


//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

/// init_nonstencil_() is called twice: first by
/// init_graph_nonstencil_() with phase == phase_graph to set nonstencil
/// graph entries, and again by init_matrix_nonstencil_() with phase
/// == phase_matrix to set nonstencil matrix entries.
void AMRsolve_FLD::init_nonstencil_(AMRsolve_Grid& grid, phase_enum phase)
{
  // Input parameter check
  if ( !(phase == phase_graph || phase == phase_matrix) ) {
    char error_message[80];
    sprintf(error_message, "init_matrix_nonstencil_ called with phase = %d", int(phase));
    ERROR(error_message);
  } // if phase unexpected

  // global grid index limits
  int index_global[3][2];
  grid.indices(index_global);
  
  // Determine the discretization method (constant, linear, quadratic)
  // ***Currently only constant is implemented***
  enum discret_type_enum {discret_type_unknown, discret_type_const};
  discret_type_enum discret_type;

  if (parameters_->value("discret") == "constant") 
    discret_type = discret_type_const;

  // Loop over each face zone in the grid, adding non-stencil graph
  // entries wherever a zone is adjacent to a coarse zone.  Both
  // fine-to-coarse and coarse-to-fine entries are added.
  int level_fine   = grid.level();
  int level_coarse = grid.level() - 1;
  assert (level_coarse >= 0);

  for (int axis=0; axis<3; axis++) {

    // axis0:       axis normal to face
    // axis1,axis2: axes within face
    int axis0 = axis;
    int axis1 = (axis+1)%3;
    int axis2 = (axis+2)%3;

    // n0     grid size normal to face
    // n1,n2: grid size within face
    int n0 = index_global[axis0][1] - index_global[axis0][0];
    int n1 = index_global[axis1][1] - index_global[axis1][0];
    int n2 = index_global[axis2][1] - index_global[axis2][0];

    // index_global[][] should be divisible by r_factor_**level.  Just
    // test r_factor_ here.
    bool l0 = (index_global[axis1][0]/r_factor_)*r_factor_ == index_global[axis1][0];
    bool l1 = (index_global[axis1][1]/r_factor_)*r_factor_ == index_global[axis1][1];
    if (!l0) printf("grid %i,  index_global[%d][0] = %d,  r_factor = %i\n",
		    grid.id(),axis1,index_global[axis1][0], r_factor_);
    assert(l0);
    if (!l1) printf("grid %i,  index_global[%d][1] = %d,  r_factor = %i\n",
		    grid.id(),axis1,index_global[axis1][1], r_factor_);
    assert(l1);

    // loop over faces orthogonal to this axis
    for (int face=0; face<2; face++) {

      // Loop over face zones that are aligned with coarse zones (hence "+= r")
      for (int index1=0; index1<n1; index1 += r_factor_) {
	for (int index2=0; index2<n2; index2 += r_factor_) {

	  AMRsolve_Grid* adjacent = grid.faces().adjacent(axis,face,index1,index2);

	  AMRsolve_Faces::Label& fz = grid.faces().label(axis,face,index1,index2);

	  // Add graph entries iff grid or adjacent grid is local, and if
	  // adjacent grid (if it exists) is in the next-coarser level
	  bool is_local =
	    ((adjacent != NULL) && (adjacent->is_local() || grid.is_local()));
	  bool is_coarse = (fz == AMRsolve_Faces::_coarse_);

	  if (is_local && is_coarse) {

	    // (fine) grid global indices
	    int index_fine[3]; 
	    index_fine[axis0] = index_global[axis0][0] + face*(n0 - r_factor_);
	    index_fine[axis1] = index_global[axis1][0] + index1;
	    index_fine[axis2] = index_global[axis2][0] + index2;

	    // (coarse) adjacent global indices
	    int index_coarse[3]; 
	    index_coarse[axis0] = (index_fine[axis0]) / r_factor_ + (face*r_factor_-1);
	    index_coarse[axis1] = (index_fine[axis1]) / r_factor_;
	    index_coarse[axis2] = (index_fine[axis2]) / r_factor_;

	    // adjust for periodicity
	    if (hierarchy_->is_periodic(axis0)) {
	      int period = hierarchy_->period_index(axis0,level_coarse);
 	      index_coarse[axis0] = (index_coarse[axis0] + period) % period;
 	    }

	    //-------------------------------------------------------
	    // GRAPH ENTRY: FINE-TO-COARSE (ADJUSTS FINE GRID MATRIX)
	    //-------------------------------------------------------

	    if (discret_type == discret_type_const) {

	      update_fine_coarse_const_(face,grid,axis0,phase,
					level_fine,level_coarse,
					index_fine,index_coarse);

	      //-------------------------------------------------------
	      // GRAPH ENTRY: COARSE-TO-FINE (ADJUSTS FINE GRID MATRIX)
	      //-------------------------------------------------------
	      
	      if (adjacent->is_local()) {

		update_coarse_fine_const_(face,*adjacent,axis0,phase,
					  level_fine,level_coarse,
					  index_fine,index_coarse);

	      }

	    } else {
	      char error_message[80];
	      strcpy(error_message,"Unknown parameter discret = ");
	      strcat(error_message,parameters_->value("discret").c_str());
	      ERROR(error_message);
	    } // if discret unexpected
	  } // if is_local && fz == AMRsolve_Faces::_coarse_
	} // for index2
      } // for index1
    } // for face
  } // for axis
} // AMRsolve_FLD::init_nonstencil_()

//------------------------------------------------------------------------

/// Initialize matrix stencil and graph entries
void AMRsolve_FLD::init_elements_matrix_()
{

  ItHierarchyLevels itl (*hierarchy_);
  while (AMRsolve_Level* level = itl++) {

    int part = level->index();

    // 1. Set stencil values within level
    ItLevelGridsLocal itlg (*level);
    while (AMRsolve_Grid* grid = itlg++)  init_matrix_stencil_(*grid);

    if (part > 0) {
      // *** WARNING: POSSIBLE SCALING ISSUE.  Below we loop over all
      // *** grids; however, we only need to loop over "parent-child
      // *** pairs such that either child or parent is local to this
      // *** MPI process."
 
      // Set matrix values between levels
      ItLevelGridsAll itag (*level);
      while (AMRsolve_Grid* grid = itag++)  init_matrix_nonstencil_(*grid);

    } // while level > 0
  } // while level = itl++

  // Clean up stencil connections between levels
  for (int part=1; part<hierarchy_->num_levels(); part++) 
    init_matrix_clear_(part);

} // init_elements_matrix_()

//------------------------------------------------------------------------

/// Utility function to handle FLD limiter
Scalar AMRsolve_FLD::limiter_(Scalar E1, Scalar E2, Scalar k1, Scalar k2, 
				    Scalar nUn, Scalar lUn, Scalar dxi) 
{
  Scalar c = 2.99792458e10;
  Scalar Rmin = 1.0e-2 / lUn_;
  //  Rmin = MIN(Rmin, 1.e-20);   // 1st is astro/cosmo, 2nd is lab frame
  Scalar Emin = 1.0e-30;
  Scalar Dmax = 2.0539e-3 * c * lUn_;
  // Dmax = MAX(Dmax, 1.e20);     // 1st is astro/cosmo, 2nd is lab frame

  Scalar Eavg = MAX((E1 + E2)*0.5, Emin);
  Scalar kap = 2.0*k1*k2/(k1+k2)*nUn_;
  Scalar R = MAX(dxi*ABS(E1 - E2)/Eavg, Rmin);
  return MIN(c/sqrt(9.0*kap*kap + R*R), Dmax);

} // limiter_()

//------------------------------------------------------------------------

/// Set right-hand-side elements based on values in grids
void AMRsolve_FLD::init_elements_rhs_()
{
  // declare shortcut variables
  Scalar Ed_zl, Ed_yl, Ed_xl, Ed_xr, Ed_yr, Ed_zr, R, R0, kap, kap0;
  Scalar D_zl, D_yl, D_xl, D_xr, D_yr, D_zr;
  Scalar D0_zl, D0_yl, D0_xl, D0_xr, D0_yr, D0_zr;
  Scalar afac  = adot_  / aval_;
  Scalar afac0 = adot0_ / aval0_;
  Scalar dtfac  = dt_ * theta_;
  Scalar dtfac0 = dt_ * (1.0 - theta_);
  Scalar c = 2.99792458e10;
  int i0, i1, i2, i, ierr;

  _TRACE_;
  // iterate over all grids local to this processor
  ItHierarchyGridsLocal itg(*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // for each local grid, use our "f" array to store RHS entries 
    // for the linear system.
    int n0,n1,n2;
    Scalar* values = grid->get_f(&n0,&n1,&n2);
    for (i=0; i<n0*n1*n2; i++)  values[i] = 0.0;

    // access relevant arrays from this grid to compute RHS
    Scalar* E     = grid->get_E();
    Scalar* eta   = grid->get_eta();
    Scalar* kappa = grid->get_kap();

    // get buffering information on relating amrsolve grid to Enzo data
    int ghosts[3][2]; 
    grid->get_Ghosts(ghosts);
    int en0 = n0 + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
    int en1 = n1 + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
    int en2 = n2 + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions

    // access this grid's mesh spacing, and set relevant shortcuts
    Scalar hx,hy,hz;
    grid->h(hx,hy,hz);
    Scalar dxi  = 1.0 / hx / lUn_;
    Scalar dyi  = 1.0 / hy / lUn_;
    Scalar dzi  = 1.0 / hz / lUn_;
    Scalar dxi0 = 1.0 / hx / lUn0_;
    Scalar dyi0 = 1.0 / hy / lUn0_;
    Scalar dzi0 = 1.0 / hz / lUn0_;
    Scalar dxfac = dtfac * dxi * dxi;
    Scalar dyfac = dtfac * dyi * dyi;
    Scalar dzfac = dtfac * dzi * dzi;
    Scalar dxfac0 = dtfac0 * dxi0 * dxi0;
    Scalar dyfac0 = dtfac0 * dyi0 * dyi0;
    Scalar dzfac0 = dtfac0 * dzi0 * dzi0;

    // iterate over the domain, computing the RHS array
    int k2, k1, k0, k_l00, k_r00, k_0l0, k_0r0, k_00l, k_00r, k_000;
    for (i2=0; i2<n2; i2++) {
      k2 = ghosts[2][0] + i2;

      for (i1=0; i1<n1; i1++) {
	k1 = ghosts[1][0] + i1;

	for (i0=0; i0<n0; i0++) {
	  k0 = ghosts[0][0] + i0;

	  // compute indices of neighboring Enzo cells
	  k_l00 = k0-1 + en0*(k1   + en1*k2);
	  k_0l0 = k0   + en0*(k1-1 + en1*k2);
	  k_00l = k0   + en0*(k1   + en1*(k2-1));
	  k_000 = k0   + en0*(k1   + en1*k2);
	  k_r00 = k0+1 + en0*(k1   + en1*k2);
	  k_0r0 = k0   + en0*(k1+1 + en1*k2);
	  k_00r = k0   + en0*(k1   + en1*(k2+1));

	  //--------------
	  // z-directional limiter, lower face
	  Ed_zl = E[k_000] - E[k_00l];
	  D_zl  = limiter_(E[k_000], E[k_00l], kappa[k_000], 
			   kappa[k_00l], nUn_, lUn_, dzi);
	  D0_zl = limiter_(E[k_000], E[k_00l], kappa[k_000], 
			   kappa[k_00l], nUn0_, lUn0_, dzi0);

	  //--------------
	  // y-directional limiter, lower face
	  Ed_yl = E[k_000] - E[k_0l0];
	  D_yl  = limiter_(E[k_000], E[k_0l0], kappa[k_000], 
			   kappa[k_0l0], nUn_, lUn_, dyi);
	  D0_yl = limiter_(E[k_000], E[k_0l0], kappa[k_000], 
			   kappa[k_0l0], nUn0_, lUn0_, dyi0);

	  //--------------
	  // x-directional limiter, lower face
	  Ed_xl = E[k_000] - E[k_l00];
	  D_xl  = limiter_(E[k_000], E[k_l00], kappa[k_000],
			   kappa[k_l00], nUn_, lUn_, dxi);
	  D0_xl = limiter_(E[k_000], E[k_l00], kappa[k_000],
			   kappa[k_l00], nUn0_, lUn0_, dxi0);

	  //--------------
	  // x-directional limiter, upper face
	  Ed_xr = E[k_r00] - E[k_000];
	  D_xr  = limiter_(E[k_r00], E[k_000], kappa[k_000],
			   kappa[k_r00], nUn_, lUn_, dxi);
	  D0_xr = limiter_(E[k_r00], E[k_000], kappa[k_000],
			   kappa[k_r00], nUn0_, lUn0_, dxi0);

	  //--------------
	  // y-directional limiter, upper face
	  Ed_yr = E[k_0r0] - E[k_000];
	  D_yr  = limiter_(E[k_0r0], E[k_000], kappa[k_000],
			   kappa[k_0r0], nUn_, lUn_, dyi);
	  D0_yr = limiter_(E[k_0r0], E[k_000], kappa[k_000],
			   kappa[k_0r0], nUn0_, lUn0_, dyi0);

	  //--------------
	  // z-directional limiter, upper face
	  Ed_zr = E[k_00r] - E[k_000];
	  D_zr  = limiter_(E[k_00r], E[k_000], kappa[k_000],
			   kappa[k_00r], nUn_, lUn_, dzi);
	  D0_zr = limiter_(E[k_00r], E[k_000], kappa[k_000],
			   kappa[k_00r], nUn0_, lUn0_, dzi0);

	  // opacity values in this cell
	  kap = kappa[k_000];
	  kap0 = kap*nUn0_;
	  kap *= nUn_;

	  // set rhs in this cell
	  i = i0 + n0*(i1 + n1*i2);
	  values[i] = ( (dtfac/rUn_ + dtfac0/rUn0_)*eta[k_000]
		      + (1.0 - dtfac0*(afac0+c*kap0))*E[k_000]
		      + dxfac0*(D0_xr*Ed_xr - D0_xl*Ed_xl)
		      + dyfac0*(D0_yr*Ed_yr - D0_yl*Ed_yl)
		      + dzfac0*(D0_zr*Ed_zr - D0_zl*Ed_zl)
		      - (1.0 + dtfac*(afac+c*kap))*E[k_000]
		      + dxfac*(D_xr*Ed_xr - D_xl*Ed_xl)
		      + dyfac*(D_yr*Ed_yr - D_yl*Ed_yl)
		      + dzfac*(D_zr*Ed_zr - D_zl*Ed_zl) );

	  // check that value is legal, otherwise output an error message
	  if (isinf(values[i]) || isnan(values[i])) {
	    fprintf(stderr,"init_elements_rhs_ ERROR: encountered illegal value (%g)\n   eta = %g, E = %g, dtfac = %g, dtfac0 = %g, afac = %g, kap = %g\n   D* = %g, %g, %g, %g, %g, %g\n   Ed* = %g %g %g %g %g %g\n\n",
		    values[i], eta[k_000], E[k_000], dtfac, dtfac0, afac, kap, D_xl, D_xr, D_yl, D_yr, D_zl, D_zr, Ed_xl, Ed_xr, Ed_yl, Ed_yr, Ed_zl, Ed_zr);
	    ERROR("NaN error in init_elements_rhs_\n");
	  }

	}
      }
    }

    // Set Hypre B_ vector to RHS values
    int part = grid->level();
    int lower[3] = { grid->index_lower(0), 
		     grid->index_lower(1), 
		     grid->index_lower(2) };
    int upper[3] = { grid->index_upper(0), 
		     grid->index_upper(1), 
		     grid->index_upper(2) };
    ierr = HYPRE_SStructVectorAddToBoxValues(B_, part, lower, 
					     upper, 0, values);
    if (ierr != 0) ERROR("could not AddToBoxValues in B_\n");

  }  // while grid = itg++

  // call HYPRE to zero out overlapped cell values
  int nlevels = hierarchy_->num_levels();
  if (nlevels > 1) {
    int plevels[nlevels];
    int rfactors[nlevels][3];
    for (int part=0; part<nlevels; part++) {
      plevels[part] = part;
      rfactors[part][0] = r_factor_;
      rfactors[part][1] = r_factor_;
      rfactors[part][2] = r_factor_;
    }
    ierr = HYPRE_SStructFACZeroAMRVectorData(B_, plevels, rfactors);
    if (ierr != 0) ERROR("could not ZeroAMRVectorData in B_\n");
  }
  

} // init_elements_rhs_()

//------------------------------------------------------------------------

/// Compute a relative difference norm between E and E0, using the pnorm
/// and absolute tolerances given as input.  This routine ensures that 
/// for (pnorm != 0), the vector values are scaled by the cell volume, 
/// and that overlapped cells are not double-counted.
double AMRsolve_FLD::rdiff_norm(double pnorm, double atol)
{

  // local variables
  int n0, n1, n2, en0, en1, en2, ghosts[3][2], part, lower[3], upper[3];
  int i, i0, i1, i2, k, ierr;
  Scalar *E, *E0, *values;
  Scalar hx, hy, hz, dV, weight, diff;
  double tmp, loc_diff, proc_diff=0.0, glob_diff=0.0;

  _TRACE_;
  // for each local grid, fill the "b" array to store difference values 
  //    iterate over all grids local to this processor
  ItHierarchyGridsLocal itg(*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // access AMRsolve array
    values = grid->get_f(&n0,&n1,&n2);

    // access relevant arrays from this grid to compute RHS
    E  = grid->get_E();
    E0 = grid->get_E0();

    // get buffering information on relating amrsolve grid to Enzo data
    grid->get_Ghosts(ghosts);
    en0 = n0 + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
    en1 = n1 + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
    en2 = n2 + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions

    // iterate over the domain, computing the difference array
    for (i2=0; i2<n2; i2++)
      for (i1=0; i1<n1; i1++)
	for (i0=0; i0<n0; i0++) {

	  // compute Enzo cell index (k), AMRsolve grid index (i)
	  k = (ghosts[0][0] + i0) + en0*((ghosts[1][0] + i1) + en1*(ghosts[2][0] + i2));
	  i = i0 + n0*(i1 + n1*i2);

	  // compute local weight factor
	  weight = sqrt(fabs(E[k]*E0[k])) + atol;

	  // compute local weighted difference, and store inside b
	  values[i] = fabs(E[k] - E0[k])/weight;
	}

    // Set Hypre B_ vector to RHS values
    part = grid->level();
    grid->get_limits(lower,upper);
    ierr = HYPRE_SStructVectorSetBoxValues(B_, part, lower, 
    					   upper, 0, values);
    if (ierr != 0) ERROR("could not SetBoxValues in B_\n");
    
  }  // while grid = itg++

  // call HYPRE to zero out overlapped cell values
  int nlevels = hierarchy_->num_levels();
  if (nlevels > 1) {
    int plevels[nlevels];
    int rfactors[nlevels][3];
    for (int part=0; part<nlevels; part++) {
      plevels[part] = part;
      rfactors[part][0] = r_factor_;
      rfactors[part][1] = r_factor_;
      rfactors[part][2] = r_factor_;
    }
    ierr = HYPRE_SStructFACZeroAMRVectorData(B_, plevels, rfactors);
    if (ierr != 0) ERROR("could not ZeroAMRVectorData in B_\n");
  }

  
  // iterate back over local grids, accumulating the difference norm
  ItHierarchyGridsLocal itg2(*hierarchy_);
  while (AMRsolve_Grid* grid = itg2++) {

    // get level, grid information
    int level = grid->level();
    grid->get_limits(lower,upper);

    // extract Enzo solution
    values = grid->get_f(&n0,&n1,&n2);
    ierr = HYPRE_SStructVectorGather(B_);
    if (ierr != 0) ERROR("could not gather B_\n");
    ierr = HYPRE_SStructVectorGetBoxValues(B_, level, lower, 
    					   upper, 0, values);
    if (ierr != 0) ERROR("could not GetBoxValues in B_\n");

    // access this grid's mesh spacing, and set relevant shortcuts
    grid->h(hx,hy,hz);
    dV = (pnorm > 0.0) ? hx*hy*hz : 1.0;

    // iterate over the domain, updating local p-norm 
    loc_diff = 0.0;
    if (pnorm > 0.0) { 
      for (i2=0; i2<n2; i2++)
	for (i1=0; i1<n1; i1++)
	  for (i0=0; i0<n0; i0++) {
	    tmp = values[i0 + n0*(i1 + n1*i2)];
	    loc_diff += pow(tmp,pnorm)*dV;
	  }
      proc_diff += loc_diff;
    // iterate over the domain, updating local max norm 
    } else { 
      for (i2=0; i2<n2; i2++)
	for (i1=0; i1<n1; i1++)
	  for (i0=0; i0<n0; i0++) {
	    tmp = values[i0 + n0*(i1 + n1*i2)];
	    loc_diff = (loc_diff > tmp) ? loc_diff : tmp;
	  }
      proc_diff = (loc_diff > proc_diff) ? loc_diff : proc_diff;
    }
   
  } // grid = itg2++

  // communicate to obtain overall sum/max
  if (pnorm > 0.0) {
    if (MPI_Allreduce(&proc_diff,&glob_diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD) != 0) {
      char error_message[80];
      strcpy(error_message,"rdiff_norm: Error in MPI_Allreduce!\n");
      ERROR(error_message);
    }
    glob_diff = pow(glob_diff, 1.0/pnorm);
  } else {
    if (MPI_Allreduce(&proc_diff,&glob_diff,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) != 0) {
      char error_message[80];
      strcpy(error_message,"rdiff_norm: Error in MPI_Allreduce!\n");
      ERROR(error_message);
    }
  }

  // return with overall difference value
  return glob_diff;

} // rdiff_norm()

//------------------------------------------------------------------------

/// Set matrix stencil values for this grid's interior
void AMRsolve_FLD::init_matrix_stencil_(AMRsolve_Grid& grid)
{
  int n          = grid.num_unknowns();
  int entries[7] = { 0,1,2,3,4,5,6 };
  int n3[3]      = {grid.n(0),grid.n(1),grid.n(2)};
  int ierr;

  double* v0;         // Diagonal elements
  double* v1[3][2];   // Off-diagonal elements
  double  vtmp;
  int     badvalue;

  // declare shortcut variables
  double Ed_zl, Ed_yl, Ed_xl, Ed_xr, Ed_yr, Ed_zr, kap;
  double D_zl, D_yl, D_xl, D_xr, D_yr, D_zr;
  double afac = adot_ / aval_;
  double dtfac = dt_ * theta_;
  double c = 2.99792458e10;

  // access relevant arrays from this grid to compute RHS
  Scalar* E     = grid.get_E();
  Scalar* kappa = grid.get_kap();
  
  // get buffering information on relating amrsolve grid to Enzo data
  int ghosts[3][2]; 
  grid.get_Ghosts(ghosts);
  int en0 = n3[0] + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
  int en1 = n3[1] + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
  int en2 = n3[2] + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions

  // access this grid's mesh spacing, and set relevant shortcuts
  double dxi = 1.0 / grid.h(0) / lUn_;
  double dyi = 1.0 / grid.h(1) / lUn_;
  double dzi = 1.0 / grid.h(2) / lUn_;
  double dxfac = dtfac * dxi * dxi;
  double dyfac = dtfac * dyi * dyi;
  double dzfac = dtfac * dzi * dzi;

  // determine whether grid is on coarse level
  bool coarsegrid = (grid.level() == 0);
  
  // Allocate storage
  v0 = new double[n];
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      v1[axis][face] = new double[n];
    } // for face
  } // for axis

  //-----------------------------------------------------------
  // Set stencil for all unknowns, ignoring boundary conditions
  //-----------------------------------------------------------

  int k2, k1, k0, k_l00, k_r00, k_0l0, k_0r0, k_00l, k_00r, k_000;
  int i0, i1, i2, i;
  for (i2=0; i2<n3[2]; i2++) {
    k2 = ghosts[2][0] + i2;

    for (i1=0; i1<n3[1]; i1++) {
      k1 = ghosts[1][0] + i1;

      for (i0=0; i0<n3[0]; i0++) {
	k0 = ghosts[0][0] + i0;
	  
	// compute indices of neighboring Enzo cells
	k_l00 = k0-1 + en0*(k1   + en1*k2);
	k_0l0 = k0   + en0*(k1-1 + en1*k2);
	k_00l = k0   + en0*(k1   + en1*(k2-1));
	k_000 = k0   + en0*(k1   + en1*k2);
	k_r00 = k0+1 + en0*(k1   + en1*k2);
	k_0r0 = k0   + en0*(k1+1 + en1*k2);
	k_00r = k0   + en0*(k1   + en1*(k2+1));

	//--------------
	// z-directional limiter, lower face
	Ed_zl = E[k_000] - E[k_00l];
	D_zl = limiter_(E[k_000], E[k_00l], kappa[k_000],
			kappa[k_00l], nUn_, lUn_, dzi);
	
	//--------------
	// y-directional limiter, lower face
	Ed_yl = E[k_000] - E[k_0l0];
	D_yl = limiter_(E[k_000], E[k_0l0], kappa[k_000],
			kappa[k_0l0], nUn_, lUn_, dyi);
	
	//--------------
	// x-directional limiter, lower face
	Ed_xl = E[k_000] - E[k_l00];
	D_xl = limiter_(E[k_000], E[k_l00], kappa[k_000],
			kappa[k_l00], nUn_, lUn_, dxi);
	
	//--------------
	// x-directional limiter, upper face
	Ed_xr = E[k_r00] - E[k_000];
	D_xr = limiter_(E[k_000], E[k_r00], kappa[k_000],
			kappa[k_r00], nUn_, lUn_, dxi);
	
	//--------------
	// y-directional limiter, upper face
	Ed_yr = E[k_0r0] - E[k_000];
	D_yr = limiter_(E[k_000], E[k_0r0], kappa[k_000],
			kappa[k_0r0], nUn_, lUn_, dyi);
	
	//--------------
	// z-directional limiter, upper face
	Ed_zr = E[k_00r] - E[k_000];
	D_zr = limiter_(E[k_000], E[k_00r], kappa[k_000],
			kappa[k_00r], nUn_, lUn_, dzi);
	
	// opacity values in this cell
	kap = kappa[k_000]*nUn_;
	
	// get the linear grid index
	i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);

	// set the matrix entries
	v1[2][0][i] = -dzfac*D_zl;    // z-left
	v1[1][0][i] = -dyfac*D_yl;    // y-left
	v1[0][0][i] = -dxfac*D_xl;    // x-left
	v1[0][1][i] = -dxfac*D_xr;    // x-right
	v1[1][1][i] = -dyfac*D_yr;    // y-right
	v1[2][1][i] = -dzfac*D_zr;    // z-right
	v0[i] = 1.0 + dtfac*(afac + c*kap) + dxfac*(D_xl+D_xr) 
		    + dyfac*(D_yl+D_yr) + dzfac*(D_zl+D_zr);  // self

	// check that value is legal, if not issue an error message
	badvalue = 0;
	vtmp = v1[2][0][i];
	if (isinf(vtmp) || isnan(vtmp)) {
	  fprintf(stderr,"init_matrix_stencil_ ERROR: illegal value (%g)\n   dtfac = %g, dzi = %g, D_zl = %g, i* = %i %i %i %i\n\n", vtmp, dtfac, dzi, D_zl, i, i0, i1, i2);
	  badvalue = 1;
	}
	vtmp = v1[1][0][i];
	if (isinf(vtmp) || isnan(vtmp)) {
	  fprintf(stderr,"init_matrix_stencil_ ERROR: illegal value (%g)\n   dtfac = %g, dyi = %g, D_yl = %g, i* = %i %i %i %i\n\n", vtmp, dtfac, dyi, D_yl, i, i0, i1, i2);
	  badvalue = 1;
	}
	vtmp = v1[0][0][i];
	if (isinf(vtmp) || isnan(vtmp)) {
	  fprintf(stderr,"init_matrix_stencil_ ERROR: illegal value (%g)\n   dtfac = %g, dxi = %g, D_xl = %g, i* = %i %i %i %i\n\n", vtmp, dtfac, dxi, D_xl, i, i0, i1, i2);
	  badvalue = 1;
	}
	vtmp = v1[0][1][i];
	if (isinf(vtmp) || isnan(vtmp)) {
	  fprintf(stderr,"init_matrix_stencil_ ERROR: illegal value (%g)\n   dtfac = %g, dxi = %g, D_xr = %g, i* = %i %i %i %i\n\n", vtmp, dtfac, dxi, D_xr, i, i0, i1, i2);
	  badvalue = 1;
	}
	vtmp = v1[1][1][i];
	if (isinf(vtmp) || isnan(vtmp)) {
	  fprintf(stderr,"init_matrix_stencil_ ERROR: illegal value (%g)\n   dtfac = %g, dyi = %g, D_yr = %g, i* = %i %i %i %i\n\n", vtmp, dtfac, dyi, D_yr, i, i0, i1, i2);
	  badvalue = 1;
	}
	vtmp = v1[2][1][i];
	if (isinf(vtmp) || isnan(vtmp)) {
	  fprintf(stderr,"init_matrix_stencil_ ERROR: illegal value (%g)\n   dtfac = %g, dzi = %g, D_zr = %g, i* = %i %i %i %i\n\n", vtmp, dtfac, dzi, D_zr, i, i0, i1, i2);
	  badvalue = 1;
	}
	vtmp = v0[i];
	if (isinf(vtmp) || isnan(vtmp)) {
	  fprintf(stderr,"init_matrix_stencil_ ERROR: illegal value (%g)\n   dtfac = %g, afac = %g, kap = %g, d*i = %g %g %g, D* = %g %g %g %g %g %g, i* = %i %i %i %i \n\n", vtmp, dtfac, afac, kap, dxi, dyi, dzi, D_xl, D_xr, D_yl, D_yr, D_zl, D_zr, i, i0, i1, i2);
	  badvalue = 1;
	}
	if (badvalue == 1)
	  ERROR("init_matrix_stencil failure!\n")

      } // for i0
    } // for i1
  } // for i2


  //-----------------------------------------------------------
  // Adjust stencil at grid boundaries
  //-----------------------------------------------------------

  // update matrix/rhs based on boundary conditions/location
  //    z-left face
  if (grid.faces().label(2,0,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[2][0] == 1) {           // Dirichlet
      i2 = 0;
      for (i1=0; i1<n3[1]; i1++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[2][0][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[2][0] == 2) {    // Neumann
      i2 = 0;
      for (i1=0; i1<n3[1]; i1++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[2][0][i];
	  v1[2][0][i] = 0.0;
	} 
      }
    }  // BdryType_ == 2
  }  // label == boundary

  //    y-left face
  if (grid.faces().label(1,0,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[1][0] == 1) {           // Dirichlet
      i1 = 0;
      for (i2=0; i2<n3[2]; i2++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[1][0][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[1][0] == 2) {    // Neumann
      i1 = 0;
      for (i2=0; i2<n3[2]; i2++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[1][0][i];
	  v1[1][0][i] = 0.0;
	}
      }
    }  // BdryType_ == 2
  }  // label == boundary

  //    x-left face
  if (grid.faces().label(0,0,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[0][0] == 1) {           // Dirichlet
      i0 = 0;
      for (i2=0; i2<n3[2]; i2++) {
	for (i1=0; i1<n3[1]; i1++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[0][0][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[0][0] == 2) {    // Neumann
      i0 = 0;
      for (i2=0; i2<n3[2]; i2++) {
	for (i1=0; i1<n3[1]; i1++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[0][0][i];
	  v1[0][0][i] = 0.0;
	}
      }
    }  // BdryType_ == 2
  }  // label == boundary

  //    x-right face
  if (grid.faces().label(0,1,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[0][1] == 1) {           // Dirichlet
      i0 = n3[0]-1;
      for (i2=0; i2<n3[2]; i2++) {
	for (i1=0; i1<n3[1]; i1++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[0][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[0][1] == 2) {    // Neumann
      i0 = n3[0]-1;
      for (i2=0; i2<n3[2]; i2++) {
	for (i1=0; i1<n3[1]; i1++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[0][1][i];
	  v1[0][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 2
  }  // label == boundary

  //    y-right face
  if (grid.faces().label(1,1,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[1][1] == 1) {           // Dirichlet
      i1 = n3[1]-1;
      for (i2=0; i2<n3[2]; i2++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[1][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[1][1] == 2) {    // Neumann
      i1 = n3[1]-1;
      for (i2=0; i2<n3[2]; i2++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[1][1][i];
	  v1[1][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 2
  }  // label == boundary

  //    z-right face
  if (grid.faces().label(2,1,0,0) == AMRsolve_Faces::_boundary_) {
    if (BdryType_[2][1] == 1) {           // Dirichlet
      i2 = n3[2]-1;
      for (i1=0; i1<n3[1]; i1++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[2][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 1
    else if (BdryType_[2][1] == 2) {    // Neumann
      i2 = n3[2]-1;
      for (i1=0; i1<n3[1]; i1++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[2][1][i];
	  v1[2][1][i] = 0.0;
	}
      }
    }  // BdryType_ == 2
  }  // label == boundary
  

  //-----------------------------------------------------------
  // insert matrix entries into Hypre matrices A_ and Ac_
  //-----------------------------------------------------------

  //   AMRsolve_Faces& faces = grid.faces();
  int level = grid.level();
  int index_lower[3] = { grid.index_lower(0), 
			 grid.index_lower(1), 
			 grid.index_lower(2) };
  int index_upper[3] = { grid.index_upper(0), 
			 grid.index_upper(1), 
			 grid.index_upper(2) };

  // Update matrix stencil values
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper,
					 0, 1, &entries[0], v0);
  if (ierr != 0) ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[1], v1[0][1]);
  if (ierr != 0) ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[2], v1[0][0]);
  if (ierr != 0) ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[3], v1[1][1]);
  if (ierr != 0) ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[4], v1[1][0]);
  if (ierr != 0) ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[5], v1[2][1]);
  if (ierr != 0) ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[6], v1[2][0]);
  if (ierr != 0) ERROR("could not SetBoxValues in A_\n");

  if (coarsegrid && use_prec) {
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 1, &entries[0], v0);
    if (ierr != 0) ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 1, &entries[1], v1[0][1]);
    if (ierr != 0) ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 1, &entries[2], v1[0][0]);
    if (ierr != 0) ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 1, &entries[3], v1[1][1]);
    if (ierr != 0) ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 1, &entries[4], v1[1][0]);
    if (ierr != 0) ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 1, &entries[5], v1[2][1]);
    if (ierr != 0) ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 1, &entries[6], v1[2][0]);
    if (ierr != 0) ERROR("could not SetBoxValues in Ac_\n");
  }

  // Deallocate arrays
  delete [] v0;
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      delete [] v1[axis][face];
    } // for face
  } // for axis

} // AMRsolve_FLD::init_matrix_stencil_()

//------------------------------------------------------------------------

/// Clean up stencil connections between parts for FAC solver
void AMRsolve_FLD::init_matrix_clear_(int part)
{
  if (part > 0) {
    int r_factors[3] = {r_factor_,r_factor_,r_factor_}; 
    int ierr = HYPRE_SStructFACZeroAMRMatrixData(A_, part-1, r_factors);
    if (ierr != 0) ERROR("could not ZeroAMRMatrixData in A_\n");
  }

} // AMRsolve_FLD::init_matrix_clear_()

//------------------------------------------------------------------------

/// Solve using BiCGStab
void AMRsolve_FLD::solve_bicgstab_(int itmax, double restol)
{
  _TRACE_;
  int ierr;

  // local variables
  scalar beta_n, beta_d, beta, omega_d, omega_n, omega, alpha_n, alpha_d, alpha;
  int flag = 1;
  iter_ = 0;
  resid_ = 0.0;

  // set up the preconditioner (if requested)
  AMRsolve_HG_prec *precond;
  if (use_prec) {

    // set default solver parameters for HG preconditioner (if unset)
    if (parameters_->value("prec_itmax") == "")
      parameters_->add_parameter("prec_itmax","1"); 
    if (parameters_->value("prec_restol") == "")  
      parameters_->add_parameter("prec_restol","0.0");
    if (parameters_->value("prec_rlxtype") == "")  
      parameters_->add_parameter("prec_rlxtype","2");
    if (parameters_->value("prec_npre") == "")  
      parameters_->add_parameter("prec_npre","3"); 
    if (parameters_->value("prec_npost") == "")  
      parameters_->add_parameter("prec_npost","3");
    if (parameters_->value("prec_printl") == "")  
      parameters_->add_parameter("prec_printl","1");
    if (parameters_->value("prec_log") == "")  
      parameters_->add_parameter("prec_log","1");
    if (parameters_->value("prec_Jaciters") == "")  
      parameters_->add_parameter("prec_Jaciters","5");

    precond = new AMRsolve_HG_prec(*hierarchy_, BdryType_);
    ierr = precond->Initialize_(parameters_, &Ac_, &Xc_, &Bc_, &Y_);
    if (ierr != 0) ERROR("could not initialize HG preconditioner\n");
    // precond->Setup_();
  }

  // Solve the linear system
  //    compute initial residual, and prepare to begin iterations
  matvec_(ix, iv);                                   // v = A*x
  linear_combination(1.0, ib, -1.0, iv, 0.0, irs);   // rs = b-v
  copy_(irs, ir);                                    // r = rs
  beta_d = dot_(irs, irs);                           // beta = <rs,rs>
  resid_ = sqrt(beta_d);                             // resid = ||rs||
  if (resid_ > restol) {
    copy_(irs, ip);                                  // p = rs
    for (iter_=0; iter<itmax; iter_++) {             // iteration loop

      if (use_prec) {                                // v = M\p
	//precond->Solve_();
	copy_(ip, iv);
      } else {
	copy_(ip, iv);
      }
      matvec_(iv, is);                               // s = A*v
      ddot_(alpha_n, ir, irs, alpha_d, is, irs);     // alpha_n=<r,rs>, alpha_d=<s,rs>
      alpha = alpha_n / alpha_s;
      linear_combination_(1.0, ir, -alpha, is, 0.0, iq);   // q = r - alpha*s
      if (use_prec) {                                // v = M\q
	//precond->Solve_();
	copy_(iq, iv);
      } else {
	copy_(iq, iv);
      }
      matvec_(iv, ir);                               // r = A*v
      ddot_(omega_n, ir, iq, omega_d, ir, ir);       // omega_n=<r,q>, omega_d=<r,r>
      omega = (omega_d == 0.0) ? omega_n : omega_n/omega_d;
      if (omega == 0.0) {                            // check for breakdown
	flag = -1;
	break;
      }
      linear_combination_(alpha, ip, omega, iq, 1.0, ix);    // x = x + alpha*p + omega*q
      linear_combination_(1.0, iq, 0.0, iNULL, -omega, ir);  // r = q - omega*r
      ddot_(resid_, ir, ir, beta_n, ir, irs);        // resid=<r,r>, beta_n=<r,rs>
      resid_ = sqrt(resid_);                         // resid = ||r||
      if (resid_ < restol) {
	flag = 0;
	break;
      }
      beta = (beta_n / beta_d) * (alpha / omega);
      beta_d = beta_n;
      if (beta == 0.0) {                            // check for breakdown
	flag = -2;
	break;
      }

    } // end for

    // if successful, compute final solution
    if (flag == 0) {
      if (use_prec) {                                // v = M\x
	//precond->Solve_();
	//copy_(iv, ix);
      }
    // if unsuccessful, report on why and fail
    } else {
      if (flag == 1)
	fprintf(stderr, "BiCGStab could not converge in allowed iterations:\n");
      if (flag == -1)
	fprintf(stderr, "BiCGStab breakdown (omega = 0)\n");
      if (flag == -2)
	fprintf(stderr, "BiCGStab breakdown (beta = 0)\n");
      fprintf(stderr, "  iterations = %d\n  final residual = %g\n",iter_,resid_);
      hierarchy_->print();
      ERROR("could not solve with BiCGStab\n");
    }
  }  // end check whether initial guess was good enough


  // Write out some diagnostic info about the solve
  if (debug && pmpi->is_root()) {
    printf("hypre BiCGSTAB num iterations: %d\n",iter_);
    printf("hypre BiCGSTAB final relative residual norm: %g\n",resid_);
  }

  // Delete the preconditioner
  if (use_prec)  delete precond;

  _TRACE_;
} // AMRsolve_FLD::solve_bicgstab_()

//------------------------------------------------------------------------

/// matvec_ operation:  b = A*x
void AMRsolve_FLD::matvec_(int x, int b)
{
  _TRACE_;

  _TRACE_;
} // AMRsolve_FLD::matvec_()

//------------------------------------------------------------------------

/// copy_ operation:  y=x
void AMRsolve_FLD::copy_(int x, int y)
{
  _TRACE_;

  _TRACE_;
} // AMRsolve_FLD::copy_()

//------------------------------------------------------------------------

/// dot_ operation:  <x,y>
double AMRsolve_FLD::dot_(int x, int y)
{
  _TRACE_;

  _TRACE_;
} // AMRsolve_FLD::dot_()

//------------------------------------------------------------------------

/// ddot_ operation:  dot1=<x,y>, dot2=<p,q>
void AMRsolve_FLD::ddot_(double& dot1, int x, int y,
			 double& dot2, int p, int q)
{
  _TRACE_;

  _TRACE_;
} // AMRsolve_FLD::ddot_()

//------------------------------------------------------------------------

/// linear_combination_ operation:  z = a*x + b*y * c*z
void AMRsolve_FLD::linear_combination_(scalar a, int x, scalar b, int y, scalar c, int z)
{
  _TRACE_;

  _TRACE_;
} // AMRsolve_FLD::linear_combination_()

//------------------------------------------------------------------------


/// Update the matrix at a coarse/fine interface; this routine changes the 
/// fine-grid matrix to account for the coarse grid neighbor.  Called in 
/// two phases, with phase == phase_graph (via init_graph_nonstencil_()) 
/// for the nonzero structure, and with phase == phase_matrix
/// (via init_matrix_nonstencil_()) for the matrix nonzeros.
void AMRsolve_FLD::update_fine_coarse_const_(int face, 
						   AMRsolve_Grid& grid_fine, 
						   int axis0, 
						   phase_enum phase,
						   int level_fine, 
						   int level_coarse,
						   int index_fine[3], 
						   int index_coarse[3])
{
  int ierr;
  int axis1 = (axis0+1)%3;
  int axis2 = (axis0+2)%3;

  // declare shortcut variables
  Scalar Ed, D;
  Scalar afac = adot_ / aval_;
  Scalar dtfac = dt_ * theta_;

  // get active enzo grid size
  int n3[3];
  grid_fine.get_size(n3);

  // get buffering information on relating amrsolve grid to Enzo data
  int ghosts[3][2]; 
  grid_fine.get_Ghosts(ghosts);
  int en0 = n3[0] + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
  int en1 = n3[1] + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
  int en2 = n3[2] + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions

  // indices for Enzo fine grid cells on both sides of interface
  int k_row, k_col, k0, k1, k2;
  int adj0=0, adj1=0, adj2=0;
  if (axis0 == 0)  adj0 = (face == 0) ? -1 : 1;  
  if (axis0 == 1)  adj1 = (face == 0) ? -1 : 1;  
  if (axis0 == 2)  adj2 = (face == 0) ? -1 : 1;  

  // global grid index limits for fine grid
  int index_global[3][2];
  grid_fine.indices(index_global);
  
  // set this grid's mesh spacing in this direction
  Scalar dxi = 1.0 / grid_fine.h(axis0) / lUn_;
  Scalar dxfac = dtfac * dxi * dxi;

  // get location of fine grid
  double xl0, xl1, xl2, xu0, xu1, xu2;
  grid_fine.x_lower(xl0, xl1, xl2);
  grid_fine.x_upper(xu0, xu1, xu2);
  
  //--------------------------------------------------
  // (*) CONSTANT
  //     Scale        = 2/3 (distance between coarse/fine cell 
  //                         centers, as a fraction of h_coarse; 
  //                         since we have 1/h^2, we use 4/9)
  //     Coefficients = determined on the fly
  //--------------------------------------------------

  int index_increment[][3] = {{face*(r_factor_-1),0,0},
			      {0,1,0},
			      {0,0,1},
			      {0,-1,0},
			      {-face*(r_factor_-1),0,-1}};

  if (grid_fine.is_local()) {

    if (phase == phase_graph) {

      int k=0;
      index_fine[axis0] += index_increment[k][0];
      index_fine[axis1] += index_increment[k][1];
      index_fine[axis2] += index_increment[k][2];

      for (k=1; k<5; k++) {
	ierr = HYPRE_SStructGraphAddEntries(graph_, level_fine, index_fine, 
					    0, level_coarse, index_coarse, 0);
	if (ierr != 0) ERROR("could not add graph entries\n");
	index_fine[axis0] += index_increment[k][0];
	index_fine[axis1] += index_increment[k][1];
	index_fine[axis2] += index_increment[k][2];
      } // for k = 1:4

    } else if (phase == phase_matrix) {

      // access relevant arrays from this grid to compute RHS
      Scalar* E     = grid_fine.get_E();
      Scalar* kappa = grid_fine.get_kap();

      // fine->coarse off-diagonal scaling
      //      double val_s = 2.0 / 3.0;
      double val_s = 4.0 / 9.0;
      int entry;
      double val, value;

      int k=0;
      index_fine[axis0] += index_increment[k][0];
      index_fine[axis1] += index_increment[k][1];
      index_fine[axis2] += index_increment[k][2];

      for (k=1; k<5; k++) {

	// // Query and zero out existing matrix entry
	// double val, val2 = 0.0;
	// int entry = 2*axis0 + 1 + face;
	// ierr = HYPRE_SStructMatrixGetValues(A_, level_fine, index_fine, 
	// 				    0, 1, &entry, &val);
	// if (ierr != 0) ERROR("could not GetValues from A_\n");
	// ierr = HYPRE_SStructMatrixSetValues(A_, level_fine, index_fine, 
	// 				    0, 1, &entry, &val2);
	// if (ierr != 0) ERROR("could not SetValues in A_\n");
    
	// // Set new values across interface, scaling to depend equally on each fine neighbor
	// val2 = val;
	// entry = grid_fine.counter(index_fine)++;
	// ierr = HYPRE_SStructMatrixSetValues(A_, level_fine, index_fine,
	// 				    0, 1, &entry, &val2);
	// if (ierr != 0) ERROR("could not SetValues in A_\n");
	
	// set indices for Enzo fine cells on both sides of interface
	k0 = index_fine[0] - index_global[0][0] + ghosts[0][0];
	k1 = index_fine[1] - index_global[1][0] + ghosts[1][0];
	k2 = index_fine[2] - index_global[2][0] + ghosts[2][0];
	k_row = k0 + en0*(k1 + en1*k2);
	k_col = k0 + adj0 + en0*(k1 + adj1 + en1*(k2 + adj2));

	// Compute limiter at this face
	Ed = E[k_row] - E[k_col];
	D = limiter_(E[k_row], E[k_col], kappa[k_row],
		     kappa[k_col], nUn_, lUn_, dxi);

	// Set matrix values across coarse/fine face
	val = -val_s * dxfac * D;
	if (isinf(val) || isnan(val)) {
	  fprintf(stderr,"update_fine_coarse_const ERROR: encountered illegal value (%g)\n   val_s = %g, dtfac = %g, dxi = %g, D = %g, k_row = %i, k_col = %i\n\n", val, val_s, dtfac, dxi, D, k_row, k_col);
	  ERROR("NaN or Inf value in setting matrix entries\n");
	}
	
	//   Update off-diagonal
	entry = grid_fine.counter(index_fine)++;
	value = val;
	ierr = HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
					      0, 1, &entry, &value);
	if (ierr != 0) ERROR("could not AddToValues in A_\n");

	//   Update diagonal
	entry = 0;
	value = -val;
	ierr = HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
					      0, 1, &entry, &value);
	if (ierr != 0) ERROR("could not AddToValues in A_\n");


	// Clear original matrix values from stencil
	val = -dxfac * D;
	if (isinf(val) || isnan(val)) {
	  fprintf(stderr,"update_fine_coarse_const ERROR: encountered illegal value (%g)\n   dtfac = %g, dxi = %g, D = %g, k_row = %i, k_col = %i\n\n", val, dtfac, dxi, D, k_row, k_col);
	  ERROR("NaN or Inf value in setting matrix entries\n");
	}
	
	//   Update off-diagonal, stencil xp=1,xm,yp,ym,zp,zm=6
	entry = 2*axis0 + 1 + (1-face);
	value = -val;
	ierr = HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
					      0, 1, &entry, &value);
	if (ierr != 0) ERROR("could not AddToValues in A_\n");

	//   Update diagonal
	entry = 0;
	value = val;
	ierr = HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
					      0, 1, &entry, &value);
	if (ierr != 0) ERROR("could not AddToValues in A_\n");

	// Update indices
	index_fine[axis0] += index_increment[k][0];
	index_fine[axis1] += index_increment[k][1];
	index_fine[axis2] += index_increment[k][2];

      } // for k = 1,4
     } // if phase == phase_matrix
  } // if grid_fine.is_local()
} // AMRsolve_FLD::update_fine_coarse_const_

//------------------------------------------------------------------------

/// Update the matrix at a coarse/fine interface; this routine changes the 
/// coarse-grid matrix to account for the fine grid neighbors.  Called in 
/// two phases, with phase == phase_graph (via init_graph_nonstencil_()) 
/// for the nonzero structure, and with phase == phase_matrix
/// (via init_matrix_nonstencil_()) for the matrix nonzeros.
void AMRsolve_FLD::update_coarse_fine_const_(int face, 
						   AMRsolve_Grid& grid_coarse, 
						   int axis0, 
						   phase_enum phase,
						   int level_fine, 
						   int level_coarse,
						   int index_fine[3], 
						   int index_coarse[3])
{
  int ierr;
  // set graph entry
  if (phase == phase_graph) {
    int index_increment[][3] = {{1,0,0},
				{0,1,0},
				{-1,0,0},
				{0,0,1},
				{1,0,0},
				{0,-1,0},
				{-1,0,0},
				{0,0,-1}};
    for (int k=0; k<8; k++) {
      ierr = HYPRE_SStructGraphAddEntries(graph_, level_coarse, index_coarse, 
					  0, level_fine, index_fine, 0);
      if (ierr != 0) ERROR("could not add entries to graph_\n");
      index_fine[0] += index_increment[k][0];
      index_fine[1] += index_increment[k][1];
      index_fine[2] += index_increment[k][2];
    } // for k=0,7

  // set matrix entry
  } else if (phase == phase_matrix) {

    // Get existing matrix entry to neighboring cell (and remove from matrix)
    double val, val2 = 0.0;
    int entry = 2*axis0 + 1 + face;
    ierr = HYPRE_SStructMatrixGetValues(A_, level_coarse, index_coarse, 
					0, 1, &entry, &val);
    if (ierr != 0) ERROR("could not GetValues from A_\n");
    ierr = HYPRE_SStructMatrixSetValues(A_, level_coarse, index_coarse, 
					0, 1, &entry, &val2);
    if (ierr != 0) ERROR("could not SetValues in A_\n");
    
    // Set new values across interface, scaling to depend equally on each fine neighbor
    val2 = val / 8.0;
    for (int i=0; i<8; i++) {
      entry = grid_coarse.counter(index_coarse)++;
      ierr = HYPRE_SStructMatrixSetValues(A_, level_coarse, index_coarse, 
					  0, 1, &entry, &val2);
      if (ierr != 0) ERROR("could not SetValues in A_\n");
    }  // for i=0,7

    // // declare shortcut variables
    // double Eavg, Ed, R, kap, D;
    // double afac = adot_ / aval_;
    // double dtfac = dt_ * theta_;
    // // double Rmin = 1.0e-20 / lUn_;
    // double Rmin = 1.0e-20;
    // double c = 2.99792458e10;
    
    // // access relevant arrays from this grid to compute RHS
    // Scalar* E     = grid_coarse.get_E();
    // Scalar* kappa = grid_coarse.get_kap();
    
    // // get active enzo grid size
    // int n3[3];
    // grid_coarse.get_size(n3);

    // // get buffering information on relating amrsolve grid to Enzo data
    // int ghosts[3][2]; 
    // grid_coarse.get_Ghosts(ghosts);
    // int en0 = n3[0] + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
    // int en1 = n3[1] + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
    // int en2 = n3[2] + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions
    
    // // global grid index limits for coarse grid
    // int index_global[3][2];
    // grid_coarse.indices(index_global);
  
    // // set this grid's mesh spacing in this direction
    // Scalar dxi = 1.0 / grid_coarse.h(axis0) / lUn_;
    // Scalar dxfac = dtfac * dxi * dxi;
    
    // // set indices for Enzo coarse grid cells on both sides of interface
    // int adj0=0, adj1=0, adj2=0;
    // if (axis0 == 0)  adj0 = (face == 0) ? -1 : 1;  
    // if (axis0 == 1)  adj1 = (face == 0) ? -1 : 1;  
    // if (axis0 == 2)  adj2 = (face == 0) ? -1 : 1;  
    
    // int k0 = index_coarse[0] - index_global[0][0] + ghosts[0][0];
    // int k1 = index_coarse[1] - index_global[1][0] + ghosts[1][0];
    // int k2 = index_coarse[2] - index_global[2][0] + ghosts[2][0];
    // int k_row = k0 + en0*(k1 + en1*k2);
    // int k_col = k0 + adj0 + en0*(k1 + adj1 + en1*(k2 + adj2));
    
    // // Compute limiter at this face
    // Ed = E[k_row] - E[k_col];
    // Eavg = (E[k_row] + E[k_col])*0.5;
    // if (Eavg != 0.0)
    //   R = MAX(dxi*fabs(Ed)/Eavg, Rmin);
    // else
    //   R = Rmin;
    // kap = 2.0 * kappa[k_row] * kappa[k_col] / (kappa[k_row] + kappa[k_col]);
    // // kap = 0.5 * (kappa[k_row] + kappa[k_col]);
    // D = c/sqrt(9.0*kap*kap*nUn_*nUn_ + R*R);

    // // set matrix entry across coarse/fine face
    // double val_s = 1.0 / 8.0;
    // double val = -val_s * dxfac * D;
    // if (isnan(val) || isinf(val)) {
    //   fprintf(stderr,"update_coarse_fine_const ERROR: encountered NaN/Inf value (%g)\n   val_s = %g, dtfac = %g, dxi = %g, D = %g, k_row = %i, k_col = %i\n\n", val, val_s, dtfac, dxi, D, k_row, k_col);
    //   ERROR("NaN or Inf value in setting matrix entries\n");
    // }
    // int    entry;
    // double value;

    // // Adjust coarse-fine nonstencil values
    // for (int i=0; i<8; i++) {
    //   // Set new nonstencil coarse-fine entry
    //   entry = grid_coarse.counter(index_coarse)++;
    //   value = val;
    //   ierr = HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
    // 					    0, 1, &entry, &value);
    //   if (ierr != 0) ERROR("could not AddToValues in A_\n");
    //   // Adjust stencil diagonal
    //   entry = 0;
    //   value = -val;
    //   ierr = HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
    // 					    0, 1, &entry, &value);
    //   if (ierr != 0) ERROR("could not AddToValues in A_\n");
    // } // for i=0,7

    // // Clear original matrix values from stencil
    // val = -dxfac * D;

    // //   Update off-diagonal, stencil xp=1,xm,yp,ym,zp,zm=6
    // //   (note: "face" is for fine grid, but we want coarse)
    // entry = 2*axis0 + 1 + face;
    // value = -val;
    // ierr = HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
    // 					  0, 1, &entry, &value);
    // if (ierr != 0) ERROR("could not AddToValues in A_\n");

    // //   Update diagonal
    // entry = 0;
    // value = val;
    // ierr = HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
    // 					  0, 1, &entry, &value);
    // if (ierr != 0) ERROR("could not AddToValues in A_\n");

  } // if phase == phase_matrix
} // AMRsolve_FLD::update_coarse_fine_const_

//------------------------------------------------------------------------

/// Stub routine for inserting unit tests; should be called at same point 
/// in code as solve() routine. 
void AMRsolve_FLD::tester()
{
#define TEST7

#ifdef TEST1
  ////////  test 1 -- see if AMRsolve_Grid::overlap_indices works ////////
  // iterate over all grids on this processor
  ItHierarchyGridsLocal itgl(*hierarchy_);
  while (AMRsolve_Grid* grid1 = itgl++) {

    // output this grid's dimensions
    Scalar xl0, xu0, xl1, xu1, xl2, xu2;
    grid1->x_lower(xl0, xl1, xl2);
    grid1->x_upper(xu0, xu1, xu2);
    printf("Grid %i dims = %i %i %i, bounds = (%g:%g, %g:%g, %g:%g)\n",
  	   grid1->id(),grid1->n(0),grid1->n(1),grid1->n(2),xl0,xu0,xl1,xu1,xl2,xu2);

    // iterate over all other grids
    ItHierarchyGridsAll itga(*hierarchy_);
    while (AMRsolve_Grid* grid2 = itga++) {

      // get grid IDs
      int g1_ID = grid1->id();
      int g2_ID = grid2->id();

      // see if grids overlap, if so output the overlap region
      int g1_ilo[3], g1_ihi[3], g2_ilo[3], g2_ihi[3];
      if (grid1->overlap_indices(grid2, g1_ilo, g1_ihi, g2_ilo, g2_ihi)) {
  	printf("Grid overlap %i with %i: regions (%i:%i,%i:%i,%i:%i) and (%i:%i,%i:%i,%i:%i)\n\n",
  	       g1_ID, g2_ID, g1_ilo[0], g1_ihi[0], g1_ilo[1], g1_ihi[1], g1_ilo[2], g1_ihi[2],
  	       g2_ilo[0], g2_ihi[0], g2_ilo[1], g2_ihi[1], g2_ilo[2], g2_ihi[2]);
      } else {
  	printf("Grids %i and %i do not overlap\n\n", g1_ID, g2_ID);
      }

    } // while itga
  } // while itgl
#endif


#ifdef TEST2
  ////////  test 2 -- see if AMRsolve_Hierarchy::restrict works ////////
  // set a value of "1" into a single cell somewhere in finest grid, zeroing out all other cells
  int numlevels = hierarchy_->num_levels();
  ItHierarchyGridsLocal itgl(*hierarchy_);
  while (AMRsolve_Grid* grid = itgl++) {
    int n3[3];
    Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
    for (int id=0; id<n3[0]*n3[1]*n3[2]; id++)  u[id] = 0.0;

    if (grid->level() == numlevels-1) {
      int i = n3[0]/2;
      int j = n3[1]/2;
      int k = n3[2]/2;
      printf("grid %i on finest level, setting center zone (%i,%i,%i) to 1\n",
  	     grid->id(),i,j,k);
      u[(k*n3[1] + j)*n3[0] + i] = 1.0;
      break;
    }
  }

  // create AMRsolve_HG_prec object
  AMRsolve_HG_prec *precond = new AMRsolve_HG_prec(*hierarchy_, BdryType_);

  // perform restriction
  if (precond->restrict(numlevels-1, 0) != 0)
    printf("ERROR: restriction failed\n");

  // delete HG_prec object
  delete precond;

  // output nonzero values/locations from all grids
  // iterate over all grids on this processor
  ItHierarchyGridsLocal itgl2(*hierarchy_);
  while (AMRsolve_Grid* grid = itgl2++) {
    int n3[3];
    Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
    Scalar usum=0.0;
    for (int id=0; id<n3[0]*n3[1]*n3[2]; id++)  usum += fabs(u[id]);
    if (usum > 0.0) {
      Scalar xl0, xu0, xl1, xu1, xl2, xu2;
      grid->x_lower(xl0, xl1, xl2);
      grid->x_upper(xu0, xu1, xu2);
      printf("Grid %i has nonzeros. level %i, dims %i %i %i, extents (%g:%g,%g:%g,%g:%g):\n",
  	     grid->id(),grid->level(),grid->n(0),grid->n(1),grid->n(2),xl0,xu0,xl1,xu1,xl2,xu2);
      for (int k=0; k<n3[2]; k++)
  	for (int j=0; j<n3[1]; j++)
  	  for (int i=0; i<n3[0]; i++) 
  	    if (fabs(u[(k*n3[1] + j)*n3[0] + i]) > 0.0)
  	      printf("  u(%i,%i,%i) = %g\n",i,j,k,u[(k*n3[1] + j)*n3[0] + i]);
    }
  }  
#endif


#ifdef TEST3
  ////////  test 3 -- see if AMRsolve_Hierarchy::prolong works ////////
  // set a value of "1" into a single cell in center of coarsest grid, zeroing out all other cells
  int numlevels = hierarchy_->num_levels();
  ItHierarchyGridsLocal itgl(*hierarchy_);
  while (AMRsolve_Grid* grid = itgl++) {
    int n3[3];
    Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
    for (int id=0; id<n3[0]*n3[1]*n3[2]; id++)  u[id] = 0.0;

    if (grid->level() == 0) {
      int i = n3[0]/2;
      int j = n3[1]/2;
      int k = n3[2]/2;
      printf("grid %i on coarsest level, setting center zone (%i,%i,%i) to 1\n",
  	     grid->id(),i,j,k);
      u[(k*n3[1] + j)*n3[0] + i] = 1.0;
      break;
    }
  }

  // create AMRsolve_HG_prec object
  AMRsolve_HG_prec *precond = new AMRsolve_HG_prec(*hierarchy_, BdryType_);

  // perform prolongation, method 0 (piecewise constant)
  if (precond->prolong(0, numlevels-1, 0) != 0)
    printf("ERROR: prolongation failed\n");

  // delete HG_prec object
  delete precond;

  // output nonzero values/locations from all grids
  // iterate over all grids on this processor
  ItHierarchyGridsLocal itgl2(*hierarchy_);
  while (AMRsolve_Grid* grid = itgl2++) {
    int n3[3];
    Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
    Scalar usum=0.0;
    for (int id=0; id<n3[0]*n3[1]*n3[2]; id++)  usum += fabs(u[id]);
    if (usum > 0.0) {
      Scalar xl0, xu0, xl1, xu1, xl2, xu2;
      grid->x_lower(xl0, xl1, xl2);
      grid->x_upper(xu0, xu1, xu2);
      printf("Grid %i has nonzeros. level %i, dims %i %i %i, extents (%g:%g,%g:%g,%g:%g):\n",
  	     grid->id(),grid->level(),grid->n(0),grid->n(1),grid->n(2),xl0,xu0,xl1,xu1,xl2,xu2);
      for (int k=0; k<n3[2]; k++)
  	for (int j=0; j<n3[1]; j++)
  	  for (int i=0; i<n3[0]; i++) 
  	    if (fabs(u[(k*n3[1] + j)*n3[0] + i]) > 0.0)
  	      printf("  u(%i,%i,%i) = %g\n",i,j,k,u[(k*n3[1] + j)*n3[0] + i]);
    }
  }  
#endif


#ifdef TEST4
  ////////  test 4 -- check spread of restrict->prolong ////////
  // set a value of "1" into a single cell in center of finestgrid, zeroing out all other cells
  int numlevels = hierarchy_->num_levels();
  ItHierarchyGridsLocal itgl(*hierarchy_);
  while (AMRsolve_Grid* grid = itgl++) {
    int n3[3];
    Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
    for (int id=0; id<n3[0]*n3[1]*n3[2]; id++)  u[id] = 0.0;
    if (grid->level() == numlevels-1) {
      int i = n3[0]/2;
      int j = n3[1]/2;
      int k = n3[2]/2;
      printf("grid %i on finest level, setting center zone (%i,%i,%i) to 1\n",
  	     grid->id(),i,j,k);
      u[(k*n3[1] + j)*n3[0] + i] = 1.0;
      break;
    }
  }

  // create AMRsolve_HG_prec object
  AMRsolve_HG_prec *precond = new AMRsolve_HG_prec(*hierarchy_, BdryType_);

  // perform restriction (piecewise constant)
  if (precond->restrict(numlevels-1, 0) != 0)
    printf("ERROR: restriction failed\n");

  // perform prolongation, method 0 (piecewise constant)
  if (precond->prolong(0, numlevels-1, 0) != 0)
    printf("ERROR: prolongation failed\n");

  // delete HG_prec object
  delete precond;

  // output nonzero values/locations from finest grid
  // iterate over all grids on this processor
  ItHierarchyGridsLocal itgl2(*hierarchy_);
  while (AMRsolve_Grid* grid = itgl2++) {
    if (grid->level() == numlevels-1) {
      int n3[3];
      Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
      printf("Grid %i. level %i, dims %i %i %i:\n",
	     grid->id(), grid->level(), grid->n(0), grid->n(1), grid->n(2));
      for (int k=0; k<n3[2]; k++)
	for (int j=0; j<n3[1]; j++)
	  for (int i=0; i<n3[0]; i++) 
	    if (fabs(u[(k*n3[1] + j)*n3[0] + i]) > 0.0)
	      printf("  u(%i,%i,%i) = %g\n",i,j,k,u[(k*n3[1] + j)*n3[0] + i]);
    }
  }  
#endif


#ifdef TEST5
  ////////  test 5 -- see whether prolong->restrict = identity ////////
  int icell, jcell, kcell, cgrid, myid;

  // set a value of "1" into a single cell in center of coarsest grid, zeroing out all other cells
  int numlevels = hierarchy_->num_levels();
  ItHierarchyGridsLocal itgl(*hierarchy_);
  while (AMRsolve_Grid* grid = itgl++) {
    int n3[3];
    Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
    for (int id=0; id<n3[0]*n3[1]*n3[2]; id++)  u[id] = 0.0;

    if (grid->level() == 0) {
      icell = n3[0]/2;
      jcell = n3[1]/2;
      kcell = n3[2]/2;
      cgrid = grid->id();
      myid = grid->ip();
      u[(kcell*n3[1] + jcell)*n3[0] + icell] = 1.0;
      break;
    }
  }

  // create AMRsolve_HG_prec object
  AMRsolve_HG_prec *precond = new AMRsolve_HG_prec(*hierarchy_, BdryType_);
  
  // perform prolongation, method 0 (piecewise constant)
  if (precond->prolong(0, numlevels-1, 0) != 0)
    printf("ERROR: prolongation failed\n");
  
  // perform restriction (piecewise constant)
  if (precond->restrict(numlevels-1, 0) != 0)
    printf("ERROR: restriction failed\n");

  // delete HG_prec object
  delete precond;


  // compute error: (coarsest grid) - (perturbation)
  Scalar err=0.0;
  // iterate over all grids on this processor
  ItHierarchyGridsLocal itgl2(*hierarchy_);
  while (AMRsolve_Grid* grid = itgl2++) {
    if (grid->level() == 0) {
      int n3[3];
      Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
      for (int k=0; k<n3[2]; k++)
  	for (int j=0; j<n3[1]; j++)
  	  for (int i=0; i<n3[0]; i++) 
	    if (i==icell && j==jcell && k==kcell && grid->id()==cgrid) {
	      err += fabs(u[(k*n3[1] + j)*n3[0] + i] - 1.0);
	    } else {
	      err += fabs(u[(k*n3[1] + j)*n3[0] + i]);
	    }	      
    }
  }
  printf("proc %i: prolong->interp error = %g\n",myid,err);
  
#endif

#ifdef TEST6
  ////////  test 6 -- AMRsolve_HG_prec initialization and smoother ////////

  // set up the preconditioner
  AMRsolve_HG_prec *precond;
  precond = new AMRsolve_HG_prec(*hierarchy_, BdryType_);
  int ierr = precond->Initialize_(parameters_, &Ac_, &Xc_, &Bc_, &Y_);
  if (ierr)
    fprintf(stderr,"Error in AMRsolve_HG_prec::Initialize = %i\n",ierr);

  // check initial linear residual
  double resid2;
  ierr = HYPRE_SStructVectorCopy(B_, Y_);
  ierr = HYPRE_SStructMatrixMatvec(-1.0, A_, X_, 1.0, Y_);
  ierr = HYPRE_SStructInnerProd(Y_, Y_, &resid2);
  printf("tester: initial linear residual norm = %g\n",sqrt(resid2));

  // test Jacobi smoother
  int nsweeps = 8;
  for (int i=0; i<nsweeps; i++) {
    ierr = precond->Jacobi_smooth_(A_, X_, B_, Y_);
    if (ierr)
      fprintf(stderr,"Error in AMRsolve_HG_prec::Jacobi_smooth = %i\n",ierr);

    // check resulting linear residual
    ierr = HYPRE_SStructVectorCopy(B_, Y_);
    ierr = HYPRE_SStructMatrixMatvec(-1.0, A_, X_, 1.0, Y_);
    ierr = HYPRE_SStructInnerProd(Y_, Y_, &resid2);
    printf("tester: sweep %i, linear residual norm = %g\n",i,sqrt(resid2));
  }

  // Delete the preconditioner
  delete precond;

#endif

#ifdef TEST7
  ////////  test 7 -- AMRsolve_HG_prec init, setup, solve, stats ////////

  // set up the preconditioner
  AMRsolve_HG_prec *precond;
  precond = new AMRsolve_HG_prec(*hierarchy_, BdryType_);
  int ierr = precond->Initialize_(parameters_, &Ac_, &Xc_, &Bc_, &Y_);
  if (ierr)
    fprintf(stderr,"Error in AMRsolve_HG_prec::Initialize = %i\n",ierr);

  // set up the preconditioner
  ierr = precond->Setup_(A_, B_, X_);
  if (ierr)
    fprintf(stderr,"Error in AMRsolve_HG_prec::Setup = %i\n",ierr);

  // check initial linear residual
  double resid2;
  ierr = HYPRE_SStructVectorCopy(B_, Y_);
  ierr = HYPRE_SStructMatrixMatvec(-1.0, A_, X_, 1.0, Y_);
  ierr = HYPRE_SStructInnerProd(Y_, Y_, &resid2);
  printf("tester: initial linear residual norm = %g\n",sqrt(resid2));

  // solve with the preconditioner
  ierr = precond->Solve_(A_, B_, X_);
  if (ierr)
    fprintf(stderr,"Error in AMRsolve_HG_prec::Solve = %i\n",ierr);

  // extract preconditioner iterations and residual
  Scalar rnorm = precond->GetResid_();
  int    coarse_iters = precond->GetIters_();
  printf("AMRsolve_HG_Prec achieved coarse residual %g in %i iters\n",coarse_iters,rnorm);

  // check final linear residual
  ierr = HYPRE_SStructVectorCopy(B_, Y_);
  ierr = HYPRE_SStructMatrixMatvec(-1.0, A_, X_, 1.0, Y_);
  ierr = HYPRE_SStructInnerProd(Y_, Y_, &resid2);
  printf("tester: final linear residual norm = %g\n",sqrt(resid2));

  // Delete the preconditioner
  delete precond;

#endif

}  // AMRsolve_FLD::tester

//------------------------------------------------------------------------
