/// @file      AMRsolve_hypre_grav.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Implementation of the AMRsolve_Hypre_Grav class

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
//#include "AMRsolve_point.h"
#include "AMRsolve_faces.h"
#include "AMRsolve_domain.h"
#include "AMRsolve_grid.h"
#include "AMRsolve_level.h"
#include "AMRsolve_hierarchy.h"
#include "AMRsolve_parameters.h"
//#include "AMRsolve_problem.h"
#include "AMRsolve_hypre_grav.h"
#include "AMRsolve_error.h"

//======================================================================

// Coefficient for Poisson problem

inline Scalar acoef(Scalar x, Scalar y, Scalar z)
{
  return 1.0;
}

//======================================================================
// PUBLIC MEMBER FUNCTIONS
//======================================================================

/// AMRsolve_Hypre_Grav constructor
AMRsolve_Hypre_Grav::AMRsolve_Hypre_Grav(AMRsolve_Hierarchy& hierarchy, 
					 AMRsolve_Parameters& parameters)
  : grid_(0), graph_(0), stencil_(0), A_(0), B_(0), X_(0),
    solver_(0), parameters_(&parameters), hierarchy_(&hierarchy),
    resid_(-1.0), iter_(-1), r_factor_(const_r_factor), matrix_scale_(1.0)
{
  //
}

//----------------------------------------------------------------------

/// AMRsolve_Hypre_Grav destructor
AMRsolve_Hypre_Grav::~AMRsolve_Hypre_Grav()
{
  // destroy HYPRE objects that we created along the way
  char error_message[100];
  if (HYPRE_SStructVectorDestroy(B_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy B_\n");
    ERROR(error_message);
  }
  if (HYPRE_SStructVectorDestroy(X_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy X_\n");
    ERROR(error_message);
  }
  if (HYPRE_SStructMatrixDestroy(A_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy A_\n");
    ERROR(error_message);
  }
  if (HYPRE_SStructGraphDestroy(graph_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy graph_\n");
    ERROR(error_message);
  }
  if (HYPRE_SStructStencilDestroy(stencil_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy stencil_\n");
    ERROR(error_message);
  }
  if (HYPRE_SStructGridDestroy(grid_) != 0) {
    sprintf(error_message, "AMRsolve_Hypre_FLD::could not destroy grid_\n");
    ERROR(error_message);
  }

}

//----------------------------------------------------------------------

/// Initialize the Grid Hierarchy
/** Creates a hypre grid, with one part per level and one box per Grid
    patch object, for an AMR problem.  Sets grid box extents, grid
    part variables, and periodicity. */
void AMRsolve_Hypre_Grav::init_hierarchy(AMRsolve_Mpi& mpi)
{

  int dim       = hierarchy_->dimension();
  int num_parts = hierarchy_->num_levels();

  // Create the hypre grid
  _TRACE_;
  HYPRE_SStructGridCreate(MPI_COMM_WORLD, dim, num_parts, &grid_);

  ItHierarchyLevels itl (*hierarchy_);

  _TRACE_;
  while (AMRsolve_Level * level = itl++) {

    _TRACE_;
    int part = level->index();

    ItLevelGridsLocal itgl (*level);

    // Set extents for boxes that comprise the hypre grid
    while (AMRsolve_Grid * grid = itgl++) {

      int lower[3] = {grid->index_lower(0),
		      grid->index_lower(1),
		      grid->index_lower(2)};
      int upper[3] = {grid->index_upper(0),
		      grid->index_upper(1),
		      grid->index_upper(2)};
      HYPRE_SStructGridSetExtents(grid_, part, lower, upper);
      
    } // while grid = itgl++

    _TRACE_;
    // Create a single cell-centered variable for each grid part (level)
    HYPRE_SStructVariable variable_types[] = { HYPRE_SSTRUCT_VARIABLE_CELL };
    const int numvars = 1;
    HYPRE_SStructGridSetVariables(grid_, part, numvars, variable_types);

    // Set periodicity of the grid part
    int period[3] = { hierarchy_->period_index(0,part),
		      hierarchy_->period_index(1,part),
		      hierarchy_->period_index(2,part) };

    _TRACE_;
    HYPRE_SStructGridSetPeriodic(grid_, part, period);

  } // while level = itl++

  // When finished, assemble the hypre grid
  _TRACE_;

  HYPRE_SStructGridAssemble(grid_);
  _TRACE_;
  
} // AMRsolve_Hypre_Grav::init_hierarchy()

//----------------------------------------------------------------------

/// Initialize the discretization stencils.  
/** Creates and initializes a hypre stencil object. */
void AMRsolve_Hypre_Grav::init_stencil()
{

  _TRACE_;
  int dim = hierarchy_->dimension();

  _TRACE_;
  HYPRE_SStructStencilCreate(dim,dim*2+1,&stencil_);

  int entries[][3] = { {  0, 0, 0 },     // center
		       {  1, 0, 0 },     // X+
		       { -1, 0, 0 },     // X-
		       {  0, 1, 0 },     // Y+
		       {  0,-1, 0 },     // Y-
		       {  0, 0, 1 },     // Z+
		       {  0, 0,-1 } };   // Z-

  _TRACE_;
  if (dim >= 1) HYPRE_SStructStencilSetEntry(stencil_, 0, entries[0], 0);
  if (dim >= 1) HYPRE_SStructStencilSetEntry(stencil_, 1, entries[1], 0);
  if (dim >= 1) HYPRE_SStructStencilSetEntry(stencil_, 2, entries[2], 0);
  if (dim >= 2) HYPRE_SStructStencilSetEntry(stencil_, 3, entries[3], 0);
  if (dim >= 2) HYPRE_SStructStencilSetEntry(stencil_, 4, entries[4], 0);
  if (dim >= 3) HYPRE_SStructStencilSetEntry(stencil_, 5, entries[5], 0);
  if (dim >= 3) HYPRE_SStructStencilSetEntry(stencil_, 6, entries[6], 0);
  _TRACE_;

} // AMRsolve_Hypre_Grav::init_stencil()

//----------------------------------------------------------------------

/// Initialize the graph.
/** Creates a graph containing the matrix nonzero structure.  Graph
    edges include both those for nonzeros from the stencil within each
    part (level), and nonzeros for graph entries connecting linked
    parts.  The matrix nonzero structure is generally nonsymmetric.
    Only the stencil step is required for unigrid problems. */
void AMRsolve_Hypre_Grav::init_graph()
{
  // Create the hypre graph object
  _TRACE_;
  HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid_, &graph_);
  
  _TRACE_;
  HYPRE_SStructGraphSetObjectType(graph_, HYPRE_SSTRUCT);

  _TRACE_;
  ItHierarchyLevels itl (*hierarchy_);

  _TRACE_;
  while (AMRsolve_Level * level = itl++) {

    int part = level->index();

    // Define stencil connections within each level
    _TRACE_;
    HYPRE_SStructGraphSetStencil(graph_, part, 0, stencil_);

    // Define graph connections between levels
    _TRACE_;
    if (part > 0) {
      ItLevelGridsAll itag (*level);
      while (AMRsolve_Grid * grid = itag++) {
	init_graph_nonstencil_(*grid);
      } // while grid = itag++
    } // if part > 0

    ItLevelGridsAll itag (*level);

    // Initialize face counters for subsequent matrix inter-level entries
    _TRACE_;
    while (AMRsolve_Grid * grid = itag++) {
      int dim = hierarchy_->dimension();
      grid->init_counter(dim*2+1);
    } // while grid = itag++
  } // while level = itl++

  // Assemble the hypre graph
  _TRACE_;
  HYPRE_SStructGraphAssemble(graph_);
  _TRACE_;

} // AMRsolve_Hypre_Grav::init_graph()

//----------------------------------------------------------------------

/// Initialize the matrix A and right-hand-side vector b
/* Creates a matrix with a given nonzero structure, and sets nonzero
   values. */
/* void AMRsolve_Hypre_Grav::init_elements(std::vector<AMRsolve_Point *> points, 
                                           Scalar f_scale) */
void AMRsolve_Hypre_Grav::init_elements(Scalar f_scale)
{
  // Create the hypre matrix A_, solution X_, and right-hand side B_ objects
  HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph_, &A_);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid_,  &X_);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid_,  &B_);

  // Set the object types
  if (parameters_->value("solver") == "bicgstab-boomer") {
    HYPRE_SStructMatrixSetObjectType(A_,HYPRE_PARCSR);
    HYPRE_SStructVectorSetObjectType(X_,HYPRE_PARCSR);
    HYPRE_SStructVectorSetObjectType(B_,HYPRE_PARCSR);
  } else {
    HYPRE_SStructMatrixSetObjectType(A_,HYPRE_SSTRUCT);
    HYPRE_SStructVectorSetObjectType(X_,HYPRE_SSTRUCT);
    HYPRE_SStructVectorSetObjectType(B_,HYPRE_SSTRUCT);
  }

  // Initialize the hypre matrix and vector objects
  HYPRE_SStructMatrixInitialize(A_);
  HYPRE_SStructVectorInitialize(X_);
  HYPRE_SStructVectorInitialize(B_);

  //--------------------------------------------------
  // Initialize the matrix A_ elements
  //--------------------------------------------------

  init_elements_matrix_();

  //--------------------------------------------------
  // Initialize B_ elements 
  //--------------------------------------------------

  //  init_elements_rhs_(points,f_scale);
  init_elements_rhs_(f_scale);

  // Assemble the matrix and vectors
  HYPRE_SStructMatrixAssemble(A_);
  HYPRE_SStructVectorAssemble(B_);
  HYPRE_SStructVectorAssemble(X_);

  // Optionally write the matrix and vectors to a file for debugging
  if (parameters_->value("dump_a") == "true")  HYPRE_SStructMatrixPrint("A-hypre",A_,0);
  if (parameters_->value("dump_x") == "true")  HYPRE_SStructVectorPrint("X0-hypre",X_,0);
  if (parameters_->value("dump_b") == "true")  HYPRE_SStructVectorPrint("B-hypre",B_,0);

} // AMRsolve_Hypre_Grav::init_elements()

//----------------------------------------------------------------------

/// Initialize and solve the linear solver
void AMRsolve_Hypre_Grav::solve()
{
  std::string solver = parameters_->value("solver");
  int         levels = hierarchy_->num_levels();

  int    itmax  = 0;
  double restol = 0.0;

  // Check solver parameters
  std::string sitmax  = parameters_->value("solver_itmax");
  std::string srestol = parameters_->value("solver_restol");

  // If not defined, then define them
  if (sitmax == "")  parameters_->add_parameter("solver_itmax","200");
  if (srestol == "") parameters_->add_parameter("solver_restol","1e-6");

  // Set local variables
  itmax  = atoi(sitmax.c_str());
  restol = atof(srestol.c_str());

  if (solver == "pfmg" && levels == 1) {
    solve_pfmg_(itmax,restol);

  } else if (solver == "fac"  && levels > 1) {
    solve_fac_(itmax,restol);

  } else if (solver == "bicgstab") {
    solve_bicgstab_(itmax,restol);

  } else if (solver == "bicgstab-boomer") {
    solve_bicgstab_boomer_(itmax,restol);

  } else if (solver == "gmres") {
    solve_gmres_(itmax,restol);

  } else {
    char error_message[100];
    sprintf(error_message, "AMRsolve_Hypre_Grav::solve called with illegal "
	    "combination of solver %s on %d levels", solver.c_str(),levels);
    ERROR(error_message);
  }
  
  if (parameters_->value("dump_x") == "true")  
    HYPRE_SStructVectorPrint("X-hypre",X_,1);

} // AMRsolve_Hypre_Grav::solve()

//----------------------------------------------------------------------

/// Evaluate the success of the solve
int AMRsolve_Hypre_Grav::evaluate()
{
  // check whether solution/rhs was requested
  if (parameters_->value("dump_x") == "true" || 
      parameters_->value("dump_b") == "true") {

    // iterate over processor-local grids
    ItHierarchyGridsLocal itg(*hierarchy_);
    while (AMRsolve_Grid * grid = itg++) {
      
      char filename[80];

      // get level & grid information
      int level = grid->level();
      int lower[3],upper[3];
      grid->get_limits(lower,upper);

      if (parameters_->value("dump_x") == "true") {
	// extract Enzo solution
	int nx[3];
	HYPRE_SStructVectorGather(X_);
	HYPRE_SStructVectorGetBoxValues(X_, level, lower, upper, 0,
					grid->get_u(&nx[0],&nx[1],&nx[2]));  
	sprintf(filename,"X.%d",grid->id());
	grid->write("header",filename);
	grid->write("u",filename);
      }
    
      if (parameters_->value("dump_b") == "true") {
	// extract Enzo rhs
	int nb[3];
	HYPRE_SStructVectorGather(B_);
	HYPRE_SStructVectorGetBoxValues(B_, level, lower, upper, 0,
				    grid->get_f(&nb[0],&nb[1],&nb[2]));  
	sprintf(filename,"B.%d",grid->id());
	grid->write("f",filename);
      }

    } // grid = itg++
  } // if dump_x or dump_b

  // check for error flags; output info to stdout, if requested
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
      printf("Stalled: %d >= %d\n", iterations(),itmax);
    err_flag = 1;
  }

//   // Appears to have completed successfully
//   if (err_flag == 0)  
//     if (pmpi->is_root()) {
//       printf("AMRsolve_Hypre_FLD Success!\n"); 
//       fflush(stdout); 
//     }

  return err_flag;

} // AMRsolve_Hypre_Grav::evaluate()


//----------------------------------------------------------------------

/// Extracts HYPRE solution and updates potential field
void AMRsolve_Hypre_Grav::update_enzo()
{
  // iterate over grids on this processor
  ItHierarchyGridsLocal itg(*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // get level, grid information
    int level = grid->level();
    int lower[3],upper[3];
    grid->get_limits(lower,upper);

    // extract Enzo solution
    int n3[3];
    Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
    HYPRE_SStructVectorGather(X_);
    HYPRE_SStructVectorGetBoxValues(X_, level, lower, upper, 0, u);  

    // access Enzo PotentialField
    Scalar* phi = grid->get_phi();

    // get buffering information on relating amrsolve grid to Enzo data
    int ghosts[3][2]; 
    grid->get_GravGhosts(ghosts);
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
	  phi[k] += u[i];
	}
      }
    }

  } // grid = itg++

} // AMRsolve_Hypre_Grav::update_enzo()


//----------------------------------------------------------------------

/// dumps HYPRE matrix and RHS (called when aborting solve)
void AMRsolve_Hypre_Grav::abort_dump()
{

  // have HYPRE dump out everything it knows to disk
  HYPRE_SStructMatrixPrint("A.mat",A_,0);
  HYPRE_SStructVectorPrint("x.vec",X_,0);
  HYPRE_SStructVectorPrint("b.vec",B_,0);

} // AMRsolve_Hypre_Grav::abort_dump()


//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

/// init_nonstencil_() is called twice: first by
/// init_graph_nonstencil_() with phase == phase_graph to set nonstencil
/// graph entries, and again by init_matrix_nonstencil_() with phase
/// == phase_matrix to set nonstencil matrix entries.
void AMRsolve_Hypre_Grav::init_nonstencil_(AMRsolve_Grid& grid, phase_enum phase)
{
  // Input parameter check
  if ( !(phase == phase_graph || phase == phase_matrix) ) {
    char error_message[80];
    sprintf(error_message,"init_matrix_nonstencil_ called with phase = %d",
	    int(phase));
    ERROR(error_message);
  } // if phase unexpected

  // global grid index limits
  int index_global[3][2];
  grid.indices(index_global);
  
  // Determine the discretization method (constant, linear, quadratic)
  // Currently only constant is implemented
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

    if (!l0) printf("index_global[%d][0] = %d\n",axis1,index_global[axis1][0]);
    assert(l0);
    if (!l1) printf("index_global[%d][1] = %d\n",axis1,index_global[axis1][1]);
    assert(l1);

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
	    index_coarse[axis0] = (index_fine[axis0]) / r_factor_  + (face*r_factor_-1);
	    index_coarse[axis1] = (index_fine[axis1]) / r_factor_;
	    index_coarse[axis2] = (index_fine[axis2]) / r_factor_;

	    // adjust for periodicity
	    if (hierarchy_->is_periodic(axis0)) {
	      int period = hierarchy_->period_index(axis0,level_coarse);
 	      index_coarse[axis0] = (index_coarse[axis0] + period) % period;
 	    }

	    //--------------------------------------------------
	    // GRAPH ENTRY: FINE-TO-COARSE 
	    //--------------------------------------------------

	    if (discret_type == discret_type_const) {

	      update_fine_coarse_const_(face,grid,axis0,phase,
					level_fine,level_coarse,
					index_fine,index_coarse);

	      //--------------------------------------------------
	      // GRAPH ENTRY: COARSE-TO-FINE
	      //--------------------------------------------------
	      
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
} // AMRsolve_Hypre_Grav::init_nonstencil_()

//------------------------------------------------------------------------

/// Initialize matrix stencil and graph entries
void AMRsolve_Hypre_Grav::init_elements_matrix_()
{

  ItHierarchyLevels itl (*hierarchy_);
  while (AMRsolve_Level * level = itl++) {

    int part = level->index();

    // 1. Set stencil values within level
    ItLevelGridsLocal itlg (*level);
    while (AMRsolve_Grid * grid = itlg++)  init_matrix_stencil_(*grid);

    if (part > 0) {
      // *** WARNING: POSSIBLE SCALING ISSUE.  Below we loop over all
      // *** grids; however, we only need to loop over "parent-child
      // *** pairs such that either child or parent is local to this
      // *** MPI process."
 
      // Set matrix values between levels
      ItLevelGridsAll itag (*level);
      while (AMRsolve_Grid * grid = itag++)  init_matrix_nonstencil_(*grid);

    } // while level > 0
  } // while level = itl++

  // Clean up stencil connections between levels
  for (int part = 1; part < hierarchy_->num_levels(); part++) 
    init_matrix_clear_(part);

} // init_elements_matrix_()

//------------------------------------------------------------------------

/// Set right-hand-side elements
/*void AMRsolve_Hypre_Grav::init_elements_rhs_(std::vector<AMRsolve_Point *>& points,
                                             Scalar f_scale) */
void AMRsolve_Hypre_Grav::init_elements_rhs_(Scalar f_scale)
{
  // declare local variables
  Scalar local_shift_b_sum = 0.0;
  long long shift_b_count  = 0.0;

  _TRACE_;
  // iterate over all grids local to this processor
  ItHierarchyGridsLocal itg (*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // for each local grid, use our "f" array to store RHS entries 
    // for the linear system.
    int n0,n1,n2;
    Scalar* values = grid->get_f(&n0,&n1,&n2);

    // access Enzo's GravitatingMassField
    Scalar* GravMassField = grid->get_gmass();

    // get buffering information on relating amrsolve grid to Enzo data
    int ghosts[3][2]; 
    grid->get_GravGhosts(ghosts);
    int en0 = n0 + ghosts[0][0] + ghosts[0][1];  // enzo data dimensions
    int en1 = n1 + ghosts[1][0] + ghosts[1][1];  // enzo data dimensions
    int en2 = n2 + ghosts[2][0] + ghosts[2][1];  // enzo data dimensions

    // iterate over the domain, copying Enzo's GravitatingMassField into f
    int k0, k1, k2, k, i0, i1, i2, i;
    for (i2=0; i2<n2; i2++) {
      k2 = ghosts[2][0] + i2;

      for (i1=0; i1<n1; i1++) {
	k1 = ghosts[1][0] + i1;

	for (i0=0; i0<n0; i0++) {
	  k0 = ghosts[0][0] + i0;

	  // compute indices of amrsolve, enzo cells
	  k = k0 + en0*(k1 + en1*k2);
	  i = i0 + n0*(i1 + n1*i2);

	  // fill in f with Enzo RHS data
	  values[i] = GravMassField[k];
	}
      }
    }
   
    // set RHS scaling to input scaling multiplied by 1/(hx*hy*hz)
    Scalar hx,hy,hz;
    grid->h(hx,hy,hz);
    Scalar scale = hx*hy*hz * f_scale;
    for (int i=0; i<n0*n1*n2; i++)  values[i] *= scale;

    // Set Hypre B_ vector to grid f_ values
    int part = grid->level();
    int lower[3] = { grid->index_lower(0), 
		     grid->index_lower(1), 
		     grid->index_lower(2) };
    int upper[3] = { grid->index_upper(0), 
		     grid->index_upper(1), 
		     grid->index_upper(2) };
    HYPRE_SStructVectorAddToBoxValues(B_,part,lower,upper,0,values);

  } //  while grid = itg++

  // for periodic BCs, need to shift RHS to have zero average
  // value to deflate the null space
  if ( parameters_->value("boundary") == "periodic" ) {

    // Clear under overlapped areas
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
      HYPRE_SStructFACZeroAMRVectorData(B_, plevels, rfactors);
    }

    // Accumulate local sums
    local_shift_b_sum = 0.0;

    // Compute total number of variables (excluding overlap)
    shift_b_count = 0; 
    const int r_factor3 = r_factor_*r_factor_*r_factor_;

    ItHierarchyLevels itl (*hierarchy_);
    while (AMRsolve_Level * level = itl++) {
      
      int part = level->index();
      ItLevelGridsAll itgl (*level);
      while (AMRsolve_Grid * grid = itgl++) {

	// Adjust shift_b_count: add grid; subtract overlap
	shift_b_count += grid->num_unknowns();
	if (part > 0)  shift_b_count -= grid->num_unknowns() / r_factor3;

	if (grid->is_local()) {

 	  // Get grid info
 	  int lower[3],upper[3];
	  grid->get_limits(lower, upper);
	  int n = grid->num_unknowns();

	  // Create space for the patch
 	  double* tmpvals = new double[n];
 	  for (int i=0; i<n; i++)  tmpvals[i] = 0.0;

 	  // Copy vector values to the array
 	  HYPRE_SStructVectorGetBoxValues(B_, part, lower, upper, 0, tmpvals);

	  // Accumulate the sum
	  for (int i=0; i<n; i++)  local_shift_b_sum += tmpvals[i];

	  // Delete the zeroed values
 	  delete [] tmpvals;
	}
      }
    }

    // Get global sum from local sums
    Scalar shift_b_sum = 0.0;
    MPI_Allreduce(&local_shift_b_sum, &shift_b_sum, 1, 
		  MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD);

    // Compute the shift given the global sum and global count
    Scalar shift_b_amount = -shift_b_sum/shift_b_count;

    // Perform the shift of B 
    ItHierarchyLevels itl2 (*hierarchy_);
    while (AMRsolve_Level* level = itl2++) {

      int part = level->index();

      ItLevelGridsLocal itg2 (*level);
      while (AMRsolve_Grid* grid = itg2++) {

	int lower[3],upper[3];
	grid->get_limits(lower, upper);
	Scalar* tmpvals = new Scalar[grid->num_unknowns()];

	for (int i=0; i<grid->num_unknowns(); i++) tmpvals[i] = shift_b_amount;

	HYPRE_SStructVectorAddToBoxValues(B_,part,lower,upper,0,tmpvals);

	delete [] tmpvals;

      } // while grid = itg2++
    } // while level = itl2++

    // Re-clear under overlapped area that got shifted
    nlevels = hierarchy_->num_levels();
    if (nlevels > 1) {
      int plevels[nlevels];
      int rfactors[nlevels][3];
      for (int part=0; part<nlevels; part++) {
	plevels[part] = part;
	rfactors[part][0] = r_factor_;
	rfactors[part][1] = r_factor_;
	rfactors[part][2] = r_factor_;
      }
      HYPRE_SStructFACZeroAMRVectorData(B_, plevels, rfactors);
    }

  } // if periodic

} // init_elements_rhs_()

//------------------------------------------------------------------------

/// Set matrix stencil values for the grid interior
void AMRsolve_Hypre_Grav::init_matrix_stencil_(AMRsolve_Grid& grid)
{
  int n          = grid.num_unknowns();
  int entries[7] = { 0,1,2,3,4,5,6 };
  double h3[3]   = {grid.h(0),grid.h(1),grid.h(2)};
  int    n3[3]   = {grid.n(0),grid.n(1),grid.n(2)};

  double h120 = h3[1]*h3[2] / h3[0];
  double h201 = h3[2]*h3[0] / h3[1];
  double h012 = h3[0]*h3[1] / h3[2];

  double* v0;         // Diagonal elements
  double* v1[3][2];   // Off-diagonal elements

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

//  WHERE; printf("h120  h201  h012 = %g %g %g\n",h120, h201, h012);

  int i0,i1,i2,i;
  for (i2=0; i2<n3[2]; i2++) {
    for (i1=0; i1<n3[1]; i1++) {
      for (i0=0; i0<n3[0]; i0++) {

	// get the linear grid index
	i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);

	// set the matrix entries
	v1[0][0][i] = matrix_scale_ * h120;
	v1[0][1][i] = matrix_scale_ * h120;
	v1[1][0][i] = matrix_scale_ * h201;
	v1[1][1][i] = matrix_scale_ * h201;
	v1[2][0][i] = matrix_scale_ * h012;
	v1[2][1][i] = matrix_scale_ * h012;

	v0[i] = -( v1[0][0][i] + v1[0][1][i] +
		   v1[1][0][i] + v1[1][1][i] + 
		   v1[2][0][i] + v1[2][1][i] );

      } // for i0
    } // for i1
  } // for i2

//  WHERE; printf("v0[0]=%g\n",v0[0]);

  //-----------------------------------------------------------
  // Adjust stencil at grid boundaries (not implemented)
  //  -> if we want to allow Neumann/Dirichlet conditions, 
  //     insert code here
  //-----------------------------------------------------------


  //-----------------------------------------------------------
  // insert matrix entries into Hypre matrix A_
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
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper,
				  0, 1, &entries[0], v0);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[1], v1[0][1]);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[2], v1[0][0]);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[3], v1[1][1]);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[4], v1[1][0]);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[5], v1[2][1]);
  HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
				  0, 1, &entries[6], v1[2][0]);

  // Deallocate arrays
  delete [] v0;
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      delete [] v1[axis][face];
    } // for face
  } // for axis

} // AMRsolve_Hypre_Grav::init_matrix_stencil_()

//------------------------------------------------------------------------

/// Clean up stencil connections between parts for FAC solver
void AMRsolve_Hypre_Grav::init_matrix_clear_(int part)
{
  if (part > 0) {
    int r_factors[3] = {r_factor_,r_factor_,r_factor_}; 
    HYPRE_SStructFACZeroAMRMatrixData(A_, part-1, r_factors);
  }

} // AMRsolve_Hypre_Grav::init_matrix_clear_()

//------------------------------------------------------------------------

// /// Add contributions from point sources to right-hand side B
// Scalar AMRsolve_Hypre_Grav::init_vector_points_(std::vector<AMRsolve_Point *>& points)
// {
//   const Scalar scaling0 = -4.0*const_G*const_pi;

//   Scalar shift_b_sum = 0.0;

//   int i;
//   for (i=0; i<int(points.size()); i++) {

//     AMRsolve_Point& point = *points[i];
//     AMRsolve_Grid& grid = hierarchy_->return_grid(point.igrid());
//     if (grid.is_local()) {

//       Scalar cell_volume = grid.h(0) * grid.h(1) * grid.h(2);
//       Scalar density     = point.mass() / cell_volume;
//       Scalar value       = scaling0 * density;

//       // Add contribution of the point to the right-hand side vector
//       int index[3];
//       Scalar lower[3],upper[3];
//       grid.x_lower(lower[0],lower[1],lower[2]);
//       grid.x_upper(upper[0],upper[1],upper[2]);
//       for (int k=0; k<3; k++) {
// 	Scalar ap = point.x(k)      - lower[k];
// 	Scalar ag = upper[k] - lower[k];
// 	int    ig = grid.num_unknowns(k);
// 	int    i0 = grid.index_lower(k);
// 	index[k] = int (ap/ag*ig) + i0;
//       } // for k=0,2
//       if (index[0] < grid.index_lower(0) || grid.index_upper(0) < index[0] ||
// 	  index[1] < grid.index_lower(1) || grid.index_upper(1) < index[1] ||
// 	  index[2] < grid.index_lower(2) || grid.index_upper(2) < index[2]) {
// 	printf("WARNING: Point apparently not in grid: \n");
// 	printf("WARNING:    Point: (%g,%g,%g)\n",
// 	       point.x(0),point.x(1),point.x(2));
// 	printf("WARNING:    Grid:  (%g,%g,%g) - (%g,%g,%g)\n",
// 	       lower[0],lower[1],lower[2],upper[0],upper[1],upper[2]);
//       } // if index not in grid
//       if (debug) {
// 	point.print();
// 	grid.print();
// 	printf("Point index  = %d %d %d)\n",index[0],index[1],index[2]);
// 	printf("Cell size    = %g %g %g\n",grid.h(0),grid.h(1),grid.h(2));
// 	printf("Cell volume  = %g\n",cell_volume);
// 	printf("Cell density = %g\n",density);
// 	printf("RHS contribution = %g\n",value);
//       } // if debug
    
//       shift_b_sum += value;
//       HYPRE_SStructVectorAddToValues(B_, grid.level(), index, 0, &value);

//     } // if grid.is_local()
//   } // for i=0 to # points

//   return shift_b_sum;

// } // AMRsolve_Hypre_Grav::init_vector_points_()

// //------------------------------------------------------------------------

// /// Add contributions from Density in enzo HDF5 files to right-hand side B
// Scalar AMRsolve_Hypre_Grav::init_vector_file_(std::string file_prefix,
// 					      bool enzo_packed)
// {
//   ItHierarchyGridsLocal itg (*hierarchy_);
//   char error_message[80];

//   herr_t status;
//   hid_t  file_id;
//   hid_t  dataset_id;

//   Scalar shift_b_sum = 0.0;

//   while (AMRsolve_Grid* grid = itg++) {

//     // Open the HDF5 grid file
//     char grid_num_str[10];
//     sprintf(grid_num_str,"%04d",grid->id() + 1);
//     std::string file_name = file_prefix + ".grid" + grid_num_str;
//     file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
//     if (file_id < 0) {
//       strcpy(error_message,"H5Fopen cannot open file ");
//       strcat(error_message,file_name.c_str());
//       ERROR(error_message);
//     } else {
//       printf("DEBUG %s:%d %s opened successfully\n",
// 	     __FILE__,__LINE__,file_name.c_str());
//     }

//     // Open the dataset Density
//     dataset_id = H5Dopen(file_id, "Density");
//     if (dataset_id < 0) {
//       strcpy(error_message,"H5Dopen cannot open dataset Density");
//       ERROR(error_message);
//     } else {
//       printf("DEBUG %s:%d Density opened successfully\n",
// 	     __FILE__,__LINE__);
//     }

//     // Read the dataset
//     double* values = new double[grid->n()];
//     status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
// 		     H5P_DEFAULT, values);
//     if (status < 0) {
//       strcpy(error_message,"H5Dread exited with status ");
//       char status_str[10];
//       sprintf(status_str,"%d",status);
//       strcat(error_message,status_str);
//       ERROR(error_message);
//     } else {
//       printf("DEBUG %s:%d Density read successfully\n",
// 	     __FILE__,__LINE__);
//       printf("%s:%d %d %d %d  %d  [%g %g]\n",
// 	     __FILE__,__LINE__,
// 	     grid->n(0),grid->n(1),grid->n(2),grid->n(),
// 	     values[0],values[grid->n()-1]);
//     }

//     // Copy the values to the hypre vector
//     int part = grid->level();
//     int lower[3] = { grid->index_lower(0), 
// 		     grid->index_lower(1), 
// 		     grid->index_lower(2) };
//     int upper[3] = { grid->index_upper(0), 
// 		     grid->index_upper(1), 
// 		     grid->index_upper(2) };

//     // Set Hypre B_ vector to values from the file
//     HYPRE_SStructVectorAddToBoxValues(B_,part,lower,upper,0,values);

//     // Compute sum to return in case we need to shift rhs for periodic problems
//     int n0 = grid->num_unknowns(0);
//     int n1 = grid->num_unknowns(1);
//     int n2 = grid->num_unknowns(2);

//     for (int i0=0; i0<n0; i0++) {
//       for (int i1=0; i1<n1; i1++) {
// 	for (int i2=0; i2<n2; i2++) {
// 	  int k = i0 + n0*(i1 + n1*i2);
// 	  shift_b_sum += values[k];
// 	}
//       }
//     }

//     delete [] values;

//     // Close the HDF5 grid file
//     status = H5Fclose(file_id);
//     if (status < 0) {
//       strcpy(error_message,"H5Fclose exited with status ");
//       char status_str[10];
//       sprintf(status_str,"%d",status);
//       strcat(error_message,status_str);
//       ERROR(error_message);
//     }
    
//   }

//   return 0.0;
// }

// //------------------------------------------------------------------------

// /// Use existing f_ for right-hand side B
// Scalar AMRsolve_Hypre_Grav::init_vector_attach_(Scalar f_scale)
// {
//   _TRACE_;
//   Scalar shift_b_sum = 0.0;
//   ItHierarchyGridsLocal itg (*hierarchy_);
//   while (AMRsolve_Grid* grid = itg++) {
//     int n0,n1,n2;
//     Scalar* values = grid->get_f(&n0,&n1,&n2);

//     // Scale by 1/(hx*hy*hz)
//     Scalar hx,hy,hz;
//     grid->h(hx,hy,hz);
//     Scalar scale = hx*hy*hz * f_scale;
//     for (int i0=0; i0<n0; i0++) {
//       for (int i1=0; i1<n1; i1++) {
// 	for (int i2=0; i2<n2; i2++) {
// 	  int k = i0 + n0*(i1 + n1*i2);
// 	  values[k] *= scale;
// 	}
//       }
//     }

//     // Set Hypre B_ vector to grid f_ values
//     int part = grid->level();
//     int lower[3] = { grid->index_lower(0), 
// 		     grid->index_lower(1), 
// 		     grid->index_lower(2) };
//     int upper[3] = { grid->index_upper(0), 
// 		     grid->index_upper(1), 
// 		     grid->index_upper(2) };
//     HYPRE_SStructVectorAddToBoxValues(B_,part,lower,upper,0,values);

//     // Compute sum to return in case we need to shift rhs for periodic problems
//     for (int i0=0; i0<n0; i0++) {
//       for (int i1=0; i1<n1; i1++) {
// 	for (int i2=0; i2<n2; i2++) {
// 	  int k = i0 + n0*(i1 + n1*i2);
// 	  shift_b_sum += values[k];
// 	}
//       }
//     }
//   }

//   return shift_b_sum;
// }

//------------------------------------------------------------------------

/// Initialize the PFMG hypre solver
void AMRsolve_Hypre_Grav::solve_pfmg_(int itmax, double restol)
{
  _TRACE_;
  // Create and initialize the solver
  HYPRE_SStructSysPFMGCreate(MPI_COMM_WORLD, &solver_);

  // extract some additional solver parameters
  std::string srlxtype = parameters_->value("solver_rlxtype");
  std::string snpre    = parameters_->value("solver_npre");
  std::string snpost   = parameters_->value("solver_npost");
  std::string sprintl  = parameters_->value("solver_printl");
  std::string slog     = parameters_->value("solver_log");

  //   if not defined, then define them
  if (srlxtype == "")  srlxtype = "1";
  if (snpre == "")     snpre = "1";
  if (snpost == "")    snpost = "1";
  if (sprintl == "")   sprintl = "1";
  if (slog == "")      slog = "1";

  //   set local variables
  int rlxtype = atoi(srlxtype.c_str());
  int npre    = atoi(snpre.c_str());
  int npost   = atoi(snpost.c_str());
  int printl  = atoi(sprintl.c_str());
  int log     = atoi(slog.c_str());

  // set solver options
  if (itmax != 0 )   HYPRE_SStructSysPFMGSetMaxIter(solver_,itmax);
  if (restol != 0.0) HYPRE_SStructSysPFMGSetTol(solver_,    restol);
  HYPRE_SStructSysPFMGSetRelaxType(solver_,    rlxtype);
  HYPRE_SStructSysPFMGSetNumPreRelax(solver_,  npre);
  HYPRE_SStructSysPFMGSetNumPostRelax(solver_, npost);
  HYPRE_SStructSysPFMGSetPrintLevel(solver_,   printl);
  HYPRE_SStructSysPFMGSetLogging(solver_,      log);

  // setup solver 
  HYPRE_SStructSysPFMGSetup(solver_,A_,B_,X_);

  // Solve the linear system
  HYPRE_SStructSysPFMGSolve(solver_,A_,B_,X_);

  // Write out some diagnostic info about the solve
  HYPRE_SStructSysPFMGGetNumIterations(solver_,&iter_);
  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm(solver_,&resid_);
  if (debug && pmpi->is_root()) {
    printf("hypre PFMG num iterations: %d\n",iter_);
    printf("hypre PFMG final relative residual norm: %g\n",resid_);
  }

  // Delete the solver
  HYPRE_SStructSysPFMGDestroy(solver_);

  _TRACE_;
} // AMRsolve_Hypre_Grav::solve_pfmg_()

//------------------------------------------------------------------------

/// Initialize the FAC hypre solver
void AMRsolve_Hypre_Grav::solve_fac_(int itmax, double restol)
{
  _TRACE_;
  int i;

  // Create the solver
  HYPRE_SStructFACCreate(MPI_COMM_WORLD, &solver_);

  // Initialize parts
  int num_parts = hierarchy_->num_levels();
  HYPRE_SStructFACSetMaxLevels(solver_, num_parts);
  int *parts = new int[num_parts];
  for (i=0; i<num_parts; i++) parts[i] = i;
  HYPRE_SStructFACSetPLevels(solver_, num_parts, parts);

  // Initialize refinement factors
  int3 *refinements = new int3[num_parts];
  for (i=0; i<num_parts; i++) {
    refinements[i][0] = r_factor_;
    refinements[i][1] = r_factor_;
    refinements[i][2] = r_factor_;
  }
  HYPRE_SStructFACSetPRefinements(solver_, num_parts, refinements);

  // extract some additional solver parameters
  std::string srlxtype = parameters_->value("solver_rlxtype");
  std::string snpre    = parameters_->value("solver_npre");
  std::string snpost   = parameters_->value("solver_npost");
  std::string scsolve  = parameters_->value("solver_csolve");
  std::string sprintl  = parameters_->value("solver_printl");
  std::string slog     = parameters_->value("solver_log");

  //   if not defined, then define them
  if (srlxtype == "")  srlxtype = "2";
  if (snpre == "")     snpre = "2";
  if (snpost == "")    snpost = "2";
  if (scsolve == "")   scsolve = "1";
  if (sprintl == "")   sprintl = "1";
  if (slog == "")      slog = "1";

  //   set local variables
  int relax   = atoi(srlxtype.c_str());
  int npre    = atoi(snpre.c_str());
  int npost   = atoi(snpost.c_str());
  int printl  = atoi(sprintl.c_str());
  int log     = atoi(slog.c_str());
  int csolve  = atoi(scsolve.c_str());

  // solver parameters
  HYPRE_SStructFACSetNumPreRelax(solver_,      npre);
  HYPRE_SStructFACSetNumPostRelax(solver_,     npost);
  HYPRE_SStructFACSetCoarseSolverType(solver_, csolve);
  HYPRE_SStructFACSetRelaxType(solver_,        relax);

  // stopping criteria
  if (itmax != 0 )    HYPRE_SStructFACSetMaxIter(solver_, itmax);
  if (restol != 0.0)  HYPRE_SStructFACSetTol(solver_,     restol);

  // output amount
  HYPRE_SStructFACSetLogging(solver_, 1);

  // prepare for solve
  HYPRE_SStructFACSetup2(solver_, A_, B_, X_);

  // Solve the linear system
  HYPRE_SStructFACSolve3(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve
  HYPRE_SStructFACGetNumIterations(solver_, &iter_);
  HYPRE_SStructFACGetFinalRelativeResidualNorm(solver_, &resid_);
  if (debug && pmpi->is_root()) {
    printf("hypre FAC num iterations: %d\n",iter_);
    printf("hypre FAC final relative residual norm: %g\n",resid_);
  }

  // Delete the solver
  HYPRE_SStructFACDestroy2(solver_);

  // Delete local dynamic storage
  delete [] parts;
  delete [] refinements;

  _TRACE_;
} // AMRsolve_Hypre_Grav::solve_fac_()

//------------------------------------------------------------------------

/// Initialize the BICGSTAB hypre solver
void AMRsolve_Hypre_Grav::solve_bicgstab_(int itmax, double restol)
{
  _TRACE_;
  // Create the solver
  HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver_);

  // extract some additional solver parameters
  std::string slog = parameters_->value("solver_log");
  if (slog == "")      slog = "1";
  int log = atoi(slog.c_str());

  // stopping criteria
  if (itmax != 0 )   HYPRE_SStructBiCGSTABSetMaxIter(solver_,itmax);
  if (restol != 0.0) HYPRE_SStructBiCGSTABSetTol(solver_,    restol);

  // output amount
  HYPRE_SStructBiCGSTABSetLogging(solver_, 1);

  // Initialize the solver
  HYPRE_SStructBiCGSTABSetup(solver_, A_, B_, X_);

  // Solve the linear system
  HYPRE_SStructBiCGSTABSolve(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve
  HYPRE_SStructBiCGSTABGetNumIterations(solver_, &iter_);
  HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver_, &resid_);
  if (debug && pmpi->is_root()) {
    printf("hypre BiCGSTAB num iterations: %d\n",iter_);
    printf("hypre BiCGSTAB final relative residual norm: %g\n",resid_);
  }

  // Delete the solver
  HYPRE_SStructBiCGSTABDestroy(solver_);

  _TRACE_;
} // AMRsolve_Hypre_Grav::solve_bicgstab_()

//------------------------------------------------------------------------

/// Initialize the BiCGSTAB hypre solver with BoomerAMG preconditioning
void AMRsolve_Hypre_Grav::solve_bicgstab_boomer_(int itmax, double restol)
{
  _TRACE_;
  HYPRE_Solver solver;

  // Create the solver
  HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &solver);

  // stopping criteria
  if (itmax != 0 )   HYPRE_BiCGSTABSetMaxIter(solver,itmax);
  if (restol != 0.0) HYPRE_BiCGSTABSetTol(solver,    restol);

  // extract some additional solver parameters
  std::string slog    = parameters_->value("solver_log");
  std::string sprintl = parameters_->value("solver_printl");
  if (slog == "")      slog = "1";
  if (sprintl == "")   sprintl = "1";
  int log    = atoi(slog.c_str());
  int printl = atoi(sprintl.c_str());

  // output amount
  HYPRE_BiCGSTABSetLogging(solver, log);

  // Set BoomerAMG preconditioner
  HYPRE_Solver par_precond;
  HYPRE_BoomerAMGCreate(&par_precond);
  HYPRE_BoomerAMGSetCoarsenType(par_precond, 6);
  HYPRE_BoomerAMGSetStrongThreshold(par_precond, 0.25);
  HYPRE_BoomerAMGSetTol(par_precond, 0.0);
  HYPRE_BoomerAMGSetPrintLevel(par_precond, 1);
  HYPRE_BoomerAMGSetPrintFileName(par_precond, "sstruct.out.log");
  HYPRE_BoomerAMGSetMaxIter(par_precond, 1);
  HYPRE_BiCGSTABSetPrecond(solver,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
			   par_precond);

  // Initialize the solver
  HYPRE_BiCGSTABSetup(solver, (HYPRE_Matrix) A_, 
		      (HYPRE_Vector) B_, (HYPRE_Vector) X_);

  // Solve the linear system
  HYPRE_BiCGSTABSolve(solver, (HYPRE_Matrix) A_, 
		      (HYPRE_Vector)B_, (HYPRE_Vector) X_);

  // Write out some diagnostic info about the solve
  HYPRE_BiCGSTABGetNumIterations(solver, &iter_);
  HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver, &resid_);
  if (debug &&pmpi->is_root()) {
    printf("hypre BiCGSTAB num iterations: %d\n",iter_);
    printf("hypre BiCGSTAB final relative residual norm: %g\n",resid_);
  }


  // Delete the solver and preconditiner
  HYPRE_BoomerAMGDestroy(par_precond);
  HYPRE_BiCGSTABDestroy(solver);

  _TRACE_;
} // AMRsolve_Hypre_Grav::solve_bicgstab_boomer_()

//------------------------------------------------------------------------

/// Initialize the GMRES hypre solver
void AMRsolve_Hypre_Grav::solve_gmres_(int itmax, double restol)
{
  _TRACE_;
  // Create the solver
  HYPRE_SStructGMRESCreate(MPI_COMM_WORLD, &solver_);

  // extract some additional solver parameters
  std::string slog = parameters_->value("solver_log");
  if (slog == "")      slog = "1";
  int log = atoi(slog.c_str());

  // stopping criteria
  if (itmax != 0 )   HYPRE_SStructGMRESSetMaxIter(solver_, itmax);
  if (restol != 0.0) HYPRE_SStructGMRESSetTol(solver_,     restol);

  // output amount
  HYPRE_SStructGMRESSetLogging(solver_, log);

  // Initialize the solver
  HYPRE_SStructGMRESSetup(solver_, A_, B_, X_);

  // Solve the linear system
  HYPRE_SStructGMRESSolve(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve
  HYPRE_SStructGMRESGetNumIterations(solver_, &iter_);
  HYPRE_SStructGMRESGetFinalRelativeResidualNorm(solver_, &resid_);
  if (debug && pmpi->is_root()) {
    printf("hypre GMRES num iterations: %d\n",iter_);
    printf("hypre GMRES final relative residual norm: %g\n",resid_);
  }

  // Delete the solver
  HYPRE_SStructGMRESDestroy(solver_);

  _TRACE_;
} // AMRsolve_Hypre_Grav::solve_gmres_()

//------------------------------------------------------------------------

/// Update the matrix at a coarse/fine interface; this routine changes the
/// fine-grid matrix to account for teh coarse grid neighbor.  Called in
/// two phases, with phase == phase_graph (via init_graph_nonstencil_()) 
/// for the nonzero structure, and with phase == phase_matrix 
/// (via init_matrix_nonstencil_()) for the matrix nonzeros.
void AMRsolve_Hypre_Grav::update_fine_coarse_const_(int face, 
						    AMRsolve_Grid& grid_fine, 
						    int axis0, 
						    phase_enum phase,
						    int level_fine, 
						    int level_coarse,
						    int index_fine[3], 
						    int index_coarse[3])
{
  int axis1 = (axis0+1)%3;
  int axis2 = (axis0+2)%3;

  //--------------------------------------------------
  // (*) CONSTANT
  //     Scale        = 2/3
  //     Coefficients = 1
  //--------------------------------------------------

  Scalar val_h_fine = grid_fine.h(axis1) * grid_fine.h(axis2) / grid_fine.h(axis0);

  int index_increment[][3] = {{face*(r_factor_-1),0,0},
			      {0,1,0},
			      {0,0,1},
			      {0,-1,0},
			      {-face*(r_factor_-1),0,-1}};

  if (grid_fine.is_local()) {

    if (phase == phase_graph) {

      int k = 0;
      index_fine[axis0] += index_increment[k][0];
      index_fine[axis1] += index_increment[k][1];
      index_fine[axis2] += index_increment[k][2];

      for (k=1; k<5; k++) {
	HYPRE_SStructGraphAddEntries(graph_, level_fine, index_fine, 
				     0, level_coarse, index_coarse, 0);
	index_fine[axis0] += index_increment[k][0];
	index_fine[axis1] += index_increment[k][1];
	index_fine[axis2] += index_increment[k][2];
      } // for k = 1:4

    } else if (phase == phase_matrix) {

      // fine->coarse off-diagonal
      double val_s = 2.0 / 3.0;
      int entry;
      double val,value;

      int k=0;
      index_fine[axis0] += index_increment[k][0];
      index_fine[axis1] += index_increment[k][1];
      index_fine[axis2] += index_increment[k][2];

      for (k=1; k<5; k++) {

	// Set matrix entry
	val = matrix_scale_ * val_h_fine * val_s;

	// Update off-diagonal
	entry = grid_fine.counter(index_fine)++;
	value = val;
	HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
				       0, 1, &entry, &value);

	// Update diagonal
	entry = 0;
	value = -val;
	HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
				       0, 1, &entry, &value);

	// Clear coarse-fine stencil values
	val = matrix_scale_ * val_h_fine;

	//   Update off-diagonal, stencil xp=1,xm,yp,ym,zp,zm=6
	entry = 2*axis0 + 1 + (1-face);
	value = -val;
	HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
				       0, 1, &entry, &value);

	//   Update diagonal
	entry = 0;
	value = val;
	HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
				       0, 1, &entry, &value);

	// Update indices
	index_fine[axis0] += index_increment[k][0];
	index_fine[axis1] += index_increment[k][1];
	index_fine[axis2] += index_increment[k][2];

      } // for k = 1,4
    } // if phase == phase_matrix
  } // if grid_fine.is_local()
} // if discret_type_const

//------------------------------------------------------------------------

/// Update the matrix at a coarse/fine interface; this routine changes the 
/// coarse-grid matrix to account for the fine grid neighbors.  Called in 
/// two phases, with phase == phase_graph (via init_graph_nonstencil_()) 
/// for the nonzero structure, and with phase == phase_matrix
/// (via init_matrix_nonstencil_()) for the matrix nonzeros.
void AMRsolve_Hypre_Grav::update_coarse_fine_const_(int face, 
						    AMRsolve_Grid& grid_coarse, 
						    int axis0, 
						    phase_enum phase,
						    int level_fine, 
						    int level_coarse,
						    int index_fine[3], 
						    int index_coarse[3])
{
  int axis1 = (axis0+1)%3;
  int axis2 = (axis0+2)%3;

  Scalar val_h_coarse = grid_coarse.h(axis1) * grid_coarse.h(axis2) / grid_coarse.h(axis0);

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
      HYPRE_SStructGraphAddEntries(graph_, level_coarse, index_coarse, 
				   0, level_fine, index_fine, 0);
      index_fine[0] += index_increment[k][0];
      index_fine[1] += index_increment[k][1];
      index_fine[2] += index_increment[k][2];
    } // for k=0,7

  // set matrix entry
  } else if (phase == phase_matrix) {

    // Set matrix entry using diffusion coefficient 
    // (x,y,z) == (x3[0],x3[1],x3[2])
    double val_s = 1.0 / 8.0;
    double val   = matrix_scale_ * val_h_coarse * val_s;
    int    entry;
    double value;

    // Adjust coarse-fine nonstencil values
    for (int i=0; i<8; i++) {
      // Set new nonstencil coarse-fine entry
      entry = grid_coarse.counter(index_coarse)++;
      value = val;
      HYPRE_SStructMatrixAddToValues(A_, level_coarse,index_coarse, 
				     0, 1, &entry, &value);
      // Adjust stencil diagonal
      entry = 0;
      value = -val;
      HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
				     0, 1, &entry, &value);
    } // for i=0,7

    // Clear original matrix values from stencil
    val = matrix_scale_ * val_h_coarse;

    //   Update off-diagonal, stencil xp=1,xm,yp,ym,zp,zm=6
    //   (note: "face" is for fine grid, but we want coarse)
    entry = 2*axis0 + 1 + face;
    value = -val;
    HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
				   0, 1, &entry, &value);

    //   Update diagonal
    entry = 0;
    value = val;
    HYPRE_SStructMatrixAddToValues(A_, level_coarse, index_coarse, 
				   0, 1, &entry, &value);

  } // if phase == phase_matrix
}

//------------------------------------------------------------------------
