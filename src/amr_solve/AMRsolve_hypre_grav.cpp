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
#include "AMRsolve_HG_prec.h"
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
					 AMRsolve_Parameters& parameters,
					 int precflag, int zeroguess)
  : grid_(0), graph_(0), stencil_(0), A_(0), B_(0), X_(0), Y_(0),
    solver_(0), Ac_(0), Bc_(0), Xc_(0), cgrid_(0), cstencil_(0), citer_(-1), 
    parameters_(&parameters), hierarchy_(&hierarchy),
    resid_(-1.0), iter_(-1), r_factor_(const_r_factor), matrix_scale_(1.0)
{
  // set preconditioner and initial guess flags
  use_prec_   = (precflag != 0);
  zero_guess_ = (zeroguess != 0);

  // set array-valued items
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      BdryType_[i][j] = -1;
}

//----------------------------------------------------------------------

/// AMRsolve_Hypre_Grav destructor
AMRsolve_Hypre_Grav::~AMRsolve_Hypre_Grav()
{
  // destroy HYPRE objects that we created along the way
  if (HYPRE_SStructVectorDestroy(B_) != 0) 
    ERROR("AMRsolve_Hypre_FLD::could not destroy B_\n");
  if (HYPRE_SStructVectorDestroy(X_) != 0) 
    ERROR("AMRsolve_Hypre_FLD::could not destroy X_\n");
  if (HYPRE_SStructMatrixDestroy(A_) != 0) 
    ERROR("AMRsolve_Hypre_FLD::could not destroy A_\n");
  if (HYPRE_SStructGraphDestroy(graph_) != 0) 
    ERROR("AMRsolve_Hypre_FLD::could not destroy graph_\n");
  if (HYPRE_SStructStencilDestroy(stencil_) != 0) 
    ERROR("AMRsolve_Hypre_FLD::could not destroy stencil_\n");
  if (HYPRE_SStructGridDestroy(grid_) != 0) 
    ERROR("AMRsolve_Hypre_FLD::could not destroy grid_\n");

  if (use_prec_) {
    if (HYPRE_SStructVectorDestroy(Y_) != 0)  
      ERROR("could not destroy Y_\n");
    if (HYPRE_StructVectorDestroy(Bc_) != 0)
      ERROR("could not destroy Bc_\n");
    if (HYPRE_StructVectorDestroy(Xc_) != 0)
      ERROR("could not destroy Xc_\n");
    if (HYPRE_StructMatrixDestroy(Ac_) != 0)
      ERROR("could not destroy Ac_\n");
    if (HYPRE_StructStencilDestroy(cstencil_) != 0)
      ERROR("could not destroy cstencil_\n");
    if (HYPRE_StructGridDestroy(cgrid_) != 0)
      ERROR("could not destroy cgrid_\n");
  }
}

//----------------------------------------------------------------------

/// Initialize the Grid Hierarchy
/** Creates a hypre grid, with one part per level and one box per Grid
    patch object, for an AMR problem.  Sets grid box extents, grid
    part variables, and periodicity. */
void AMRsolve_Hypre_Grav::init_hierarchy()
{

  int ierr;
  int dim       = hierarchy_->dimension();
  int num_parts = hierarchy_->num_levels();

  // Create the hypre grid
  _TRACE_;
  ierr = HYPRE_SStructGridCreate(MPI_COMM_WORLD, dim, num_parts, &grid_);
  if (ierr != 0)  ERROR("could not create grid_\n");
  if (use_prec_) {
    ierr = HYPRE_StructGridCreate(MPI_COMM_WORLD, dim, &cgrid_);
    if (ierr != 0)  ERROR("could not create cgrid_\n");
  }

  _TRACE_;
  ItHierarchyLevels itl (*hierarchy_);
  while (AMRsolve_Level * level = itl++) {

    _TRACE_;
    int part = level->index();

    // Set extents for boxes that comprise the hypre grid
    ItLevelGridsLocal itgl (*level);
    while (AMRsolve_Grid * grid = itgl++) {
      int lower[3] = {grid->index_lower(0),
		      grid->index_lower(1),
		      grid->index_lower(2)};
      int upper[3] = {grid->index_upper(0),
		      grid->index_upper(1),
		      grid->index_upper(2)};
      ierr = HYPRE_SStructGridSetExtents(grid_, part, lower, upper);
      if (ierr != 0)  ERROR("could not set grid_ extents\n");
      
      // set up cgrid for coarse level
      if (use_prec_ && part==0) {
	ierr = HYPRE_StructGridSetExtents(cgrid_, lower, upper);
	if (ierr != 0)  ERROR("could not set cgrid_ extents\n");
      }
    } // while grid = itgl++

    _TRACE_;
    // Create a single cell-centered variable for each grid part (level)
    HYPRE_SStructVariable variable_types[] = { HYPRE_SSTRUCT_VARIABLE_CELL };
    const int numvars = 1;
    ierr = HYPRE_SStructGridSetVariables(grid_, part, numvars, variable_types);
    if (ierr != 0)  ERROR("could not set grid_ variables\n");

    // Set periodicity of the grid part
    int period[3] = { hierarchy_->period_index(0,part),
		      hierarchy_->period_index(1,part),
		      hierarchy_->period_index(2,part) };

    _TRACE_;
    ierr = HYPRE_SStructGridSetPeriodic(grid_, part, period);
    if (ierr != 0)  ERROR("could not set grid_ periodicity\n");

    if (use_prec_ && part==0) {
      ierr = HYPRE_StructGridSetPeriodic(cgrid_, period);
      if (ierr != 0)  ERROR("could not set cgrid_ periodicity\n");
    }

  } // while level = itl++

  // When finished, assemble the hypre grid
  _TRACE_;

  ierr = HYPRE_SStructGridAssemble(grid_);
  if (ierr != 0)  ERROR("could not assemble grid_\n");
  if (use_prec_) {
    ierr = HYPRE_StructGridAssemble(cgrid_);
    if (ierr != 0)  ERROR("could not assemble cgrid_\n");
  }
  _TRACE_;
} // AMRsolve_Hypre_Grav::init_hierarchy()

//----------------------------------------------------------------------

/// Initialize the discretization stencils.  
/** Creates and initializes a hypre stencil object. */
void AMRsolve_Hypre_Grav::init_stencil()
{
  int ierr;

  _TRACE_;
  int dim = hierarchy_->dimension();

  _TRACE_;
  ierr = HYPRE_SStructStencilCreate(dim,dim*2+1,&stencil_);
  if (ierr != 0)  ERROR("could not initialize stencil_\n");
  if (use_prec_) {
    ierr = HYPRE_StructStencilCreate(dim, dim*2+1, &cstencil_);
    if (ierr != 0)  ERROR("could not initialize stencil_\n");
  }

  int entries[][3] = { {  0, 0, 0 },     // center
		       {  1, 0, 0 },     // X+
		       { -1, 0, 0 },     // X-
		       {  0, 1, 0 },     // Y+
		       {  0,-1, 0 },     // Y-
		       {  0, 0, 1 },     // Z+
		       {  0, 0,-1 } };   // Z-

  _TRACE_;
  if (dim >= 1) {
    ierr = HYPRE_SStructStencilSetEntry(stencil_, 0, entries[0], 0);
    if (ierr != 0)  ERROR("could not set stencil_ entry\n");
    if (use_prec_) {
      ierr = HYPRE_StructStencilSetElement(cstencil_, 0, entries[0]);
      if (ierr != 0)  ERROR("could not set cstencil_ entry\n");
    }
  }
  if (dim >= 1) {
    ierr = HYPRE_SStructStencilSetEntry(stencil_, 1, entries[1], 0);
    if (ierr != 0)  ERROR("could not set stencil_ entry\n");
    if (use_prec_) {
      ierr = HYPRE_StructStencilSetElement(cstencil_, 1, entries[1]);
      if (ierr != 0)  ERROR("could not set cstencil_ entry\n");
    }
  }
  if (dim >= 1) {
    ierr = HYPRE_SStructStencilSetEntry(stencil_, 2, entries[2], 0);
    if (ierr != 0)  ERROR("could not set stencil_ entry\n");
    if (use_prec_) {
      ierr = HYPRE_StructStencilSetElement(cstencil_, 2, entries[2]);
      if (ierr != 0)  ERROR("could not set cstencil_ entry\n");
    }
  }
  if (dim >= 2) {
    ierr = HYPRE_SStructStencilSetEntry(stencil_, 3, entries[3], 0);
    if (ierr != 0)  ERROR("could not set stencil_ entry\n");
    if (use_prec_) {
      ierr = HYPRE_StructStencilSetElement(cstencil_, 3, entries[3]);
      if (ierr != 0)  ERROR("could not set cstencil_ entry\n");
    }
  }
  if (dim >= 2) {
    ierr = HYPRE_SStructStencilSetEntry(stencil_, 4, entries[4], 0);
    if (ierr != 0)  ERROR("could not set stencil_ entry\n");
    if (use_prec_) {
      ierr = HYPRE_StructStencilSetElement(cstencil_, 4, entries[4]);
      if (ierr != 0)  ERROR("could not set cstencil_ entry\n");
    }
  }
  if (dim >= 3) {
    ierr = HYPRE_SStructStencilSetEntry(stencil_, 5, entries[5], 0);
    if (ierr != 0)  ERROR("could not set stencil_ entry\n");
    if (use_prec_) {
      ierr = HYPRE_StructStencilSetElement(cstencil_, 5, entries[5]);
      if (ierr != 0)  ERROR("could not set cstencil_ entry\n");
    }
  }
  if (dim >= 3) {
    ierr = HYPRE_SStructStencilSetEntry(stencil_, 6, entries[6], 0);
    if (ierr != 0)  ERROR("could not set stencil_ entry\n");
    if (use_prec_) {
      ierr = HYPRE_StructStencilSetElement(cstencil_, 6, entries[6]);
      if (ierr != 0)  ERROR("could not set cstencil_ entry\n");
    }
  }
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
  int ierr;
  // Create the hypre graph object
  _TRACE_;
  ierr = HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid_, &graph_);
  if (ierr != 0)  ERROR("could not create graph_\n");
  
  _TRACE_;
  ierr = HYPRE_SStructGraphSetObjectType(graph_, HYPRE_SSTRUCT);
  if (ierr != 0)  ERROR("could not set graph_ type\n");

  _TRACE_;
  ItHierarchyLevels itl (*hierarchy_);
  while (AMRsolve_Level * level = itl++) {

    int part = level->index();

    // Define stencil connections within each level
    _TRACE_;
    ierr = HYPRE_SStructGraphSetStencil(graph_, part, 0, stencil_);
    if (ierr != 0)  ERROR("could not set stencil_ into graph_\n");

    // Define graph connections between levels
    _TRACE_;
    if (part > 0) {
      ItLevelGridsAll itag (*level);
      while (AMRsolve_Grid * grid = itag++) init_graph_nonstencil_(*grid);
    } // if part > 0

    // Initialize face counters for subsequent matrix inter-level entries
    _TRACE_;
    ItLevelGridsAll itag (*level);
    while (AMRsolve_Grid * grid = itag++) {
      int dim = hierarchy_->dimension();
      grid->init_counter(dim*2+1);
    } // while grid = itag++

  } // while level = itl++

  // Assemble the hypre graph
  _TRACE_;
  ierr = HYPRE_SStructGraphAssemble(graph_);
  if (ierr != 0)  ERROR("could not assemble graph_\n");
  _TRACE_;
} // AMRsolve_Hypre_Grav::init_graph()

//----------------------------------------------------------------------

/// Initialize the matrix A and right-hand-side vector b
/* Creates a matrix with a given nonzero structure, and sets nonzero
   values. */
/* void AMRsolve_Hypre_Grav::init_elements(std::vector<AMRsolve_Point *> points, 
                                           Scalar f_scale) */
void AMRsolve_Hypre_Grav::init_elements(int BdryType, Scalar f_scale)
{
  // set input arguments into AMRsolve_Hypre_Grav object
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      BdryType_[i][j] = BdryType;
 
  // Create the hypre matrix A_, solution X_, and right-hand side B_ objects
  int ierr;
  ierr = HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph_, &A_);
  if (ierr != 0)  ERROR("could not create A_\n");
  ierr = HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid_,  &X_);
  if (ierr != 0)  ERROR("could not create X_\n");
  ierr = HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid_,  &B_);
  if (ierr != 0)  ERROR("could not create B_\n");

  if (use_prec_) {
    ierr = HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid_,  &Y_);
    if (ierr != 0)  ERROR("could not create Y_\n");
    ierr = HYPRE_StructMatrixCreate(MPI_COMM_WORLD, cgrid_, cstencil_,  &Ac_);
    if (ierr != 0)  ERROR("could not create Ac_\n");
    ierr = HYPRE_StructVectorCreate(MPI_COMM_WORLD, cgrid_, &Xc_);
    if (ierr != 0)  ERROR("could not create Xc_\n");
    ierr = HYPRE_StructVectorCreate(MPI_COMM_WORLD, cgrid_, &Bc_);
    if (ierr != 0)  ERROR("could not create Bc_\n");
  }

  // Set the object types
  if (parameters_->value("solver") == "bicgstab-boomer") {
    ierr = HYPRE_SStructMatrixSetObjectType(A_,HYPRE_PARCSR);
    if (ierr != 0)  ERROR("could not set A_ type\n");
    ierr = HYPRE_SStructVectorSetObjectType(X_,HYPRE_PARCSR);
    if (ierr != 0)  ERROR("could not set X_ type\n");
    ierr = HYPRE_SStructVectorSetObjectType(B_,HYPRE_PARCSR);
    if (ierr != 0)  ERROR("could not set B_ type\n");
    if (use_prec_) {
      ierr = HYPRE_SStructVectorSetObjectType(Y_,HYPRE_PARCSR);
      if (ierr != 0)  ERROR("could not set Y_ type\n");
    }
  } else {
    ierr = HYPRE_SStructMatrixSetObjectType(A_,HYPRE_SSTRUCT);
    if (ierr != 0)  ERROR("could not set A_ type\n");
    ierr = HYPRE_SStructVectorSetObjectType(X_,HYPRE_SSTRUCT);
    if (ierr != 0)  ERROR("could not set X_ type\n");
    ierr = HYPRE_SStructVectorSetObjectType(B_,HYPRE_SSTRUCT);
    if (ierr != 0)  ERROR("could not set B_ type\n");
    if (use_prec_) {
      ierr = HYPRE_SStructVectorSetObjectType(Y_,HYPRE_SSTRUCT);
      if (ierr != 0)  ERROR("could not set Y_ type\n");
    }
  }

  // Initialize the hypre matrix and vector objects
  ierr = HYPRE_SStructMatrixInitialize(A_);
  if (ierr != 0)  ERROR("could not initialize A_\n");
  ierr = HYPRE_SStructVectorInitialize(X_);
  if (ierr != 0)  ERROR("could not initialize X_\n");
  ierr = HYPRE_SStructVectorInitialize(B_);
  if (ierr != 0)  ERROR("could not initialize B_\n");

  if (use_prec_) {
    ierr = HYPRE_SStructVectorInitialize(Y_);
    if (ierr != 0)  ERROR("could not initialize Y_\n");
    ierr = HYPRE_StructMatrixInitialize(Ac_);
    if (ierr != 0)  ERROR("could not initialize Ac_\n");
    ierr = HYPRE_StructVectorInitialize(Xc_);
    if (ierr != 0)  ERROR("could not initialize Xc_\n");
    ierr = HYPRE_StructVectorInitialize(Bc_);
    if (ierr != 0)  ERROR("could not initialize Bc_\n");
  }

  //--------------------------------------------------
  // Initialize the matrix A_ and Ac_ elements
  //--------------------------------------------------
  init_elements_matrix_();

  //--------------------------------------------------
  // Initialize B_ elements 
  //--------------------------------------------------
  init_elements_rhs_(f_scale);

  // Assemble the matrix and vectors
  ierr = HYPRE_SStructMatrixAssemble(A_);
  if (ierr != 0)  ERROR("could not assemble A_\n");
  ierr = HYPRE_SStructVectorAssemble(B_);
  if (ierr != 0)  ERROR("could not assemble B_\n");
  ierr = HYPRE_SStructVectorAssemble(X_);
  if (ierr != 0)  ERROR("could not assemble X_\n");

  if (use_prec_) {
    ierr = HYPRE_SStructVectorSetConstantValues(Y_, 0.0);
    if (ierr != 0)  ERROR("could not initialize Y_ to 0.0\n");
    ierr = HYPRE_SStructVectorAssemble(Y_);
    if (ierr != 0)  ERROR("could not assemble Y_\n");
    ierr = HYPRE_StructMatrixAssemble(Ac_);
    if (ierr != 0)  ERROR("could not assemble Ac_\n");
    ierr = HYPRE_StructVectorAssemble(Bc_);
    if (ierr != 0)  ERROR("could not assemble Ac_\n");
    ierr = HYPRE_StructVectorAssemble(Xc_);
    if (ierr != 0)  ERROR("could not assemble Ac_\n");
  }

  // Optionally write the matrix and vectors to a file for debugging
  if (parameters_->value("dump_a") == "true") {
    ierr = HYPRE_SStructMatrixPrint("A-hypre",A_,0);
    if (ierr != 0)  ERROR("could not print A_\n");
  }
  if (parameters_->value("dump_b") == "true") {
    ierr = HYPRE_SStructVectorPrint("B-hypre",B_,0);
    if (ierr != 0)  ERROR("could not print A_\n");
  }
  if (parameters_->value("dump_ac") == "true" && use_prec_) {
    ierr = HYPRE_StructMatrixPrint("Ac-hypre",Ac_,0);
    if (ierr != 0)  ERROR("could not print Ac_\n");
  }

} // AMRsolve_Hypre_Grav::init_elements()

//----------------------------------------------------------------------

/// Initialize and solve the linear solver
void AMRsolve_Hypre_Grav::solve()
{
  int ierr;
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

  // recheck solver parameters
  sitmax  = parameters_->value("solver_itmax");
  srestol = parameters_->value("solver_restol");

  // Set local variables
  itmax  = atoi(sitmax.c_str());
  restol = atof(srestol.c_str());

  // call solver
  if (solver == "pfmg" && levels == 1) {
    solve_pfmg_(itmax,restol);

  } else if (solver == "fac"  && levels == 1) {
    solve_pfmg_(itmax,restol);     // if only one level, FAC defaults to standard MG

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
    sprintf(error_message, "AMRsolve_Hypre_Grav::solve called with illegal combination of solver %s on %d levels", solver.c_str(),levels);
    ERROR(error_message);
  }
  
  if (parameters_->value("dump_x") == "true") {
    ierr = HYPRE_SStructVectorPrint("X-hypre",X_,1);
    if (ierr != 0)  ERROR("could not print X_\n");
  }

} // AMRsolve_Hypre_Grav::solve()

//----------------------------------------------------------------------

/// Evaluate the success of the solve
int AMRsolve_Hypre_Grav::evaluate()
{
//   // check whether solution/rhs was requested
//   if (parameters_->value("dump_x") == "true" || 
//       parameters_->value("dump_b") == "true") {

//     // iterate over processor-local grids
//     ItHierarchyGridsLocal itg(*hierarchy_);
//     while (AMRsolve_Grid * grid = itg++) {
      
//       int ierr;
//       char filename[80];

//       // get level & grid information
//       int level = grid->level();
//       int lower[3],upper[3];
//       grid->get_limits(lower,upper);

//       if (parameters_->value("dump_x") == "true") {
// 	// extract Enzo solution
// 	int nx[3];
// 	ierr = HYPRE_SStructVectorGather(X_);
// 	if (ierr != 0)  ERROR("could not gather X_\n");
// 	ierr = HYPRE_SStructVectorGetBoxValues(X_, level, lower, upper, 0,
// 					       grid->get_u(&nx[0],&nx[1],&nx[2]));  
// 	if (ierr != 0)  ERROR("could not GetBoxValues from X_\n");
// 	sprintf(filename,"X.%d",grid->id());
// 	grid->write("header",filename);
// 	grid->write("u",filename);
//       }
    
//       if (parameters_->value("dump_b") == "true") {
// 	// extract Enzo rhs
// 	int nb[3];
// 	ierr = HYPRE_SStructVectorGather(B_);
// 	if (ierr != 0)  ERROR("could not gather B_\n");
// 	ierr = HYPRE_SStructVectorGetBoxValues(B_, level, lower, upper, 0,
// 					       grid->get_f(&nb[0],&nb[1],&nb[2]));  
// 	if (ierr != 0)  ERROR("could not GetBoxValues from B_\n");
// 	sprintf(filename,"B.%d",grid->id());
// 	grid->write("f",filename);
//       }

//     } // grid = itg++
//   } // if dump_x or dump_b

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

  return err_flag;

} // AMRsolve_Hypre_Grav::evaluate()


//----------------------------------------------------------------------

/// Extracts HYPRE solution and updates potential field
void AMRsolve_Hypre_Grav::update_enzo()
{
  int ierr;

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
    ierr = HYPRE_SStructVectorGather(X_);
    if (ierr != 0)  ERROR("could not gather X_\n");
    ierr = HYPRE_SStructVectorGetBoxValues(X_, level, lower, upper, 0, u);  
    if (ierr != 0)  ERROR("could not GetBoxValues from X_\n");

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
	  phi[k] = u[i];
	}
      }
    }

  } // grid = itg++

} // AMRsolve_Hypre_Grav::update_enzo()


//----------------------------------------------------------------------

/// Writes Enzo potential field to disk
void AMRsolve_Hypre_Grav::write_potential()
{
  // iterate over grids on this processor
  ItHierarchyGridsLocal itgl(*hierarchy_);
  while (AMRsolve_Grid* grid = itgl++) {

    // get level, grid information
    int level = grid->level();
    int lower[3],upper[3];
    grid->get_limits(lower,upper);

    // access AMRsolve_hypre_grav and Enzo solutions
    int n3[3];
    Scalar* u = grid->get_u(&n3[0],&n3[1],&n3[2]);
    Scalar* phi = grid->get_phi();

    // get buffering information on relating amrsolve grid to Enzo data
    int ghosts[3][2]; 
    grid->get_GravGhosts(ghosts);
    int en0 = n3[0] + ghosts[0][0] + ghosts[0][1];
    int en1 = n3[1] + ghosts[1][0] + ghosts[1][1];
    int en2 = n3[2] + ghosts[2][0] + ghosts[2][1];

    // open output file
    FILE *fptr=NULL;
    char filename[80];
    sprintf(filename, "phi_g%d.txt", grid->id());
    fptr = fopen(filename,"w");

    // output Enzo solution
    int k0, k1, k2, i0, i1, i2;
    for (k2=0; k2<n3[2]+ghosts[2][0]+ghosts[2][1]; k2++) 
      for (k1=0; k1<n3[1]+ghosts[1][0]+ghosts[1][1]; k1++) 
	for (k0=0; k0<n3[0]+ghosts[0][0]+ghosts[0][1]; k0++) 
	  fprintf(fptr,"%i   %i   %i   %.16g\n",
		  k0-ghosts[0][0],
		  k1-ghosts[1][1],
		  k2-ghosts[2][1],
		  phi[k0 + en0*(k1 + en1*k2)]);

    // close output file
    fclose(fptr);

  } // grid = itgl++

} // AMRsolve_Hypre_Grav::write_potential()


//----------------------------------------------------------------------

/// dumps HYPRE matrix and RHS (called when aborting solve)
void AMRsolve_Hypre_Grav::abort_dump()
{
  int ierr;

  // have HYPRE dump out everything it knows to disk
  ierr = HYPRE_SStructMatrixPrint("A.mat",A_,0);
  if (ierr != 0)  ERROR("could not print A_\n");
  ierr = HYPRE_SStructVectorPrint("x.vec",X_,0);
  if (ierr != 0)  ERROR("could not print X_\n");
  ierr = HYPRE_SStructVectorPrint("b.vec",B_,0);
  if (ierr != 0)  ERROR("could not print B_\n");
  if (use_prec_) {
    ierr = HYPRE_StructMatrixPrint("Ac.mat",Ac_,0);
    if (ierr != 0)  ERROR("could not print Ac_\n");
  }

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
    sprintf(error_message,"init_nonstencil_ called with phase = %d", int(phase));
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
		    grid.id(), axis1, index_global[axis1][0], r_factor_);
    assert(l0);
    if (!l1) printf("grid %i,  index_global[%d][1] = %d,  r_factor = %i\n",
		    grid.id(), axis1, index_global[axis1][1], r_factor_);
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
  int i, ierr;

  _TRACE_;
  // iterate over all grids local to this processor
  ItHierarchyGridsLocal itg (*hierarchy_);
  while (AMRsolve_Grid* grid = itg++) {

    // for each local grid, use our "f" array to store vector entries 
    // for the linear system.
    int n0,n1,n2;
    Scalar* values = grid->get_f(&n0,&n1,&n2);
    for (i=0; i<n0*n1*n2; i++)  values[i] = 0.0;

    // access Enzo's GravitatingMassField
    Scalar* GravMassField = grid->get_gmass();

    // get buffering information on relating amrsolve grid to Enzo data
    int ghosts[3][2]; 
    grid->get_GravGhosts(ghosts);
    int en0 = n0 + ghosts[0][0] + ghosts[0][1];  // enzo gravity data dimensions
    int en1 = n1 + ghosts[1][0] + ghosts[1][1];  // enzo gravity data dimensions
    int en2 = n2 + ghosts[2][0] + ghosts[2][1];  // enzo gravity data dimensions

    // iterate over the domain, copying Enzo's GravitatingMassField into f
    int k0, k1, k2, i0, i1, i2;
    for (i2=0; i2<n2; i2++) {
      k2 = ghosts[2][0] + i2;
      for (i1=0; i1<n1; i1++) {
	k1 = ghosts[1][0] + i1;
	for (i0=0; i0<n0; i0++) {
	  k0 = ghosts[0][0] + i0;
	  values[i0 + n0*(i1 + n1*i2)] = GravMassField[k0 + en0*(k1 + en1*k2)];
	}
      }
    }
   
    // scale RHS by input scaling multiplied by cell volume
    // [input contains GravitationalConstant/a]
    Scalar hx,hy,hz;
    grid->h(hx,hy,hz);
    Scalar scale = f_scale * hx * hy * hz;
    for (int i=0; i<n0*n1*n2; i++)  values[i] *= scale;

    // Set Hypre B_ vector to grid f_ values
    int part = grid->level();
    int lower[3] = { grid->index_lower(0), 
		     grid->index_lower(1), 
		     grid->index_lower(2) };
    int upper[3] = { grid->index_upper(0), 
		     grid->index_upper(1), 
		     grid->index_upper(2) };
    ierr = HYPRE_SStructVectorAddToBoxValues(B_,part,lower,upper,0,values);
    if (ierr != 0)  ERROR("could not AddToBoxValues in B_\n");

    // insert current values for PotentialField into HYPRE vector X_
    if (!zero_guess_) {
      Scalar* phi = grid->get_phi();
      for (i0=0; i0<n0*n1*n2; i0++)  values[i0] = 0.0;
      for (i2=0; i2<n2; i2++) {
	k2 = ghosts[2][0] + i2;
	for (i1=0; i1<n1; i1++) {
	  k1 = ghosts[1][0] + i1;
	  for (i0=0; i0<n0; i0++) {
	    k0 = ghosts[0][0] + i0;
	    values[i0 + n0*(i1 + n1*i2)] = phi[k0 + en0*(k1 + en1*k2)];
	  }
	}
      }
      ierr = HYPRE_SStructVectorSetBoxValues(X_,part,lower,upper,0,values);
      if (ierr != 0)  ERROR("could not SetBoxValues in X_\n");
    }

  } //  while grid = itg++
  
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
    if (ierr != 0)  ERROR("could not ZeroAMRVectorData in B_\n");
    if (!zero_guess_) {
      ierr = HYPRE_SStructFACZeroAMRVectorData(X_, plevels, rfactors);
      if (ierr != 0)  ERROR("could not ZeroAMRVectorData in X_\n");
    }
  }


  // for periodic BCs, need to shift RHS to have zero average
  // value to deflate the null space
  if ( hierarchy_->is_periodic(0) && 
       hierarchy_->is_periodic(1) && 
       hierarchy_->is_periodic(2) ) {

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
 	  ierr = HYPRE_SStructVectorGetBoxValues(B_, part, lower, upper, 0, tmpvals);
	  if (ierr != 0)  ERROR("could not GetBoxValues in B_\n");

	  // Accumulate the sum
	  for (int i=0; i<n; i++)  local_shift_b_sum += tmpvals[i];

	  // Delete the patch
 	  delete [] tmpvals;
	}
      }
    }

    // Accumulate global sum from local sums
    Scalar shift_b_sum = 0.0;
    ierr = MPI_Allreduce(&local_shift_b_sum, &shift_b_sum, 1, 
			 MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD);
    if (ierr != 0)  ERROR("error in MPI_Allreduce\n");

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

	ierr = HYPRE_SStructVectorAddToBoxValues(B_,part,lower,upper,0,tmpvals);
	if (ierr != 0)  ERROR("could not AddToBoxValues in B_\n");

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
      ierr = HYPRE_SStructFACZeroAMRVectorData(B_, plevels, rfactors);
      if (ierr != 0)  ERROR("could not ZeroAMRVectorData in B_\n");
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
  int ierr;

  double h120 = h3[1]*h3[2] / h3[0];
  double h201 = h3[2]*h3[0] / h3[1];
  double h012 = h3[0]*h3[1] / h3[2];

  double* v0;         // Diagonal elements
  double* v1[3][2];   // Off-diagonal elements

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

  //-----------------------------------------------------------
  // Adjust stencil at grid boundaries:
  //-----------------------------------------------------------

  // If purely periodic, we must deflate the null space by setting
  // a single Dirichlet condition in a corner of the root grid.
  if (hierarchy_->is_periodic(0) && 
      hierarchy_->is_periodic(1) && 
      hierarchy_->is_periodic(2) && coarsegrid &&
      (grid.index_lower(0)+grid.index_lower(1)+grid.index_lower(2) == 0)) {

    for (i0=0; i0<3; i0++)
      for (i1=0; i1<2; i1++)
	v1[i0][i1][0] = 0.0;

  // Otherwise, Enzo is using isolated boundary conditions
  } else {

    //    z-left face
    if (grid.faces().label(2,0,0,0) == AMRsolve_Faces::_boundary_) {
      i2 = 0;
      for (i1=0; i1<n3[1]; i1++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[2][0][i];
	  v1[2][0][i] = 0.0;
	} 
      }
    }  // label == boundary

    //    y-left face
    if (grid.faces().label(1,0,0,0) == AMRsolve_Faces::_boundary_) {
      i1 = 0;
      for (i2=0; i2<n3[2]; i2++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[1][0][i];
	  v1[1][0][i] = 0.0;
	}
      }
    }  // label == boundary

    //    x-left face
    if (grid.faces().label(0,0,0,0) == AMRsolve_Faces::_boundary_) {
      i0 = 0;
      for (i2=0; i2<n3[2]; i2++) {
	for (i1=0; i1<n3[1]; i1++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[0][0][i];
	  v1[0][0][i] = 0.0;
	}
      }
    }  // label == boundary

    //    x-right face
    if (grid.faces().label(0,1,0,0) == AMRsolve_Faces::_boundary_) {
      i0 = n3[0]-1;
      for (i2=0; i2<n3[2]; i2++) {
	for (i1=0; i1<n3[1]; i1++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[0][1][i];
	  v1[0][1][i] = 0.0;
	}
      }
    }  // label == boundary

    //    y-right face
    if (grid.faces().label(1,1,0,0) == AMRsolve_Faces::_boundary_) {
      i1 = n3[1]-1;
      for (i2=0; i2<n3[2]; i2++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[1][1][i];
	  v1[1][1][i] = 0.0;
	}
      }
    }  // label == boundary

    //    z-right face
    if (grid.faces().label(2,1,0,0) == AMRsolve_Faces::_boundary_) {
      i2 = n3[2]-1;
      for (i1=0; i1<n3[1]; i1++) {
	for (i0=0; i0<n3[0]; i0++) {
	  i = AMRsolve_Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v0[i] += v1[2][1][i];
	  v1[2][1][i] = 0.0;
	}
      }
    }  // label == boundary
  }  // isolated boundary conditions


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
  if (ierr != 0)  ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[1], v1[0][1]);
  if (ierr != 0)  ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[2], v1[0][0]);
  if (ierr != 0)  ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[3], v1[1][1]);
  if (ierr != 0)  ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[4], v1[1][0]);
  if (ierr != 0)  ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[5], v1[2][1]);
  if (ierr != 0)  ERROR("could not SetBoxValues in A_\n");
  ierr = HYPRE_SStructMatrixSetBoxValues(A_, level, index_lower, index_upper, 
					 0, 1, &entries[6], v1[2][0]);
  if (ierr != 0)  ERROR("could not SetBoxValues in A_\n");

  if (coarsegrid && use_prec_) {
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 
					  1, &entries[0], v0);
    if (ierr != 0)  ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 
					  1, &entries[1], v1[0][1]);
    if (ierr != 0)  ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 
					  1, &entries[2], v1[0][0]);
    if (ierr != 0)  ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 
					  1, &entries[3], v1[1][1]);
    if (ierr != 0)  ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 
					  1, &entries[4], v1[1][0]);
    if (ierr != 0)  ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 
					  1, &entries[5], v1[2][1]);
    if (ierr != 0)  ERROR("could not SetBoxValues in Ac_\n");
    ierr = HYPRE_StructMatrixSetBoxValues(Ac_, index_lower, index_upper, 
					  1, &entries[6], v1[2][0]);
    if (ierr != 0)  ERROR("could not SetBoxValues in Ac_\n");
  }

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
    int ierr = HYPRE_SStructFACZeroAMRMatrixData(A_, part-1, r_factors);
    if (ierr != 0)  ERROR("could not ZeroAMRMatrixData in A_\n");
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
  int ierr;
  _TRACE_;
  // Create and initialize the solver
  ierr = HYPRE_SStructSysPFMGCreate(MPI_COMM_WORLD, &solver_);
  if (ierr != 0)  ERROR("could not create SysPFMG solver_\n");

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
  if (itmax != 0 ) {
    ierr = HYPRE_SStructSysPFMGSetMaxIter(solver_,itmax);
    if (ierr != 0)  ERROR("could not set itmax\n");
  }
  if (restol != 0.0) {
    ierr = HYPRE_SStructSysPFMGSetTol(solver_,restol);
    if (ierr != 0)  ERROR("could not set restol\n");
  }
  ierr = HYPRE_SStructSysPFMGSetRelaxType(solver_,rlxtype);
  if (ierr != 0)  ERROR("could not set rlxtype\n");
  ierr = HYPRE_SStructSysPFMGSetNumPreRelax(solver_,npre);
  if (ierr != 0)  ERROR("could not set npre\n");
  ierr = HYPRE_SStructSysPFMGSetNumPostRelax(solver_,npost);
  if (ierr != 0)  ERROR("could not set npost\n");
  ierr = HYPRE_SStructSysPFMGSetPrintLevel(solver_,printl);
  if (ierr != 0)  ERROR("could not set printl\n");
  ierr = HYPRE_SStructSysPFMGSetLogging(solver_,log);
  if (ierr != 0)  ERROR("could not set log\n");

  // setup solver 
  ierr = HYPRE_SStructSysPFMGSetup(solver_,A_,B_,X_);
  if (ierr != 0)  ERROR("could not setup solver_\n");

  // Solve the linear system
  ierr = HYPRE_SStructSysPFMGSolve(solver_,A_,B_,X_);
  if (ierr != 0) {
    fprintf(stderr, "Error %i detected, printing current hierarchy:\n",ierr);
    hierarchy_->print();
    ERROR("could not solve with PFMG\n");
  }

  // Write out some diagnostic info about the solve
  ierr = HYPRE_SStructSysPFMGGetNumIterations(solver_,&iter_);
  if (ierr != 0)  ERROR("could not get iter_\n");
  ierr = HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm(solver_,&resid_);
  if (ierr != 0)  ERROR("could not get resid_\n");
  if (debug && pmpi->is_root()) {
    printf("hypre PFMG num iterations: %d\n",iter_);
    printf("hypre PFMG final relative residual norm: %g\n",resid_);
  }

  // Delete the solver
  ierr = HYPRE_SStructSysPFMGDestroy(solver_);
  if (ierr != 0)  ERROR("could not destroy solver_\n");
  _TRACE_;
} // AMRsolve_Hypre_Grav::solve_pfmg_()

//------------------------------------------------------------------------

/// Initialize the FAC hypre solver
void AMRsolve_Hypre_Grav::solve_fac_(int itmax, double restol)
{
  _TRACE_;
  int i, ierr;

  // Create the solver
  ierr = HYPRE_SStructFACCreate(MPI_COMM_WORLD, &solver_);
  if (ierr != 0)  ERROR("could not create FAC solver_\n");

  // Initialize parts
  int num_parts = hierarchy_->num_levels();
  ierr = HYPRE_SStructFACSetMaxLevels(solver_, num_parts);
  if (ierr != 0)  ERROR("could not set max FAC levels\n");

  int *parts = new int[num_parts];
  for (i=0; i<num_parts; i++) parts[i] = i;
  ierr = HYPRE_SStructFACSetPLevels(solver_, num_parts, parts);
  if (ierr != 0)  ERROR("could not set FAC PLevels\n");

  // Initialize refinement factors
  int3 *refinements = new int3[num_parts];
  for (i=0; i<num_parts; i++) {
    refinements[i][0] = r_factor_;
    refinements[i][1] = r_factor_;
    refinements[i][2] = r_factor_;
  }
  ierr = HYPRE_SStructFACSetPRefinements(solver_, num_parts, refinements);
  if (ierr != 0)  ERROR("could not set FAC refinements\n");

  // extract some additional solver parameters
  std::string srlxtype = parameters_->value("solver_rlxtype");
  std::string snpre    = parameters_->value("solver_npre");
  std::string snpost   = parameters_->value("solver_npost");
  std::string scsolve  = parameters_->value("solver_csolve");
  std::string sprintl  = parameters_->value("solver_printl");
  std::string slog     = parameters_->value("solver_log");

  //   if not defined, then define them
  if (srlxtype == "")  srlxtype = "2";
  if (snpre    == "")  snpre    = "2";
  if (snpost   == "")  snpost   = "2";
  if (scsolve  == "")  scsolve  = "1";
  if (sprintl  == "")  sprintl  = "1";
  if (slog     == "")  slog     = "1";

  //   set local variables
  int relax   = atoi(srlxtype.c_str());
  int npre    = atoi(snpre.c_str());
  int npost   = atoi(snpost.c_str());
  int printl  = atoi(sprintl.c_str());
  int log     = atoi(slog.c_str());
  int csolve  = atoi(scsolve.c_str());

  // solver parameters
  ierr = HYPRE_SStructFACSetNumPreRelax(solver_,npre);
  if (ierr != 0)  ERROR("could not set npre\n");
  ierr = HYPRE_SStructFACSetNumPostRelax(solver_,npost);
  if (ierr != 0)  ERROR("could not set npost\n");
  ierr = HYPRE_SStructFACSetCoarseSolverType(solver_,csolve);
  if (ierr != 0)  ERROR("could not set csolve\n");
  ierr = HYPRE_SStructFACSetRelaxType(solver_,relax);
  if (ierr != 0)  ERROR("could not set relax\n");

  // stopping criteria
  if (itmax != 0 ) {
    ierr = HYPRE_SStructFACSetMaxIter(solver_,itmax);
    if (ierr != 0)  ERROR("could not set itmax\n");
  }
  if (restol != 0.0) {
    ierr = HYPRE_SStructFACSetTol(solver_,restol);
    if (ierr != 0)  ERROR("could not set restol\n");
  }

  // output amount
  ierr = HYPRE_SStructFACSetLogging(solver_, 1);
  if (ierr != 0)  ERROR("could not set logging\n");

  // prepare for solve
  ierr = HYPRE_SStructFACSetup2(solver_, A_, B_, X_);
  if (ierr != 0)  ERROR("could not setup FAC solver_\n");

  // Solve the linear system
  ierr = HYPRE_SStructFACSolve3(solver_, A_, B_, X_);
  if (ierr != 0) {
    fprintf(stderr, "Error %i detected, printing current hierarchy:\n",ierr);
    hierarchy_->print();
    ERROR("could not solve with FAC\n");
  }

  // Write out some diagnostic info about the solve
  ierr = HYPRE_SStructFACGetNumIterations(solver_, &iter_);
  if (ierr != 0)  ERROR("could not get iter_\n");
  ierr = HYPRE_SStructFACGetFinalRelativeResidualNorm(solver_, &resid_);
  if (ierr != 0)  ERROR("could not get resid_\n");
  if (debug && pmpi->is_root()) {
    printf("hypre FAC num iterations: %d\n",iter_);
    printf("hypre FAC final relative residual norm: %g\n",resid_);
  }

  // Delete the solver
  ierr = HYPRE_SStructFACDestroy2(solver_);
  if (ierr != 0)  ERROR("could not destroy solver_\n");

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
  int ierr;

  // Create the solver
  ierr = HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver_);
  if (ierr != 0)  ERROR("could not create solver_\n");

  // extract some additional solver parameters
  std::string slog = parameters_->value("solver_log");
  if (slog == "")      slog = "1";
  int log = atoi(slog.c_str());

  // stopping criteria
  if (itmax != 0 ) {
    ierr = HYPRE_SStructBiCGSTABSetMaxIter(solver_,itmax);
    if (ierr != 0)  ERROR("could not set itmax\n");
  }
  if (restol != 0.0) {
    ierr = HYPRE_SStructBiCGSTABSetTol(solver_,restol);
    if (ierr != 0)  ERROR("could not set restol\n");
  }

  // output amount
  ierr = HYPRE_SStructBiCGSTABSetLogging(solver_, 1);
  if (ierr != 0)  ERROR("could not set logging\n");

  // set up the preconditioner (if requested)
  AMRsolve_HG_prec *precond;
  if (use_prec_) {

    // set default solver parameters for HG preconditioner (if unset)
    if (parameters_->value("prec_itmax") == "")
      parameters_->add_parameter("prec_itmax","1"); 
    if (parameters_->value("prec_restol") == "")  
      parameters_->add_parameter("prec_restol","0.0");
    if (parameters_->value("prec_rlxtype") == "")  
      parameters_->add_parameter("prec_rlxtype","2");
    if (parameters_->value("prec_npre") == "")  
      parameters_->add_parameter("prec_npre","2"); 
    if (parameters_->value("prec_npost") == "")  
      parameters_->add_parameter("prec_npost","2");
    if (parameters_->value("prec_printl") == "")  
      parameters_->add_parameter("prec_printl","1");
    if (parameters_->value("prec_log") == "")  
      parameters_->add_parameter("prec_log","1");
    if (parameters_->value("prec_Jaciters") == "")  
      parameters_->add_parameter("prec_Jaciters","3");

    precond = new AMRsolve_HG_prec(*hierarchy_, BdryType_);
    ierr = precond->Initialize_(parameters_, &Ac_, &Xc_, &Bc_, &Y_);
    if (ierr != 0)  ERROR("could not initialize HG preconditioner\n");
    ierr = HYPRE_SStructBiCGSTABSetPrecond(solver_,
					   (HYPRE_PtrToSStructSolverFcn) HG_prec_solve,
					   (HYPRE_PtrToSStructSolverFcn) HG_prec_setup,
					   (HYPRE_SStructSolver) precond);
    if (ierr != 0)  ERROR("could not set preconditioner\n");
  }

  // Initialize the solver
  ierr = HYPRE_SStructBiCGSTABSetup(solver_, A_, B_, X_);
  if (ierr != 0)  ERROR("could not setup solver_\n");

  // Solve the linear system
  ierr = HYPRE_SStructBiCGSTABSolve(solver_, A_, B_, X_);
  if (ierr != 0) {
    HYPRE_SStructBiCGSTABGetNumIterations(solver_, &iter_);
    HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver_, &resid_);
    fprintf(stderr, "Error %i detected, BiCGSTAB solver stats:\n  iterations = %d\n  final relative residual = %g\n",ierr,iter_,resid_);
    fprintf(stderr, "Error %i detected, writing B_ to 'B-hypre.*':\n",ierr);
    HYPRE_SStructVectorPrint("B-hypre",B_,0);
    fprintf(stderr, "Error %i detected, writing X_ to 'X-hypre.*':\n",ierr);
    HYPRE_SStructVectorPrint("X-hypre",X_,0);
    fprintf(stderr, "Error %i detected, writing A_ to 'A-hypre.*':\n",ierr);
    HYPRE_SStructMatrixPrint("A-hypre",A_,0);
    if (use_prec_) {
      fprintf(stderr, "Error %i detected, writing Bc_ to 'Bc-hypre.*':\n",ierr);
      HYPRE_StructVectorPrint("Bc-hypre",Bc_,0);
      fprintf(stderr, "Error %i detected, writing Xc_ to 'Xc-hypre.*':\n",ierr);
      HYPRE_StructVectorPrint("Xc-hypre",Xc_,0);
      fprintf(stderr, "Error %i detected, writing Ac_ to 'Ac-hypre.*':\n",ierr);
      HYPRE_StructMatrixPrint("Ac-hypre",Ac_,0);
    }
    fprintf(stderr, "Error %i detected, printing current hierarchy:\n",ierr);
    hierarchy_->print();
    ERROR("could not solve with BiCGStab\n");
  }

  // Write out some diagnostic info about the solve
  ierr = HYPRE_SStructBiCGSTABGetNumIterations(solver_, &iter_);
  if (ierr != 0)  ERROR("could not get iter_\n");
  ierr = HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver_, &resid_);
  if (ierr != 0)  ERROR("could not get resid_\n");
  if (use_prec_) {
    Scalar presid = precond->GetResid_();
    int piters = precond->GetIters_();
    if (debug && pmpi->is_root()) {
      printf("hypre BiCGSTAB num iterations: %d,  prec: %d\n",iter_,piters);
      printf("hypre BiCGSTAB final relative residual norm: %g,  prec: %g\n",resid_,presid);
    }
  } else {
    if (debug && pmpi->is_root()) {
      printf("hypre BiCGSTAB num iterations: %d\n",iter_);
      printf("hypre BiCGSTAB final relative residual norm: %g\n",resid_);
    }
  }


  // project solutions from finest through coarsest grids in overlapped regions
  // (currently requires that preconditioner has been set up)
  if (use_prec_) {
    ierr = precond->HYPRE_to_AMRsolve_(&X_, 1); // copy X_ to AMRsolve's u_
    if (ierr != 0)  ERROR("could not copy X_ to u_\n");
    ierr = precond->restrict_(hierarchy_->num_levels()-1, 0);  // restriction
    if (ierr != 0)  ERROR("could not restrict u_ to coarse levels\n");
    ierr = precond->AMRsolve_to_HYPRE_(&X_, 1); // copy u_ to X_
    if (ierr != 0)  ERROR("could not copy u_ to X_\n");
  }


  // Delete the solver
  ierr = HYPRE_SStructBiCGSTABDestroy(solver_);
  if (ierr != 0)  ERROR("could not destroy solver_\n");
  if (use_prec_)  delete precond;

  _TRACE_;
} // AMRsolve_Hypre_Grav::solve_bicgstab_()

//------------------------------------------------------------------------

/// Initialize the BiCGSTAB hypre solver with BoomerAMG preconditioning
void AMRsolve_Hypre_Grav::solve_bicgstab_boomer_(int itmax, double restol)
{
  _TRACE_;
  int ierr;
  HYPRE_Solver solver;

  // Create the solver
  ierr = HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &solver);
  if (ierr != 0)  ERROR("could not create solver\n");

  // stopping criteria
  if (itmax != 0 ) {
    ierr = HYPRE_BiCGSTABSetMaxIter(solver,itmax);
    if (ierr != 0)  ERROR("could not set itmax\n");
  }
  if (restol != 0.0) {
    ierr = HYPRE_BiCGSTABSetTol(solver,restol);
    if (ierr != 0)  ERROR("could not set restol\n");
  }

  // extract some additional solver parameters
  std::string slog    = parameters_->value("solver_log");
  std::string sprintl = parameters_->value("solver_printl");
  if (slog    == "")   slog    = "1";
  if (sprintl == "")   sprintl = "1";
  int log    = atoi(slog.c_str());
  int printl = atoi(sprintl.c_str());

  // output amount
  ierr = HYPRE_BiCGSTABSetLogging(solver, log);
  if (ierr != 0)  ERROR("could not set logging\n");

  // Set BoomerAMG preconditioner
  HYPRE_Solver par_precond;
  ierr = HYPRE_BoomerAMGCreate(&par_precond);
  if (ierr != 0)  ERROR("could not create AMG precond\n");
  ierr = HYPRE_BoomerAMGSetCoarsenType(par_precond, 6);
  if (ierr != 0)  ERROR("could not set coarsening type\n");
  ierr = HYPRE_BoomerAMGSetStrongThreshold(par_precond, 0.25);
  if (ierr != 0)  ERROR("could not set threshold\n");
  ierr = HYPRE_BoomerAMGSetTol(par_precond, 0.0);
  if (ierr != 0)  ERROR("could not set tolerance\n");
  ierr = HYPRE_BoomerAMGSetPrintLevel(par_precond, 1);
  if (ierr != 0)  ERROR("could not set print level\n");
  ierr = HYPRE_BoomerAMGSetPrintFileName(par_precond, "sstruct.out.log");
  if (ierr != 0)  ERROR("could not SetPrintFileName\n");
  ierr = HYPRE_BoomerAMGSetMaxIter(par_precond, 1);
  if (ierr != 0)  ERROR("could not set maxiter\n");
  ierr = HYPRE_BiCGSTABSetPrecond(solver,
				  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
				  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
				  par_precond);
  if (ierr != 0)  ERROR("could not set preconditioner\n");

  // Initialize the solver
  ierr = HYPRE_BiCGSTABSetup(solver, (HYPRE_Matrix) A_, 
			     (HYPRE_Vector) B_, (HYPRE_Vector) X_);
  if (ierr != 0)  ERROR("could not setup solver\n");

  // Solve the linear system
  ierr = HYPRE_BiCGSTABSolve(solver, (HYPRE_Matrix) A_, 
			     (HYPRE_Vector)B_, (HYPRE_Vector) X_);
  if (ierr != 0) {
    fprintf(stderr, "Error %i detected, printing current hierarchy:\n",ierr);
    hierarchy_->print();
    ERROR("could not solve with BiCGStab\n");
  }

  // Write out some diagnostic info about the solve
  ierr = HYPRE_BiCGSTABGetNumIterations(solver, &iter_);
  if (ierr != 0)  ERROR("could not get iter_\n");
  ierr = HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver, &resid_);
  if (ierr != 0)  ERROR("could not get resid_\n");
  if (debug &&pmpi->is_root()) {
    printf("hypre BiCGSTAB num iterations: %d\n",iter_);
    printf("hypre BiCGSTAB final relative residual norm: %g\n",resid_);
  }


  // Delete the solver and preconditiner
  ierr = HYPRE_BoomerAMGDestroy(par_precond);
  if (ierr != 0)  ERROR("could not destroy precond\n");
  ierr = HYPRE_BiCGSTABDestroy(solver);
  if (ierr != 0)  ERROR("could not destroy solver\n");
  _TRACE_;
} // AMRsolve_Hypre_Grav::solve_bicgstab_boomer_()

//------------------------------------------------------------------------

/// Initialize the GMRES hypre solver
void AMRsolve_Hypre_Grav::solve_gmres_(int itmax, double restol)
{
  _TRACE_;
  int ierr;

  // Create the solver
  ierr = HYPRE_SStructGMRESCreate(MPI_COMM_WORLD, &solver_);
  if (ierr != 0)  ERROR("could not create solver_\n");

  // extract some additional solver parameters
  std::string slog = parameters_->value("solver_log");
  if (slog == "")      slog = "1";
  int log = atoi(slog.c_str());

  // stopping criteria
  if (itmax != 0 ) {
    ierr = HYPRE_SStructGMRESSetMaxIter(solver_, itmax);
    if (ierr != 0)  ERROR("could not set itmax\n");
  }
  if (restol != 0.0) {
    ierr = HYPRE_SStructGMRESSetTol(solver_, restol);
    if (ierr != 0)  ERROR("could not set restol\n");
  }

  // output amount
  ierr = HYPRE_SStructGMRESSetLogging(solver_, log);
  if (ierr != 0)  ERROR("could not set logging\n");

  // set up the preconditioner (if requested)
  AMRsolve_HG_prec *precond;
  if (use_prec_) {

    // set default solver parameters for HG preconditioner (if unset)
    if (parameters_->value("prec_itmax") == "")
      parameters_->add_parameter("prec_itmax","1"); 
    if (parameters_->value("prec_restol") == "")  
      parameters_->add_parameter("prec_restol","0.0");
    if (parameters_->value("prec_rlxtype") == "")  
      parameters_->add_parameter("prec_rlxtype","2");
    if (parameters_->value("prec_npre") == "")  
      parameters_->add_parameter("prec_npre","2"); 
    if (parameters_->value("prec_npost") == "")  
      parameters_->add_parameter("prec_npost","2");
    if (parameters_->value("prec_printl") == "")  
      parameters_->add_parameter("prec_printl","1");
    if (parameters_->value("prec_log") == "")  
      parameters_->add_parameter("prec_log","1");
    if (parameters_->value("prec_Jaciters") == "")  
      parameters_->add_parameter("prec_Jaciters","3");

    precond = new AMRsolve_HG_prec(*hierarchy_, BdryType_);
    ierr = precond->Initialize_(parameters_, &Ac_, &Xc_, &Bc_, &Y_);
    if (ierr != 0)  ERROR("could not initialize HG preconditioner\n");
    ierr = HYPRE_SStructGMRESSetPrecond(solver_,
					(HYPRE_PtrToSStructSolverFcn) HG_prec_solve,
					(HYPRE_PtrToSStructSolverFcn) HG_prec_setup,
					(HYPRE_SStructSolver) precond);
    if (ierr != 0)  ERROR("could not set preconditioner\n");
  }

  // Initialize the solver
  ierr = HYPRE_SStructGMRESSetup(solver_, A_, B_, X_);
  if (ierr != 0)  ERROR("could not setup solver_\n");

  // Solve the linear system
  ierr = HYPRE_SStructGMRESSolve(solver_, A_, B_, X_);
  if (ierr != 0) {
    HYPRE_SStructGMRESGetNumIterations(solver_, &iter_);
    HYPRE_SStructGMRESGetFinalRelativeResidualNorm(solver_, &resid_);
    fprintf(stderr, "Error %i detected, GMRES solver stats:\n  iterations = %d\n  final relative residual = %g\n",ierr,iter_,resid_);
    fprintf(stderr, "Error %i detected, writing B_ to 'B-hypre.*':\n",ierr);
    HYPRE_SStructVectorPrint("B-hypre",B_,0);
    fprintf(stderr, "Error %i detected, writing X_ to 'X-hypre.*':\n",ierr);
    HYPRE_SStructVectorPrint("X-hypre",X_,0);
    fprintf(stderr, "Error %i detected, writing A_ to 'A-hypre.*':\n",ierr);
    HYPRE_SStructMatrixPrint("A-hypre",A_,0);
    if (use_prec_) {
      fprintf(stderr, "Error %i detected, writing Bc_ to 'Bc-hypre.*':\n",ierr);
      HYPRE_StructVectorPrint("Bc-hypre",Bc_,0);
      fprintf(stderr, "Error %i detected, writing Xc_ to 'Xc-hypre.*':\n",ierr);
      HYPRE_StructVectorPrint("Xc-hypre",Xc_,0);
      fprintf(stderr, "Error %i detected, writing Ac_ to 'Ac-hypre.*':\n",ierr);
      HYPRE_StructMatrixPrint("Ac-hypre",Ac_,0);
    }
    fprintf(stderr, "Error %i detected, printing current hierarchy:\n",ierr);
    hierarchy_->print();
    ERROR("could not solve with GMRES\n");
  }

  // Write out some diagnostic info about the solve
  ierr = HYPRE_SStructGMRESGetNumIterations(solver_, &iter_);
  if (ierr != 0)  ERROR("could not get iter_\n");
  ierr = HYPRE_SStructGMRESGetFinalRelativeResidualNorm(solver_, &resid_);
  if (ierr != 0)  ERROR("could not get resid_\n");
  if (use_prec_) {
    Scalar presid = precond->GetResid_();
    int piters = precond->GetIters_();
    if (debug && pmpi->is_root()) {
      printf("hypre GMRES num iterations: %d,  prec: %d\n",iter_,piters);
      printf("hypre GMRES final relative residual norm: %g,  prec: %g\n",resid_,presid);
    }
  } else {
    if (debug && pmpi->is_root()) {
      printf("hypre GMRES num iterations: %d\n",iter_);
      printf("hypre GMRES final relative residual norm: %g\n",resid_);
    }
  }

  // project solutions from finest through coarsest grids in overlapped regions
  // (currently requires that preconditioner has been set up)
  if (use_prec_) {
    ierr = precond->HYPRE_to_AMRsolve_(&X_, 1); // copy X_ to AMRsolve's u_
    if (ierr != 0)  ERROR("could not copy X_ to u_\n");
    ierr = precond->restrict_(hierarchy_->num_levels()-1, 0);  // restriction
    if (ierr != 0)  ERROR("could not restrict u_ to coarse levels\n");
    ierr = precond->AMRsolve_to_HYPRE_(&X_, 1); // copy u_ to X_
    if (ierr != 0)  ERROR("could not copy u_ to X_\n");
  }

  // Delete the solver
  ierr = HYPRE_SStructGMRESDestroy(solver_);
  if (ierr != 0)  ERROR("could not destroy GMRES solver\n");
  if (use_prec_)  delete precond;

  _TRACE_;
} // AMRsolve_Hypre_Grav::solve_gmres_()

//------------------------------------------------------------------------

/// Update the matrix at a coarse/fine interface; this routine changes the
/// fine-grid matrix to account for the coarse grid neighbor.  Called in
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
  int ierr;
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
	ierr = HYPRE_SStructGraphAddEntries(graph_, level_fine, index_fine, 
					    0, level_coarse, index_coarse, 0);
	if (ierr != 0)  ERROR("could not add graph entries\n");
	index_fine[axis0] += index_increment[k][0];
	index_fine[axis1] += index_increment[k][1];
	index_fine[axis2] += index_increment[k][2];
      } // for k = 1:4

    } else if (phase == phase_matrix) {

      // fine->coarse off-diagonal
      double val_s = 2.0 / 3.0;
      int entry;
      double val, value;

      int k=0;
      index_fine[axis0] += index_increment[k][0];
      index_fine[axis1] += index_increment[k][1];
      index_fine[axis2] += index_increment[k][2];

      for (k=1; k<5; k++) {

	// remove old off-diagonal, stencil xp=1,xm,yp,ym,zp,zm=6
	entry = 2*axis0 + 1 + (1-face);
	value = 0.0;
	ierr = HYPRE_SStructMatrixSetValues(A_, level_fine, index_fine, 
					    0, 1, &entry, &value);
	if (ierr != 0)  ERROR("could not AddToValues in A_\n");

	// set new off-diagonal (2/3 of the normal weight)
	entry = grid_fine.counter(index_fine)++;
	value = matrix_scale_ * val_h_fine * val_s;
	ierr = HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
					      0, 1, &entry, &value);
	if (ierr != 0)  ERROR("could not AddToValues in A_\n");

	// update diagonal (add 1/3 of neighbor weight due to linear interp.)
	entry = 0;
	value = matrix_scale_ * val_h_fine * (1.0 - val_s);
	ierr = HYPRE_SStructMatrixAddToValues(A_, level_fine, index_fine, 
					      0, 1, &entry, &value);
	if (ierr != 0)  ERROR("could not AddToValues in A_\n");

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
      if (ierr != 0)  ERROR("could not add entries to graph_\n");
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
    if (ierr != 0)  ERROR("could not GetValues from A_\n");
    ierr = HYPRE_SStructMatrixSetValues(A_, level_coarse, index_coarse, 
					0, 1, &entry, &val2);
    if (ierr != 0)  ERROR("could not SetValues in A_\n");
    
    // Set new values across interface, scaling to depend equally on each fine neighbor
    val2 = val / 8.0;
    for (int i=0; i<8; i++) {
      entry = grid_coarse.counter(index_coarse)++;
      ierr = HYPRE_SStructMatrixSetValues(A_, level_coarse, index_coarse, 
					  0, 1, &entry, &val2);
      if (ierr != 0)  ERROR("could not SetValues in A_\n");
    }  // for i=0,7

  } // if phase == phase_matrix
}

//------------------------------------------------------------------------
