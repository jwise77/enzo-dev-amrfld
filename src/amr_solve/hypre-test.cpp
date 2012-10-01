// $Id: hypre-test.cpp 14 2010-03-19 19:49:39Z jbordner $

/// @file      hypre-test.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @brief     Minimalistic driver for hypre AMR problem 1: two nested grids on two processors

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>

#include "HYPRE_sstruct_ls.h"

//----------------------------------------------------------------------
// DEBUG
//----------------------------------------------------------------------

const int    trace       = 0;
const int    trace_hypre = 0;

#define TRACE \
  if (trace) printf ("TRACE %d %s:%d\n",mpi_rank,__FILE__,__LINE__); \
  fflush(stdout);

//----------------------------------------------------------------------
// PARAMETERS
//----------------------------------------------------------------------

const double MASS        = 1e43;
const double BOX         = 8e9;
const int    itmax       = 50;
const double restol      = 1e-6;

//----------------------------------------------------------------------
// CONSTANTS
//----------------------------------------------------------------------

const double PI             = M_PI;
const double G              = 6.67428e-8;
const int    num_parts      = 2; // Number of hypre "parts"
const int    part_coarse    = 0; // Coarse part id
const int    part_fine      = 1; // Fine part id
const int    r              = 2; // refinement factor TESTED FOR 2 ONLY

//----------------------------------------------------------------------
// FUNCTIONS
//----------------------------------------------------------------------

int min (int i, int j) { return i<j ? i : j;};
int index(int i0, int i1, int i2, int N) {
  return i0 + N*(i1 + N*i2);
}
void usage(int mpi_rank,char ** argv) {
  if (mpi_rank==0) {
    fprintf (stderr, "Usage: %s <int size> [\"cg\"|\"mg\"]\n",argv[0]);
  }
  MPI_Finalize();

  exit(1);
}

//----------------------------------------------------------------------
// MACROS
//----------------------------------------------------------------------

#define ASSERT_BOUND(A,B,C) \
  { \
    if ((B) < (A) || (C) < (B) ) { \
      printf ("%s:%d  ASSERT_BOUND: %d %d %d\n", \
	      __FILE__,__LINE__,A,B,C); \
      assert (0); \
    } \
  }

#define TEMPORARY(X) { \
   printf ("%s:%d TEMPORARY: " X " is temporary code that should be removed\n",__FILE__,__LINE__); \
   fflush(stdout); \
}

//======================================================================
// GLOBALS!
//======================================================================

FILE *mpi_fp;
char mpi_file[20];

//======================================================================
// MAIN
//======================================================================

int main(int argc, char * argv[])
{

  int i,i0,i1,i2;

  //------------------------------------------------------------
  // INITIALIZE MPI
  //    Out: mpi_size      (number of MPI processors)
  //    Out: mpi_rank      (rank of this MPI process)
  //    Out: is_mpi_coarse (whether MPI process owns coarse grid)
  //    Out: is_mpi_fine   (whether MPI process owns fine grid)
  //------------------------------------------------------------

  int mpi_rank,mpi_size;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);

  // MPI rank for grids

  const int mpi_rank_coarse   = 0;            
  const int mpi_rank_fine     = min(1,mpi_size-1);
  assert (mpi_rank_fine == 0 || mpi_rank_fine == 1);

  // Boolean variables for whether MPI process owns coarse or fine grids

  const bool is_mpi_coarse = (mpi_rank == mpi_rank_coarse);
  const bool is_mpi_fine   = (mpi_rank == mpi_rank_fine);

  //------------------------------------------------------------
  // PARSE AND CHECK ARGUMENTS
  //    Out: int N        (problem size)
  //    Out: int use_fac  (use FAC solver?)
  //------------------------------------------------------------

  int N;
  char * arg_solver;
  int use_fac = 0;

  if (argc == 3) {

    // Argument 1

    N       = atoi (argv[1]);

    // Argument 2

    arg_solver = argv[2];
    if (strcmp(arg_solver,"cg") == 0) {
      use_fac = 0;
    } else if (strcmp(arg_solver,"mg") == 0) {
      use_fac = 1;
    } else {
      usage(mpi_rank,argv);
    }

  } else {

    usage(mpi_rank,argv);

  }

  //------------------------------------------------------------
  // OPEN HYPRE TRACE OUTPUT FILE
  //------------------------------------------------------------

  if (trace_hypre) {
    sprintf (mpi_file,"hypre-test.out.%d",mpi_rank);
    mpi_fp = fopen (mpi_file,"w");
  } // trace_hypre

  //------------------------------------------------------------
  // CREATE AND INITIALIZE HYPRE 'GRID'
  //    Out: grid  (a two-level hypre grid)
  //------------------------------------------------------------

  HYPRE_SStructGrid grid;

  // Create the grid

  HYPRE_SStructGridCreate (MPI_COMM_WORLD, 3, num_parts, &grid);
  if (trace_hypre) {
    fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGridCreate ()\n",
	     __FILE__,__LINE__,mpi_rank);
    fflush(mpi_fp);
  } // trace_hypre


  // Initialize coarse grid extents on coarse grid MPI process

  HYPRE_SStructVariable variable_type[] = { HYPRE_SSTRUCT_VARIABLE_CELL };

  int lower_coarse[3] = {    0,      0,      0};
  int upper_coarse[3] = {  N-1,    N-1,    N-1};

  if (is_mpi_coarse) {

    HYPRE_SStructGridSetExtents(grid, part_coarse, lower_coarse, upper_coarse);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGridSetExtents()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  }

  // Initialize fine grid extents on fine grid MPI process

  int lower_fine[3]   = {  N/r,    N/r,    N/r};
  int upper_fine[3]   = {3*N/r-1,3*N/r-1,3*N/r-1};

  if (is_mpi_fine) {

    HYPRE_SStructGridSetExtents(grid, part_fine,   lower_fine,   upper_fine);

    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGridSetExtents()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  }

  HYPRE_SStructGridSetVariables(grid, part_fine,   1, variable_type);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGridSetVariables()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructGridSetVariables(grid, part_coarse, 1, variable_type);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGridSetVariables()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  // Assemble the grid
  
  HYPRE_SStructGridAssemble (grid);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGridAssemble()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  //------------------------------------------------------------
  // CREATE HYPRE STENCIL
  //    Out: stencil  (structured connections)
  //------------------------------------------------------------

  HYPRE_SStructStencil stencil;
  
  int entries[][3] = { {  0, 0, 0 },
		       {  1, 0, 0 },
		       { -1, 0, 0 },
		       {  0, 1, 0 },
		       {  0,-1, 0 },
		       {  0, 0, 1 },
		       {  0, 0,-1 } };

  HYPRE_SStructStencilCreate (3,7,&stencil);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructStencilCreate()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  HYPRE_SStructStencilSetEntry (stencil, 0, entries[0], 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructStencilSetEntry()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructStencilSetEntry (stencil, 1, entries[1], 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructStencilSetEntry()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructStencilSetEntry (stencil, 2, entries[2], 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructStencilSetEntry()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructStencilSetEntry (stencil, 3, entries[3], 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructStencilSetEntry()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructStencilSetEntry (stencil, 4, entries[4], 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructStencilSetEntry()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructStencilSetEntry (stencil, 5, entries[5], 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructStencilSetEntry()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructStencilSetEntry (stencil, 6, entries[6], 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructStencilSetEntry()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  //------------------------------------------------------------
  // CREATE HYPRE GRAPH
  //    Out: graph  (unstructured connections)
  //------------------------------------------------------------

  HYPRE_SStructGraph   graph;

  // Create the graph

  HYPRE_SStructGraphCreate (MPI_COMM_WORLD, grid, &graph);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphCreate()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  // Initialize stencil part

  HYPRE_SStructGraphSetStencil (graph, part_coarse, 0, stencil);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphSetStencil()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructGraphSetStencil (graph, part_fine,   0, stencil);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphSetStencil()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  // Initialize nonstencil part

  int axis,face;     // 0 <= axis <= 2;  0 <= face <= 1
  int ic0,ic1;       // coarse face indices
  int ind_fine[3];   // fine grid indices
  int ind_coarse[3]; // coarse grid indices

  //----------------------------------------
  // GRAPH ENTRIES: FINE-TO-COARSE
  //----------------------------------------

  if (is_mpi_fine) {

    for (axis=0; axis<3; axis++) {

      int j0 = axis;
      int j1 = (axis+1)%3;
      int j2 = (axis+2)%3;

      for (face=0; face<2; face++) {

	// loop over coarse face zones

	for (ic0=0; ic0<N/r; ic0++) {
	  for (ic1=0; ic1<N/r; ic1++) {

	    // coarse zone index

	    ind_coarse[j0] =    (1-face) * (N/4-1) + (face) * (3*N/4);
	    ind_coarse[j1] = N/4 + ic0;
	    ind_coarse[j2] = N/4 + ic1;

	    // fine zone index 000

	    ind_fine[j0]   = r*((1-face) * (N/4)   + (face) * (3*N/4-1));
	    ind_fine[j1]   = r*ind_coarse[j1];
	    ind_fine[j2]   = r*ind_coarse[j2];

	    // ind_fine computed above corresponds to coarse zone
	    if (face == 1)  ++ ind_fine[j0];

	    // 000 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries 
	      (graph, part_fine, ind_fine, 0, part_coarse, ind_coarse, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    ++ ind_fine[j1]; 	  
	    // 010 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries
	      (graph, part_fine, ind_fine, 0, part_coarse, ind_coarse, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    ++ ind_fine[j2];	  
	    // 011 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries
	      (graph, part_fine, ind_fine, 0, part_coarse, ind_coarse, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    -- ind_fine[j1];	  
	    // 001 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries
	      (graph, part_fine, ind_fine, 0, part_coarse, ind_coarse, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    -- ind_fine[j2];	  
	    // 000 ---------------------------------------------
	  }
	}
      }
    }
  }

  //----------------------------------------
  // GRAPH ENTRIES: COARSE-TO-FINE
  //----------------------------------------

  if (is_mpi_coarse) {

    for (axis=0; axis<3; axis++) {

      int j0 = axis;
      int j1 = (axis+1)%3;
      int j2 = (axis+2)%3;

      for (face=0; face<2; face++) {

	// loop over coarse face zones

	for (ic0=0; ic0<N/r; ic0++) {
	  for (ic1=0; ic1<N/r; ic1++) {

	    // coarse zone index

	    ind_coarse[j0] =    (1-face) * (N/4-1) + (face) * (3*N/4);
	    ind_coarse[j1] = N/4 + ic0;
	    ind_coarse[j2] = N/4 + ic1;

	    // fine zone index 000

	    ind_fine[j0]   = r*((1-face) * (N/4)   + (face) * (3*N/4-1));
	    ind_fine[j1]   = r*ind_coarse[j1];
	    ind_fine[j2]   = r*ind_coarse[j2];

	    HYPRE_SStructGraphAddEntries (graph, 
					  part_coarse, ind_coarse, 0, 
					  part_fine, ind_fine, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    ++ ind_fine[j0]; 	  

	    // 100 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries 
	      (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    ++ ind_fine[j1];	  
	    // 110 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries 
	      (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    -- ind_fine[j0];	  
	    // 010 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries 
	      (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    ++ ind_fine[j2];	  
	    // 011 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries 
	      (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    ++ ind_fine[j0]; 	  
	    // 111 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries 
	      (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    -- ind_fine[j1];	  
	    // 101 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries 
	      (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    -- ind_fine[j0];	  
	    // 001 ---------------------------------------------
	    HYPRE_SStructGraphAddEntries 
	      (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAddEntries()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    -- ind_fine[j2];	  
	    // 000 ---------------------------------------------
	  }
	}
      }
    }
  }

  HYPRE_SStructGraphAssemble (graph);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGraphAssemble()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  //============================================================
  // CREATE HYPRE MATRIX A AND VECTORS X,B
  //============================================================

  //--------------------------------------------------
  // Create the hypre vector right-hand side B
  //--------------------------------------------------

  HYPRE_SStructVector  B;       
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid,  &B);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorCreate()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructVectorSetObjectType (B,HYPRE_SSTRUCT);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SSTRUCT()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructVectorInitialize (B);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorInitialize()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  double fine_h   = 0.5 * BOX / N;  // Fine mesh width

  int ind3[8][3] = {
    {N-1, N-1, N-1},
    {N-1, N-1, N},
    {N-1, N,   N-1},
    {N-1, N,   N},
    {N,   N-1, N-1},
    {N,   N-1, N},
    {N,   N,   N-1},
    {N,   N,   N}};

  if (is_mpi_fine) {

    double bval = -4.0 * G * PI * MASS / (fine_h*fine_h*fine_h);

    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[0], 0, &bval);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[1], 0, &bval);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[2], 0, &bval);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[3], 0, &bval);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[4], 0, &bval);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[5], 0, &bval);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[6], 0, &bval);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[7], 0, &bval);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  }

  HYPRE_SStructVectorAssemble (B);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorAssemble()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  // Create the hypre matrix A

  HYPRE_SStructMatrix  A;
  HYPRE_SStructMatrixCreate (MPI_COMM_WORLD, graph, &A);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixCreate()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructMatrixSetObjectType (A,HYPRE_SSTRUCT);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SSTRUCT()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructMatrixInitialize (A);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixInitialize()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  // Initialize A

  // Stencil entries

  int nums[7] = {0,1,2,3,4,5,6};

  double * v0  = new double [N*N*N];
  double * vxp = new double [N*N*N];
  double * vxm = new double [N*N*N];
  double * vyp = new double [N*N*N];
  double * vym = new double [N*N*N];
  double * vzp = new double [N*N*N];
  double * vzm = new double [N*N*N];

  double coarse_h =       BOX / N;  // Coarse mesh width

  if (is_mpi_coarse) {

    //--------------------------------------------------
    // MATRIX STENCIL VALUES: GRID INTERIOR
    //--------------------------------------------------

    for (i=0; i<N*N*N; i++) {
      v0[i]  = -6*coarse_h;
      vxp[i] =  1*coarse_h;
      vxm[i] =  1*coarse_h;
      vyp[i] =  1*coarse_h;
      vym[i] =  1*coarse_h;
      vzp[i] =  1*coarse_h;
      vzm[i] =  1*coarse_h;
    }

    //--------------------------------------------------
    // SET HYPRE MATRIX STENCIL VALUES
    //--------------------------------------------------

    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[0],v0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[1],vxp);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[2],vxm);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[3],vyp);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[4],vym);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[5],vzp);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[6],vzm);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  }

  if (is_mpi_fine) {

    //--------------------------------------------------
    // MATRIX STENCIL VALUES: GRID INTERIOR
    //--------------------------------------------------

    for (i=0; i<N*N*N; i++) {
      v0[i]  = -6*fine_h;
      vxp[i] =  1*fine_h;
      vxm[i] =  1*fine_h;
      vyp[i] =  1*fine_h;
      vym[i] =  1*fine_h;
      vzp[i] =  1*fine_h;
      vzm[i] =  1*fine_h;
    }

    //--------------------------------------------------
    // MATRIX STENCIL VALUES: GRID BOUNDARY
    //--------------------------------------------------

    for (i0=0; i0<N; i0++) {
      for (i1=0; i1<N; i1++) {
	vxp[index(N-1,i0,i1,N)] -= fine_h;
	v0 [index(N-1,i0,i1,N)] += fine_h;
	vxm[index(  0,i0,i1,N)] -= fine_h;
	v0 [index(  0,i0,i1,N)] += fine_h;
	vyp[index(i0,N-1,i1,N)] -= fine_h;
	v0 [index(i0,N-1,i1,N)] += fine_h;
	vym[index(i0,  0,i1,N)] -= fine_h;
	v0 [index(i0,  0,i1,N)] += fine_h;
	vzp[index(i0,i1,N-1,N)] -= fine_h;
	v0 [index(i0,i1,N-1,N)] += fine_h;
	vzm[index(i0,i1,  0,N)] -= fine_h;
	v0 [index(i0,i1,  0,N)] += fine_h;
      }
    }

    //--------------------------------------------------
    // SET HYPRE MATRIX STENCIL VALUES
    //--------------------------------------------------

    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[0],v0);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[1],vxp);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[2],vxm);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[3],vyp);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[4],vym);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[5],vzp);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[6],vzm);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixSetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  }

  //--------------------------------------------------
  // MATRIX CLEAR STENCIL OVERLAP
  //--------------------------------------------------

  int refinements[3] = {2,2,2};

  // removing breaks MG

  HYPRE_SStructFACZeroCFSten (A,grid, part_fine, refinements);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACZeroCFSten()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  // removing apparently does nothing

  HYPRE_SStructFACZeroFCSten (A,grid, part_fine);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACZeroFCSten()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  // removing breaks CG

  HYPRE_SStructFACZeroAMRMatrixData (A, part_coarse, refinements);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACZeroAMRMatrixData()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
      

  delete [] v0;
  delete [] vxp;
  delete [] vxm;
  delete [] vyp;
  delete [] vym;
  delete [] vzp;
  delete [] vzm;

  //--------------------------------------------------
  // MATRIX NONSTENCIL ENTRIES
  //--------------------------------------------------

  // Declare counts for zones
  int * count_coarse = new int [N*N*N];
  int * count_fine   = new int [N*N*N];

  for (i=0; i<N*N*N; i++) {
    count_coarse[i] = 7; // skip the stencil elements:
    count_fine[i]   = 7; // graph entries start at 7
  }

  //--------------------------------------------------
  // MATRIX ENTRIES: FINE-TO-CORSE 
  //--------------------------------------------------

  if (is_mpi_fine) {

    for (axis=0; axis<3; axis++) {

      int j0 = axis;
      int j1 = (axis+1)%3;
      int j2 = (axis+2)%3;

      for (face=0; face<2; face++) {

	// loop over coarse face zones

	for (ic0=0; ic0<N/r; ic0++) {
	  for (ic1=0; ic1<N/r; ic1++) {

	    // coarse zone index

	    ind_coarse[j0] =    (1-face) * (N/4-1) + (face) * (3*N/4);
	    ind_coarse[j1] = N/4 + ic0;
	    ind_coarse[j2] = N/4 + ic1;

	    // fine zone index 000

	    ind_fine[j0]   = r*((1-face) * (N/4)   + (face) * (3*N/4-1));
	    ind_fine[j1]   = r*ind_coarse[j1];
	    ind_fine[j2]   = r*ind_coarse[j2];
  
	    // ind_fine computed above corresponds to coarse zone
	    if (face == 1)  ++ ind_fine[j0];

	    int    entry;
	    double value;

	    int icount;
	    int ishift = -index(N/r,N/r,N/r,N);
	    int num_entries = 1;
	  
	    double a = 2.0/3.0;

	    // 000 ---------------------------------------------
	    icount = ishift + index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	    ASSERT_BOUND(0,icount,N*N*N);
	    // update diagonal
	    entry = count_fine[icount]++;
	    value = a * fine_h;
	    HYPRE_SStructMatrixAddToValues 
	      (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    // update off-diagonal
	    entry = 0;
	    value = -value;
	    HYPRE_SStructMatrixAddToValues 
	      (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    ++ ind_fine[j1]; 	  
	    // 010 ---------------------------------------------
	    icount = ishift + index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	    ASSERT_BOUND(0,icount,N*N*N);
	    // update diagonal
	    entry = count_fine[icount]++;
	    value = a * fine_h;
	    HYPRE_SStructMatrixAddToValues 
	      (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    // update off-diagonal
	    entry = 0;
	    value = -value;
	    HYPRE_SStructMatrixAddToValues 
	      (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    ++ ind_fine[j2];	  
	    // 011 ---------------------------------------------
	    icount = ishift + index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	    ASSERT_BOUND(0,icount,N*N*N);
	    // update diagonal
	    entry = count_fine[icount]++;
	    value = a * fine_h;
	    HYPRE_SStructMatrixAddToValues 
	      (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    // update off-diagonal
	    entry = 0;
	    value = -value;
	    HYPRE_SStructMatrixAddToValues 
	      (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    -- ind_fine[j1];	  
	    // 001 ---------------------------------------------
	    icount = ishift + index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	    ASSERT_BOUND(0,icount,N*N*N);
	    // update diagonal
	    entry = count_fine[icount]++;
	    value = a * fine_h;
	    HYPRE_SStructMatrixAddToValues 
	      (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    // update off-diagonal
	    entry = 0;
	    value = -value;
	    HYPRE_SStructMatrixAddToValues 
	      (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    -- ind_fine[j2];	  
	    // 000 ---------------------------------------------

	  }
	}
      }
    }
  }

  //--------------------------------------------------
  // MATRIX ENTRIES: COARSE-TO-FINE
  //--------------------------------------------------

  if (is_mpi_coarse) {

    for (axis=0; axis<3; axis++) {

      int j0 = axis;
      int j1 = (axis+1)%3;
      int j2 = (axis+2)%3;

      for (face=0; face<2; face++) {

	// loop over coarse face zones

	for (ic0=0; ic0<N/r; ic0++) {
	  for (ic1=0; ic1<N/r; ic1++) {

	    // coarse zone index

	    ind_coarse[j0] =    (1-face) * (N/4-1) + (face) * (3*N/4);
	    ind_coarse[j1] = N/4 + ic0;
	    ind_coarse[j2] = N/4 + ic1;

	    // fine zone index 000

	    ind_fine[j0]   = r*((1-face) * (N/4)   + (face) * (3*N/4-1));
	    ind_fine[j1]   = r*ind_coarse[j1];
	    ind_fine[j2]   = r*ind_coarse[j2];

	    int    entry;
	    double value;

	    int icount = index(ind_coarse[0],ind_coarse[1],ind_coarse[2],N);
	    int num_entries = 1;

	    for (int k=0; k<8; k++) {
	      ASSERT_BOUND(0,icount,N*N*N);
	      entry = count_coarse[icount]++;
	      value = (1./8.) * coarse_h;
	      HYPRE_SStructMatrixAddToValues 
		(A, part_coarse, ind_coarse, 0, num_entries, &entry, &value);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixAddToValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
	    }
	  }
	}
      }
    }
  }

  delete [] count_coarse;
  delete [] count_fine;

  HYPRE_SStructMatrixAssemble (A);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructMatrixAssemble()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  // Create the hypre vector solution X

  HYPRE_SStructVector  X;
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid,  &X);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorCreate()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructVectorSetObjectType (X,HYPRE_SSTRUCT);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SSTRUCT()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructVectorInitialize (X);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorInitialize()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  HYPRE_SStructVectorAssemble (X);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorAssemble()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  //------------------------------------------------------------
  // Solve the linear system
  //------------------------------------------------------------

  HYPRE_SStructSolver  solver;  // hypre solver
  int iter;
  double resid;

  if (use_fac) {
    
    //--------------------------------------------------
    // FAC SOLVER
    //--------------------------------------------------

    HYPRE_SStructFACCreate(MPI_COMM_WORLD, &solver);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACCreate()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

    int num_parts = 2;
    HYPRE_SStructFACSetMaxLevels(solver,  num_parts);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetMaxLevels()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

    HYPRE_SStructFACSetMaxIter(solver,itmax);
    //    HYPRE_SStructFACSetMaxIter(solver,20);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetMaxIter()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

    HYPRE_SStructFACSetTol(solver,    restol);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetTol()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre


    int parts[2] = {0,1};
    HYPRE_SStructFACSetPLevels(solver, num_parts, parts);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetPLevels()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

    int refinements[2][3] = {{1,1,1},{2,2,2}};
    HYPRE_SStructFACSetPRefinements(solver, num_parts, refinements);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetPRefinements()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

    TEMPORARY("Adding SStructFACSetRelChange");
    //    HYPRE_SStructFACSetRelChange(solver,0);

    HYPRE_SStructFACSetRelaxType(solver,        2);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetRelaxType()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

    HYPRE_SStructFACSetNumPreRelax(solver,      1);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetNumPreRelax()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructFACSetNumPostRelax(solver,     1);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetNumPostRelax()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

    HYPRE_SStructFACSetCoarseSolverType(solver, 2);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetCoarseSolverType()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

    HYPRE_SStructFACSetLogging(solver, 1);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetLogging()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

    TRACE;
    HYPRE_SStructFACSetup2(solver, A, B, X);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSetup2()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    TRACE;

    // SOLVE
    HYPRE_SStructFACSolve3(solver, A, B, X);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACSolve3()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

    HYPRE_SStructFACGetNumIterations(solver, &iter);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACGetNumIterations()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructFACGetFinalRelativeResidualNorm(solver, &resid);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACGetFinalRelativeResidualNorm()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre

  } else {

    //--------------------------------------------------
    // BiCGSTAB SOLVER
    //--------------------------------------------------

    HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructBiCGSTABCreate()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructBiCGSTABSetLogging(solver, 1);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructBiCGSTABSetLogging()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructBiCGSTABSetup(solver, A, B, X);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructBiCGSTABSetup()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    // SOLVE
    HYPRE_SStructBiCGSTABSolve(solver, A, B, X);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructBiCGSTABSolve()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructBiCGSTABGetNumIterations(solver, &iter);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructBiCGSTABGetNumIterations()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver, &resid);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  }

  printf ("solver:     %s\n", use_fac ? "FAC" : "BiCGSTAB");
  printf ("iterations: %d\n",iter);
  printf ("residual:   %g\n",resid);

  if (use_fac) {
    HYPRE_SStructFACDestroy2(solver);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructFACDestroy2()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  } else {
    HYPRE_SStructBiCGSTABDestroy(solver);
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructBiCGSTABDestroy()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
  }

  //------------------------------------------------------------
  // WRITE MATRIX
  //------------------------------------------------------------

  HYPRE_SStructMatrixPrint ("A",A,0);

  //------------------------------------------------------------
  // WRITE SOLUTION
  //------------------------------------------------------------

  HYPRE_SStructVectorPrint ("X",X,0);
  HYPRE_SStructVectorPrint ("B",B,0);

  FILE *fp;

  //-------------------------------------------------------------
  // WRITE COARSE GRID
  //-------------------------------------------------------------

  if (is_mpi_coarse) {

    fp = fopen("X.0","w");
    fprintf (fp,"Grid\n");
    fprintf (fp,"   id             0\n");
    fprintf (fp,"   parent id      -1\n");
    fprintf (fp,"   processor      %d\n",mpi_rank);
    fprintf (fp,"   lower position %d %d %d\n",
	     lower_coarse[0],lower_coarse[1],lower_coarse[2]);
    fprintf (fp,"   upper position %d %d %d\n",
	     upper_coarse[0],upper_coarse[1],upper_coarse[2]);
    fprintf (fp,"   lower index    0 0 0\n");
    fprintf (fp,"   zones          %d %d %d\n",N,N,N);
    fprintf (fp,"   level          0\n");

    double * x_coarse = new double [N*N*N];
    HYPRE_SStructVectorGetBoxValues (X,part_coarse,lower_coarse,upper_coarse,0,x_coarse);  
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorGetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    for (i0=0; i0<N; i0++) {
      for (i1=0; i1<N; i1++) {
	for (i2=0; i2<N; i2++) {
	  i = i0 + N*(i1 + N*i2);
	  fprintf (fp,"%d %d %d %g\n",i0,i1,i2,x_coarse[i]);
	}
      }
    }
    delete [] x_coarse;

    fclose(fp);

  }

  //-------------------------------------------------------------
  // WRITE FINE GRID
  //-------------------------------------------------------------

  if (is_mpi_fine) {

    fp = fopen("X.1","w");
    fprintf (fp,"Grid\n");
    fprintf (fp,"   id             1\n");
    fprintf (fp,"   parent id      0\n");
    fprintf (fp,"   processor      %d\n",mpi_rank);
    fprintf (fp,"   lower position %d %d %d\n",
	     lower_fine[0],lower_fine[1],lower_fine[2]);
    fprintf (fp,"   upper position %d %d %d\n",
	     upper_fine[0],upper_fine[1],upper_fine[2]);
    fprintf (fp,"   lower index    %d %d %d\n",
	     lower_fine[0],lower_fine[1],lower_fine[2]);
    fprintf (fp,"   zones          %d %d %d\n",N,N,N);
    fprintf (fp,"   level          1\n");

    double * x_fine = new double [N*N*N];
    HYPRE_SStructVectorGetBoxValues (X,part_fine,lower_fine,upper_fine,0,x_fine);  
    if (trace_hypre) {
      fprintf (mpi_fp,"%s:%d %d HYPRE_SStructVectorGetBoxValues()\n",
	       __FILE__,__LINE__,mpi_rank);
      fflush(mpi_fp);
    } // trace_hypre
    for (int i0=0; i0<N; i0++) {
      for (int i1=0; i1<N; i1++) {
	for (int i2=0; i2<N; i2++) {
	  i = i0 + N*(i1 + N*i2);
	  fprintf (fp,"%d %d %d %g\n",i0,i1,i2,x_fine[i]);
	}
      }
    }
    delete [] x_fine;

    fclose(fp);

  }

  //------------------------------------------------------------
  // EXIT SUCCESSFULLY
  //------------------------------------------------------------


  if (mpi_rank==0) {
    printf ("Finished!\n");
  }
  MPI_Finalize();

  exit (0);

}

