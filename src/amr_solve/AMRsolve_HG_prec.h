/// @file      AMRsolve_HG_prec.h
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Declaration of the AMRsolve_HG_prec class

#ifndef AMRSOLVE_HG_PREC_H
#define AMRSOLVE_HG_PREC_H

class AMRsolve_HG_prec {

  /// @class    AMRsolve_HG_prec
  /// @brief    AMRsolve_HG_prec class for preconditioning a block-structured AMR linear solver

private:

  HYPRE_SStructVector *Y_;             // hypre vector temporary
  HYPRE_StructMatrix  *Ac_;            // coarse preconditioning matrix
  HYPRE_StructVector  *Bc_;            // hypre coarse vector right-hand side
  HYPRE_StructVector  *Xc_;            // hypre coarse vector solution
  HYPRE_StructSolver   csolver_;       // hypre coarse solver

  AMRsolve_Hierarchy  *hierarchy_;     // Pointer to the hierarchy

  double               resid_;         // Solver residual
  int                  citer_;         // Coarse solver iterations
  
  int                  Jacobi_iters;   // number of Jacobi iterations to apply on full grid

  int                  BdryType_[3][2];   // boundary condition types (global problem)

  bool                 initialized;    // flag denoting whether solver has been initialized

public:
  
  // Main solver interface routines
  AMRsolve_HG_prec(AMRsolve_Hierarchy& hierarchy, int BdryType[3][2]);
  ~AMRsolve_HG_prec();
  int Initialize_(AMRsolve_Parameters *parameters, 
		  HYPRE_StructMatrix *A,
		  HYPRE_StructVector *X,
		  HYPRE_StructVector *B,
		  HYPRE_SStructVector *Y);
  int Setup_(HYPRE_SStructMatrix A, 
	     HYPRE_SStructVector b, 
	     HYPRE_SStructVector x);
  int Solve_(HYPRE_SStructMatrix A, 
	     HYPRE_SStructVector b, 
	     HYPRE_SStructVector x);
  Scalar GetResid_() { return(resid_); }
  int    GetIters_() { return(citer_); }

  // Utility routines to copy data between AMRsolve and HYPRE formats
  int AMRsolve_to_HYPRE_(HYPRE_SStructVector *V, int u_vs_f);
  int HYPRE_to_AMRsolve_(HYPRE_SStructVector *V, int u_vs_f);
  int AMRsolve_to_HYPRE_coarse_(HYPRE_StructVector *V, int u_vs_f);
  int HYPRE_to_AMRsolve_coarse_(HYPRE_StructVector *V, int u_vs_f);

  // Restriction/Prolongation operators over hierarchy
  int restrict(int level_fine, int level_coarse) throw();
  void do_restrict(Scalar *mydata, Scalar *overlap, int *mysize, int *g1_ilo, 
		   int *g1_ihi, int *ovsize) throw();
  int prolong(int level_coarse, int level_fine, int method) throw();
  void do_prolong(Scalar *mydata, Scalar *overlap, int *mysize, int *g1_ilo, 
		  int *g1_ihi, int *ovsize, int method) throw();

  // Jacobi smoother on full hierarchy
  int Jacobi_smooth_(HYPRE_SStructMatrix A,
		     HYPRE_SStructVector v, 
		     HYPRE_SStructVector b,
		     HYPRE_SStructVector y);

};


// wrapper routines for passing HG_prec_setup_ and HG_prec_solve_ to HYPRE solvers
HYPRE_Int HG_prec_setup(HYPRE_SStructSolver solver, HYPRE_SStructMatrix A, 
			HYPRE_SStructVector b, HYPRE_SStructVector x);
HYPRE_Int HG_prec_solve(HYPRE_SStructSolver solver, HYPRE_SStructMatrix A, 
			HYPRE_SStructVector b, HYPRE_SStructVector x);

#endif /* AMRSOLVE_HG_PREC_H */
