/// @file      AMRsolve_hypre_fld.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Declaration of the AMRsolve_Hypre_FLD class

#ifndef AMRSOLVE_HYPRE_FLD_H
#define AMRSOLVE_HYPRE_FLD_H

class AMRsolve_Hypre_FLD {

  /// @class    AMRsolve_Hypre_FLD
  /// @brief    AMRsolve_Hypre_FLD class for interfacing to the LLNL hypre solver

private:

  HYPRE_SStructGrid    grid_;          // hypre grid
  HYPRE_SStructGraph   graph_;         // hypre graph
  HYPRE_SStructStencil stencil_;       // hypre stencil
  HYPRE_SStructMatrix  A_;             // hypre matrix
  HYPRE_SStructVector  B_;             // hypre vector right-hand side
  HYPRE_SStructVector  X_;             // hypre vector solution
  HYPRE_SStructVector  Y_;             // hypre vector temporary
  HYPRE_SStructSolver  solver_;        // hypre solver

  bool                 use_prec;       // flag whether to setup/use preconditioner (0 -> no)
  HYPRE_StructGrid     cgrid_;         // hypre coarse grid (structured)
  HYPRE_StructStencil  cstencil_;      // hypre coarse stencil (structured)
  HYPRE_StructMatrix   Ac_;            // coarse preconditioning matrix
  HYPRE_StructVector   Bc_;            // hypre coarse vector right-hand side
  HYPRE_StructVector   Xc_;            // hypre coarse vector solution

  AMRsolve_Parameters *parameters_;    // Pointer to parameters
  AMRsolve_Hierarchy  *hierarchy_;     // Pointer to the hierarchy

  double               resid_;         // Solver residual
  int                  iter_;          // Solver iterations
  int                  citer_;         // Coarse solver iterations

  const int            r_factor_;      // Refinement factor
  int                  bin_;           // current radiation bin
  int                  Nchem_;         // number of chemical species
  double               theta_;         // time discretization parameter
  double               dt_;            // time step size
  double               aval_;          // cosmological expansion constant
  double               aval0_;         // cosmological expansion constant (old time)
  double               adot_;          // rate of cosmological expansion
  double               adot0_;         // rate of cosmological expansion (old time)
  double               nUn_;           // number density units
  double               nUn0_;          // number density units (old time)
  double               lUn_;           // length units
  double               lUn0_;          // length units (old time)
  double               rUn_;           // radiation units
  double               rUn0_;          // radiation units (old time)
  int                  BdryType_[3][2];   // boundary condition types (global problem)

public:

  AMRsolve_Hypre_FLD(AMRsolve_Hierarchy& hierarchy, 
		     AMRsolve_Parameters& parameters,
		     int bin, int precflag);
  ~AMRsolve_Hypre_FLD();

  void init_hierarchy();
  void init_stencil();
  void init_graph();
  void init_elements(double dt, int Nchem, double theta, double aval, 
		     double aval0, double adot, double adot0, 		     
		     double nUn, double nUn0, double lUn, double lUn0, 
		     double rUn, double rUn0, int BdryType[3][2]);
  double rdiff_norm(double pnorm, double atol);
  void solve();
  int  evaluate();
  void update_enzo();
  void abort_dump();
  void tester();

  int    iterations() { return iter_; };
  double residual() { return resid_; };

private:

  enum phase_enum  {phase_unknown,phase_graph,phase_matrix};

  // init_graph() functions
  void init_graph_nonstencil_(AMRsolve_Grid& grid)
  { init_nonstencil_(grid, phase_graph); };

  // init_elements() functions
  void init_elements_matrix_();
  void init_elements_rhs_();

  // init_matrix() functions
  void init_matrix_stencil_(AMRsolve_Grid& grid);
  void init_matrix_clear_(int part);
  void init_matrix_nonstencil_(AMRsolve_Grid& grid)
  { init_nonstencil_(grid, phase_matrix); };

  void init_nonstencil_(AMRsolve_Grid& grid, phase_enum phase);

  Scalar limiter_(Scalar E1, Scalar E2, Scalar k1, Scalar k2, 
		  Scalar nUn, Scalar lUn, Scalar dxi);
  
  // solve() functions
  void solve_fac_(int itmax, double restol);
  void solve_bicgstab_(int itmax, double restol);
  void solve_bicgstab_boomer_(int itmax, double restol);
  void solve_gmres_(int itmax, double restol);
  void solve_pfmg_(int itmax, double restol);

  // matrix graph update functions
  void update_fine_coarse_const_(int face, AMRsolve_Grid& grid, 
				 int axis0, phase_enum phase,
				 int level_fine, int level_coarse,
				 int igg3[3], int ign3[3]);

  void update_coarse_fine_const_(int face, AMRsolve_Grid& grid, 
				 int axis0, phase_enum phase,
				 int level_fine, int level_coarse,
				 int igg3[3], int ign3[3]);
};

#endif /* AMRSOLVE_HYPRE_FLD_H */
