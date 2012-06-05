/// @file      AMRsolve_hypre_grav.h
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Declaration of the AMRsolve_Hypre_Grav class

#ifndef AMRSOLVE_HYPRE_GRAV_H
#define AMRSOLVE_HYPRE_GRAV_H

class AMRsolve_Hypre_Grav {

  /// @class    AMRsolve_Hypre_Grav
  /// @brief    AMRsolve_Hypre_Grav class for interfacing to the LLNL hypre solver

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
  Scalar               matrix_scale_;  // 1.0:  1 1 1 -6 1 1 1

public:

  AMRsolve_Hypre_Grav(AMRsolve_Hierarchy& hierarchy, 
		      AMRsolve_Parameters& parameters,
		      int precflag);

  ~AMRsolve_Hypre_Grav();

  void init_hierarchy(AMRsolve_Mpi& mpi);
  void init_stencil();
  void init_graph();
  /*   void init_elements(std::vector<AMRsolve_Point *> points,
                          Scalar f_scale=1.0); */
  void init_elements(Scalar f_scale=1.0);  // can't we just get our RHS from Enzo directly?
  void solve();
  int  evaluate();
  void update_enzo();
  void abort_dump();

  int    iterations() { return iter_; };
  double residual() { return resid_; };

private:

  enum phase_enum  {phase_unknown,phase_graph,phase_matrix};

  // init_graph() functions
  void init_graph_nonstencil_(AMRsolve_Grid& grid)
  { init_nonstencil_(grid, phase_graph); };

  // init_elements() functions
  void init_elements_matrix_();
  /*  void init_elements_rhs_(std::vector<AMRsolve_Point *>& points, 
      Scalar f_scale=1.0); */
  void init_elements_rhs_(Scalar f_scale=1.0);

  // init_matrix() functions
  void init_matrix_stencil_(AMRsolve_Grid& grid);
  void init_matrix_clear_(int part);
  void init_matrix_nonstencil_(AMRsolve_Grid& grid)
  { init_nonstencil_(grid, phase_matrix); };

  void init_nonstencil_(AMRsolve_Grid& grid, phase_enum phase);
  
  // init_vector() functions
  /*  Scalar init_vector_points_(std::vector<AMRsolve_Point *>& points);
  Scalar init_vector_file_(std::string file_prefix, bool is_packed);
  Scalar init_vector_attach_(Scalar f_scale=1.0); */

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

#endif /* AMRSOLVE_HYPRE_GRAV_H */
