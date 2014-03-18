/// @file      AMRsolve_fld.h
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Declaration of the AMRsolve_FLD class
#ifndef AMRSOLVE_HYPRE_H
#define AMRSOLVE_HYPRE_H

class AMRsolve_FLD {

  /// @class    AMRsolve_FLD
  /// @brief    AMRsolve_FLD class for solving hierarchical FLD systems
  /// 
  ///           This solver differs from AMRsolve_Hypre_FLD in that most
  ///           of the algorithm occurs internally in a matrix-free fashion
  ///           (without creating the large/multilevel matrices), although 
  ///           the uniform coarse grid is still preconditioned using the 
  ///           LLNL HYPRE solver

private:

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

  int                  level_start;    // coarsest level to use in solve
  int                  level_end;      // finest level to use in solve

  int                  iNULL;          // flag indicating no vector argument
  int                  ix;             // flag indicating to use Grid's u_ vector
  int                  ib;             // flag indicating to use Grid's f_ vector
  int                  iv;             // index of BiCGStab "v" vector in temp array
  int                  ir;             // index of BiCGStab "r" vector in temp array
  int                  irs;            // index of BiCGStab "rs" vector in temp array
  int                  ip;             // index of BiCGStab "p" vector in temp array
  int                  is;             // index of BiCGStab "s" vector in temp array
  int                  iq;             // index of BiCGStab "q" vector in temp array

  const int            r_factor_;      // Refinement factor
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
		     int precflag);
  ~AMRsolve_Hypre_FLD();

  void init_hierarchy();
  void init_stencil();
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

  // BiCGStab internal solver functions
  void matvec_(int ix, int ib);
  void copy_(int ix, int iy);
  double dot_(int ix, int iy);
  void ddot_(double& dot1, int ix, int iy, 
	     double& dot2, int ip, int iq);
  void linear_combination_(Scalar a, int ix, Scalar b, int iy, 
			   Scalar c, int iz);
  
  // solve() functions
  void solve_bicgstab_(int itmax, double restol);

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
