/*****************************************************************************
 *                                                                           *
 * Copyright 2010 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Group/Frequency, AMR, Flux-Limited Diffusion Solver
/  Split Implicit Problem Class
/
/  written by: Daniel Reynolds
/  date:       December 2010
/
/  PURPOSE: This class defines problem-specific functions for an 
/           implicit gray flux-limited diffusion solve.
/
/           The variables are stored in the following order: 
/              0 -> radiation energy density
/              1 -> fluid energy correction
/              2:Nspecies+1 -> chemical species (Nspecies may be 0)
/
************************************************************************/

#ifdef TRANSFER
#ifndef AMRFLD_SPLIT_PROBLEM_DEFINED__
#define AMRFLD_SPLIT_PROBLEM_DEFINED__

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
/* #include "EnzoVector.h" */
#include "ImplicitProblemABC.h"


// increase the maximum number of FLD sources as needed
#define MAX_FLD_SOURCES       100


// AMRFLDSplit is hard-coded with two source types:
//     0 -> monochromatic at hnu = 13.6 eV
//     1 -> blackbody with temperature 1e5 K
#define NUM_FLD_SOURCE_TYPES  2

// AMRFLDSplit is hard-coded with three boundary condition types:
//     0 -> periodic
//     1 -> dirichlet
//     2 -> neumann
#define NUM_FLD_BDRY_TYPES  3

// AMRFLDSplit is hard-coded with five linear solver types:
//     0 -> FAC
//     1 -> BiCGStab
//     2 -> BiCGStab-BoomerAMG
//     3 -> GMRES
//     4 -> PFMG (unigrid only)
#define NUM_FLD_SOL_TYPES  5

// AMRFLDSplit is hard-coded with four time step adaptivity algorithms:
//     -1 -> original time controller
//      0 -> I controller
//      1 -> PI controller
//      2 -> PID controller
#define NUM_FLD_DT_CONTROLLERS  4

// HYPRE is hard-coded with four relaxation algorithms:
//    0 -> weighted Jacobi
//    1 -> weighted Jacobi
//    2 -> red-black Gauss-Seidel
//    3 -> red-black Gauss-Seidel
#define NUM_HYPRE_RLX_TYPES  4


class AMRFLDSplit : public virtual ImplicitProblemABC {

 private:
  
  // AMRsolve parameter interface structure
#ifdef AMR_SOLVE
  AMRsolve_Parameters* amrsolve_params;
#endif

  // solver-specific parameters
  float  sol_tolerance;          // desired solver tolerance
  int    sol_maxit;              // maximum number of iterations
  int    sol_type;               // HYPRE solver type
  int    sol_printl;             // print output level
  int    sol_log;                // amount of logging
  int    sol_prec;               // enable HG preconditioner within AMRsolve
  int    sol_precmaxit;          // # PFMG iterations within preconditioner
  int    sol_precnpre;           // # pre-relaxation sweeps within precond.
  int    sol_precnpost;          // # post-relaxation sweeps within precond.
  int    sol_precJacit;          // # Jacobi iterations within precond.
  int    sol_precrelax;          // smoother type for PFMG iteration
  float  sol_precrestol;         // tolerance for coarse PFMG solve

  // PFMG-specific solver parameters
  int    sol_rlxtype;            // relaxation type (PFMG only)
  int    sol_npre;               // num. pre-relaxation sweeps
  int    sol_npost;              // num. post-relaxation sweeps

  // limiter-specific parameters
  float LimiterRmin;             // Rmin constant to use in limiter
  float LimiterDmax;             // Dmax constant to use in limiter

  // solver diagnostics
  bool  diags;                   // flag on whether to compute/output diagnostics
  float RTtime;                  // total AMRFLDSplit module time
  float AMRSolTime;              // total time in AMRsolve routines
  int   totIters[MAX_FLD_FIELDS];  // cumulative iterations for solves

  // grid information
  int    rank;                   // Rank of problem
  int    LocDims[3];             // top grid local dims (no ghost or bdry cells)
  bool   OnBdry[3][2];           // denotes if proc owns piece of boundary
  float *BdryVals[MAX_FLD_FIELDS][3][2];   // boundary values for radiation BCs

  // time-stepping related data
  float initdt;        // initial radiation time step size
  float maxdt;         // maximum radiation time step size
  float mindt;         // minimum radiation time step size
  float maxsubcycles;  // max subcycle factor for rad time step within hydro step
  float timeAccuracy;  // desired relative change in radiation fields per step
  int   dt_control;    // time step controller algorithm:
  float Err_cur;       // storage for error estimates in adaptive time stepper
  float Err_new;
  float dtnorm;        // norm choice for computing relative change (default=2.0):
                       //    0 -> max pointwise norm
                       //   >0 -> rms p-norm over entire domain
  float dtgrowth;      // maximum time step growth factor between steps
  float tnew;          // new time
  float told;          // old time
  float dt;            // time step size
  float dtrad;         // radiation time step size (if subcycled)
  float theta;         // implicitness parameter (1->BE, 0.5->CN, 0->FE)
  
  // radiation and chemistry problem-defining data
  int   NumRadiationFields;                 // # radiation frequencies/groups to use
  float FrequencyBand[MAX_FLD_FIELDS][2];   // frequencies or bin boundaries
  bool  FieldMonochromatic[MAX_FLD_FIELDS]; // flag if the field is monochromatic
  bool  FieldNeighbors[MAX_FLD_FIELDS][2];  // flag if field has neighbors (for redshifting): 
                                            // ignored for monochromatic fields,
                                            // groups may have lower/upper neighbors
                                            // if their frequency bands touch
  int   Isothermal;                         // flag denoting temperature-dependence of run

  // ionization source parameters
  int   NumSources;                          // number of input ionization sources
  float SourceLocation[MAX_FLD_SOURCES][3];  // ionization source location
  float SourceGroupEnergy[MAX_FLD_SOURCES][MAX_FLD_FIELDS];   // energy per source/field
  int   WeakScaling;                         // whether to replicate inputs on each process
  float OriginalSourceLocation[MAX_FLD_SOURCES][3];  // necessary to restart WeakScaling runs

  // cosmology and scaling constants
  FLOAT a;             // cosmology expansion coefficient
  FLOAT a0;            // cosmology expansion coefficient (old time)
  FLOAT adot;          // time-derivative of a
  FLOAT adot0;         // time-derivative of a (old time)
  float aUnits;        // expansion parameter scaling
  bool  autoScale[MAX_FLD_FIELDS];   // flag to enable/disable automatic scaling factors
  bool  StartAutoScale[MAX_FLD_FIELDS];  // flag to turn begin automatic scaling in a run
  float ErScale[MAX_FLD_FIELDS];     // radiation energy density scaling factors
  float ErUnits[MAX_FLD_FIELDS];     // radiation energy density unit conversion factor
  float ErUnits0[MAX_FLD_FIELDS];    // radiation energy density unit conversion factor
  float NiUnits;       // species density unit conversion factor
  float NiUnits0;      // species density unit conversion factor
  float DenUnits;      // density scaling factor
  float DenUnits0;     // density scaling factor
  float LenUnits;      // length scaling factor
  float LenUnits0;     // length scaling factor
  float TimeUnits;     // time scaling factor
  float VelUnits;      // velocity scaling factor

  // integrals over frequency space (set during initialization)
  float intOpacity_HI[MAX_FLD_FIELDS];     // 1/|binwidth| * int_{bin} sigmaHI d nu
  float intOpacity_HeI[MAX_FLD_FIELDS];    // 1/|binwidth| * int_{bin} sigmaHeI d nu
  float intOpacity_HeII[MAX_FLD_FIELDS];   // 1/|binwidth| * int_{bin} sigmaHeII d nu
  float intIonizing_HI[MAX_FLD_FIELDS];    // 1/|binwidth| * int_{bin} sigmaHI/nu d nu
  float intIonizing_HeI[MAX_FLD_FIELDS];   // 1/|binwidth| * int_{bin} sigmaHeI/nu d nu
  float intIonizing_HeII[MAX_FLD_FIELDS];  // 1/|binwidth| * int_{bin} sigmaHeII/nu d nu
  float intHeating_HI[MAX_FLD_FIELDS];     // 1/|binwidth| * int_{bin} sigmaHI*(1-nuHI/nu) d nu
  float intHeating_HeI[MAX_FLD_FIELDS];    // 1/|binwidth| * int_{bin} sigmaHeI*(1-nuHeI/nu) d nu
  float intHeating_HeII[MAX_FLD_FIELDS];   // 1/|binwidth| * int_{bin} sigmaHeII*(1-nuHeII/nu) d nu

  // private computation routines
  int EnforceBoundary(int Bin, LevelHierarchyEntry *LevelArray[]);
  int RadiationSource(LevelHierarchyEntry *LevelArray[], int level, float time);
  int Opacity(int Bin, LevelHierarchyEntry *LevelArray[], int level, float time);
  int ComputeRadiationIntegrals();
  int Redshifting(LevelHierarchyEntry *LevelArray[], int level);
  int FillRates(LevelHierarchyEntry *LevelArray[], int level);
#ifdef AMR_SOLVE
  int RadStep(int Bin, LevelHierarchyEntry *LevelArray[], int level, 
	      AMRsolve_Hierarchy *hierarchy, float Etyp, 
	      float Emax, Eflt64 *Echange);
#endif

 public:

  // boundary type in each dimension, face:
  Eint32 BdryType[3][2];

  ///////////////////////////////////////
  // FLD-Specific Routines

  // Constructor
  AMRFLDSplit();
  
  // Destructor
  ~AMRFLDSplit();

  // Problem Initializer
  int Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Problem Evolver
  //  int Evolve(LevelHierarchyEntry *LevelArray[], int level, float deltat);
  int Evolve(LevelHierarchyEntry *LevelArray[], int level, 
	     HierarchyEntry *Grids[], int NumberOfGrids,
	     TopGridData *MetaData, ExternalBoundary *Exterior, 
#ifdef FAST_SIB
	     SiblingGridList SiblingList[],
#endif
	     float deltat);
  
  // Write module parameters to file
  int WriteParameters(FILE *fptr);

  // Problem Boundary Condition setup (called once or at each time step, 
  //    must be called for each locally-owned external face separately)
  int SetupBoundary(int Bin, int Dimension, int Face, int BdryConst, float *BdryData);

  // Return the maximum rad-hydro time step size
  float ComputeTimeStep(Eflt64 CurError);

 private:

  float* AccessRadiationField(int ibin, HierarchyEntry *ThisGrid);
  float* AccessEmissivityField(int ibin, HierarchyEntry *ThisGrid);

};


#endif
#endif
