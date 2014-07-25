/***********************************************************************
/
/  INITIALIZE RADIATION-HYDRODYNAMICS TEST -- IONIZATION TEST IN A 
/  R^{-2} DENSITY PROFILE
/
/  written by: Daniel Reynolds
/  date:       December 2007
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"


/* default constants */
#define DEFAULT_MU 0.6           // mean molecular mass
#define MIN_TEMP 1.0             // minimum temperature [K]
#define MAX_INITIAL_PATCHES 100  // max number of parameter file defined subgrids


// function prototypes
int InitializeRateData(FLOAT Time);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);



int RHIonizationSteepInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local)
{
#ifdef TRANSFER
//   if (MyProcessorNumber == ROOT_PROCESSOR)
//     fprintf(stdout,"Entering RHIonizationSteepInitialize routine\n");

  char *kphHIName    = "HI_kph";
  char *kphHeIName   = "HeI_kph";
  char *kphHeIIName  = "HeII_kph";
  char *gammaName    = "PhotoGamma";
  char *kdissH2IName = "H2I_kdiss";
  char *DensName  = "Density";
  char *TEName    = "TotalEnergy";
  char *IEName    = "GasEnergy";
  char *Vel0Name  = "x-velocity";
  char *Vel1Name  = "y-velocity";
  char *Vel2Name  = "z-velocity";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *DeName    = "Electron_Density";
  char *RadName   = "Grey_Radiation_Energy";
  char *RadName0  = "Radiation0";
  char *RadName1  = "Radiation1";
  char *RadName2  = "Radiation2";
  char *RadName3  = "Radiation3";
  char *RadName4  = "Radiation4";
  char *RadName5  = "Radiation5";
  char *RadName6  = "Radiation6";
  char *RadName7  = "Radiation7";
  char *RadName8  = "Radiation8";
  char *RadName9  = "Radiation9";
  char *EtaName   = "Emissivity";
  char *EtaName0  = "Emissivity0";
  char *EtaName1  = "Emissivity1";
  char *EtaName2  = "Emissivity2";
  char *EtaName3  = "Emissivity3";
  char *EtaName4  = "Emissivity4";
  char *EtaName5  = "Emissivity5";
  char *EtaName6  = "Emissivity6";
  char *EtaName7  = "Emissivity7";
  char *EtaName8  = "Emissivity8";
  char *EtaName9  = "Emissivity9";

  // local declarations
  char line[MAX_LINE_LENGTH];
  int  i, j, dim, gridnum, ret, level, patch;

  // Setup and parameters:
  float RadHydroX0Velocity           = 0.0;
  float RadHydroX1Velocity           = 0.0;
  float RadHydroX2Velocity           = 0.0;
  float RadHydroNumDensity           = 3.2;           // [cm^{-3}]
  float RadHydroDensityRadius        = 1.14375e-1;    // [code units]
  float DensityCenter0               = 0.0;
  float DensityCenter1               = 0.0;
  float DensityCenter2               = 0.0;
  float RadHydroTemperature          = 100.0;         // [K]
  float RadHydroRadiationEnergy      = 1.0e-20;
  float RadHydroHydrogenMassFraction = 1.0;
  float RadHydroInitialFractionHII   = 0.0;
  float RadHydroInitialFractionHeII  = 0.0;
  float RadHydroInitialFractionHeIII = 0.0;
  int   RadHydroChemistry            = 1;
  int   AMRFLDNumRadiationFields     = 0;    // grey solver

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	ret = 0;
	// read relevant problem parameters
	ret += sscanf(line, "RadHydroVelocity = %"FSYM" %"FSYM" %"FSYM,
		      &RadHydroX0Velocity, &RadHydroX1Velocity, 
		      &RadHydroX2Velocity);
	ret += sscanf(line, "RadHydroChemistry = %"ISYM, 
		      &RadHydroChemistry);
	ret += sscanf(line, "AMRFLDNumRadiationFields = %"ISYM, 
		      &AMRFLDNumRadiationFields);
	ret += sscanf(line, "RadHydroNumDensity = %"FSYM, 
		      &RadHydroNumDensity);
	ret += sscanf(line, "RadHydroDensityRadius = %"FSYM, 
		      &RadHydroDensityRadius);
	ret += sscanf(line, "RadHydroTemperature = %"FSYM, 
		      &RadHydroTemperature);
	ret += sscanf(line, "RadHydroRadiationEnergy = %"FSYM, 
		      &RadHydroRadiationEnergy);
	ret += sscanf(line, "RadHydroInitialFractionHII = %"FSYM, 
		      &RadHydroInitialFractionHII);
	ret += sscanf(line, "RadHydroHFraction = %"FSYM, 
		      &RadHydroHydrogenMassFraction);
	ret += sscanf(line, "RadHydroInitialFractionHeII = %"FSYM, 
		      &RadHydroInitialFractionHeII);
	ret += sscanf(line, "RadHydroInitialFractionHeIII = %"FSYM, 
		      &RadHydroInitialFractionHeIII);
	ret += sscanf(line, "EtaCenter = %"FSYM" %"FSYM" %"FSYM, 
		      &DensityCenter0, &DensityCenter1, &DensityCenter2);

      } // end input from parameter file
      fclose(RHfptr);
    }
  }

  /* error checking */
  if (Mu != DEFAULT_MU) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: mu =%f assumed in initialization; setting Mu = %f for consistency.\n", DEFAULT_MU);
    Mu = DEFAULT_MU;
  }

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }

  // since AMRFLD solver no longer relies on RadHydroChemistry input, deduce the value here
  if (ImplicitProblem == 6) 
    RadHydroChemistry = (RadiativeTransferHydrogenOnly) ? 1 : 3;

  // convert input temperature to internal energy
  RadHydroTemperature = max(RadHydroTemperature,MIN_TEMP); // enforce minimum
  float mp = 1.67262171e-24;    // proton mass [g]
  float kb = 1.3806504e-16;     // boltzmann constant [erg/K]
  float nH, HI, HII, nHe, HeI, HeII, HeIII, ne, num_dens, mu;
  if (RadHydroChemistry == 1) {
    HI = 1.0 - RadHydroInitialFractionHII;
    HII = RadHydroInitialFractionHII;
    ne = HII;
    num_dens = HI + HII + ne;
    mu = 1.0/num_dens;
  }
  else if (RadHydroChemistry == 3) {
    nH = RadHydroHydrogenMassFraction;
    nHe = (1.0 - RadHydroHydrogenMassFraction);
    HI = nH*(1.0 - RadHydroInitialFractionHII);
    HII = nH*RadHydroInitialFractionHII;
    HeII = nHe*RadHydroInitialFractionHeII;
    HeIII = nHe*RadHydroInitialFractionHeIII;
    HeI = nHe - HeII - HeIII;
    ne = HII + HeII/4.0 + HeIII/2.0;
    num_dens = 0.25*(HeI + HeII + HeIII) + HI + HII + ne;
    mu = 1.0/num_dens;
  }
  // compute the internal energy
  float RadHydroIEnergy = kb*RadHydroTemperature/mu/mp/(Gamma-1.0);	

  /////////////////
  // Set up the TopGrid as usual
  HierarchyEntry *TempGrid = &TopGrid;
  while (TempGrid != NULL) {
    if (TempGrid->GridData->RHIonizationSteepInitializeGrid(
                        RadHydroChemistry, AMRFLDNumRadiationFields, RadHydroNumDensity, 
			RadHydroDensityRadius, DensityCenter0, 
			DensityCenter1, DensityCenter2, RadHydroX0Velocity, 
			RadHydroX1Velocity, RadHydroX2Velocity, 
			RadHydroIEnergy, RadHydroRadiationEnergy, 
			RadHydroHydrogenMassFraction, 
			RadHydroInitialFractionHII, 
			RadHydroInitialFractionHeII, 
			RadHydroInitialFractionHeIII, local) == FAIL) {
      fprintf(stderr, "Error in RHIonizationSteepInitializeGrid.\n");
      return FAIL;
    }
    TempGrid = TempGrid->NextGridThisLevel;
  }


  // set up field names and units
  // note: we must set up He species fields as well since Enzo 
  //       requires them for H chemistry (initialized to zero)
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  if (DualEnergyFormalism) 
    DataLabel[BaryonField++] = IEName;
  DataLabel[BaryonField++] = Vel0Name;
  DataLabel[BaryonField++] = Vel1Name;
  DataLabel[BaryonField++] = Vel2Name;
  if (AMRFLDNumRadiationFields == 0)
    DataLabel[BaryonField++] = RadName;
  if (AMRFLDNumRadiationFields > 0)
    DataLabel[BaryonField++] = RadName0;
  if (AMRFLDNumRadiationFields > 1)
    DataLabel[BaryonField++] = RadName1;
  if (AMRFLDNumRadiationFields > 2)
    DataLabel[BaryonField++] = RadName2;
  if (AMRFLDNumRadiationFields > 3)
    DataLabel[BaryonField++] = RadName3;
  if (AMRFLDNumRadiationFields > 4)
    DataLabel[BaryonField++] = RadName4;
  if (AMRFLDNumRadiationFields > 5)
    DataLabel[BaryonField++] = RadName5;
  if (AMRFLDNumRadiationFields > 6)
    DataLabel[BaryonField++] = RadName6;
  if (AMRFLDNumRadiationFields > 7)
    DataLabel[BaryonField++] = RadName7;
  if (AMRFLDNumRadiationFields > 8)
    DataLabel[BaryonField++] = RadName8;
  if (AMRFLDNumRadiationFields > 9)
    DataLabel[BaryonField++] = RadName9;
  DataLabel[BaryonField++] = DeName;
  DataLabel[BaryonField++] = HIName;
  DataLabel[BaryonField++] = HIIName;
  if ((RadHydroChemistry == 3) || (MultiSpecies > 0)) {
    DataLabel[BaryonField++] = HeIName;
    DataLabel[BaryonField++] = HeIIName;
    DataLabel[BaryonField++] = HeIIIName;
  }

  // if using external chemistry/cooling, set rate labels and update params
  if (RadiativeCooling) {
    DataLabel[BaryonField++] = kphHIName;
    DataLabel[BaryonField++] = gammaName;
    if (RadiativeTransferHydrogenOnly == FALSE) {
      DataLabel[BaryonField++] = kphHeIName;
      DataLabel[BaryonField++] = kphHeIIName;
    }
    if (MultiSpecies > 1)
      DataLabel[BaryonField++] = kdissH2IName;
  }

  // if using the AMRFLDSplit solver, set fields for the emissivity
  if (ImplicitProblem == 6) {
    if (AMRFLDNumRadiationFields > 0)
      DataLabel[BaryonField++] = EtaName0;
    if (AMRFLDNumRadiationFields > 1)
      DataLabel[BaryonField++] = EtaName1;
    if (AMRFLDNumRadiationFields > 2)
      DataLabel[BaryonField++] = EtaName2;
    if (AMRFLDNumRadiationFields > 3)
      DataLabel[BaryonField++] = EtaName3;
    if (AMRFLDNumRadiationFields > 4)
      DataLabel[BaryonField++] = EtaName4;
    if (AMRFLDNumRadiationFields > 5)
      DataLabel[BaryonField++] = EtaName5;
    if (AMRFLDNumRadiationFields > 6)
      DataLabel[BaryonField++] = EtaName6;
    if (AMRFLDNumRadiationFields > 7)
      DataLabel[BaryonField++] = EtaName7;
    if (AMRFLDNumRadiationFields > 8)
      DataLabel[BaryonField++] = EtaName8;
    if (AMRFLDNumRadiationFields > 9)
      DataLabel[BaryonField++] = EtaName9;
  }

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;
 
#endif

}
