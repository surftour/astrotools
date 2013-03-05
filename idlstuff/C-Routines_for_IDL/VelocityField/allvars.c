#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


/* system of units  */

 double UnitTime_in_s,
              UnitMass_in_g,
              UnitLength_in_cm,
              UnitPressure_in_cgs,
              UnitDensity_in_cgs,
              UnitCoolingRate_in_cgs,
              UnitEnergy_in_cgs,
              G,
              FOUR_PI_G; 



/* variables for input/output */

 char   ParameterFile[100],
              InitCondFile[100],
              OutputDir[100],
              SnapshotFileBase[100],
              EnergyFile[100],
              NgbFile[100],
              InfoFile[100],
              TimeBinFile[100],
              TimeCritFile[100],
              StarsFile[100],
              SfrFile[100],
              RestartFile[100];

   int    SnapshotFileCount;
   double TimeLastSnapshot;
   double TimeBetRestartFile,
          TimeLastRestartFile;

       FILE   *FdInfo,
              *FdEnergy,
              *FdNgb,
              *FdTimeBin,
              *FdTimeCrit,
              *FdStars,
              *FdSfr;

 int    RestartFlag;



/* Current time of the simulation */

 int     NumCurrentTiStep;

 double  Time, 
               TimeMax;





/* cosmological quantities */

 double Z_Current,    
              Z_AtEnd, 
              Z_AtStart,
              A_Scale,
              
              OmegaMatter,
              OmegaBary,
              OmegaDark,
              OmegaLambda,
              RhoUniverse,
              RhoBary,
              RhoDark,
              RhoCriticalToday,
              Hubble_h,
              HubbleToday;

/* parameters determining output frequency */

 double TimeBetSnapshot;



/* parameters of SPH algorithm */

 double CourantFac,
              ErrTolForce,
              ErrTolDensity,
              ErrTolEnergy,
              ErrTolTheta,
              ParHsmlDet,
              ArtBulkViscConst,
              MassFrac_H,
              MassFrac_He;

 int    DesNumNgb;

 int    CoolingOn,
              StarFormationOn,
              FeedbackOn;

 double FactorSFR,
              StarMass,
              ThermalizeTimescale,
              FeedbackEnergy;


/* Cut-off values */

 double MinTemperature,
              MinEnergySpecific;



/* data for massbins & zones */

 int     BegMassZone[MASSBINS+1][ZONES],
               EndMassZone[MASSBINS+1][ZONES];

 int     *PartZone,
               *IndOrdMassZone;

 int     MinMassBin,
               MaxMassBin,
               MinZone,
               MaxZone;
             
 double  MassOfBin[MASSBINS+1],
               SofteningOfBin[MASSBINS+1];


/* Length scales and softening */

 double  MinLengthScale,
               MaxLengthScale,
               ZoneFactor,
               SofteningGas,
               SofteningHalo,
               SofteningDisk,
               SofteningBulge,
               SofteningStars,
               MaxSoftening;


/* GRAPE hardware specifications */

 int     GrNumChips, 
               GrNumChipsPerBoard,
               GrNumBoards,
               GrMaxPart,
               GrMaxNbr;

 int     UseGrape,              /* If set, GRAPE board is used, otherwise */
               GrapeAquired,          /* Tree-Code */
               FreeGrape;



/* tabulated smoothing kernel */

 double  Kernel[KERNEL_TABLE+1],
               KernelDer[KERNEL_TABLE+1],
               KernelRad[KERNEL_TABLE+1];


/* state of total system */

 struct state_of_system 
   SysState,SysStateAtStart,SysStateAtEnd;




/* Quantities for hierachy of timesteps */

 double  Min_dt,
               MinSizeTimestep,
               MaxSizeTimestep;

 double  Lev_dt[MAX_TIMELEVELS+1];

 int     LevelFac[MAX_TIMELEVELS+1],
               LevelSync[MAX_TIMELEVELS+1], 
               TimeLevels, 
               CountTiLevFac, 
               MinOccTimeLevel, 
               MaxOccTimeLevel;

               

/*  particle numbers */

 int     NumPart, 
               NumPartBary, 
               N_gas,
               N_halo,
               N_disk,
               N_bulge,
               N_stars;

/* Quantities for all particles */

 struct particle_data 
 **P;


/* Quantities for SPH particles only */

 struct sph_particle_data
 **SphP;










