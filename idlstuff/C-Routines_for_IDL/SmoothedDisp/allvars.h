
#define  DEBUG

#define  MADAU 0             /* the two options mutually exclude each */
#define  ADIABATIC 1         /* other !!! Either Madau or adiabatic!  */



/* ... often used constants (cgs units) */

#define  MAX_REAL_NUMBER  1e200
#define  MIN_REAL_NUMBER  1e-200

#define  THIRD            (1.0/3.0)
#define  PI               3.1415927
#define  PI_INV           (1/PI)
#define  LN2              0.69314718

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define  GAMMA         (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)







/* define array sizes */



#define MAX_TIMELEVELS  28

#define MAX_GRCHIPS  16 
#define MAX_GRNBR    1024                 

#define MAX_NGB      500

#define KERNEL_TABLE 10000

#define MASSBINS     10
#define ZONES        6


/* system of units  */

extern double UnitTime_in_s,
              UnitMass_in_g,
              UnitLength_in_cm,
              UnitPressure_in_cgs,
              UnitDensity_in_cgs,
              UnitCoolingRate_in_cgs,
              UnitEnergy_in_cgs,
              G,
              FOUR_PI_G; 



/* variables for input/output */

extern char   ParameterFile[100],
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

extern int    SnapshotFileCount;
extern double TimeLastSnapshot;
extern double TimeBetRestartFile,
              TimeLastRestartFile;

extern FILE   *FdInfo,
              *FdEnergy,
              *FdNgb,
              *FdTimeBin,
              *FdTimeCrit,
              *FdStars,
              *FdSfr;

extern int    RestartFlag;



/* Current time of the simulation */

extern int     NumCurrentTiStep;

extern double  Time, 
               TimeMax;





/* cosmological quantities */

extern double Z_Current,    
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

extern double TimeBetSnapshot;



/* parameters of SPH algorithm */

extern double CourantFac,
              ErrTolForce,
              ErrTolDensity,
              ErrTolEnergy,
              ErrTolTheta,
              ParHsmlDet,
              ArtBulkViscConst,
              MassFrac_H,
              MassFrac_He;

extern int    DesNumNgb;

extern int    CoolingOn,
              StarFormationOn,
              FeedbackOn;

extern double FactorSFR,
              StarMass,
              ThermalizeTimescale,
              FeedbackEnergy;


/* Cut-off values */

extern double MinTemperature,
              MinEnergySpecific;



/* data for massbins & zones */

extern int     BegMassZone[MASSBINS+1][ZONES],
               EndMassZone[MASSBINS+1][ZONES];

extern int     *PartZone,
               *IndOrdMassZone;

extern int     MinMassBin,
               MaxMassBin,
               MinZone,
               MaxZone;
             
extern double  MassOfBin[MASSBINS+1],
               SofteningOfBin[MASSBINS+1];


/* Length scales and softening */

extern double  MinLengthScale,
               MaxLengthScale,
               ZoneFactor,
               SofteningGas,
               SofteningHalo,
               SofteningDisk,
               SofteningBulge,
               SofteningStars,
               MaxSoftening;


/* GRAPE hardware specifications */

extern int     GrNumChips, 
               GrNumChipsPerBoard,
               GrNumBoards,
               GrMaxPart,
               GrMaxNbr;

extern int     UseGrape,              /* If set, GRAPE board is used, otherwise */
               GrapeAquired,          /* Tree-Code */
               FreeGrape;



/* tabulated smoothing kernel */

extern double  Kernel[KERNEL_TABLE+1],
               KernelDer[KERNEL_TABLE+1],
               KernelRad[KERNEL_TABLE+1];


/* state of total system */

extern struct state_of_system 
{
  double  Mass,
          EnergyKin,
          EnergyPot,
          EnergyInt,
          EnergyTot,
          Momentum[4],
          AngMomentum[4],
          CenterOfMass[4],

          MassComp[2],
          EnergyKinComp[2],
          EnergyPotComp[2],
          EnergyIntComp[2],
          EnergyTotComp[2],
          MomentumComp[2][4],
          AngMomentumComp[2][4],
          CenterOfMassComp[2][4];

}   SysState,SysStateAtStart,SysStateAtEnd;




/* Quantities for hierachy of timesteps */

extern double  Min_dt,
               MinSizeTimestep,
               MaxSizeTimestep;

extern double  Lev_dt[MAX_TIMELEVELS+1];

extern int     LevelFac[MAX_TIMELEVELS+1],
               LevelSync[MAX_TIMELEVELS+1], 
               TimeLevels, 
               CountTiLevFac, 
               MinOccTimeLevel, 
               MaxOccTimeLevel;

               

/*  particle numbers */

extern int     NumPart, 
               NumPartBary, 
               N_gas,
               N_halo,
               N_disk,
               N_bulge,
               N_stars;

/* Quantities for all particles */

extern struct particle_data 
{
  double Pos[3], 
         Mass, 
         Accel[3],
         Potential;
  int    MassBin;
  double Vel[3],
         VelPred[3]; 

  int    TimeLevel,
         TimeLevelCause;

} **P;


/* Quantities for SPH particles only */

extern struct sph_particle_data
{
  double Density,
         DtDensity,
         Pressure,
         DivVel, 
         CurlVel, 
         Hsml,  
         DtHsml,
         EgySpec, 
         EgySpecPred,
         DtEgySpec,

         EgySpecTp,
         EgySpecPredTp,
         DtEgySpecTp,

         MaxUij,
         FormedStellarMass,
         RateOfSF,
         AvailableStellarMass;


  int    NumNgb,
         NgbList[MAX_NGB];



} **SphP;










