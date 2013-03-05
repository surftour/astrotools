
/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h" 
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
//#include <gsl/gsl_rng.h>
/* #include "tags.h" */

#define  GADGETVERSION   "2.0"   /*!< code version string */

#define  GENERATIONS     2       /*!< Number of star particles that may be created per gas particle
                                  */
 

#define  TIMEBASE        (1<<28) /*!< The simulated timespan is mapped onto the integer interval [0,TIMEPAN],
                                  *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds
                                  *   to 2^29
                                  */

#define MAXTOPNODES    10000


typedef  long long  peanokey;

#define  BITS_PER_DIMENSION 18	/* for Peano-Hilbert order. Note: Maximum is 10 to fit in 32-bit integer ! */
#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))




#define  GAMMA         (5.0/3)    /*!< adiabatic index of simulated gas */
#define  GAMMA_MINUS1  (GAMMA-1)

#define  HYDROGEN_MASSFRAC 0.76   /*!< mass fraction of hydrogen, relevant only for radiative cooling */

#define  METAL_YIELD       0.02   /*!< effective metal yield for star formation */

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  RNDTABLE 128

/* ... often used physical constants (cgs units) */

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
#define  ELECTRONMASS 9.10953e-28
#define  THOMPSON     6.65245e-25
#define  ELECTRONCHARGE  4.8032e-10
#define  HUBBLE          3.2407789e-18	/* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#ifndef ASMTH
/*! ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#define ASMTH 1.25		
#endif
#ifndef RCUT
/*! RCUT gives the maximum distance (in units of the scale used for the force split) out to which short-range
 * forces are evaluated in the short-range tree walk.
 */
#define RCUT  4.5		
#endif

#ifndef HPM_SMTH
#define HPM_SMTH ASMTH
#endif

#define COND_TIMESTEP_PARAMETER 0.25

#define MAX_NGB  20000		/*!< defines maximum length of neighbour list */

#define MAXLEN_OUTPUTLIST 350	/*!< maxmimum number of entries in output list */

#define DRIFT_TABLE_LENGTH  1000  /*!< length of the lookup table used to hold the drift and kick factors */ 


#define MAXITER 500

#define LINKLENGTH 0.2
#define GROUP_MIN_LEN 32

#define MINRESTFAC 0.05

#ifdef DOUBLEPRECISION
#define FLOAT double
#else
#define FLOAT float
#endif

#ifndef  TWODIMS
#define  NUMDIMS 3                                      /*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470                 /*!< Coefficients for SPH spline kernel and its derivative */ 
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786                 /*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */ 
#else
#define  NUMDIMS 2                                      /*!< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)         /*!< Coefficients for SPH spline kernel and its derivative */ 
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI                           /*!< Coefficient for kernel normalization. */
#endif


#if defined(SFR_METALS) || defined (BLACK_HOLES)
#define PPP P
#ifdef SFR_FEEDBACK
#define EgySNcgs 1e51  
#endif
extern int Flag_phase; 
#else
#define PPP SphP
#endif


#define DMAX(a,b) (dmax1=(a),dmax2=(b),(dmax1>dmax2)?dmax1:dmax2)
#define DMIN(a,b) (dmin1=(a),dmin2=(b),(dmin1<dmin2)?dmin1:dmin2)
#define IMAX(a,b) (imax1=(a),imax2=(b),(imax1>imax2)?imax1:imax2)
#define IMIN(a,b) (imin1=(a),imin2=(b),(imin1<imin2)?imin1:imin2)




/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!<  total particle numbers (global value) */
  long long TotN_gas;		/*!<  total gas particle number (global value) */

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one
				     processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one
				     processor. */

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;      /*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel; /*!< maximum number of files that may be written simultaneously when
                                   writing/reading restart-files, or when writing snapshot files */ 

  int BufferSize;		/*!< size of communication buffer in MB */
  int BunchSizeForce;		/*!< number of particles fitting into the buffer in the parallel tree-force
				   algorithm  */
  int BunchSizeDensity;         /*!< number of particles fitting into the communication buffer in the density
                                  computation */
  int BunchSizeFoF;

  int BunchSizeMetal;
  int BunchSizeHotNgbs;

  int BunchSizeBlackhole;

  int BunchSizeKinetic;

  int BunchSizeHydro;           /*!< number of particles fitting into the communication buffer in the SPH
                                  hydrodynamical force computation */
  int BunchSizeDomain;          /*!< number of particles fitting into the communication buffer in the domain
                                  decomposition */

  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				   NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				   the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  /* some SPH parameters */

  int DesNumNgb;                /*!< Desired number of SPH neighbours */
  int MaxNumNgbDeviation;       /*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double InitGasU;		/*!< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */



  /* some force counters  */

  long long TotNumOfForces;	/*!< counts total number of force computations  */

  long long NumForcesSinceLastDomainDecomp;  /*!< count particle updates since last domain decomposition */

  /* system of units  */

  double UnitTime_in_s,   	/*!< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,             	/*!< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,   /*!< factor to convert intqernal velocity unit to cm/sec */
    UnitLength_in_cm,           /*!< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,        /*!< factor to convert internal pressure unit to cgs units (little 'h' still
                                     around!) */
    UnitDensity_in_cgs,         /*!< factor to convert internal length unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,     /*!< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,          /*!< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,      /*!< factor to convert internal time to megayears/h */
    GravityConstantInternal,    /*!< If set to zero in the parameterfile, the internal value of the
                                  gravitational constant is set to the Newtonian value based on the system of
                                  units specified. Otherwise the value provided is taken as internal gravity
                                  constant G. */
    G;                          /*!< Gravity-constant in internal units */

  /* Cosmology */

  double Hubble;  /*!< Hubble-constant in internal units */
  double Omega0,  /*!< matter density in units of the critical density (at z=0)*/
    OmegaLambda,  /*!< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,  /*!< baryon density in units of the critical density (at z=0)*/
    HubbleParam;  /*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
		   * physical values for cooling physics
                   */

  double BoxSize; /*!< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;     /*!< flags that periodic boundaries are enabled */
  int ResubmitOn;               /*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;   /*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
                                  criterion */
  int TypeOfTimestepCriterion;  /*!< gives type of timestep criterion (only 0 supported right now - unlike
                                  gadget-1.1) */
  int OutputListOn;             /*!< flags that output times are listed in a specified file */
  int CoolingOn;                /*!< flags that cooling is enabled */
  int StarformationOn;          /*!< flags that star formation is enabled */


  /* parameters determining output frequency */

  int SnapshotFileCount;     /*!< number of snapshot that is written next */
  double TimeBetSnapshot,    /*!< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,     /*!< simulation time of first snapshot files */
    CpuTimeBetRestartFile,   /*!< cpu-time between regularly generated restart files */
    TimeLastRestartFile,     /*!< cpu-time when last restart-file was written */
    TimeBetStatistics,       /*!< simulation time interval between computations of energy statistics */
    TimeLastStatistics;      /*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;      /*!< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,  /*!< current time of the simulation */
    TimeBegin,  /*!< time of initial conditions of the simulation */
    TimeStep,   /*!< difference between current times of previous and current timestep */
    TimeMax;	/*!< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval; /*!< factor to convert from floating point time interval to integer timeline */
  int Ti_Current;           /*!< current time on integer timeline */ 
  int Ti_nextoutput;        /*!< next output time on integer timeline */

   int PresentMinStep;


  int PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], UpperCorner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];


  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_TreeConstruction;
  double CPU_TreeWalk;
  double CPU_Gravity;
  double CPU_Potential;
  double CPU_Domain;
  double CPU_Snapshot;
  double CPU_Total;
  double CPU_CommSum;
  double CPU_Imbalance;
  double CPU_HydCompWalk;
  double CPU_HydCommSumm;
  double CPU_HydImbalance;
  double CPU_Hydro;
  double CPU_EnsureNgb;
  double CPU_Predict;
  double CPU_TimeLine;
  double CPU_PM;
  double CPU_Peano;
  double CPU_SfrCool;

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
                                  timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep,       /*!< minimum allowed timestep. Normally, the simulation terminates if the
                                  timestep determined by the timestep criteria falls below this limit. */ 
         MaxSizeTimestep;       /*!< maximum allowed timestep */

  double MaxRMSDisplacementFac; /*!< this determines a global timestep criterion for cosmological simulations
                                     in comoving coordinates.  To this end, the code computes the rms velocity
                                     of all particles, and limits the timestep such that the rms displacement
                                     is a fraction of the mean particle separation (determined from the
                                     particle mass and the cosmological parameters). This parameter specifies
                                     this fraction. */



  double CourantFac;		/*!< SPH-Courant factor */


  /* frequency of tree reconstruction/domain decomposition */


  double TreeDomainUpdateFrequency; /*!< controls frequency of domain decompositions  */


  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   * five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars
   */
  double MinGasHsmlFractional, /*!< minimum allowed SPH smoothing length in units of SPH gravitational
                                  softening length */
    MinGasHsml;                /*!< minimum allowed SPH smoothing length */


  double SofteningGas,    /*!< for type 0 */ 
    SofteningHalo,        /*!< for type 1 */ 
    SofteningDisk,        /*!< for type 2 */ 
    SofteningBulge,       /*!< for type 3 */ 
    SofteningStars,       /*!< for type 4 */ 
    SofteningBndry;       /*!< for type 5 */ 

  double SofteningGasMaxPhys,   /*!< for type 0 */ 
    SofteningHaloMaxPhys,       /*!< for type 1 */ 
    SofteningDiskMaxPhys,       /*!< for type 2 */ 
    SofteningBulgeMaxPhys,      /*!< for type 3 */ 
    SofteningStarsMaxPhys,      /*!< for type 4 */ 
    SofteningBndryMaxPhys;      /*!< for type 5 */ 

  double SofteningTable[6];  /*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];  /*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[6];


  /* some filenames */
  char InitCondFile[100],
    OutputDir[100],
    SnapshotFileBase[100],
    EnergyFile[100],
    CpuFile[100],
    InfoFile[100], 
    TimingsFile[100], 
    RestartFile[100], 
    ResubmitCommand[100], 
    OutputListFilename[100];

  /*! table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];

  int OutputListLength; /*!< number of times stored in table of desired output times */



  double OrigGasMass;
  double EgySpecCold;
  double EgySpecSN;
  double OverDensThresh;
  double PhysDensThresh;
  double FeedbackEnergy;
  double TempSupernova;
  double TempClouds;
  double CritOverDensity;
  double CritPhysDensity;
  double FactorSN;
  double FactorEVP;
  double MaxSfrTimescale;
  double WindEfficiency;
  double WindEnergyFraction;
  double WindFreeTravelLength;
  double WindFreeTravelDensFac;
  double FactorForSofterEQS;

  double FactorSFR;
  double MinTlifeSNI;
  double MaxTlifeSNI;
  double TlifeSNII;
  double RateSNI;
  double FactorSN_Phase;
  double Tcrit_Phase;
  double DensFrac_Phase;
  double FracEnergySN_Phase;
  double DensityTailThreshold; 

  double DarkEnergyParam;	/*!< fixed w for equation of state */
  char DarkEnergyFile[100];	/*!< tabelized w for equation of state */

  double VelIniScale;		/*!< Scale the initial velocities by this amount */

  double ViscSource0;		/*!< Given sourceterm in viscosity evolution */
  double DecayLength;		/*!< Number of h for the viscosity decay */
  double ViscSource;		/*!< Reduced sourceterm in viscosity evolution*/
  double DecayTime;		/*!< Calculated decaytimescale */
  double AlphaMin;		/*!< Minimum of allowed viscosity parameter */

  double ConductionCoeff;         /*!< Thermal Conductivity */
  double ElectronFreePathFactor;  /*!< Factor to get electron mean free path */

  double BiniX,BiniY,BiniZ;       /*!< Initial values for B */

  int BSmoothInt;
  double BSmoothFrac;
  int MainTimestepCounts;

  double TimeNextBlackHoleCheck;
  double TimeBetBlackHoleSearch;
  double BlackHoleAccretionFactor;  /*!< Fraction of BH bondi accretion rate */
  double BlackHoleFeedbackFactor;   /*!< Fraction of the black luminosity feed into thermal feedback */
  double SeedBlackHoleMass;         /*!< Seed black hole mass */
  double MinFoFMassForNewSeed;      /*!< Halo mass required before new seed is put in */
  double BlackHoleNgbFactor;        /*!< Factor by which the normal SPH neighbour should be increased/decreased */
  double BlackHoleActiveTime;
  double BlackHoleRefDensity;
  double BlackHoleRefSoundspeed;

  double CR_Alpha;               /*!< Cosmic ray spectral index [2..3]*/
  double CR_SNEff;               /*!< SN injection efficiency [0..1] */
  double CR_SNAlpha;             /*!< SN injection spectral index [2..3] */
  int bDebugFlag;                /*!< enables debug outputs after triggered */

  double CR_Diffusion_Gamma;     /*!< CR diffusion coefficient d_gamma */
  double CR_Diffusion_Proton;    /*!< CR diffusion coefficient d_p     */

  double CR_DiffusionCoeff;      /*!< (temporary) fixed value for CR diffusivity */

  double CR_ShockAlpha;          /*!< spectral index to be used in shock injection */
  double CR_ShockEfficiency;     /*!< energy fraction of shock energy fed into CR */


  double HPM_entr0, HPM_entr1;
  double HPM_ne0, HPM_ne1;
  double HPM_rho0, HPM_rho1;
  double HPM_P0, HPM_P1, HPM_alpha;


  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
  int    BiggestGroupLen;
  float  BiggestGroupCM[3];
  double BiggestGroupMass;   
}
All;

#endif
