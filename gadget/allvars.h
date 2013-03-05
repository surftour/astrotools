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

#include <mpi.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#ifdef MPISENDRECV_CHECKSUM
#define MPI_Sendrecv MPI_Check_Sendrecv
#endif

#ifdef MPISENDRECV_SIZELIMIT
#define MPI_Sendrecv MPI_Sizelimited_Sendrecv
#endif

#include "tags.h"
#ifdef CHEMISTRY
#include "chemistry.h"
#endif
#include "assert.h"


#define  GADGETVERSION   "3.0"	/*!< code version string */

#ifdef GASRETURN
#define  GENERATIONS     1	/*!< Number of star particles that may be created per gas particle */
#else
#define  GENERATIONS     2	/*!< Number of star particles that may be created per gas particle */
#endif


#define  TIMEBINS        29

#define  TIMEBASE        (1<<TIMEBINS)	/*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
					 *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds
					 *   to 2^29
					 */
#ifndef  MULTIPLEDOMAINS
#define  MULTIPLEDOMAINS     1
#endif

#ifndef  TOPNODEFACTOR
#define  TOPNODEFACTOR       2.5
#endif

#define  NODELISTLENGTH      8

typedef unsigned long long peanokey;


#define  BITS_PER_DIMENSION 21	/* for Peano-Hilbert order. Note: Maximum is 10 to fit in 32-bit integer ! */
#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))



#ifndef  GAMMA
#define  GAMMA         (5.0/3)	/*!< adiabatic index of simulated gas */
#endif

#define  GAMMA_MINUS1  (GAMMA-1)

#define  HYDROGEN_MASSFRAC 0.76	/*!< mass fraction of hydrogen, relevant only for radiative cooling */

#define  METAL_YIELD       0.02	/*!< effective metal yield for star formation */

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  RNDTABLE 1048576

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
#define  LYMAN_ALPHA      1215.6e-8	/* 1215.6 Angstroem */
#define  LYMAN_ALPHA_HeII  303.8e-8	/* 303.8 Angstroem */
#define  OSCILLATOR_STRENGTH       0.41615
#define  OSCILLATOR_STRENGTH_HeII  0.41615

#ifdef NAVIERSTOKES
#define  LOG_LAMBDA      37.8	/* logarithmic Coulomb factor */
#endif

#ifdef CHEMISTRY
#define  T_CMB0      2.728	/* present-day CMB temperature */
#endif

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

/*Determines the maximum size of arrays related to the number of CR populations */ 
#ifndef NUMCRPOP   /*!< Number of CR populations pressent in parameter file */
#define NUMCRPOP 1 
#endif



#ifndef FOF_PRIMARY_LINK_TYPES
#define FOF_PRIMARY_LINK_TYPES 2
#endif

#ifndef FOF_SECONDARY_LINK_TYPES
#define FOF_SECONDARY_LINK_TYPES 0
#endif


/* some flags for the field "flag_ic_info" in the file header */
#define FLAG_ZELDOVICH_ICS     1
#define FLAG_SECOND_ORDER_ICS  2
#define FLAG_EVOLVED_ZELDOVICH 3 
#define FLAG_EVOLVED_2LPT      4
#define FLAG_NORMALICS_2LPT    5


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

#define COND_TIMESTEP_PARAMETER 0.25
#define VISC_TIMESTEP_PARAMETER 0.25

#define MAXLEN_OUTPUTLIST 350	/*!< maxmimum number of entries in output list */

#define DRIFT_TABLE_LENGTH  1000	/*!< length of the lookup table used to hold the drift and kick factors */


#define MAXITER 150

#ifndef LINKLENGTH
#define LINKLENGTH 0.2
#endif

#ifndef FOF_GROUP_MIN_LEN
#define FOF_GROUP_MIN_LEN 32
#endif  

#define MINRESTFAC 0.05


#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif


#ifndef DOUBLEPRECISION     /* default is single-precision */
typedef float  MyFloat;
typedef float  MyDouble;
#else
#if (DOUBLEPRECISION == 2)   /* mixed precision */
typedef float   MyFloat;
typedef double  MyDouble;
#else                        /* everything double-precision */
typedef double  MyFloat;
typedef double  MyDouble;
#endif
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
#else
typedef float MyOutputFloat;
#endif

#ifdef INPUT_IN_DOUBLEPRECISION
typedef double MyInputFloat;
#else
typedef float MyInputFloat;
#endif

struct unbind_data
{
  int index;
};


#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
extern MPI_Status mpistat;
#undef MPI_STATUS_IGNORE
#define MPI_STATUS_IGNORE &mpistat
#endif

#ifdef FLTROUNDOFFREDUCTION  
#define FLT(x) ((MyFloat)(x))
#ifdef SOFTDOUBLEDOUBLE      /* this requires a C++ compilation */
#include "dd.h"
typedef dd MyLongDouble;
#else
typedef long double MyLongDouble;
#endif
#else  /* not enabled */
#define FLT(x) (x)
typedef MyFloat MyLongDouble;
#endif  /* end FLTROUNDOFFREDUCTION */


#define CPU_ALL            0
#define CPU_TREEWALK1      1
#define CPU_TREEWALK2      2
#define CPU_TREEWAIT1      3
#define CPU_TREEWAIT2      4
#define CPU_TREESEND       5
#define CPU_TREERECV       6
#define CPU_TREEMISC       7
#define CPU_TREEBUILD      8
#define CPU_TREEUPDATE     9
#define CPU_TREEHMAXUPDATE 10
#define CPU_DOMAIN         11
#define CPU_DENSCOMPUTE    12
#define CPU_DENSWAIT       13
#define CPU_DENSCOMM       14
#define CPU_DENSMISC       15
#define CPU_HYDCOMPUTE     16
#define CPU_HYDWAIT        17
#define CPU_HYDCOMM        18
#define CPU_HYDMISC        19
#define CPU_DRIFT          20
#define CPU_TIMELINE       21
#define CPU_POTENTIAL      22
#define CPU_MESH           23
#define CPU_PEANO          24
#define CPU_COOLINGSFR     25        
#define CPU_SNAPSHOT       26
#define CPU_FOF            27
#define CPU_BLACKHOLES     28
#define CPU_MISC           29
#define CPU_SMTHCOMPUTE    30
#define CPU_SMTHWAIT       31
#define CPU_SMTHCOMM       32
#define CPU_SMTHMISC       33
#define CPU_HOTNGBS        34
#define CPU_WEIGHTS_HOT    35
#define CPU_ENRICH_HOT     36
#define CPU_WEIGHTS_COLD   37
#define CPU_ENRICH_COLD    38
#define CPU_CSMISC         39
#define CPU_HYDNETWORK     40

#define CPU_PARTS          41  /* this gives the number of parts above (must be last) */

#define CPU_STRING_LEN 120


#ifndef  TWODIMS
#define  NUMDIMS 3		/*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470	/*!< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786	/*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */
#else
#define  NUMDIMS 2		/*!< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)	/*!< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI	/*!< Coefficient for kernel normalization. */
#endif



#if defined (BLACK_HOLES) || defined(CS_MODEL) || defined(RADTRANSFER) || defined(SNIA_HEATING)
#define PPP P
#else
#define PPP SphP
#endif


#define DMAX(a,b) (dmax1=(a),dmax2=(b),(dmax1>dmax2)?dmax1:dmax2)
#define DMIN(a,b) (dmin1=(a),dmin2=(b),(dmin1<dmin2)?dmin1:dmin2)
#define IMAX(a,b) (imax1=(a),imax2=(b),(imax1>imax2)?imax1:imax2)
#define IMIN(a,b) (imin1=(a),imin2=(b),(imin1<imin2)?imin1:imin2)

#ifdef PERIODIC
extern MyDouble boxSize, boxHalf;
#ifdef LONG_X
extern MyDouble boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
extern MyDouble boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
extern MyDouble boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif

#ifdef PERIODIC
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Y(x) (xtmp=fabs(x),(xtmp>boxHalf_Y)?(boxSize_Y-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Z(x) (xtmp=fabs(x),(xtmp>boxHalf_Z)?(boxSize_Z-xtmp):xtmp)
#else
#define NGB_PERIODIC_LONG_X(x) fabs(x)
#define NGB_PERIODIC_LONG_Y(x) fabs(x)
#define NGB_PERIODIC_LONG_Z(x) fabs(x)
#endif

#define FACT1 0.366025403785	/* FACT1 = 0.5 * (sqrt(3)-1) */



/*********************************************************/
/*  Global variables                                     */
/*********************************************************/


extern int FirstActiveParticle;
extern int *NextActiveParticle;
extern unsigned char *ProcessedFlag;

extern int TimeBinCount[TIMEBINS];
extern int TimeBinCountSph[TIMEBINS];
extern int TimeBinActive[TIMEBINS];

extern int FirstInTimeBin[TIMEBINS];
extern int LastInTimeBin[TIMEBINS];
extern int *NextInTimeBin;
extern int *PrevInTimeBin;

#ifdef SFR
extern double TimeBinSfr[TIMEBINS];
#endif

#ifdef GASRETURN
extern double TimeBinGasReturn[TIMEBINS];
#endif

#ifdef BLACK_HOLES
extern double TimeBin_BH_mass[TIMEBINS];
extern double TimeBin_BH_dynamicalmass[TIMEBINS];
extern double TimeBin_BH_Mdot[TIMEBINS];
extern double TimeBin_BH_Medd[TIMEBINS];
#endif

extern int ThisTask;		/*!< the number of the local processor  */
extern int NTask;		/*!< number of processors */
extern int PTask;		/*!< note: NTask = 2^PTask */

#ifdef INVARIANCETEST
extern int World_ThisTask;
extern int World_NTask;
extern int Color;
extern MPI_Comm MPI_CommLocal;
#ifndef DO_NOT_REDEFINE_MPI_COMM_WORLD
#undef  MPI_COMM_WORLD
#define MPI_COMM_WORLD MPI_CommLocal
#endif
#endif


extern double CPUThisRun;	/*!< Sums CPU time of current process */

extern int NumForceUpdate;	/*!< number of active particles on local processor in current timestep  */
extern long long GlobNumForceUpdate;

extern int NumSphUpdate;	/*!< number of active SPH particles on local processor in current timestep  */

extern int MaxTopNodes;	        /*!< Maximum number of nodes in the top-level tree used for domain decomposition */

extern int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */
extern int RestartSnapNum;

extern int *Exportflag;	        /*!< Buffer used for flagging whether a particle needs to be exported to another process */
extern int *Exportnodecount;
extern int *Exportindex;

extern int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount_matrix;

extern size_t AllocatedBytes;
extern size_t FreeBytes;

extern double CPU_Step[CPU_PARTS];
extern char CPU_Symbol[CPU_PARTS];
extern char CPU_SymbolImbalance[CPU_PARTS];
extern char CPU_String[CPU_STRING_LEN + 1];

extern double WallclockTime;    /*!< This holds the last wallclock time measurement for timings measurements */

extern int Flag_FullStep;	/*!< Flag used to signal that the current step involves all particles */


extern int TreeReconstructFlag;
extern int GlobFlag;

extern char DumpFlag;

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
extern long long Ntype[6];	/*!< total number of particles of each type */
extern int NtypeLocal[6];	/*!< local number of particles of each type */

extern gsl_rng *random_generator;	/*!< the random number generator used */


#ifdef SFR
extern int Stars_converted;	/*!< current number of star particles in gas particle block */
#endif

#ifdef GASRETURN
extern int Gas_converted;	/*!< current number of gas particles in the star particle block */
#endif


extern double TimeOfLastTreeConstruction;	/*!< holds what it says */

extern int *Ngblist;		/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */

extern double *R2ngblist;

extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern int *DomainStartList, *DomainEndList;

extern double *DomainWork;
extern int *DomainCount;
extern int *DomainCountSph;
extern int *DomainTask;
extern int *DomainNodeIndex;
extern int *DomainList, DomainNumChanged;

extern peanokey *Key, *KeySorted;

extern struct topnode_data
{
  peanokey Size;
  peanokey StartKey;
  long long Count;
  MyFloat GravCost;
  int Daughter;
  int Pstart;
  int Blocks;
  int Leaf;
} *TopNodes;

extern int NTopnodes, NTopleaves;

extern double RndTable[RNDTABLE];


#ifdef SUBFIND
extern int GrNr;
extern int NumPartGroup;
#endif


/* variables for input/output , usually only used on process 0 */


extern char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

extern FILE *FdInfo,		/*!< file handle for info.txt log-file. */
 *FdEnergy,			/*!< file handle for energy.txt log-file. */
 *FdTimings,			/*!< file handle for timings.txt log-file. */
 *FdBalance,			/*!< file handle for balance.txt log-file. */
 *FdCPU;			/*!< file handle for cpu.txt log-file. */

#ifdef SFR
extern FILE *FdSfr;		/*!< file handle for sfr.txt log-file. */
#endif

#ifdef GASRETURN
extern FILE *FdGasReturn;	/*!< file handle for GasReturn.txt log-file */
#endif

#ifdef RADTRANSFER
extern FILE *FdEddington;     /*!< file handle for eddington.txt log-file. */
extern FILE *FdRadtransfer;		/*!< file handle for radtransfer.txt log-file. */
#endif

#ifdef DISTORTIONTENSORPS
#ifdef PMGRID
extern FILE *FdTidaltensor;     /*!< file handle for tidaltensor.txt log-file. */
#endif
extern FILE *FdCaustics;	/*!< file handle for Caustics.txt log-file. */
#endif

#ifdef BLACK_HOLES
extern FILE *FdBlackHoles;	/*!< file handle for blackholes.txt log-file. */
extern FILE *FdBlackHolesDetails;
#endif


#ifdef FORCETEST
extern FILE *FdForceTest;	/*!< file handle for forcetest.txt log-file. */
#endif

#ifdef DARKENERGY
extern FILE *FdDE;  /*!< file handle for darkenergy.txt log-file. */
#endif

#ifdef XXLINFO
extern FILE *FdXXL;		/*!< file handle for xxl.txt log-file. */

#ifdef MAGNETIC
extern double MeanB;

#ifdef TRACEDIVB
extern double MaxDivB;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
extern double MeanAlpha;
#endif
#endif

/*! table for the cosmological drift factors */
extern double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
extern double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
extern double HydroKickTable[DRIFT_TABLE_LENGTH];



extern void *CommBuffer;	/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!<  total particle numbers (global value) */
  long long TotN_gas;		/*!<  total gas particle number (global value) */

#ifdef BLACK_HOLES
  int TotBHs;
#endif

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one
				   processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one
				   processor. */

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int DoDynamicUpdate;

  int NumFilesPerSnapshot;	/*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;	/*!< maximum number of files that may be written simultaneously when
					   writing/reading restart-files, or when writing snapshot files */

  int BufferSize;		/*!< size of communication buffer in MB */
  int BunchSize;     	        /*!< number of particles fitting into the buffer in the parallel tree algorithm  */


  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				   NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				   the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  double TopNodeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

#ifdef SCALARFIELD
  double ScalarBeta;
  double ScalarScreeningLength;
#endif

  /* some SPH parameters */

  int DesNumNgb;		/*!< Desired number of SPH neighbours */
#ifdef SUBFIND
  int DesLinkNgb;
  double ErrTolThetaSubfind;
#endif

  double MaxNumNgbDeviation;	/*!< Maximum allowed deviation neighbour number */
#ifdef START_WITH_EXTRA_NGBDEV
  double MaxNumNgbDeviationStart;    /*!< Maximum allowed deviation neighbour number to start with*/
#endif

  double ArtBulkViscConst;	/*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double InitGasU;		/*!< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;		/*!< the minimum allowed temperature expressed as energy per unit mass */

#ifdef GASRETURN
  double GasReturnFraction;
#endif


  /* some force counters  */

  long long TotNumOfForces;	/*!< counts total number of force computations  */

  long long NumForcesSinceLastDomainDecomp;	/*!< count particle updates since last domain decomposition */

  /* some variable for dynamic work-load adjustment based on CPU measurements */

  double Cadj_Cost;
  double Cadj_Cpu;

  /* system of units  */

  double UnitTime_in_s,		/*!< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,		/*!< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,	/*!< factor to convert intqernal velocity unit to cm/sec */
    UnitLength_in_cm,		/*!< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,	/*!< factor to convert internal pressure unit to cgs units (little 'h' still
				   around!) */
    UnitDensity_in_cgs,		/*!< factor to convert internal length unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,	/*!< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,		/*!< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,	/*!< factor to convert internal time to megayears/h */
    GravityConstantInternal,	/*!< If set to zero in the parameterfile, the internal value of the
				   gravitational constant is set to the Newtonian value based on the system of
				   units specified. Otherwise the value provided is taken as internal gravity
				   constant G. */
    G;				/*!< Gravity-constant in internal units */
  double UnitDensity_in_Gev_per_cm3; /*!< factor to convert internal density unit to GeV/c^2 / cm^3 */
  /* Cosmology */

  double Hubble;		/*!< Hubble-constant in internal units */
  double Omega0,		/*!< matter density in units of the critical density (at z=0) */
    OmegaLambda,		/*!< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,		/*!< baryon density in units of the critical density (at z=0) */
    HubbleParam;		/*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
				 * physical values for cooling physics
				 */

  double BoxSize;		/*!< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;	/*!< flags that periodic boundaries are enabled */
  int ResubmitOn;		/*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;	/*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
				   criterion */
  int TypeOfTimestepCriterion;	/*!< gives type of timestep criterion (only 0 supported right now - unlike
				   gadget-1.1) */
  int OutputListOn;		/*!< flags that output times are listed in a specified file */
  int CoolingOn;		/*!< flags that cooling is enabled */
  int StarformationOn;		/*!< flags that star formation is enabled */


  /* parameters determining output frequency */

  int SnapshotFileCount;	/*!< number of snapshot that is written next */
  double TimeBetSnapshot,	/*!< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,	/*!< simulation time of first snapshot files */
    CpuTimeBetRestartFile,	/*!< cpu-time between regularly generated restart files */
    TimeLastRestartFile,	/*!< cpu-time when last restart-file was written */
    TimeBetStatistics,		/*!< simulation time interval between computations of energy statistics */
    TimeLastStatistics;		/*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;		/*!< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,			/*!< current time of the simulation */
    TimeBegin,			/*!< time of initial conditions of the simulation */
    TimeStep,			/*!< difference between current times of previous and current timestep */
    TimeMax;			/*!< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval;	/*!< factor to convert from floating point time interval to integer timeline */
  int Ti_Current;		/*!< current time on integer timeline */
  int Ti_nextoutput;		/*!< next output time on integer timeline */

#ifdef PMGRID
  int PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], UpperCorner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];
#endif

#ifdef CHEMISTRY
  double Epsilon;
#endif

  int Ti_nextlineofsight;
#ifdef OUTPUTLINEOFSIGHT
  double TimeFirstLineOfSight;
#endif

  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_Sum[CPU_PARTS];    /*!< sums wallclock time/CPU consumption in whole run */

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
				   timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep,	/*!< minimum allowed timestep. Normally, the simulation terminates if the
				   timestep determined by the timestep criteria falls below this limit. */
    MaxSizeTimestep;		/*!< maximum allowed timestep */

  double MaxRMSDisplacementFac;	/*!< this determines a global timestep criterion for cosmological simulations
				   in comoving coordinates.  To this end, the code computes the rms velocity
				   of all particles, and limits the timestep such that the rms displacement
				   is a fraction of the mean particle separation (determined from the
				   particle mass and the cosmological parameters). This parameter specifies
				   this fraction. */



  double CourantFac;		/*!< SPH-Courant factor */


  /* frequency of tree reconstruction/domain decomposition */


  double TreeDomainUpdateFrequency;	/*!< controls frequency of domain decompositions  */


  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   * five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars
   */
  double MinGasHsmlFractional,	/*!< minimum allowed SPH smoothing length in units of SPH gravitational
				   softening length */
    MinGasHsml;			/*!< minimum allowed SPH smoothing length */


  double SofteningGas,		/*!< for type 0 */
    SofteningHalo,		/*!< for type 1 */
    SofteningDisk,		/*!< for type 2 */
    SofteningBulge,		/*!< for type 3 */
    SofteningStars,		/*!< for type 4 */
    SofteningBndry;		/*!< for type 5 */

  double SofteningGasMaxPhys,	/*!< for type 0 */
    SofteningHaloMaxPhys,	/*!< for type 1 */
    SofteningDiskMaxPhys,	/*!< for type 2 */
    SofteningBulgeMaxPhys,	/*!< for type 3 */
    SofteningStarsMaxPhys,	/*!< for type 4 */
    SofteningBndryMaxPhys;	/*!< for type 5 */

  double SofteningTable[6];	/*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];	/*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


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
    InfoFile[100], TimingsFile[100], RestartFile[100], ResubmitCommand[100], OutputListFilename[100];

  /*! table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];
  char OutputListFlag[MAXLEN_OUTPUTLIST];
  int OutputListLength;		/*!< number of times stored in table of desired output times */



#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
  double ReferenceGasMass;
#endif

#ifdef VORONOI_MESHRELAX
  double MeanMass;
  double MeanPressure;
#endif

#ifdef RADTRANSFER
  double IonizingLumPerSolarMass;
  double IonizingLumPerSFR;
  int Radiation_Ti_begstep;
  int Radiation_Ti_endstep;
#endif

#if defined(SIM_ADAPTIVE_SOFT) || defined(REINIT_AT_TURNAROUND)
  double CurrentTurnaroundRadius;
  double InitialTurnaroundRadius;
  double SIM_epsilon; 
#endif
#ifdef DISTORTIONTENSORPS
  /* present day velocity dispersion of DM particle in cm/s (e.g. Neutralino = 0.03 cm/s) */
  double DM_velocity_dispersion;
#endif

#ifdef SFR		/* star formation and feedback sector */
  double CritOverDensity;
  double CritPhysDensity;
  double OverDensThresh;
  double PhysDensThresh;
  double EgySpecSN;
  double FactorSN;
  double EgySpecCold;
  double FactorEVP;
  double FeedbackEnergy;
  double TempSupernova;
  double TempClouds;
  double MaxSfrTimescale;
  double WindEfficiency;
  double WindEnergyFraction;
  double WindFreeTravelLength;
  double WindFreeTravelDensFac;
  double FactorForSofterEQS;
#endif

#ifdef CS_MODEL
  double FactorSFR;
  double DecouplingParam;
  double MinTlifeSNI;
  double MaxTlifeSNI;
  double TlifeSNII;
  int    Raiteri_TlifeSNII;
  double RateSNI;
  double SN_Energy_cgs;
  double Tcrit_Phase;
  double DensFrac_Phase;
  double SN_Energy_frac_cold;
  double MaxHotHsmlParam;
  double InitialHotHsmlFactor;
  double DensityTailThreshold;
  double MaxNumHotNgbDeviation;	/*!< Maximum allowed deviation HOT neighbour number */
#endif

#ifdef DARKENERGY
  double DarkEnergyParam;	/*!< fixed w for equation of state */
#ifdef TIMEDEPDE
  char DarkEnergyFile[100];	/*!< tabelized w for equation of state */
#ifdef TIMEDEPGRAV
  double Gini;
#endif
#endif
#endif

#ifdef RESCALEVINI
  double VelIniScale;		/*!< Scale the initial velocities by this amount */
#endif

#if defined(SNIA_HEATING)
  double SnIaHeatingRate;
#endif


#ifdef TIME_DEP_ART_VISC
  double ViscSource0;		/*!< Given sourceterm in viscosity evolution */
  double DecayLength;		/*!< Number of h for the viscosity decay */
  double ViscSource;		/*!< Reduced sourceterm in viscosity evolution */
  double DecayTime;		/*!< Calculated decaytimescale */
  double AlphaMin;		/*!< Minimum of allowed viscosity parameter */
#endif

#ifdef CONDUCTION
  double ConductionCoeff;	/*!< Thermal Conductivity */
#ifdef CONDUCTION_SATURATION
  double ElectronFreePathFactor;	/*!< Factor to get electron mean free path */
#endif

  int Conduction_Ti_begstep, Conduction_Ti_endstep;
  double MaxSizeConductionStep;
#endif

#ifdef MAGNETIC
#ifdef BINISET
  double BiniX, BiniY, BiniZ;	/*!< Initial values for B */
#endif

#ifdef BSMOOTH
  int BSmoothInt;
  double BSmoothFrac;
  int MainTimestepCounts;
#ifdef SETMAINTIMESTEPCOUNT
  int MainTimestepCountIni;
#endif
#endif

#ifdef MAGNETIC_DISSIPATION
  double ArtMagDispConst;	/*!< Sets the parameter \f$\alpha\f$ of the artificial magnetic disipation */
#ifdef TIME_DEP_MAGN_DISP
  double ArtMagDispMin;
  double ArtMagDispSource;
  double ArtMagDispTime;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
  double DivBcleanParabolicSigma;
  double DivBcleanHyperbolicSigma;
#endif

#ifdef HEALPIX
  //change this to read in the Parameterfile
  int Nside;
#define NSIDE2NPIX(nside)  (12*nside*nside)
  float *healpixmap;
#endif

#ifdef MAGNETIC_DIFFUSION
  double MagneticEta;
#endif
#endif

#ifdef BLACK_HOLES
  double TimeNextBlackHoleCheck;
  double TimeBetBlackHoleSearch;
  double BlackHoleAccretionFactor;	/*!< Fraction of BH bondi accretion rate */
  double BlackHoleFeedbackFactor;	/*!< Fraction of the black luminosity feed into thermal feedback */
  double SeedBlackHoleMass;	/*!< Seed black hole mass */
  double MinFoFMassForNewSeed;	/*!< Halo mass required before new seed is put in */
  double BlackHoleNgbFactor;	/*!< Factor by which the normal SPH neighbour should be increased/decreased */
  double BlackHoleMaxAccretionRadius;
  double BlackHoleEddingtonFactor;	/*! Factor above Eddington */
#ifdef FOF
  double massDMpart;
#endif
#ifdef MODIFIEDBONDI
  double BlackHoleRefDensity;
  double BlackHoleRefSoundspeed;
#endif
#endif

#ifdef COSMIC_RAYS
  double CR_Alpha[NUMCRPOP];	/*!< Cosmic ray spectral index [2..3] */
  double CR_SNEff;		/*!< SN injection efficiency [0..1] */
  double CR_SNAlpha;		/*!< SN injection spectral index [2..3] */
  int bDebugFlag;		/*!< enables debug outputs after triggered */

#if defined(CR_DIFFUSION)
  double CR_DiffusionCoeff;	/*!< (temporary) fixed value for CR diffusivity */

  double CR_DiffusionDensScaling;	/*!< grade of density dependence of diffusivity */
  double CR_DiffusionDensZero;	/*!< Reference point density for diffusivity */

  double CR_DiffusionEntropyScaling;	/*!< grade of specific energy dependence of diffusivity */

  double CR_DiffusionEntropyZero;	/*!< Reference Entropic function for diffusivity */

  double CR_DiffusionMaxSizeTimestep;
  int CR_Diffusion_Ti_begstep, CR_Diffusion_Ti_endstep;
#endif				/* CR_DIFFUSION */

#if defined(CR_SHOCK)
#if (CR_SHOCK == 1)
  double CR_ShockAlpha;		/*!< spectral index to be used in shock injection */
#else
  double CR_ShockCutoff;	/*!< Cutoff factor x_inj for CR accel */
#endif
  double CR_ShockEfficiency;	/*!< energy fraction of shock energy fed into CR */
#endif				/* CR_SHOCK */

#ifdef FIX_QINJ
  double Shock_Fix_Qinj;	/*!< inject only CRps with threshold cutoff Shock_Fix_Qinj */
#endif

#ifdef CR_BUBBLES
  double CR_AGNEff;               /*!< AGN injection efficiency [0..1] */
#endif
#endif				/* COSMIC_RAYS */

#ifdef MACHNUM
  double Shock_Length;		/*!< length scale on which the shock is smoothed out */
  double Shock_DeltaDecayTimeMax;	/*!< maximum time interval (Dloga) for which the 
					   Mach number is kept at its maximum */
#endif

#ifdef REIONIZATION
  int not_yet_reionized;	/*!< flag that makes sure that there is only one reionization */
#endif



#ifdef BUBBLES
  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
  double FirstBubbleRedshift;
#ifdef FOF
  int BiggestGroupLen;
  float BiggestGroupCM[3];
  double BiggestGroupMass;
#endif
#endif

#ifdef BH_BUBBLES
  double BubbleDistance;
  double BubbleRadius;
  double BubbleEnergy;
  double BlackHoleRadioTriggeringFactor;
  double DefaultICMDensity;
  double RadioFeedbackFactor;
#ifdef UNIFIED_FEEDBACK
  double RadioThreshold;
#endif
#endif

#if defined(MULTI_BUBBLES) && defined(FOF)
#ifndef BLACK_HOLES
  double MinFoFMassForNewSeed;	/*!< Halo mass required before new seed is put in */
  double massDMpart;
#endif
  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
  double ClusterMass200;
  double FirstBubbleRedshift;
#endif

#ifdef NAVIERSTOKES
  double NavierStokes_ShearViscosity;
  double FractionSpitzerViscosity;
  double ShearViscosityTemperature;
#endif
#ifdef NAVIERSTOKES_BULK
  double NavierStokes_BulkViscosity;
#endif
#ifdef VISCOSITY_SATURATION
  double IonMeanFreePath;
#endif

#ifdef EOS_DEGENERATE
  char EosTable[100];
  char EosSpecies[100];
#endif

#ifdef NUCLEAR_NETWORK
  char NetworkData[100];
#endif

#ifdef GASRETURN
  double AGBGasTemp;
  double AGBGasDensity;
  double AGBGasEnergy;
  double AGBGasEntropy;
#endif
}
All;



/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data
{
  MyDouble Pos[3];   /*!< particle position at its current time */
  MyDouble Vel[3];   /*!< particle velocity at its current time */
  MyDouble Mass;     /*!< particle mass */
  MyIDType ID;

  union
  {
    MyFloat       GravAccel[3];		/*!< particle acceleration due to gravity */
    MyLongDouble dGravAccel[3];
  } g;
#ifdef PMGRID
  MyFloat GravPM[3];		/*!< particle acceleration due to long-range PM gravity force */
#endif
#ifdef FORCETEST
  MyFloat GravAccelDirect[3];	/*!< particle acceleration calculated by direct summation */
#endif
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
  union
  {
    MyFloat       Potential;		/*!< gravitational potential */
    MyLongDouble dPotential;
  } p;
#endif

#ifdef DISTORTIONTENSORPS 
  MyLongDouble distortion_tensorps[6][6];          /*!< Phase Space Distortion tensor */
  MyLongDouble tidal_tensorps[3][3];               /*!< tidal tensor (=second derivatives of grav. potential) */ 
  MyLongDouble V_matrix[3][3];                     /*!< initial orientation of CDM sheet the particle is embedded in */ 
  MyDouble init_density;                           /*!< initial stream density */ 
  MyFloat caustic_counter;                         /*!< caustic counter */  
  MyDouble last_stream_determinant;                /*!< last stream density determinant, needed to identify caustics */
#ifdef REINIT_AT_TURNAROUND 
  int turnaround_flag;                             /*!< mark when a particle turned around */
#endif
  MyDouble annihilation;                            /*!< integrated annihilation rate */
  MyDouble rho_normed_cutoff_current;               /*!< current and last normed_cutoff density in rho_max/rho_init * sqrt(sigma) */
  MyDouble rho_normed_cutoff_last; 
  MyDouble s_1_current, s_2_current, s_3_current;   /*! < current and last stretching factor */
  MyDouble s_1_last, s_2_last, s_3_last;
  MyDouble second_deriv_current;                    /*! < current and last second derivative */
  MyDouble second_deriv_last; 
  MyDouble stream_density;                          /*!< physical stream density that is going to be integrated (in terms of rho_crit) */ 
  MyFloat analytic_caustics;                        /*!< number of caustics that were integrated analyticall, 
                                                         i.e. where the physical caustic density was higher 
                                                         than the numerical GDE density */
#ifdef OUTPUT_LAST_CAUSTIC
  MyDouble lc_Time;                          /*!< time of caustic passage */ 
  MyDouble lc_Pos[3];                        /*!< position of caustic */
  MyDouble lc_Vel[3];                        /*!< particle velocity when passing through caustic */
  MyDouble lc_rho_normed_cutoff;             /*!< normed_cutoff density at caustic */
  MyDouble lc_Dir_x[3];                      /*!< principal axis frame of smear out */
  MyDouble lc_Dir_y[3];
  MyDouble lc_Dir_z[3];  
  MyDouble lc_smear_x;                       /*!< smear out length */  
  MyDouble lc_smear_y;
  MyDouble lc_smear_z;          
#endif
#ifdef PMGRID 
  MyLongDouble tidal_tensorpsPM[3][3];	    /*!< for TreePM simulations, long range tidal field */	
#endif
#endif

  MyFloat OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening
                                          criterion */
#if defined(EVALPOTENTIAL) && defined(PMGRID)
  MyFloat PM_Potential;
#endif

#ifdef STELLARAGE
  MyFloat StellarAge;		/*!< formation time of star particle */
#endif
#ifdef GASRETURN
  MyFloat GasReturnRate;
  MyFloat GasAge;
#endif
#ifdef METALS
  MyFloat Metallicity;		/*!< metallicity of gas or star particle */
#endif				/* closes METALS */

#if defined (BLACK_HOLES) || defined(CS_MODEL) || defined(RADTRANSFER) || defined(SNIA_HEATING)
  MyFloat Hsml;

  union
  {
    MyFloat       NumNgb;
    MyLongDouble dNumNgb;
  } n;
#if defined(RADTRANSFER) || defined(SNIA_HEATING)
  MyFloat DensAroundStar;
#endif
#endif


#ifdef BLACK_HOLES
  int SwallowID;
#ifdef BH_COUNTPROGS
  int BH_CountProgs;
#endif
  MyFloat BH_Mass;
  MyFloat BH_Mdot;
#ifdef BH_BUBBLES
  MyFloat BH_Mass_bubbles;
  MyFloat BH_Mass_ini;
#ifdef UNIFIED_FEEDBACK
  MyFloat BH_Mass_radio;
#endif
#endif
  union
  {
    MyFloat BH_Density;
    MyLongDouble dBH_Density;
  } b1;
  union
  {
    MyFloat BH_Entropy;
    MyLongDouble dBH_Entropy;
  } b2;
  union
  {
    MyFloat BH_SurroundingGasVel[3];
    MyLongDouble dBH_SurroundingGasVel[3];
  } b3;
  union
  {
    MyFloat BH_accreted_Mass;
    MyLongDouble dBH_accreted_Mass;
  } b4;
  union
  {
    MyFloat BH_accreted_BHMass;
    MyLongDouble dBH_accreted_BHMass;
  } b5;
  union
  {
    MyFloat BH_accreted_momentum[3];
    MyLongDouble dBH_accreted_momentum[3];
  } b6;
#ifdef BH_BUBBLES
  union
  {
    MyFloat BH_accreted_BHMass_bubbles;
    MyLongDouble dBH_accreted_BHMass_bubbles;
  } b7;
#ifdef UNIFIED_FEEDBACK
  union
  {
    MyFloat BH_accreted_BHMass_radio;
    MyLongDouble dBH_accreted_BHMass_radio;
  } b8;
#endif
#endif
#ifdef REPOSITION_ON_POTMIN
  MyFloat BH_MinPotPos[3];
  MyFloat BH_MinPot;
#endif
#ifdef BH_KINETICFEEDBACK
  MyFloat ActiveTime;
  MyFloat ActiveEnergy;
#endif
#endif

#ifdef SUBFIND
  unsigned int GrNr;
  int SubNr;
  int DM_NumNgb;
  unsigned short targettask, origintask2; 
  int origintask, submark, origindex;
  MyFloat DM_Hsml;
  union
  {
    MyFloat DM_Density;	
    MyFloat DM_Potential;
  } u;
  union
  {
    MyFloat DM_VelDisp;
    MyFloat DM_BindingEnergy;
  } v;
#ifdef DENSITY_SPLIT_BY_TYPE
  union
  {
    MyFloat int_energy;
    MyFloat density_sum;
  } w;
#endif

#ifdef SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI
  MyFloat DM_Hsml_V;	
  MyFloat DM_Density_V;	
#endif

#ifdef SAVE_HSML_IN_IC_ORDER
  MyIDType ID_ic_order;
#endif
#ifdef SUBFIND_ALTERNATIVE_COLLECTIVE
  peanokey Key;
#endif
#endif

#if defined(ORDER_SNAPSHOTS_BY_ID) && !defined(SUBFIND) 
  int     GrNr;
  int     SubNr;
#endif

#ifdef SHELL_CODE
  MyDouble radius;
  MyDouble enclosed_mass;
  MyDouble dMdr;
#endif

#ifdef CS_MODEL
  MyFloat Zm[12];
  MyFloat ZmReservoir[12];
#ifdef CS_FEEDBACK
  MyFloat EnergySN;
  MyFloat EnergySNCold;    
#endif 
#endif

  float GravCost;		/*!< weight factor used for balancing the work-load */

  int Ti_begstep;		/*!< marks start of current timestep of particle on integer timeline */
  int Ti_current;		/*!< current time of the particle */

  short int Type;		/*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  short int TimeBin;

#ifdef WAKEUP
  int dt_step;
#endif
}
 *P,				/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */



/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
extern struct sph_particle_data
{
  MyDouble Entropy;		/*!< current value of entropy (actually entropic function) of particle */
  MyFloat  Pressure;		/*!< current pressure */
  MyFloat  VelPred[3];		/*!< predicted SPH particle velocity at the current time */
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  MyFloat MinViscousDt;
#else
  MyFloat MaxSignalVel;           /*!< maximum signal velocity */
#endif
#ifdef VOLUME_CORRECTION
  MyFloat DensityOld;
  MyFloat DensityStd;
#endif

#ifdef VORONOI
  MyFloat MaxDelaunayRadius;
  MyFloat Volume;
  MyFloat Center[3];
#endif

#ifdef GASRETURN
  MyFloat GasAge;
#endif

  union
  {
    MyFloat       Density;		/*!< current baryonic mass density of particle */
    MyLongDouble dDensity;
  } d;
  union
  {
    MyFloat       DtEntropy;		/*!< rate of change of entropy */
    MyLongDouble dDtEntropy;
  } e;
  union
  {
    MyFloat       HydroAccel[3];	/*!< acceleration due to hydrodynamical force */
    MyLongDouble dHydroAccel[3];
  } a;
  union
  {
    MyFloat       DhsmlDensityFactor;	/*!< correction factor needed in entropy formulation of SPH */
    MyLongDouble dDhsmlDensityFactor;
  } h;
  union
  {
    MyFloat       DivVel;		/*!< local velocity divergence */
    MyLongDouble dDivVel;
  } v;
#ifndef NAVIERSTOKES
  union
  {
    MyFloat CurlVel;     	        /*!< local velocity curl */
    MyFloat       Rot[3];		/*!< local velocity curl */
    MyLongDouble dRot[3];
  } r;
#else
  union
  {
    MyFloat DV[3][3];
    struct
    {
      MyFloat DivVel;
      MyFloat CurlVel;
      MyFloat StressDiag[3];
      MyFloat StressOffDiag[3];
#ifdef NAVIERSTOKES_BULK
      MyFloat StressBulk;
#endif
    } s;
  } u;
#endif

#if !(defined(BLACK_HOLES) || defined(CS_MODEL) || defined(RADTRANSFER) || defined(SNIA_HEATING))
  MyFloat Hsml;			/*!< current smoothing length */
  union
  {
    MyFloat       NumNgb;
    MyLongDouble dNumNgb;
  } n;
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
  union
  {
    MyFloat       Injected_BH_Energy;
    MyLongDouble dInjected_BH_Energy;
  } i;
#endif

#ifdef COOLING
  MyFloat Ne;  /*!< electron fraction, expressed as local electron number
		    density normalized to the hydrogen number density. Gives
		    indirectly ionization state and mean molecular weight. */
#endif
#ifdef SFR
  MyFloat Sfr;
#endif
#ifdef WINDS
  MyFloat DelayTime;		/*!< remaining maximum decoupling time of wind particle */
#endif


#ifdef MAGNETIC
  MyFloat BPred[3];
#ifdef EULERPOTENTIALS
  MyFloat EulerA,EulerB;
  MyFloat dEulerA[3],dEulerB[3];
#else
  MyFloat B[3];
  MyFloat DtB[3];
#endif
#if defined(TRACEDIVB) || defined(TIME_DEP_MAGN_DISP)
  MyFloat divB;
#endif
#ifdef VECT_PRO_CLEAN
  MyFloat BPredVec[3];
#endif
#if defined(BSMOOTH) || defined(BFROMROTA)
  MyFloat BSmooth[3];
#endif
#ifdef TIME_DEP_MAGN_DISP
  MyFloat Balpha, DtBalpha;
#endif
#ifdef DIVBCLEANING_DEDNER
  MyFloat Phi, PhiPred, DtPhi;
#ifdef SMOOTH_PHI
  MyFloat SmoothPhi;
#endif
#endif
#if defined(SMOOTH_PHI)
  MyFloat SmoothDivB;
#endif
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS) || defined(VECT_PRO_CLEAN)
  MyFloat RotB[3];
#ifdef SMOOTH_ROTB
  MyFloat SmoothedRotB[3];
#endif
#endif
#if defined(BSMOOTH) || defined(SMOOTH_ROTB) || defined(SMOOTH_PHI)
  MyFloat DensityNorm;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
  MyFloat alpha, Dtalpha;
#endif
#ifdef NS_TIMESTEP
  MyFloat ViscEntropyChange;
#endif
#ifdef CONDUCTION_SATURATION
  MyFloat GradEntr[3];
#endif

#ifdef MHM
  MyFloat FeedbackEnergy;
#endif

#ifdef COSMIC_RAYS
  MyFloat CR_C0[NUMCRPOP];			/*!< Cosmic ray amplitude adiabatic invariable */
  MyFloat CR_q0[NUMCRPOP];			/*!< Cosmic ray cutoff adiabatic invariable */
  MyFloat CR_E0[NUMCRPOP];			/*!< Specific Energy at Rho0 */
  MyFloat CR_n0[NUMCRPOP];			/*!< baryon fraction in cosmic rays */

  MyFloat CR_DeltaE[NUMCRPOP];		/*!< Specific Energy growth during timestep */
  MyFloat CR_DeltaN[NUMCRPOP];		/*!< baryon fraction growth during timestep */
#ifdef MACHNUM
  MyFloat CR_Gamma0[NUMCRPOP];
#endif

#ifdef CR_OUTPUT_INJECTION
  MyFloat CR_Specific_SupernovaHeatingRate;
#endif
#endif				/* COSMIC_RAYS */

#ifdef MACHNUM
  MyFloat Shock_MachNumber;	/*!< Mach number */
  MyFloat Shock_DecayTime;	/*!< Shock decay time */
#ifdef COSMIC_RAYS
  MyFloat Shock_DensityJump;	/*!< Density jump at the shock */
  MyFloat Shock_EnergyJump;	/*!< Energy jump at the shock */
  MyFloat PreShock_PhysicalDensity;	/*!< Specific energy in the preshock regime */
  MyFloat PreShock_PhysicalEnergy;	/*!< Density in the preshock regime */
  MyFloat PreShock_XCR;		/*!< XCR = PCR / Pth in the preshock regime */
#endif
#ifdef MACHSTATISTIC
  MyFloat Shock_DtEnergy;		/*!< Change of thermal specific energy at Shocks */
#endif
#ifdef OUTPUT_PRESHOCK_CSND
  MyFloat PreShock_PhysicalSoundSpeed;	/*!< Sound speed in the preshock regime */
  MyFloat PreShock_PhysicalDensity;	/*!< Specific energy in the preshock regime */
#endif
#endif				/* Mach number estimate */


#ifdef CHEMISTRY
  MyFloat elec;
  MyFloat HI;
  MyFloat HII;

  MyFloat HeI;
  MyFloat HeII;
  MyFloat HeIII;

  MyFloat H2I;
  MyFloat H2II;

  MyFloat HM;

  MyFloat Gamma;
  MyFloat t_elec, t_cool;
#endif

#ifdef RADTRANSFER
  MyFloat ET[6];                /* eddington tensor - symmetric -> only 6 elements needed */
  MyFloat Je;                   /* emissivity */
  MyFloat nHI;                  /* HI number density */
  MyFloat nHII;                 /* HII number density */
  MyFloat n_elec;               /* electron number density */
  MyFloat n_gamma;              /* photon number density */
#ifdef RADTRANSFER_FLUXLIMITER
  MyFloat Grad_ngamma[3];
#endif
#ifndef CG
  MyFloat n_gamma_old;
#endif
#ifdef RT_VAR_DT
  int rt_flag;
#endif
#endif

#if defined CS_MODEL
  MyFloat DensityOld;
#ifdef CS_FEEDBACK
  union
  {
    MyFloat       DensityAvg;		/*!< current baryonic mass density of particle */
    MyLongDouble dDensityAvg;
  } da;
  union
  {
    MyFloat       EntropyAvg;		/*!< current baryonic mass density of particle */
    MyLongDouble dEntropyAvg;
  } ea;
  MyFloat HotHsml;
  int     HotNgbNum;
  MyFloat DensPromotion;
  MyFloat TempPromotion;
#endif
#endif

#ifdef EOS_DEGENERATE
  MyFloat u;                            /* internal energy density */
  MyFloat temp;                         /* temperature */
  MyFloat dpdr;							/* derivative of pressure with respect to density at constant entropy */
  MyFloat xnuc[EOS_NSPECIES];           /* nuclear mass fractions */
  MyFloat dxnuc[EOS_NSPECIES];          /* change of nuclear mass fractions */
#endif

#ifdef WAKEUP
  short int wakeup;             /*!< flag to wake up particle */
#endif
}
  *SphP,				/*!< holds SPH particle data on local processor */
  *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */


extern peanokey *DomainKeyBuf;

/* global state of system 
*/
extern struct state_of_system
{
  double Mass,
    EnergyKin,
    EnergyPot,
    EnergyInt,
    EnergyTot,
    Momentum[4],
    AngMomentum[4],
    CenterOfMass[4],
    MassComp[6],
    EnergyKinComp[6],
    EnergyPotComp[6],
    EnergyIntComp[6], EnergyTotComp[6], MomentumComp[6][4], AngMomentumComp[6][4], CenterOfMassComp[6][4];
}
SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */

extern struct data_index
{
  int Task;
  int Index;
  int IndexGet;
}
 *DataIndexTable;		/*!< the particles to be exported are grouped
                                  by task-number. This table allows the
                                  results to be disentangled again and to be
                                  assigned to the correct particle */

extern struct data_nodelist
{
  int NodeList[NODELISTLENGTH];
}
*DataNodeList;		

extern struct gravdata_in
{
  MyFloat Pos[3];
#ifdef UNEQUALSOFTENINGS
  int Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  MyFloat Soft;
#endif
#endif
  MyFloat OldAcc;
  int NodeList[NODELISTLENGTH];
}
 *GravDataIn,			/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


extern struct gravdata_out
{
  MyLongDouble Acc[3];
#ifdef EVALPOTENTIAL
  MyLongDouble Potential;
#endif
#ifdef DISTORTIONTENSORPS
  MyLongDouble tidal_tensorps[3][3];
#endif
  int Ninteractions;
}
 *GravDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


extern struct potdata_out
{
  MyLongDouble Potential;
}
 *PotDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


/*! Header for the standard file format.
 */
extern struct io_header
{
  int npart[6];			/*!< number of particles of each type in this file */
  double mass[6];		/*!< mass of particles of each type. If 0, then the masses are explicitly
				   stored in the mass-block of the snapshot file, otherwise they are omitted */
#ifdef COSMIC_RAYS
  double SpectralIndex_CR_Pop[NUMCRPOP]; /*!< spectral indices of cosmic ray populations */
#endif
  double time;			/*!< time of snapshot file */
  double redshift;		/*!< redshift of snapshot file */
  int flag_sfr;			/*!< flags whether the simulation was including star formation */
  int flag_feedback;		/*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];	/*!< total number of particles of each type in this snapshot. This can be
				   different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;		/*!< flags whether cooling was included  */
  int num_files;		/*!< number of files in multi-file snapshot */
  double BoxSize;		/*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;		/*!< matter density in units of critical density */
  double OmegaLambda;		/*!< cosmological constant parameter */
  double HubbleParam;		/*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;		/*!< flags whether the file contains formation times of star particles */
  int flag_metals;		/*!< flags whether the file contains metallicity values for gas and star
				   particles */
  unsigned int npartTotalHighWord[6];	/*!< High word of the total number of particles of each type */
  int flag_entropy_instead_u;	/*!< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;	/*!< flags that snapshot contains double-precision instead of single precision */

  int flag_ic_info;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                     or whether they ocontains 2nd order lagrangian perturbation theory initial conditions. 
                                     For snapshots files, the value informs whether the simulation was evolved from
                                     Zeldoch or 2lpt ICs. Encoding is as follows:
                                        FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                        FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                        FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                        FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                        FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs 
                                     All other values, including 0 are interpreted as "don't know" for backwards compatability.   
                                 */
  float lpt_scalingfactor;      /*!< scaling factor for 2lpt initial conditions */

#ifdef COSMIC_RAYS
  char fill[48-8*NUMCRPOP];	/*!< fills to 256 Bytes */
#else
  char fill[48];		/*!< fills to 256 Bytes */
#endif

}
header;				/*!< holds header for snapshot files */



enum iofields
{ IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_SECONDORDERMASS,
  IO_U,
  IO_RHO,
  IO_NE,
  IO_NH,
  IO_HSML,
  IO_SFR,
  IO_AGE,
  IO_Z,
  IO_BHMASS,
  IO_BHMDOT,
  IO_BHPROGS,
  IO_BHMBUB,
  IO_BHMINI,
  IO_BHMRAD,
  IO_POT,
  IO_ACCEL,
  IO_CR_C0,
  IO_CR_Q0,
  IO_CR_P0,
  IO_CR_E0,
  IO_CR_n0,
  IO_CR_ThermalizationTime,
  IO_CR_DissipationTime,
  IO_ELECT,
  IO_HI,
  IO_HII,
  IO_HeI,
  IO_HeII,
  IO_HeIII,
  IO_H2I,
  IO_H2II,
  IO_HM,
  IO_DTENTR,
  IO_STRESSDIAG,
  IO_STRESSOFFDIAG,
  IO_STRESSBULK,
  IO_SHEARCOEFF,
  IO_TSTP,
  IO_BFLD,
  IO_BSMTH,
  IO_DBDT,
  IO_DIVB,
  IO_ABVC,
  IO_AMDC,
  IO_PHI,
  IO_ROTB,
  IO_SROTB,
  IO_COOLRATE,
  IO_CONDRATE,
  IO_DENN,
  IO_EGYPROM,
  IO_EGYCOLD,
  IO_MACH,
  IO_DTENERGY,
  IO_PRESHOCK_CSND,
  IO_PRESHOCK_DENSITY,
  IO_PRESHOCK_ENERGY,
  IO_PRESHOCK_XCR,
  IO_DENSITY_JUMP,
  IO_ENERGY_JUMP,
  IO_CRINJECT,
  IO_TIDALTENSORPS,  
  IO_DISTORTIONTENSORPS,  
  IO_EULERA,
  IO_EULERB,
  IO_FLOW_DETERMINANT,
  IO_PHASE_SPACE_DETERMINANT,  
  IO_ANNIHILATION_RADIATION,
  IO_STREAM_DENSITY,
  IO_EOSTEMP,
  IO_EOSXNUC,
  IO_PRESSURE,
  IO_nHII,
  IO_RADGAMMA,
  IO_LAST_CAUSTIC,
  IO_SHEET_ORIENTATION,   
  IO_INIT_DENSITY,
  IO_CAUSTIC_COUNTER,  
  IO_SHELL_INFO,  
  IO_DMHSML,                    /* for 'SUBFIND_RESHUFFLE_CATALOGUE' option */
  IO_DMDENSITY,
  IO_DMVELDISP,
  IO_DMHSML_V,                 /* for 'SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI' option */
  IO_DMDENSITY_V,
  IO_VALPHA,

  IO_LASTENTRY			/* This should be kept - it signals the end of the list */
};



/*
 * Variables for Tree
 * ------------------
 */

extern struct NODE
{
  MyFloat len;			/*!< sidelength of treenode */
  MyFloat center[3];		/*!< geometrical center of node */

#ifdef RADTRANSFER
  MyFloat stellar_mass;         /*!< mass in stars in the node*/
  MyFloat stellar_s[3];         /*!< enter of mass for the stars in the node*/
#endif

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  MyFloat maxsoft;		/*!< hold the maximum gravitational softening of particles in the 
				   node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
#endif
  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      MyFloat s[3];		/*!< center of mass of node */
      MyFloat mass;		/*!< mass of node */
      unsigned int bitflags;	/*!< flags certain node properties */
      int sibling;		/*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;		/*!< this gives the next node in case the current node needs to be opened */
      int father;		/*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;
#ifdef SCALARFIELD
  MyFloat s_dm[3];
  MyFloat mass_dm;
#endif
  int Ti_current;
}
 *Nodes_base,			/*!< points to the actual memory allocted for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				   gives the first allocated node */


extern struct extNODE
{
  MyLongDouble dp[3];
#ifdef SCALARFIELD
  MyLongDouble dp_dm[3];
  MyFloat vs_dm[3];
#endif
#ifdef FLTROUNDOFFREDUCTION
  MyFloat s_base[3];
  MyFloat len_base;
#ifdef SCALARFIELD
  MyFloat s_dm_base[3];
#endif
#endif
  MyFloat vs[3];
  MyFloat vmax;
  MyFloat divVmax;
  MyFloat hmax;			/*!< maximum SPH smoothing length in node. Only used for gas particles */
  int Ti_lastkicked;
  int Flag;
}
 *Extnodes, *Extnodes_base;


extern int MaxNodes;		/*!< maximum allowed number of internal nodes */
extern int Numnodestree;	/*!< number of (internal) nodes in each tree */


extern int *Nextnode;		/*!< gives next node in tree walk  (nodes array) */
extern int *Father;		/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
extern double Rs, R200;
extern double Dc;
extern double RhoCrit, V200;
extern double fac;
#endif


#ifdef CHEMISTRY
/* ----- chemistry part ------- */

#define H_number_fraction 0.76
#define He_number_fraction 0.06

/* ----- Tables ------- */
extern double T[N_T], J0_nu[N_nu], J_nu[N_nu], nu[N_nu];
extern double k1a[N_T], k2a[N_T], k3a[N_T], k4a[N_T], k5a[N_T], k6a[N_T], k7a[N_T], k8a[N_T], k9a[N_T],
  k10a[N_T], k11a[N_T];
extern double k12a[N_T], k13a[N_T], k14a[N_T], k15a[N_T], k16a[N_T], k17a[N_T], k18a[N_T], k19a[N_T],
  k20a[N_T], k21a[N_T];
extern double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T];
extern double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
extern double sigma24[N_nu], sigma25[N_nu], sigma26[N_nu], sigma27[N_nu], sigma28[N_nu], sigma29[N_nu],
  sigma30[N_nu], sigma31[N_nu];
#endif
#endif


#endif







