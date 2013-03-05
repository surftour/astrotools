/* This file declares all global variables. Further variables should be added here, 
   and declared as 'extern'. The actual existence of these variables is provided by
   the file 'allvars.c'. To produce 'allvars.c' from 'allvars.h', do the following:

   1.) Erase all #define's
   2.) add #include "allvars.h"
   3.) delete all keywords 'extern'
   4.) delete all struct definitions enclosed in {...}, e.g.
       "extern struct global_data_all_processes {....} All;"
       becomes "struct global_data_all_processes All;"
   5.) delete extern struct xyz_data {...}; completely
*/

#include <stdio.h>


#ifdef T3E
  typedef short int int4byte;   /* Note: int has 8 Bytes on the T3E ! */
#else
  typedef int int4byte;
#endif



#ifdef  GRAPE

#define MAX_GRCHIPS  48
#define MAX_GRNBR    4096                 

#endif


/* ... often used constants (cgs units) */

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37
#define  THIRD            (1.0/3.0)
#ifndef  PI
#define  PI               3.14159265358979323846
#endif
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

#define  HUBBLE      3.2407789e-18   /* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define  GAMMA         (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)


#define  KERNEL_TABLE 1000

#define  MAX_NGB  100000  /* defines maximum length of neighbour list */

#define  MAXLEN_OUTPUTLIST 350 /* maxmimum number of entries in output list*/

#define  TIMESTEP_INCREASE_FACTOR 1.3


extern int    NumForceUpdate, NumSphUpdate, IndFirstUpdate;
extern int    TimeTreeRoot;

extern int    Num_nodeupdates, Num_nodeupdate_particles;

extern int    RestartFlag;

extern int    NumPart,        /* Note: these are the LOCAL process values */
              N_gas;
extern int    NumPart_2d, N_gas_2d;


/* variables for input/output ,  usually only used on process 0 */

extern  char   ParameterFile[100];
extern  FILE  *FdInfo,
              *FdEnergy,
              *FdTimings,
              *FdCPU;


/* tabulated smoothing kernel */

extern double  Kernel[KERNEL_TABLE+2],
               KernelDer[KERNEL_TABLE+2],
               KernelRad[KERNEL_TABLE+2];


extern double  CPUThisRun;


extern struct global_data_all_processes  /* this struct contains data which is the same for all tasks 
                                            (mostly code parameters read from the parameter file) */
{

  int   TotNumPart;         /*  particle numbers  */
  int	TotNumPart_2d;
	int 	TotN_gas;
	int 	TotN_gas_2d;
	int	TotN_halo;
	int	TotN_halo_2d;
	int	TotN_disk;
	int 	TotN_disk_2d;
	int	TotN_bulge;
	int	TotN_bulge_2d;
	int	TotN_stars;
	int	TotN_stars_2d;


  int   MaxPart,        /* These numbers give the maxmimum number of particles that can be hold */
        MaxPartSph;     /* on the current processor. (Set it to ~2-3 times the average load)    */

  int   ICFormat;


  double PartAllocFactor;  /* in order to maintain work-load balance,  
			      the particle load is usually NOT balanced 
			      each processor allocates memory for PartAllocFactor times  
			      the average number of particles */
                           
  double TreeAllocFactor;  /* similarly for the tree:
			      each processor allocates a number of nodes which is 
			      TreeAllocFactor times the maximum(!) number of particles 
			      Note, that a typical local tree for N particles needs typically
			      1.5-2 N nodes */
 
 /* some SPH parameters */

  int    DesNumNgb;
  double ArtBulkViscConst;
  double InitGasTemp;     /* may be used to set the temperature in the IC's */   
  double MinGasTemp;      /* may be used to set a floor for the gas temperature */
  double MinEgySpec;


  /* diagnostics */

  int    TotNumOfForces;   /* counts total number of force computations  */
  

  /* system of units  */

  double UnitTime_in_s,
         UnitMass_in_g,
         UnitVelocity_in_cm_per_s,
         UnitLength_in_cm,
         UnitPressure_in_cgs,
         UnitDensity_in_cgs,
         UnitCoolingRate_in_cgs,
         UnitEnergy_in_cgs,
         UnitTime_in_Megayears,
         GravityConstantInternal,
         G;

  /* Cosmology */

  double Hubble;
  double BoxSize, BoxHalf;
  double Omega0,        
         OmegaLambda,
         OmegaBaryon,
         HubbleParam;   /* little `h', i.e. Hubble constant in units of 100 km/s/Mpc. 
                         * Only needed to get absolute physical values 
			 * for cooling physics 
			 */    

  /* Code options */

  int    ComovingIntegrationOn;   /* enables comoving integration */
  int    PeriodicBoundariesOn;
  int    ResubmitOn;
  int    TypeOfOpeningCriterion;
  int    TypeOfTimestepCriterion;
  int    OutputListOn;
  int    CoolingOn;


  /* parameters determining output frequency */

  int    SnapshotFileCount;
  double TimeBetSnapshot,
         TimeOfFirstSnapshot,
         CpuTimeBetRestartFile,
         TimeLastRestartFile,
         TimeBetStatistics,
         TimeLastStatistics;

  /* Current time of the simulation */

  int     NumCurrentTiStep;

  double  Time, 
          TimeStep,
          TimeBegin,
          TimeMax;   /* marks end of the simulation */
  


  /* variables that keep track of cumulative CPU consumption */

  double  TimeLimitCPU;
  double  CPU_TreeConstruction;
  double  CPU_TreeWalk;
  double  CPU_Gravity;
  double  CPU_Potential;
  double  CPU_Snapshot;
  double  CPU_Total;
  double  CPU_Hydro;
  double  CPU_Predict;
  double  CPU_TimeLine;


  /* tree code opening criterion */
  double  ErrTolTheta;
  double  ErrTolForceAcc;


 /* adjusts accuracy of time-integration */

  double  ErrTolIntAccuracy;  /* for 1/a^{1/2} collisionless timestep criterion */
  double  ErrTolVelScale;     /* for 1/a  collisionless timestep criterion */
  double  MinSizeTimestep,
          MaxSizeTimestep;

  double  CourantFac;      /* SPH-Courant factor */
 

  /* frequency of tree reconstruction */

  double  MaxNodeMove;
  double  TreeUpdateFrequency;
  int     NumForcesSinceLastTreeConstruction;


 /* gravitational and hydrodynamical softening lengths 
   * (given in terms of an `equivalent' Plummer softening length) 
   *
   * five groups of particles are supported 
   * 0=gas,1=halo,2=disk,3=bulge,4=stars 
   */
  double  MinGasHsmlFractional, MinGasHsml;


  double  SofteningGas,
          SofteningHalo,
          SofteningDisk,
          SofteningBulge,
          SofteningStars;

  double  SofteningGasMaxPhys,
          SofteningHaloMaxPhys,
          SofteningDiskMaxPhys,
          SofteningBulgeMaxPhys,
          SofteningStarsMaxPhys;



  double  SofteningTable[6];    
  double  SofteningTableMaxPhys[6];


  /* If particle masses are all equal for one type, 
   *  the corresponding entry in MassTable is set to this value,
   * allowing the size of the snapshot files to be reduced
   */
  double  MassTable[6];   
  
  char    InitCondFile[100],
          OutputDir[100],
          SnapshotFileBase[100],
          EnergyFile[100],
          InfoFile[100],
          TimingsFile[100],
          CpuFile[100],
          RestartFile[100],
          ResubmitCommand[100],
          OutputListFilename[100];

  double  OutputListTimes[MAXLEN_OUTPUTLIST]; /* was 200 in earlier version */
  int     OutputListLength;


#ifdef GRAPE  /* GRAPE stuff */

  int     GrapeAquired,
          FreeGrape;

  /* GRAPE hardware specifications */

  int   GrNumChips, 
        GrNumChipsPerBoard,
        GrNumBoards,
        GrMaxPart,
        GrMaxNbr;
#endif

} All;




/* The following structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data 
{
  float     Pos[3];       /* particle position at its current time */  
  float     Vel[3];       /* particle velocity at its current time */  
  float     Mass;         /* particle mass */   
  int4byte  ID;           /* unique particle identifier */  
  int4byte  Type;         /* flags particle type. 0=gas, 1=halo, 2=disk, 3=bulge, 4=stars */

  float     CurrentTime;  /* current time of the particle */
  float     MaxPredTime;  /* current time plus half the particles allowed timestep */
  float     PosPred[3];   /* particle position at the global prediction time */   
  float     VelPred[3];   /* particle velocity at the global prediction time */

  float     Accel[3];     /* particle acceleration */
  float     Potential;    /* particle potential */

  float     OldAcc;       /* magnitude of old force. Used in new relative opening criterion */

  int4byte  ForceFlag;    /* points to next active particle */

#ifdef VELDISP  
  float     VelDisp;
  float     HsmlVelDisp;
  float     DensVelDisp;
#endif
} *P,*P_data, *Pn, *Pn_data;



/* the following struture holds data that is stored for each SPH particle
 * in addition to the collisionless variables.
 */
extern struct sph_particle_data
{
  float  Density;         /* particle density at its current time */ 
  float  DtDensity;       /* rate of change of density */
  float  DensityPred;     /* predicted particle density */

  float  EgySpec;         /* internal energy per unit mass */
  float  DtEgySpec;       /* rate of change of the internal energy */
  float  EgySpecPred;     /* predicted internal energy per unit mass */

  float  Pressure;        /* pressure */

  float  Hsml;            /* smoothing length */
  float  DtHsml;          /* rate of change of smoothing length */

  int    NumNgb;          /* number of SPH neighbours */

  float  DivVel;          /* local velocity divergence */
  float  CurlVel;         /* local velocity curl */

#ifdef COOLING
  float  Ne;              /* electron fraction. Gives indirectly ionization state 
                             and mean molecular weight. */
#endif

} *SphP,*SphP_data, *SphPn, *SphPn_data;



/* this structure holds nodes for the ordered binary tree of the timeline.
 */
extern struct timetree_data
{
  int4byte left,right;
} *PTimeTree;   /* for ordered binary tree of max pred. times */



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

          MassComp[5],
          EnergyKinComp[5],
          EnergyPotComp[5],
          EnergyIntComp[5],
          EnergyTotComp[5],
          MomentumComp[5][4],
          AngMomentumComp[5][4],
          CenterOfMassComp[5][4];

}   SysState,SysStateAtStart,SysStateAtEnd;




/* Header for the standard file format.
 */
extern struct io_header_1
{
  int4byte npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int4byte flag_sfr;
  int4byte flag_feedback;
  int4byte npartTotal[6];
  int4byte flag_cooling;
  int4byte num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;









