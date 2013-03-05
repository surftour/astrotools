#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define YES 1
#define NO 0 

/* GADGET Particle types: 
 * 
 *  0: gas 
 *  1: halo 
 *  2: disk 
 *  3: bulge
 *  4: star 
 *  5: black hole
 *  6: boundary
 */ 

#define PTYPE_GAS 	0
#define PTYPE_HALO	1
#define PTYPE_DISK	2
#define PTYPE_BULGE	3
#define PTYPE_STAR	4
#define PTYPE_BH	5
#define PTYPE_BOUNDARY	6

/* Various global control parameters */
#define N_CELL_DIV              2		/* number to split a tree cell in each recursion (generating an N^3 cube) */
#define MAX_DISTANCE_FROM0  100.0		/* End integration along a ray when x,y,z > this value */
//#define MAX_DISTANCE_FROM0  2.0		/* End integration along a ray when x,y,z > this value */
#define TREE_FRAC_HSML        2.0		/* Fraction of h_sml for minimum size of tree cells */
#define TREE_MIN_SIZE	     0.01		/* minimum size of a tree cell (in GADGET units) */

//#define RAY_ORIGIN_TYPE	5		/* integer ID of particle type from which rays originate */
//#define NUMBER_OF_RAYS	100		/* Number of rays from each source (will be changed to nearest perfect square) */

/* Parameters that are based on the input parameterfile of the particular simulation
(read by set_All_struct_terms - eventually should just read parameterfile of sim) */
#define SFR			YES
#define COOLING		YES
#define MOREPARAMS	YES
#define NEW_RATES   YES
#define COLD_CLOUD_TEMPERATURE		1.0e3			// Temp. of cold phase of ISM, ~10^3 K
#define SUPERNOVA_TEMPERATURE		4.0e8			// Temp. of SN, ~10^8 K
#define SUPERNOVA_FACTOR_BETA		0.1				// immediate supernova fraction (beta)
#define SUPERNOVA_EVAPORATION_FACTOR_A0 4000.0		// max. evaporation factor A0
#define STAR_FORMATION_TIMESCALE_IN_GYR 8.4			// max. star formation timescale t0*

#define OMEGA_BARYON				0.04

#define UNIT_LENGTH_IN_CM			3.085678e21		// 1.0 kpc
#define UNIT_VELOCITY_IN_CM_PER_S	1.0e5			// 1 km/sec
#define UNIT_MASS_IN_GRAMS			1.989e43		// 1.0e10 solar masses
#define MINIMUM_GAS_TEMPERATURE		1.0				// minimum gas temperature in K



/* Frequently used physical constants (in cgs) */
#define  HYDROGEN_MASSFRAC 	0.76
#define  METAL_YIELD       	0.02   /*!< effective metal yield for star formation */
#define  GAMMA         		(5.0/3)
#define  GAMMA_MINUS1  		(GAMMA-1)
#ifndef  PI
#define  PI          		3.14159265358979323846
#endif
#define  PI_INV      		(1/PI)
#define  MAX_REAL_NUMBER  	1e37
#define  MIN_REAL_NUMBER  	1e-37
#define  RNDTABLE 	 		128
#define  THIRD       		(1.0/3.0)
#define  LN2         		0.69314718
#define  GRAVITY     		6.672e-8
#define  SOLAR_MASS  		1.989e33
#define  SOLAR_LUM   		3.826e33
#define  RAD_CONST   		7.565e-15
#define  AVOGADRO    		6.0222e23
#define  BOLTZMANN   		1.3806e-16
#define  GAS_CONST   		8.31425e7
#define  C           		2.9979e10
#define  PLANCK      		6.6262e-27
#define  CM_PER_MPC  		3.085678e24
#define  PROTONMASS  		1.6726e-24
#define  ELECTRONMASS 		9.10953e-28
#define  THOMPSON     		6.65245e-25
#define  ELECTRONCHARGE 	4.8032e-10
#define  HUBBLE      		3.2407789e-18   /* in h/sec */
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7



/* Various constants of relevance to the code running */
#define KERNEL_TABLE 1000
#define MAXITER 200
/* tabulated smoothing kernel */
double  Kernel[KERNEL_TABLE+2],
        KernelDer[KERNEL_TABLE+2],
        KernelRad[KERNEL_TABLE+2];



/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes
{
  long NumPart;		/*!<  total particle numbers (global value) */
  long N_gas;		/*!<  total gas particle number (global value) */
  long N_halo, N_disk, N_bulge, N_star, N_bh, N_total;

  double Time;		/*!<   simulation time */

  /* some SPH parameters */
  int DesNumNgb;                /*!< Desired number of SPH neighbours */
  int MaxNumNgbDeviation;       /*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;			/*!< may be used to set the temperature in the IC's */
  double InitGasU;				/*!< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;			/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */


  /* system of units  */
  double UnitTime_in_s,   	/*!< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,             	/*!< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,   /*!< factor to convert intqernal velocity unit to cm/sec */
    UnitLength_in_cm,           /*!< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,        /*!< factor to convert internal pressure unit to cgs units (little 'h' still around!) */
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
		   					physical values for cooling physics */
  double BoxSize; /*!< Boxsize in case periodic boundary conditions are used */
  double BoxHalf;
  double TimeBegin;
  double Redshift;
  

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


  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[6];
  double M_halo, M_disk, M_bulge;


  /* star formation and feedback sector */
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
}
All;


/* Structures */


/* Vector for local values of the various quantities to be evaluated at 
 *  each point along the ray
 */
typedef struct LOCALVALS_s{
	float h;				    /* smoothing length h */
	float rho;				    /* density */
	float P_eff;			    /* effective pressure */
	float u;				    /* internal energy per unit mass per particle */
	float T;				    /* temperature */
	float Z;				    /* metallicity */
	float ne;				    /* electron fraction (number of electrons per H atom) */
	float neutral_frac;		    /* fraction of hydrogen in neutral atoms */
	
	float u_hot;			    /* specific energy of the hot component (2-phase ISM model) */
	float cloud_massfrac;	    /* total mass fraction in cold clouds */
	float cloud_fillingfactor;  /* volume filling factor of cold clouds */
	float rho_hot;				/* intrinsic density of the hot, diffuse component */
									/* (given as a fraction of the total density rho) */
	float rho_cold;				/* intrinsic density of the cold, cloud component */
									/* (given as a fraction of the total density rho) */
} LOCALVAL;


/* Define a structure for the lines-of-sight (rays) originating 
 *  at a black hole and extending to the 'observer'
 */
typedef struct Ray_s{
  int ORIGIN_TYPE;			/* particle type of the particle from which the ray originates */
  int ORIGIN_ID;			/* ID (particle number) for the particle which is the ray origin */
  int DEST_TYPE;			/* particle type of ray destination (if any) */ 
  int DEST_ID;				/* ID (particle number) for the particle which is the ray destination */
  float theta;				/* polar angle theta of the ray */
  float phi;				/* azimuthal angle phi of the ray */
  float n_hat[3];			/* unit vector direction of the ray */
  float pos[3];				/* current position of ray (as integrated through sim) */
  float xmax[3];			/* destination of ray */ 
  
  int N_steps;				/* number of distance steps to the edge of the box */
  float D; 				/* distance from origin destination */  
  float nh;				/* integrated gas column density along the ray */
  float nh_hot;				/* integrated gas column density of the hot phase of the ISM along the ray */
  float Z;					/* integrated mass-weighted metallicity along the ray */
  float neutral_frac;		/* integrated mass-weighted neutral fraction along the ray */
} Ray_struct;




#define N_SUBCELL   N_CELL_DIV
#define N_SUBCELL_2	N_SUBCELL*N_SUBCELL
#define N_SUBCELL_3 N_SUBCELL*N_SUBCELL*N_SUBCELL

/* Structure for the tree grid to hold local gas quantities in order to 
 *   speed up ray evaluation (and allow for scattering, radiative transport, etc.)
 */
struct CELL_STRUCT{
  float width;					/* width from one edge to another of the cell */
  float min_x[3];				/* x,y,z minima of the cell */

  int parent_ID;				/* ID number of parent cell for the given cell */
  int sub_cell_check;			/* = 1 if there are sub-cells, 0 if not */
  int sub_cell_IDs[N_SUBCELL_3];/* ID numbers of the sub-cells contained (if present) */

  float hsml_min;				/* minimum smoothing length of contained particles */
  LOCALVAL U;					/* structure to hold local variables of the cell */
  int N_particles_contained;	/* number of SPH particles contained in a cell */
  int *particle_IDs;			/* vector of particle indices of the particles contained in the cell */
} *CELL;

int TOTAL_NUMBER_OF_CELLS;
int ALL_CELL_COUNTER;




/* Header - contains general information about the snapshot file */
struct snap_header
{
  int      npart[6];		// Number of particles of each type
  double   mass[6];			// Mass of particles of each type (if = 0, then 
  							//   masses are not the same for every particle of that type)
  double   time;			// Simulation time (in Gyr)
  double   redshift;		// Redshift of cosmological sims
  int      flag_sfr;		// Flag - if star formation turned on and info. kept
  int      flag_feedback;	// Flag - if feedback of star formation turned on
  int      npartTotal[6];	// Total number of particles of each type (should be identical to npart)
  int      flag_cooling;	// Flag - if cooling physics turned on
  int      num_files;		// Indicates if snapshot split among several files (for large runs)
  double   BoxSize;			// Cosmological simulation box size (=0 for galaxy mergers)
  double   Omega0;			// (cold dark matter) 
  double   OmegaLambda;		// (vacuum energy)
  double   HubbleParam;		// little h, where H0 = h * (100 km/s/Mpc)
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes (to ensure the rest of 
  																the snapshot is read properly even if the header changes) */
} header;


struct particle_gas
{
  int 	 ID;
  float  pos[3];
  float  vel[3];
  float  mass;
  float  rho;		// density (problem units) (includes both hydrogen and helium)
  float  u;		// energy density (problem units)
  float  temp; 		// temperature (must be calculated from above)
  float  ne;		// number of electrons per hydrogen atom/nucleus
  float  nh;		// fraction of hydrogen atoms which are neutral
  float  hsml;		// local smoothing length
  float  sfr;		// star formation rate within the particle/gas
  float  z;		// metallicity abundance (mass fraction)
/* additional quantities for radiative transfer calc */ 
  int    *PTYPE; 	/* array of destination particle types [BH or Star] */ 
  int    *PID; 		/* array of particle IDs */ 
  float  *NH; 		/* NH value for each star particle */ 
  float  *D; 		/* distance to each star particle */ 
} *PG;

struct particle_halo 
{
  int 	 ID;
  float  pos[3];
  float  vel[3];
} *PH;

struct particle_disk 
{
  int 	 ID;
  float  pos[3];
  float  vel[3];
} *PD;

struct particle_bulge
{
  int 	 ID;
  float  pos[3];
  float  vel[3];
} *PB;

struct particle_star
{
  int 	 ID;
  float  pos[3];
  float  sig[3];
  float  vel[3];
  float  mass;
  float  age;		// age of stars (simulation time since creation of particle)
  float  z;			// metallicity (mass fraction)
} *PS;

struct particle_bh 
{
  int 	 ID;
  float  pos[3];
  float  vel[3];
  float  mass;
  float  bhar;		// accretion rate (changes rapidly at any given time)
} *PBH;



/*----------------------------------------------------------------------------*/
/* general purpose macros and definitions (never modified) */
#ifndef MIN
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif
#define SIGN(a) ( ((a) < 0.) ? -1. : 1. )
#define SQR(x) ( (x)*(x) )
#define CUBE(x) ( (x)*(x)*(x) )
#define STR(x) #x
#define AND &&
#define OR ||
#define SQRT2 1.4142135623730951
#define ONE_OVER_SQRT2 0.7071067811865475
#ifndef PI
#define PI       3.14159265358979323846
#endif
#define ONE_3RD  0.3333333333333333
#define TWO_3RDS 0.6666666666666667
#define TINY_NUMBER 1.0e-20
#define HUGE_NUMBER 1.0e+20
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
/*----------------------------------------------------------------------------*/
/* DMAX,DMIN,IMAX,IMIN defs */
#define DMAX(a,b) (dmax1=(a),dmax2=(b),(dmax1>dmax2)?dmax1:dmax2)
#define DMIN(a,b) (dmin1=(a),dmin2=(b),(dmin1<dmin2)?dmin1:dmin2)
#define IMAX(a,b) (imax1=(a),imax2=(b),(imax1>imax2)?imax1:imax2)
#define IMIN(a,b) (imin1=(a),imin2=(b),(imin1<imin2)?imin1:imin2)


