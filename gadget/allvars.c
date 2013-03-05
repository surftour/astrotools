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

#include "allvars.h"



#ifdef PERIODIC
MyDouble boxSize, boxHalf;

#ifdef LONG_X
MyDouble boxSize_X, boxHalf_X;
#else
#endif
#ifdef LONG_Y
MyDouble boxSize_Y, boxHalf_Y;
#else
#endif
#ifdef LONG_Z
MyDouble boxSize_Z, boxHalf_Z;
#else
#endif
#endif

#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
MPI_Status mpistat;
#endif

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/



int ThisTask;			/*!< the number of the local processor  */
int NTask;			/*!< number of processors */
int PTask;			/*!< note: NTask = 2^PTask */

#ifdef INVARIANCETEST
int World_ThisTask;
int World_NTask;
int Color;
MPI_Comm MPI_CommLocal;
#endif

double CPUThisRun;		/*!< Sums CPU time of current process */

int NumForceUpdate;		/*!< number of active particles on local processor in current timestep  */
long long GlobNumForceUpdate;
int NumSphUpdate;		/*!< number of active SPH particles on local processor in current timestep  */

int MaxTopNodes;		/*!< Maximum number of nodes in the top-level tree used for domain decomposition */

int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */
int RestartSnapNum;

int *Exportflag;		/*!< Buffer used for flagging whether a particle needs to be exported to another process */
int *Exportnodecount;
int *Exportindex;

int *Send_offset, *Send_count, *Recv_count, *Recv_offset, *Sendcount_matrix;

int FirstActiveParticle;
int *NextActiveParticle;
unsigned char *ProcessedFlag;

int TimeBinCount[TIMEBINS];
int TimeBinCountSph[TIMEBINS];
int TimeBinActive[TIMEBINS];

int FirstInTimeBin[TIMEBINS];
int LastInTimeBin[TIMEBINS];
int *NextInTimeBin;
int *PrevInTimeBin;

#ifdef SFR
double TimeBinSfr[TIMEBINS];
#endif

#ifdef GASRETURN
double TimeBinGasReturn[TIMEBINS];
#endif

#ifdef BLACK_HOLES
double TimeBin_BH_mass[TIMEBINS];
double TimeBin_BH_dynamicalmass[TIMEBINS];
double TimeBin_BH_Mdot[TIMEBINS];
double TimeBin_BH_Medd[TIMEBINS];
#endif

char DumpFlag = 1;

size_t AllocatedBytes;
size_t FreeBytes;

double CPU_Step[CPU_PARTS];
char CPU_Symbol[CPU_PARTS] =
  { '-', '*', '=', ';', '<', '[', '^', ':', '.', '~', '|', '+', '"', '/', '`', ',', '>', '@', '#', '&', '$',
  ']', '(', '?', ')', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '\\', '\%', '{', '}'
};
char CPU_SymbolImbalance[CPU_PARTS] =
  { 'a', 't', 'u', 'v', 'b', 'w', 'd', 'r', 'h', 'm', 'n', 'l', 'o', 'p', 's', 'f', 'i', 'g', 'c', 'e', 'x',
  'y', 'z', 'A', 'I', 'W', 'T', 'V', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L'
};
char CPU_String[CPU_STRING_LEN + 1];

double WallclockTime;		/*!< This holds the last wallclock time measurement for timings measurements */

int Flag_FullStep;		/*!< Flag used to signal that the current step involves all particles */


int TreeReconstructFlag;
int GlobFlag;


int NumPart;			/*!< number of particles on the LOCAL processor */
int N_gas;			/*!< number of gas particles on the LOCAL processor  */
long long Ntype[6];		/*!< total number of particles of each type */
int NtypeLocal[6];		/*!< local number of particles of each type */

gsl_rng *random_generator;	/*!< the random number generator used */


#ifdef SFR
int Stars_converted;		/*!< current number of star particles in gas particle block */
#endif

#ifdef GASRETURN
int Gas_converted;		/*!< current number of gas particles in the star particle block */
#endif


double TimeOfLastTreeConstruction;	/*!< holds what it says */

int *Ngblist;			/*!< Buffer to hold indices of neighbours retrieved by the neighbour search
				   routines */
double *R2ngblist;

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
int *DomainStartList, *DomainEndList;



double *DomainWork;
int *DomainCount;
int *DomainCountSph;
int *DomainTask;
int *DomainNodeIndex;
int *DomainList, DomainNumChanged;

peanokey *Key, *KeySorted;

struct topnode_data *TopNodes;

int NTopnodes, NTopleaves;




double RndTable[RNDTABLE];

#ifdef SUBFIND
int GrNr;
int NumPartGroup;
#endif

/* variables for input/output , usually only used on process 0 */


char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

FILE *FdInfo,			/*!< file handle for info.txt log-file. */
 *FdEnergy,			/*!< file handle for energy.txt log-file. */
 *FdTimings,			/*!< file handle for timings.txt log-file. */
 *FdBalance,			/*!< file handle for balance.txt log-file. */
 *FdCPU;			/*!< file handle for cpu.txt log-file. */

#ifdef SFR
FILE *FdSfr;			/*!< file handle for sfr.txt log-file. */
#endif

#ifdef GASRETURN
FILE *FdGasReturn;		/*!< file handle for GasReturn.txt log-file */
#endif

#ifdef RADTRANSFER
FILE *FdEddington;		/*!< file handle for eddington.txt log-file. */
FILE *FdRadtransfer;		/*!< file handle for radtransfer.txt log-file. */
#endif

#ifdef DISTORTIONTENSORPS
#ifdef PMGRID
FILE *FdTidaltensor;		/*!< file handle for Tidaltensor.txt log-file. */
#endif
FILE *FdCaustics;		/*!< file handle for Caustics.txt log-file. */
#endif

#ifdef BLACK_HOLES
FILE *FdBlackHoles;		/*!< file handle for blackholes.txt log-file. */
FILE *FdBlackHolesDetails;
#endif


#ifdef FORCETEST
FILE *FdForceTest;		/*!< file handle for forcetest.txt log-file. */
#endif

#ifdef DARKENERGY
FILE *FdDE;			/*!< file handle for darkenergy.txt log-file. */
#endif

#ifdef XXLINFO
FILE *FdXXL;			/*!< file handle for xxl.txt log-file. */

#ifdef MAGNETIC
double MeanB;

#ifdef TRACEDIVB
double MaxDivB;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
double MeanAlpha;
#endif
#endif

/*! table for the cosmological drift factors */
double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
double HydroKickTable[DRIFT_TABLE_LENGTH];


void *CommBuffer;		/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes All;




/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P,	/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */



/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data *SphP,	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */


peanokey *DomainKeyBuf;

/* global state of system 
*/
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */

struct data_index *DataIndexTable;	/*!< the particles to be exported are grouped
					   by task-number. This table allows the
					   results to be disentangled again and to be
					   assigned to the correct particle */

struct data_nodelist *DataNodeList;

struct gravdata_in *GravDataIn,	/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


struct gravdata_out *GravDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


struct potdata_out *PotDataResult,	/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */




/*! Header for the standard file format.
 */
struct io_header header;	/*!< holds header for snapshot files */





/*
 * Variables for Tree
 * ------------------
 */

struct NODE *Nodes_base,	/*!< points to the actual memory allocted for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				   gives the first allocated node */


struct extNODE *Extnodes, *Extnodes_base;


int MaxNodes;			/*!< maximum allowed number of internal nodes */
int Numnodestree;		/*!< number of (internal) nodes in each tree */


int *Nextnode;			/*!< gives next node in tree walk  (nodes array) */
int *Father;			/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
double Rs, R200;
double Dc;
double RhoCrit, V200;
double fac;
#endif


#ifdef CHEMISTRY
/* ----- chemistry part ------- */

#define H_number_fraction 0.76
#define He_number_fraction 0.06

/* ----- Tables ------- */
double T[N_T], J0_nu[N_nu], J_nu[N_nu], nu[N_nu];
double k1a[N_T], k2a[N_T], k3a[N_T], k4a[N_T], k5a[N_T], k6a[N_T], k7a[N_T], k8a[N_T], k9a[N_T],
  k10a[N_T], k11a[N_T];
double k12a[N_T], k13a[N_T], k14a[N_T], k15a[N_T], k16a[N_T], k17a[N_T], k18a[N_T], k19a[N_T],
  k20a[N_T], k21a[N_T];
double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T];
double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
double sigma24[N_nu], sigma25[N_nu], sigma26[N_nu], sigma27[N_nu], sigma28[N_nu], sigma29[N_nu],
  sigma30[N_nu], sigma31[N_nu];
#endif
#endif
