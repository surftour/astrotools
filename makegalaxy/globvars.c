
#include "globvars.h"



/***********  INPUT PARAMETERS *********/

double CC;			/* halo concentration */
double V200;			/* circular velocity v_200 */
double LAMBDA;			/* spin parameter  */
double MD;			/* disk mass fraction */
double MBH;			/* central black hole mass fraction */
double JD;			/* disk spin fraction */
double MB;			/* bulge mass fraction */
double GasFraction;
double DiskHeight;
char DiskPopMode[50];
double DiskPopAge;
double DiskPopTau;
double BulgeSize;
char BulgePopMode[50];
double BulgePopAge;
double BulgePopTau;

double REDSHIFT;		/* redshift to scale to */
double Omega_m0;		/* Omega_m(z=0) */
double Omega_L0;		/* Omega_L(z=0) */
double Hz;			/* Hubble parameter at z=REDSHIFT */

int N_HALO;			/* desired number of particles in halo */
int N_DISK;			/* desired number of collsionless particles in disk */
int N_GAS;			/* number of gas particles in stellar disk */
int N_BULGE;			/* number of gas particles in stellar disk */
int N_BLACKHOLE;


double  AnisotropyRadius; /* specifies the anisotropy radius, in units of RH */

int GasDistribution;
double GasExpAlpha;
double PowerLawGamma;
double PowerLawCutoff;
double rb;
double Zmet0;            /* gas metallicity central value (in metal fraction) */
double ZmetGrad;         /* gas metallicity gradient (in dex/kpc) */



double GravSoftening;     /* used for the force calculation */
int WriteDensity;         /* 1=Write density in mass field (used for
			     Arepo) */
int UseQrng;              /* Use quasi-random generators instead of
			     pseudorandom. */
 
/*********************************************/

double NumGasScaleLengths;


int NumPart;			/* this is the particle number used for the tree */
struct part_data *P;

double ErrTolTheta;

double *Zrho; //, *Zcumul;

double MaxGasDiskHeight;

double M200;			/* virial mass */

double RH;			/* scale radius in Hernquist profile */

double U4;

double RS;			/* scale radius for halo */
double R200;			/* virial radius */
double H;			/* disk scale length */
double Z0;			/* disk thickness */
double A;			/* bulge scale radius */


int N_TOTAL;                    /* total number of particle */

double M_HALO;			/* total dark mass */
double M_DISK;			/* mass of stellar disk (collisionless part) */
double M_GAS;			/* gas mass in disk */
double M_BULGE;			/* mass of bulge */
double M_BLACKHOLE;		/* mass of black-hole */

double halo_spinfactor;		/* computed streamin of dark matter */


double RadialDispersionFactor;

double MinGasTemp;

char OutputDir[500], OutputFile[500];

double G;			/* gravitational constant */
double H0;			/* Hubble constant */
double UnitTime_in_s;
double UnitMass_in_g;
double UnitLength_in_cm;
double UnitVelocity_in_cm_per_s;
double UnitTime_in_Megayears;

double UnitPressure_in_cgs, UnitDensity_in_cgs, UnitCoolingRate_in_cgs, UnitEnergy_in_cgs;


/* particle data */

double *vmax2_halo, *vmax2_disk, *vmax2_bulge, *vmax2_gas;

double *xp_halo, *yp_halo, *zp_halo, *mp_halo;
double *xp_disk, *yp_disk, *zp_disk, *mp_disk;
double *xp_bulge, *yp_bulge, *zp_bulge, *mp_bulge;
double *xp_gas, *yp_gas, *zp_gas, *mp_gas, *rho_gas, *u_gas;

double *vxp_halo, *vyp_halo, *vzp_halo;
double *vxp_disk, *vyp_disk, *vzp_disk;
double *vxp_bulge, *vyp_bulge, *vzp_bulge;
double *vxp_gas, *vyp_gas, *vzp_gas;





double **RhoGas, **CumulMassGas;

double rhocentral[TABSIZE];
double sigmatab[TABSIZE];

double FactorEVP;
double PhysDensThresh;
double EgySpecSN;
double MaxSfrTimescale;
double FactorSN;
double TempSupernova;
double TempClouds;
double EgySpecCold;
double EgySpecSN;
double FeedbackEnergy;
double FactorForSofterEQS;


double Baselen;


double LL, HR;			/* LL = extension of fields in R and z.
				   HR = extension of high resolution region in z */
double dR;			/* delta R */




double **phi;                           /* the total potential */
double PhiMin;                          /* minimum potential */

double **Dphi_z, **Dphi_R, **Dphi_z_dR;	/* derivatives of total potential */
double *epi_gamma2, *epi_kappa2;	/* epicycle gamma^2  */

double **VelStreamGas;

/* halo velocity fields */

double **VelDispRz_halo;
double **VelDispPhi_halo;
double **VelVc2_halo;
double **VelStreamPhi_halo;
double **VelDispRz_dR_halo;


/* bulge velocity fields */

double **VelDispRz_bulge;
double **VelDispPhi_bulge;
double **VelVc2_bulge;
double **VelStreamPhi_bulge;
double **VelDispRz_dR_bulge;


/* disk velocity fields */

double **VelDispRz_disk;
double **VelDispPhi_disk;
double **VelVc2_disk;
double **VelStreamPhi_disk;
double **VelDispRz_dR_disk;




/* auxiliary field */

double *xl, *yl, *D2yl;
double *list_z, *list_R, *list_RplusdR, *list_RminusdR;
