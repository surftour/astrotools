#include "overhead.h"
#include "proto.h"

#define DEBUG NO

void set_All_struct_terms(void)
{
/* Sets the units and constants used by the "All.xxxxx" variables throughout */
  double meanweight;
  extern struct global_data_all_processes All;

  /* actual unit conversions used throughout */
  All.UnitLength_in_cm			=	UNIT_LENGTH_IN_CM;
  All.UnitVelocity_in_cm_per_s	=	UNIT_VELOCITY_IN_CM_PER_S;
  All.UnitMass_in_g				=	UNIT_MASS_IN_GRAMS;

  All.UnitTime_in_s 			=	All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears 	= 	All.UnitTime_in_s / SEC_PER_MEGAYEAR;
  All.G 						=   GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  All.UnitDensity_in_cgs 		=   All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs 		=   All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs 	=   All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs 		=   All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */
  meanweight 		= 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */
  All.MinGasTemp 	= MINIMUM_GAS_TEMPERATURE;
  All.MinEgySpec 	= 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec   *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  /* star formation / two-phase medium parameters */
  All.TempClouds 		= COLD_CLOUD_TEMPERATURE;
  All.TempSupernova 	= SUPERNOVA_TEMPERATURE;
  All.FactorSN        	= SUPERNOVA_FACTOR_BETA;
  All.FactorEVP			= SUPERNOVA_EVAPORATION_FACTOR_A0;
  All.MaxSfrTimescale	= STAR_FORMATION_TIMESCALE_IN_GYR * (1000.0/All.UnitTime_in_Megayears);

  /* cosmology */
  All.CritOverDensity = 0.0;
  All.CritPhysDensity = 0.0;
  All.OmegaBaryon     = OMEGA_BARYON;
  All.HubbleParam     = header.HubbleParam;		// Defined by readsnap.c when reading in snapshots
  All.Hubble          = All.HubbleParam * 1.0e7 / CM_PER_MPC * All.UnitTime_in_s;
  All.ComovingIntegrationOn = 0;	// Turn off comoving integration terms
  
  /* set the rest of the units and Paul Martini's versions of the CGS units */
  set_units_sfr();
  All.PhysDensThresh = 0.0;
  init_clouds();
}

void set_units_sfr(void)
{
  double meanweight, feedbackenergyinergs;

  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  All.OverDensThresh =
    All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

  All.FeedbackEnergy = All.FactorSN / (1 - All.FactorSN) * All.EgySpecSN;

  feedbackenergyinergs = All.FeedbackEnergy / All.UnitMass_in_g * (All.UnitEnergy_in_cgs * SOLAR_MASS);

  if(DEBUG)
    {
      printf("Feedback energy per formed solar mass in stars= %g  ergs\n", feedbackenergyinergs);
      printf("OverDensThresh= %g\nPhysDensThresh= %g (internal units)\n", All.OverDensThresh,
	     All.PhysDensThresh);
    }
}


void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;
  double tsfr, y, peff, fac, neff, egyeff, factorEVP, sigma, thresholdStarburst;

  if(All.PhysDensThresh == 0)
    {
      A0 = All.FactorEVP;

      egyhot = All.EgySpecSN / A0;

      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

      u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
      u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;


      if(All.ComovingIntegrationOn)
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
      else
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
	//printf("M_PI=%g \n",M_PI); /* machine value of Pi */

      if(All.ComovingIntegrationOn)
	{
	  All.Time = 1.0;	/* to be guaranteed to get z=0 rate */
	  IonizeParams();
	}

      ne = 1.0;
      SetZeroIonization();
      tcool = GetCoolingTime(egyhot, dens, &ne);
	  //printf("A0, egyhot, u4, ne, tcool, All.FactorSN, All.EgySpecSN, All.EgySpecCold \n %e  %e  %e  %e  %e  %e %e %e \n",A0, egyhot, u4, ne, tcool, All.FactorSN, All.EgySpecSN, All.EgySpecCold);
      coolrate = egyhot / tcool / dens;

      x = (egyhot - u4) / (egyhot - All.EgySpecCold);
	  //x = (egyhot + ((A0+1.0)/A0)*(All.EgySpecCold - u4)) / egyhot;

      All.PhysDensThresh =
		x / pow(1 - x,2) * (All.FactorSN * All.EgySpecSN - 
		(1 - All.FactorSN) * All.EgySpecCold) / (All.MaxSfrTimescale * coolrate);

      if(DEBUG)
	{
	  printf("\nA0= %g  \n", A0);
	  printf("Computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n", All.PhysDensThresh,
		 All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs));
	  float od = log10(All.PhysDensThresh/(All.OmegaBaryon*dens/1.0e6));
	  printf("Corresponding Critical Overdensity= 10^(%3.2f) \n",od);
	  printf("EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g\n\n", x);
	  printf("tcool=%g dens=%g egyhot=%g\n", tcool, dens, egyhot);
	}

      dens = All.PhysDensThresh * 10;

      do
	{
	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
	  tcool = GetCoolingTime(egyhot, dens, &ne);

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  peff = GAMMA_MINUS1 * dens * egyeff;

	  fac = 1 / (log(dens * 1.025) - log(dens));
	  dens *= 1.025;

	  neff = -log(peff) * fac;

	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
	  tcool = GetCoolingTime(egyhot, dens, &ne);

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  peff = GAMMA_MINUS1 * dens * egyeff;

	  neff += log(peff) * fac;
	}
      while(neff > 4.0 / 3);

      thresholdStarburst = dens;


      if(DEBUG)
	{
	  printf("Run-away sets in for dens=%g\n", thresholdStarburst);
	  printf("Dynamic range for quiescent star formation= %g\n", thresholdStarburst / All.PhysDensThresh);
	  fflush(stdout);
	}

      //integrate_sfr();	/* (not necessary for just initializing rho_threshold) */

      if(DEBUG)
	{
	  sigma = 10.0 / All.Hubble * 1.0e-10 / pow(1.0e-3, 2);

	  printf("Isotherm sheet central density: %g   z0=%g\n",
		 M_PI * All.G * sigma * sigma / (2 * GAMMA_MINUS1) / u4,
		 GAMMA_MINUS1 * u4 / (2 * M_PI * All.G * sigma));
	  fflush(stdout);

	}

      if(All.ComovingIntegrationOn)
	{
	  All.Time = All.TimeBegin;
	  IonizeParams();
	}
    }
}




/* The following routines are all directly lifted from the gadget 
 *   codes. They are unchanged and noted with their original parent file.
 *   They're just removed and put in here to cut down on clutter and 
 *   make the set of codes more manageable (a lot of gadget routines
 *   just aren't needed at this stage)
 */

/*  FROM begrun.c  */
/*******************/
/* Here the lookup table for the kernel of the SPH calculation
 * is initialized. 
 */ 
void set_sph_kernel(void)
{
  int i;

  for(i=0;i<=KERNEL_TABLE+1;i++)
    KernelRad[i] = ((double)i)/KERNEL_TABLE;

  Kernel[KERNEL_TABLE+1] = KernelDer[KERNEL_TABLE+1]= 0;

      
  for(i=0;i<=KERNEL_TABLE;i++)
    {
      if(KernelRad[i]<=0.5)
	{
	  Kernel[i] = 8/PI *(1-6*KernelRad[i]*KernelRad[i]*(1-KernelRad[i]));
	  KernelDer[i] = 8/PI *( -12*KernelRad[i] + 18*KernelRad[i]*KernelRad[i]);
	}
      else
	{
	  Kernel[i] = 8/PI * 2*(1-KernelRad[i])*(1-KernelRad[i])*(1-KernelRad[i]);
	  KernelDer[i] = 8/PI *( -6*(1-KernelRad[i])*(1-KernelRad[i]));
	}
    }
}
