#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"allvars.h"
#include"prototypes.h"
#include"cooling.h"
#include"proto.h"


struct global_data_all_processes All;





int multiphase_info(int argc, void *argv[])
{
	int N;
	float *rho;
	float *nel;
	int i,comoving =0;
	float time  = 0.0;
	float TSN   = 3.0e8;
	float fSN   = 0.1;
	float fEVP  = 3000;
	float hubble = 0.7;
	float OmegaBaryon = 0.0;
	float tstar = 4.5;

	float ne   = 1.0;//144263;
	float x;
	float *xx;
	float *xv;
	float *xz;
	float *cfraction;
	float *vfraction;
	float *htemp;


        if(argc!=13)
        {
           fprintf(stderr,"\n\nwrong number of arguments ! (found %d) \n\n",argc);
           exit(0);
        }


	N        = *(int *)argv[0];
	tstar    = *(float *)argv[1];
	TSN      = *(float *)argv[2];
	fEVP     = *(float *)argv[3];
	fSN      = *(float *)argv[4];
	hubble   = *(float *)argv[5];
	time     = *(float *)argv[6];
	comoving = *(float *)argv[7];
	rho       = (float *)argv[8];
	nel       = (float *)argv[9];
	cfraction = (float *)argv[10];
	vfraction = (float *)argv[11];
	htemp     = (float *)argv[12]; 


	/*if(argc==2)
	{
		rho = atof(argv[1]);
	}else{
		if(argc>2)
		{
			printf("Wrong arguments\n");
			printf("./cold_fraction [rho]\n");
			exit(-1);
		}
	}*/

	if(comoving)
	{
		All.ComovingIntegrationOn = 1;
	}else{
		All.ComovingIntegrationOn = 0;
	}
	All.HubbleParam = hubble;
	All.OmegaBaryon = OmegaBaryon;
	All.Time = time;
	All.FactorEVP = fEVP;
	All.FactorSN  = fSN;
	All.TempSupernova   = TSN;
	All.MaxSfrTimescale = tstar;
	set_units();
	init_clouds();

	//printf("x %e\n",get_cold_fraction(rho/(PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs),ne));
	//printf("Mass Fraction:   %e\n",get_cold_mass_fraction(rho));
	printf("calculating multiphase properties ....");
	fflush(0);
	for(i=0;i<N;i++)
	{
		x = get_cold_mass_fraction(rho[i],nel[i]);
		xx = cfraction + i;
		*xx = x;
		x = get_cold_volume_fraction(rho[i],nel[i]);
		xv = vfraction + i;
		*xv = x;
		x = get_hotphase_temperature(rho[i],nel[i]);
		xz = htemp + i;
		*xz = x;
	}
	//printf("Volume Fraction: %e\n",get_cold_volume_fraction(rho));
	printf("done!\n");
}





/* -------------
   Set Units
   ------------- */
void set_units(void)
{
  double meanweight;

  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.UnitMass_in_g = 1.989e+43;
  All.UnitVelocity_in_cm_per_s = 100000;
  All.UnitLength_in_cm = 3.08568e+21;
  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;
  All.TempClouds    = 1000.0;
  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);


  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
  All.Hubble      = HUBBLE * All.UnitTime_in_s;

  /*All.OverDensThresh =
    All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);*/
  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;
  printf("----------------------------------\n");
  printf(" All.PhysDensThresh = %g\n",All.PhysDensThresh);

  /* stupid experiement */
  /*
  All.PhysDensThresh = 100.0 * All.PhysDensThresh;
  printf(" All.PhysDensThresh = %g (post manipulation)\n",All.PhysDensThresh);
  printf("----------------------------------\n");
  */

}


/* ------------------------------
   Determine cold phase volume 
   fraction.

   Note: this does not require knowing the
         energy at all.  it is purely a
         function of the density.
   ------------------------------ */
double get_cold_volume_fraction(double Density, double Ne)
{
  double a3inv;
  int flag;
  double tsfr;
  double factorEVP, egyhot, ne, tcool, y, x, cloudmass;
  //double Ne = 1.0;



  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;


  if((Density * a3inv) < All.PhysDensThresh)
	return 0.0;


  flag = 1;			/* default is normal cooling */


  //printf("All.PhysDensThresh %e Dens %e\n",All.PhysDensThresh,Density);
  tsfr = sqrt(All.PhysDensThresh / (Density * a3inv)) * All.MaxSfrTimescale;

  factorEVP = pow(Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

  ne = Ne;

  //printf("tsfr %e fevp %e egyh %e ne %e dens %e\n",tsfr,factorEVP,egyhot,ne,Density);
  tcool = GetCoolingTime(egyhot, Density * a3inv, &ne);

  if(tcool>0)
    y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
  else
    y = 0.0;

  if(y>0)
    x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
  else
    x = 0.0;

  /*
  if((Density * a3inv) < All.PhysDensThresh)
	  x = 0.0;
  */

  //printf("rho %e x %e ne %e uh %e uc %e\n",Density,x,ne,egyhot,All.EgySpecCold);

  return x*All.EgySpecCold/(egyhot - x*egyhot + x*All.EgySpecCold);
}



/* ---------------------------------------
   Function which returns the fraction of 
   mass in the cold phase.
   --------------------------------------- */
double get_cold_mass_fraction(double Density, double Ne)
{
  double a3inv;
  int flag;
  double tsfr;
  double factorEVP, egyhot, ne, tcool, y, x, cloudmass;
  //double Ne = 1.0;



  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;


  if((Density * a3inv) < All.PhysDensThresh)
          return 0.0;



  flag = 1;			/* default is normal cooling */


  //printf("All.PhysDensThresh %e Dens %e\n",All.PhysDensThresh,Density);
  tsfr = sqrt(All.PhysDensThresh / (Density * a3inv)) * All.MaxSfrTimescale;

  factorEVP = pow(Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

  ne = Ne;

  //printf("tsfr %e fevp %e egyh %e ne %e dens %e\n",tsfr,factorEVP,egyhot,ne,Density);
  tcool = GetCoolingTime(egyhot, Density * a3inv, &ne);

  if(tcool>0)
    y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
  else
    y = 0.0;

  if(y>0)
    x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
  else
    x = 0.0;

  /*
  if((Density * a3inv) < All.PhysDensThresh)
	  x = 0.0;
  */

  //printf("rho %e x %e ne %e\n",Density,x,ne);
  return x;
}




/* ---------------------------------------
   Function which returns the temperature
   of the hot gas.
   --------------------------------------- */
double get_hotphase_temperature(double Density, double Ne)
{
  double a3inv;
  double factorEVP, egyhot, ne;
  double meanweight, hottemp;
  

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;
  
    
  if((Density * a3inv) < All.PhysDensThresh)
        return 0.0;

    
  factorEVP = pow(Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

  ne = Ne; 

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));/* note: assuming FULL ionization */

  hottemp = egyhot * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
  hottemp = meanweight * GAMMA_MINUS1 / (BOLTZMANN / PROTONMASS) * hottemp;

/*
  if((Density * a3inv) < All.PhysDensThresh)
          hottemp = 0.0;
*/
  
  return hottemp;
} 

  
          
/* --------------------
   Initialize Stuff
   -------------------- */
void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;
  double tsfr, y, peff, fac, neff, egyeff, factorEVP, sigma, thresholdStarburst;

  InitCool();

  if(All.PhysDensThresh == 0)
    {
      A0 = All.FactorEVP;

      egyhot = All.EgySpecSN / A0;

      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));/* note: assuming FULL ionization */

      u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
      u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;


      if(All.ComovingIntegrationOn)
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
      else
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);


      if(All.ComovingIntegrationOn)
	{
	  All.Time = 1.0;	/* to be guaranteed to get z=0 rate */
	  IonizeParams();
	}

      ne = 1.0;
      SetZeroIonization();
	//printf("in init clouds Hubble %e egyhot %e dens %e ne %e\n",All.Hubble,egyhot,dens,ne);
      tcool = GetCoolingTime(egyhot, dens, &ne);
	//printf("got cooling time\n");

      coolrate = egyhot / tcool / dens;

      x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      All.PhysDensThresh =
	x / pow(1 - x,
		2) * (All.FactorSN * All.EgySpecSN - (1 -
						      All.FactorSN) * All.EgySpecCold) /
	(All.MaxSfrTimescale * coolrate);

	/*printf("\nA0= %g  \n", A0);
	printf("Computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n", All.PhysDensThresh,
		 All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs));
	printf("EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g\n\n", x);
	printf("tcool=%g dens=%g egyhot=%g\n", tcool, dens, egyhot);*/

	printf("----------------------------------\n");
	printf(" All.PhysDensThresh = %g\n",All.PhysDensThresh);

	/* stupid experiement */
	/*
	All.PhysDensThresh = 100.0 * All.PhysDensThresh;
	printf(" All.PhysDensThresh = %g (post manipulation)\n",All.PhysDensThresh);
	printf("----------------------------------\n");
	*/

    }
  //printf("tsfr %e fevp %e egyh %e ne %e dens %e\n",All.MaxSfrTimescale,A0,egyhot,ne,dens);
}

