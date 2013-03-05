#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern void emissivity_(double*, double*, double*, double*, double*, double*, double*, double*);

int   N;
float redshift;
float *zmet,*e_density,*rnh,*mass,*T,*hsml;
float *soft_lum,*hard_lum;

int raymond_smith(int argc, void *argv[])
{
	int i;
	double ebin 	= 0.0;
	double etin 	= 0.0;
	double zrin 	= 0.0;
	double rmin 	= 0.0;
	double tdegkin 	= 0.0;
	double edin 	= 0.0;
	double rnhin 	= 1.0;
	double ans	= 0.0;
	float  ans_float = 0.0;
	double conversion = 9.0561e56;
	/* conversion comes from 1 solar mass=1.989e33 g
           divided by 1.667e-24 g per proton, times 0.76
           to account for the mass fraction of H.  this, in
           conjunction with the rnh (hydrogen number den.),
           gives us the physical volume of our particle */
	double solar_lum = 3.9e33; 
	float *ans_input;

	N 		=*(int*)argv[0];
	redshift	=*(float*)argv[1];
	zmet 		=(float*)argv[2];
	e_density	=(float*)argv[3];
	rnh             =(float*)argv[4];
	mass		=(float*)argv[5];
	T		=(float*)argv[6];
	hsml		=(float*)argv[7];
	soft_lum	=(float*)argv[8];
	hard_lum	=(float*)argv[9];


	zrin = redshift;
	printf("calculating xray luminosities...\n");
	printf("N= %d\n",N);
	fflush(stdout);
	for(i=0;i<N;i++)
	{
		rmin 	= zmet[i];
		tdegkin = T[i];
		edin	= e_density[i];
		rnhin   = rnh[i];

		if(tdegkin>1.0e5)
		{
			//soft band
			;ebin = 100.0;
			ebin = 500.0;
			etin = 2000.0;
			emissivity_(&ebin,&etin,&zrin,&rmin,&tdegkin,&edin,&rnhin,&ans);
			/* ans_float = ans*mass[i]*conversion/solar_lum; */
			/* RaymondSmith returns erg s-1 cm-3, hence we need
                           the volume.  We get this from rnh, the hydrogen
                           number density, the particle mass, and the conversion 
                           factor above.  Finally, brant converted this to
                           solar luminosities. */
			ans_float = ans*mass[i]*conversion/solar_lum/rnh[i];
			ans_input = soft_lum+i;
			*ans_input = ans_float;
			if(!(i%1000))
			{
				printf("i=%d  m=%e T=%e ne=%e nh=%e z=%f L=%e  ...\n",i,mass[i],tdegkin,edin,rnhin,rmin,ans_float);
				fflush(stdout);
			}

			//hard band
			ebin = 2000.0;
			etin = 10000.0;
			emissivity_(&ebin,&etin,&zrin,&rmin,&tdegkin,&edin,&rnhin,&ans);
			/* ans_float = ans*mass[i]*conversion/solar_lum; */
			ans_float = ans*mass[i]*conversion/solar_lum/rnh[i];
			ans_input = hard_lum+i;
			*ans_input = ans_float;
		}else{
			ans_float = 0.0;
			ans_input = soft_lum+i;
			*ans_input = ans_float;
			ans_float = 0.0;
			ans_input = hard_lum+i;
			*ans_input = ans_float;
		}
	}
	printf("done!\n");

	return 0;
}
