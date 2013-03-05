#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern void emissivity_(double*, double*, double*, double*, double*, double*, double*, double*);

int   N;
float redshift;
float *zmet,*e_density,*nh_bensity,*T,*hsml;
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
	double rnhin 	= 0.0;
	double ans	= 0.0;
	float  ans_float = 0.0;
	float *ans_input;

	N 		=*(int*)argv[0];
	redshift	=*(float*)argv[1];
	zmet 		=(float*)argv[2];
	e_density	=(float*)argv[3];
	nh_density	=(float*)argv[4];
	T		=(float*)argv[5];
	hsml		=(float*)argv[6];
	soft_lum	=(float*)argv[7];
	hard_lum	=(float*)argv[8];


	zrin = redshift;
	for(i=0;i<N;i++)
	{
		rmin 	= zmet[i];
		tdegkin = T[i];
		edin	= e_density[i];
		rnhin	= nh_density[i];

		//soft band
		ebin = 200.0;
		etin = 1000.0;
		emissivity_(&ebin,&etin,&zrin,&rmin,&tdegkin,&edin,&rnhin,&ans);
		ans_float = ans;
		ans_input = soft_lum+i;
		*ans_input = ans_float;

		//hard band
		ebin = 2000.0;
		etin = 10000.0;
		emissivity_(&ebin,&etin,&zrin,&rmin,&tdegkin,&edin,&rnhin,&ans);
		ans_float = ans;
		ans_input = hard_lum+i;
		*ans_input = ans_float;
	}

	return 0;
}
