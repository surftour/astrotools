#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"colors.h"
int main(int argc, char **argv)
{
	/*Takes age  in Gyr  (no h)*/
	/*Takes Z    in solar values (assuming Z_sol = 0.02)*/
	int    i,j;
	int    N_stars;
	int    model;
	int    silent;
	float *Z, *age;
	float  *mass;
	double *mass_store;
	float  *mass_p;
	float  *colors;
	double *colors_store;
	float  *colors_p;

	/*put colors into colors_temp*/
	double *colors_temp;
	double *mass_temp;
	int    n_m = 14;

	Colors CC;

	N_stars  = *(int *)  argv[0];
	age      =  (float *)argv[1]; //Gyr  (no h)
	Z        =  (float *)argv[2]; //Solar units
	colors   =  (float *)argv[3];
	mass     =  (float *)argv[4];
	model    = *(int *)  argv[5]; //0==salpeter, 1==chabrier
	silent   = *(int *)  argv[6]; //0==loud, 1==silent

	colors_store = (double *) malloc(N_stars*n_m*sizeof(double));
	mass_store   = (double *) malloc(N_stars*2*sizeof(double));

	for(j=0;j<n_m*N_stars;j++)
		colors_store[j] = 0.0;
	for(j=0;j<2*N_stars;j++)
		mass_store[j]   = 0.0;

	CC.Initialize(model);
	for(j=0;j<N_stars;j++)
	{
		//age is in Gyr, no h
		//Z   is in solar units
		//must convert to log years and mass fraction
		colors_temp = CC.GetMags(log10(age[j]*1.0e9),Z[j]*0.02);
		mass_temp   = CC.GetMass(log10(age[j]*1.0e9),Z[j]*0.02);
		for(i=0;i<n_m;i++)
			colors_store[n_m*j+i]  = (double) colors_temp[i];
		for(i=0;i<2;i++)
			mass_store[2*j+i]      = (double) mass_temp[i];
		if(!(j%5000))
		{
			if(!silent)
				printf("i %d z %f a %f Mbol %e Mass %e\n",j,Z[j],age[j],colors_temp[0],mass_temp[0]);
		}
		free(colors_temp);
		free(mass_temp);

	}

	for(j=0;j<N_stars;j++)
	{
		for(i=0;i<n_m;i++)
		{
			colors_p  = colors + n_m*j + i;
			*colors_p = colors_store[n_m*j+i];
		}
		for(i=0;i<2;i++)
		{
			mass_p  = mass + 2*j + i;
			*mass_p = mass_store[2*j+i];
		}
	}

	free(colors_store);
	free(mass_store);
	return 0;

}
