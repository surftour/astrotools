#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"stellar_spectrum.h"
int main(int argc, char **argv)
{
	int i,j;
	int N_stars, N_lambda;
	float *Z, *age, *mass;
	float *spectrum;
	float *lambda;
	double *spectrum_store;
	float *spectrum_p;

	StellarSpectrum SS;

	N_stars  = *(int *)  argv[0];
	mass     =  (float *)argv[1]; //Msun
	Z        =  (float *)argv[2]; //Solar units
	age      =  (float *)argv[3]; //Gyr
	N_lambda = *(int *)  argv[4];
	lambda   =  (float *)argv[5];
	spectrum =  (float *)argv[6];

	spectrum_store = (double *) malloc(N_lambda*sizeof(double));

	for(j=0;j<N_lambda;j++)
		spectrum_store[j] = 0.0;
	
	for(i=0;i<N_stars;i++)
	{
		SS.Initialize(mass[i],Z[i],age[i]);

		for(j=0;j<N_lambda;j++)
			spectrum_store[j] += pow(10.0,SS.Spectrum(lambda[j]));

		if(!(i%5000))
		{
			printf("i %d mass %e z %f a %f S %e\n",i,mass[i],Z[i],age[i],spectrum_store[0]);
		}
		SS.Clear();
		SS.FreeSpectrum();
	}

	for(j=0;j<N_lambda;j++)
	{
		spectrum_p  = spectrum + j;
		*spectrum_p = log10(spectrum_store[j]);
	}

	free(spectrum_store);
	return 0;

}
