#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include"interstellar_absorption.h"


int calc_atten_lum(int argc, void *argv[])
{
	int N;
	float Band_L;
	float L_lambda = 14.0; //An example L_lambda in Lsun micron^-1
	float lambda = 0.5;//an example wavelength in microns
	float Z    = 0.1; //an example metallicity in solar units
	float NH   = 21.0; // an example column density
	float *iLum;
	float *iNH;
	float *iZ;
	float *Value;
	float *vv;
	float Ext, NewLum;
	int i;


	if(argc!=6)
	    {
	      fprintf(stderr,"\n\nwrong number of arguments ! (found %d) \n\n",argc);
	      exit(0);
	    }


	N	=*(int *)argv[0];
	Band_L  =*(float *)argv[1];
	iLum	=(float*)argv[2];
	iNH	=(float*)argv[3];
	iZ	=(float*)argv[4];
	Value	=(float*)argv[5];

	Initialize();

	lambda= Band_L/1.0e4;

	printf("Determine Attenuation, N=%d\n",N);
	printf("Band_L=%g   lambda= %g\n",Band_L,lambda);

	for(i=0;i<N;i++)
	{
		if(!(i%10000))
		    printf("%d..",i);

		L_lambda = log10(iLum[i]);
		NH = log10(iNH[i]);
		Z = iZ[i];
		Ext = Extinction(lambda,NH,Z);
		NewLum= pow(10.0, L_lambda + Ext);
/*
if(i<5){
printf("i= %d   iLum[%d]= %g  iNH[%d]= %g  Ext= %g  NewLum %g\n",i,i,iLum[i],i,iNH[i],Ext,NewLum); fflush(stdout);
}
*/
		vv = Value + i;
		*vv = NewLum;
/*
		*vv = pow(10.0, L_lambda+IA.Extinction(lambda,NH,Z));
		*vv = 10.0;
*/
	}

	printf("... done\n");
	fflush(stdout);


	FreeMemory();

	return 0;
}




