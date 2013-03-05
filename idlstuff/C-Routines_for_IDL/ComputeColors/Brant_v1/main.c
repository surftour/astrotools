#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include"luminosity.h"



int colors(int argc, void *argv[])
{
	int   N;
	float z_obs;
	float zz;
	float *a_form;
	float *Age;
	float *Zmet;
	float *Mass;
	float *Value;
	float *vv;
	float img[200];
	int nm = 8;
	int i,j;

	//float *mags;
	float mags[8];


	//mags = (float *) malloc(nm*sizeof(float));
	//Age  = (float *) malloc(N*sizeof(float));


	N    =*(int *)argv[0];
	z_obs=*(int *)argv[1];
	Mass =(float*)argv[2];
	Age  =(float*)argv[3];
	Zmet =(float*)argv[4];
	Value=(float*)argv[5];

	//printf("N %d z_obs %e mass[0] %e age[0] %e zmet[0] %e\n",N,z_obs,Mass[0],Age[0],Zmet[0]);
	printf(" ....finding colors....");
	Initialize_Luminosity();
	//for(i=0;i<1000;i++)
	for(i=0;i<N;i++)
	{
		if(!(i%10000))
			printf("%d..",i);
		L_Sloan(Mass[i],Age[i],Zmet[i],z_obs,mags);
		//L_Sloan(1.0,Age[i],Zmet[i],z_obs,mags);
		for(j=0;j<nm;j++)
		{
			//printf("mags[%d,%d] %e\n",i,j,mags[j]);
			vv  = Value + nm*i +j;
			*vv = mags[j];
		}
	}
	printf("\n");
	//free(mags);
	//for(i=0;i<1000;i++)
	for(i=0;i<N;i++)
	{
		vv = Value + nm*i;
		//printf("%e\n",*vv);
	}
	return 0;
}
