#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"pl.h"

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
	int i,j;

	int nm = 14;
	double *mags;

	N    =*(int *)argv[0];
	z_obs=*(float *)argv[1];
	Mass =(float*)argv[2];
	Age  =(float*)argv[3];
	Zmet =(float*)argv[4];
	Value=(float*)argv[5];

	/*
	printf("N=%d\n",N);
	printf("T=%g\n",z_obs);
	printf("Mass (1,2,3,4,5)=%g.%g,%g,%g,%g\n",Mass[0],Mass[1],Mass[2],Mass[3],Mass[4]);
	printf("Age  (1,2,3,4,5)=%g.%g,%g,%g,%g\n",Age[0],Age[1],Age[2],Age[3],Age[4]);
	printf("Z    (1,2,3,4,5)=%g.%g,%g,%g,%g\n",Zmet[0],Zmet[1],Zmet[2],Zmet[3],Zmet[4]);
	printf("....finding colors......");
	*/

	mags = (double *) malloc(nm*sizeof(double));

	AllocateMagnitudes();
	LoadMagnitudeData();

	printf("Get Magnitudes, N=%d\n",N);
	for(i=0;i<N;i++)
	{
		if(!(i%10000))
			printf("%d..",i);
		GetMags(log10(Age[i]*1.0e9),Zmet[i],mags);
		for(j=0;j<nm;j++)
		{
			vv = Value + nm*i +j;
			*vv = mags[j];
		}
	}
/*
	for(i=0;i<N;i++)
		vv = Value + nm*i;
*/

	printf("\n");
	FreeMagnitudes();
	free(mags);

	return 0;

}
