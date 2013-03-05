#include<stdio.h>
#define PF(x) &x[1]

extern void xspexrav_(double *, int *, int *, double *, double *);

int Ne,Nel;
float redshift;
float *agn_xray_spectrum;
float *spectrum_energies;
int marconi_agn_spectrum(int argc, void *argv[])
{
	int i,ifl = 1;
	double *ear,*spec;
	double zr;
	float *s;
	Ne                = *(int*)argv[0];
	agn_xray_spectrum = (float*)argv[1];
	spectrum_energies = (float*)argv[2];
	redshift          = *(float*)argv[3];

	zr = redshift;

	ear  = (double *) malloc((Ne+1)*sizeof(double));
	spec = (double *) malloc((Ne+1)*sizeof(double));

	for(i=0;i<(Ne+1);i++)
	{
		ear[i] = spectrum_energies[i];
		//printf("i %d ear[%d] %e\n",i,i,ear[i]);
	}
	for(i=0;i<Ne+1;i++)
		spec[i] = 0.0;

	Nel = Ne-1;
	xspexrav_(&ear[0],&Nel,&ifl,PF(spec),&zr);
	for(i=0;i<Ne;i++)
	{
		s = agn_xray_spectrum + i;
		*s = spec[i+1];
	}
	free(spec);
	free(ear);
}
#undef PF
