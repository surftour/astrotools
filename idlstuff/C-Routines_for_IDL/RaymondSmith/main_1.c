#include<stdio.h>

//extern "C"
//{
extern void emissivity_(double*, double*, double*, double*, double*, double*, double*, double*);

int main(int argc, char **argv)
{
	double ebin = 200.0;
	double etin = 1000.0;
	double zrin = 0.0;
	double rmin = 0.1;
	double tdegkin = 1.0e6;
	double edin = 0.1;
	double rnhin = 0.1;
	double ans= 0.0;
	printf("testing ebin %f \n",ebin);
	emissivity_(&ebin,&etin,&zrin,&rmin,&tdegkin,&edin,&rnhin,&ans);
	printf("answer %e\n",ans);
}
