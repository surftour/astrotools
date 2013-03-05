#include<stdio.h>
#include<stdlib.h>
#include<math.h>
/*
#include "manual_gsl_interp.h"
*/
#include "interstellar_absorption.h"



InterstellarAbsorption::InterstellarAbsorption(void)
{
	int i;
	N_lambda = 500;
	R_v = 3.1;
	NH_to_A_V = 1.79e21; //Allen's Astrophysical Quantities
	lambda   = LoadLambda();
	if(!(a   = (double *) malloc(N_lambda*sizeof(double))))
	{
		printf("Error allocating 'a' in interstellar_absorption.c.\n");
		exit(-1);
	}
	if(!(b   = (double *) malloc(N_lambda*sizeof(double))))
	{
		printf("Error allocating 'b' in interstellar_absorption.c.\n");
		exit(-1);
	}
	if(!(pl   = (double *) malloc(N_lambda*sizeof(double))))
	{
		printf("Error allocating 'pl' in interstellar_absorption.c.\n");
		exit(-1);
	}
	for(i=0;i<N_lambda;i++)
	{
		a[i] = 0.0;
		b[i] = 0.0;
		pl[i] = pow(10.0,lambda[i]);
	}
	Initialize();
}
double InterstellarAbsorption::Extinction(double wavelength, double NH, double Z)
{
	double A_V;
	double A;
	double l;

	//Z is metallicity mass fraction in gas

	A_V = pow(10.0,NH)/NH_to_A_V;
	l = log10(wavelength);
	if(l<lambda[0])
	{
		A = -0.4*A_V*A_x[0]*(Z/0.02);
	}else if(l>lambda[N_lambda-1])
	{
		A = -0.4*A_V*A_x[N_lambda-1]*(Z/0.02);
	}else{
		//interpolate between zero metallicity and metallicity z gas
		//absorption
		A = (Z/0.02)*(-0.4*A_V*gsl_interp_eval(lin,lambda,A_x,l,acc)) - (-0.4*A_V*gsl_interp_eval(lin_nm,lambda,A_x_nm,l,acc))*((Z/0.02)-1.0);
	}
	return A;
}
double * InterstellarAbsorption::LoadAbsorption(void)
{
	//From Cardelli, Clayton, and Mathis 1989
	
	int i;
	double x, y, Fa, Fb;
	double x_ryd = 10.9678;
	double x_keV;
	double dx;
	double *A;

	x_keV = (2.418e17)/(3.0e14); // 1/keV_in_microns

	if(!(A = (double *) malloc(N_lambda*sizeof(double))))
	{
		printf("Error allocating absorption array.\n");
		exit(-1);
	}
	for(i=0;i<N_lambda;i++)
	{
		x = 1.0/pow(10.0,lambda[i]);

		if( (x>=0.3)&&(x<=1.1) )
		{
			a[i] =  0.574*pow(x,1.61); //eqn 2a
			b[i] = -0.527*pow(x,1.61); //eqn 2b
			A[i] = (a[i] + b[i]/R_v);

		}else if( (x>1.1)&&(x<=3.3) )
		{
			//eqn 3a and 3b
			y = x - 1.82;
			a[i] = 1.0 + 0.17699*y - 0.50447*y*y - 0.02427*pow(y,3.0) + 0.72085*pow(y,4.0) + 0.01979*pow(y,5.0) - 0.77530*pow(y,6.0) + 0.32999*pow(y,7.0);
			b[i] = 1.41338*y + 2.28305*y*y + 1.07233*pow(y,3.0) - 5.38434*pow(y,4.0) - 0.62251*pow(y,5.0) + 5.30260*pow(y,6.0) - 2.09002*pow(y,7.0);
			A[i] = (a[i] + b[i]/R_v);

		}else if( (x>3.3)&&(x<=8.0) )
		{
			//eqn 4a and 4b
			if(x<5.9)
			{
				Fa = 0.0;
				Fb = 0.0;
			}else{
				Fa = -0.04473*(x-5.9)*(x-5.9) - 0.009779*(x-5.9)*(x-5.9)*(x-5.9);
				Fb =   0.2130*(x-5.9)*(x-5.9) +   0.1207*(x-5.9)*(x-5.9)*(x-5.9);
			}
			a[i] =  1.752 - 0.316*x - 0.104/((x-4.67)*(x-4.67) + 0.341) + Fa;
			b[i] = -3.090 + 1.825*x + 1.206/((x-4.62)*(x-4.62) + 0.263) + Fb;
			A[i] = (a[i] + b[i]/R_v);
		}else if( (x>8.0)&&(x<x_ryd) )
		{
			//eqn 5a and 5b
			a[i] = -1.073 - 0.628*(x-8.0) + 0.137*(x-8.0)*(x-8.0) - 0.070*pow(x-8.0,3.0);
			b[i] = 13.670 + 4.257*(x-8.0) - 0.420*(x-8.0)*(x-8.0) + 0.374*pow(x-8.0,3.0);
			A[i] = (a[i] + b[i]/R_v);
		}else if(x>=x_ryd)
		{
			//Morrison and McCammon 1983
			//NH_to_AV gets divided out later
			//2.5 adjusts for magnitude conversion
			//log(10) adjusts for e to 10 log base conversion
			A[i] = 2.5*Morrison1983(x/x_keV)*NH_to_A_V/log(10.0);
		}
	}
	A_flag = 1;
	return A;
}
double * InterstellarAbsorption::LoadAbsorptionNoMetals(void)
{
	//From Cardelli, Clayton, and Mathis 1989
	
	int i;
	double x, y, Fa, Fb;
	double x_ryd = 10.9678;
	double x_keV;
	double dx;
	double *A;

	x_keV = (2.418e17)/(3.0e14); // 1/keV_in_microns

	if(!(A = (double *) malloc(N_lambda*sizeof(double))))
	{
		printf("Error allocating absorption array.\n");
		exit(-1);
	}
	for(i=0;i<N_lambda;i++)
	{
		x = 1.0/pow(10.0,lambda[i]);

		if(x>x_ryd)
		{
			//Morrison and McCammon 1983
			//NH_to_AV gets divided out later
			//2.5 adjusts for magnitude conversion
			//log(10) adjusts for e to 10 log base conversion
			//essential just H+He bound-free absorption
			A[i] = 2.5*Morrison1983_no_metals(x/x_keV)*NH_to_A_V/log(10.0);
		}else{
			//no dust absorption!
			A[i]=0.0;
		}
	}
	A_flag = 1;
	return A;
}
double InterstellarAbsorption::Cruddace1974(double x)
{
	//Cruddace et al 1974
	// x is in per 1/Rydberg_in_microns
	return (6.22e-18)*pow(x,-2.43);
}
double InterstellarAbsorption::Morrison1983_no_metals(double x)
{
	//Morrison and McCammon 1983
	// x is in per 1/keV_in_microns
	double x_keV = (2.418e17)/(3.0e14); // 1/keV_in_microns
	double c0, c1, c2;

	if(x<0.03)
	{
		return Morrison1983(0.03)*pow(x/0.03,-2.43);
	}
	else if((x>=0.03)&&(x<0.1))
	{
		c0 =    17.3;
		c1 =   608.1;
		c2 = -2150.0;

	}
	else if(x>0.1)
	{
		return Morrison1983(0.1)*g_bf(x*x_keV,0.0,1)*pow(x/0.1,-3.0)/g_bf(0.1*x_keV,0.0,1);
	}else{
		return 0.0;
	}
	return (1.0e-24)*(c0+c1*x+c2*x*x)/(x*x*x); //cm^2
}
double InterstellarAbsorption::Morrison1983(double x)
{
	//Morrison and McCammon 1983
	// x is in per 1/keV_in_microns
	double x_keV = (2.418e17)/(3.0e14); // 1/keV_in_microns
	double c0, c1, c2;

	if(x<0.03)
	{
		return Morrison1983(0.03)*pow(x/0.03,-2.43);
	}
	else if((x>=0.03)&&(x<0.1))
	{
		c0 =    17.3;
		c1 =   608.1;
		c2 = -2150.0;

	}else if((x>=0.1)&&(x<0.284))
	{
		c0 =    34.6;
		c1 =   267.9;
		c2 =  -476.1;

	}else if((x>=0.284)&&(x<0.4))
	{
		c0 =    78.1;
		c1 =    18.8;
		c2 =     4.3;

	}else if((x>=0.4)&&(x<0.532))
	{
		c0 =    71.4;
		c1 =    66.8;
		c2 =   -51.4;

	}else if((x>=0.532)&&(x<0.707))
	{
		c0 =    95.5;
		c1 =   145.8;
		c2 =   -61.1;

	}else if((x>=0.707)&&(x<0.867))
	{
		c0 =   308.9;
		c1 =  -380.6;
		c2 =   294.0;

	}else if((x>=0.867)&&(x<1.303))
	{
		c0 =   120.6;
		c1 =   169.3;
		c2 =   -47.7;

	}else if((x>=1.303)&&(x<1.840))
	{
		c0 =   141.3;
		c1 =   146.8;
		c2 =   -31.5;

	}else if((x>=1.840)&&(x<2.471))
	{
		c0 =   202.7;
		c1 =   104.7;
		c2 =   -17.0;

	}else if((x>=2.471)&&(x<3.210))
	{
		c0 =   342.7;
		c1 =    18.7;
		c2 =     0.0;

	}else if((x>=3.210)&&(x<4.038))
	{
		c0 =   352.2;
		c1 =    18.7;
		c2 =     0.0;

	}else if((x>=4.038)&&(x<7.111))
	{
		c0 =   433.9;
		c1 =    -2.4;
		c2 =     0.75;

	}else if((x>=7.111)&&(x<8.331))
	{
		c0 =   629.0;
		c1 =    30.9;
		c2 =     0.0;

	}else if((x>=1.840)&&(x<2.471))
	{
		c0 =   701.2;
		c1 =    25.2;
		c2 =     0.0;
	}else{
		return 0.0;
	}
	return (1.0e-24)*(c0+c1*x+c2*x*x)/(x*x*x); //cm^2
}
double * InterstellarAbsorption::ApplyExtinction(double L_lambda[], double lambda[], int N_lambda, double NH, double Z)
{
	int i;
	double x;
	double *L;
	double l;

	if(!(L = (double *) malloc(N_lambda*sizeof(double))))
	{
		printf("Error allocating extinction array.\n");
		exit(-1);
	}
	if(A_flag)
	{
		for(i=0;i<N_lambda;i++)
		{
			l = pow(10.0,lambda[i]);
			//Apply extinction
			L[i] = L_lambda[i] + Extinction(l,NH,Z);
		}
	}
	return L;
}
InterstellarAbsorption::~InterstellarAbsorption(void)
{
	free(a);
	free(b);
	free(pl);
	Clear();
	FreeAbsorption();
	FreeWavelength();
}
void InterstellarAbsorption::FreeAbsorption(void)
{
	if(A_flag)
		free(A_x);
	A_flag = 0;
}
void InterstellarAbsorption::FreeWavelength(void)
{
	if(l_flag)
		free(lambda);
	l_flag = 0;
}
double * InterstellarAbsorption::LoadLambda(void)
{
	double *l;
	int i, ii, ia;
	if(!(l = (double *) malloc(N_lambda*sizeof(double))))
	{
		printf("Error allocating extinction array.\n");
		exit(-1);
	}
	//for(i=0;i<N_lambda;i++)
		//l[i] = log10(3.0e14/pow(10.0,-31.0*((float) i)/((float)(N_lambda-1)) +23.0));

	//use 20% to fill in shortward of 10 keV
	ii = 0;
	ia = (int) (N_lambda/5);
	for(i=ii;i<ia;i++)
		l[i] = log10(3.0e14/pow(10.0,-5.0*((float) i-ii)/((float)(ia-ii))+23.0));
		//l[i] = log10(3.0e14/pow(10.0,-5.0*((float) i-ii)/((float)(ia-1))+23.0));
	//use 60% to fill in between 10 microns and 10keV
	ii = (int) N_lambda/5;
	ia = (int) (4*N_lambda/5);
	for(i=ii;i<ia;i++)
		l[i] = log10(3.0e14/pow(10.0,-4.5*((float) i-ii)/((float)(ia-ii))+18.0));
		//l[i] = log10(3.0e14/pow(10.0,-4.5*((float) i-ii)/((float)(ia-1))+18.0));
	//use 20% to fill in longward of 10 microns
	ii = (int) (4*N_lambda/5);
	ia = N_lambda;
	for(i=ii;i<ia;i++)
		l[i] = log10(3.0e14/pow(10.0,-6.5*((float) i-ii)/((float)(ia-ii-1))+13.5));
	l_flag = 1;

	//for(i=0;i<N_lambda;i++)
	//	printf("i %d lambda[%d] %f\n",i,i,l[i]);
	return l;
}
void InterstellarAbsorption::Initialize(void)
{
	A_x = LoadAbsorption();
	A_x_nm = LoadAbsorptionNoMetals();
	acc = gsl_interp_accel_alloc();
	acc_nm = gsl_interp_accel_alloc();
	lin = gsl_interp_alloc(gsl_interp_linear,N_lambda);
	lin_nm = gsl_interp_alloc(gsl_interp_linear,N_lambda);
	gsl_interp_init(lin,lambda,A_x,N_lambda);
	gsl_interp_init(lin_nm,lambda,A_x_nm,N_lambda);
	flag = 1;
}
void InterstellarAbsorption::Clear(void)
{
	if(flag)
	{
		gsl_interp_accel_free(acc);
		gsl_interp_accel_free(acc_nm);
		gsl_interp_free(lin);
		gsl_interp_free(lin_nm);
	}
	flag = 0;
}
double InterstellarAbsorption::g_bf(double x, double Te, int n)
{
	double x_g;
	double nu_g = 3.28989e15;//Hz
	double c = 3.0e14;//microns / sec
	double u_n = 0.0;
	double n_d = 0.0;

	n_d = ((double) n);
	x_g = nu_g/c;
	u_n = n_d*n_d*(x/x_g);

	//from allen's astrophysical quantities
	return 1.0 + (0.1728*(u_n-1.0)/(n_d*pow(u_n+1.0,2.0/3.0))) - (0.0496*(u_n*u_n + (4.0/3.0)*u_n + 1.0)/(n_d*pow(u_n+1.0,4.0/3.0)));
}
