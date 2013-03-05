#include<stdio.h>
#include<stdlib.h>

class Colors
{
	public:
		int model;	//0==salpeter, 1==chabrier
		int n_m;	//number of colors
		int n_z;	//number of metallicities
		int n_age;	//number of ages


		double *age;	//age array (gyr no h)
		double *z;	//metallicity array (mass fraction)
		double *Mbol;	//Bolometric magnitude
		double *u_AB;	//Sloan u_AB
		double *g_AB;	//Sloan u_AB
		double *r_AB;	//Sloan u_AB
		double *i_AB;	//Sloan u_AB
		double *z_AB;	//Sloan u_AB
		double *U;	//U
		double *B;	//B
		double *V;	//V
		double *R;	//R
		double *I;	//I
		double *J;	//J
		double *H;	//H
		double *K;	//K
		double *Mstar;	//Mstar
		double *Mgas;	//Mgas

		//Constructor, Destructor,
		//and Initialization routines
		Colors(void);
		Colors(int model);
		~Colors(void);
		void Initialize(int model);
		void AllocateMemory(void);
		void FreeMemory(void);

		double *GetMags(double age, double z);
		double *GetMass(double age, double z);
		void   LoadSalpeter(void);
		void   LoadChabrier(void);
};
