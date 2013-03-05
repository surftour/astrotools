#include<gsl/gsl_interp.h>
#ifndef BRANT_STELLAR_SPECTRUM
#define BRANT_STELLAR_SPECTRUM

class StellarSpectrum
{
	public:
		int int_flag;	//interpolation flag
		int L_flag;	//Spectrum flag
		int l_flag;	//wavelength flag

		// number of L_lambda bins in spectrum
		int N_lambda;

		double mass; //mass of stellar population
			     //in 1e10 msun
		double Z;    //metallicity of stellar population 
		             //(mass fraction)
		double age;  //age of stellar population in gyr 

		//returns L_lambda at wavelength lambda [in mum]
		double Spectrum(double wavelength);

		//returns L_nu at frequency nu [in Hz]
		double Spectrum_nu(double nu);

		//Constructor, Destructor,
		//and Initialization routines
		StellarSpectrum(void);
		StellarSpectrum(double mass, double Z, double age);
		~StellarSpectrum(void);
		void Initialize(double mass, double Z, double age);
		void Clear(void);
		void FreeSpectrum(void);
		void FreeWavelength(void);
		double * LoadLambda(void);
		double * LoadSpectrum(double mass, double Z, double age);

		//Stores L_lambda and lambda arrays
		double *L_lambda;
		double *lambda;

		//Attenuation
		double Calzetti(double wavelength);
		void   ApplyExtinction(double EBV);
		double Rv;
		double ESBV_factor;

		//Interpolation variables
		gsl_interp_accel  *acc;
		gsl_interp *lin;

};

#endif
