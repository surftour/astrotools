#ifndef EOS_H
#define EOS_H

#ifdef EOS_DEGENERATE

#define EOS_MAXITER 40

/* basic constants in cgs units */
#define EOS_PI 3.14159265359 /* Pi */
#define EOS_EPS 1.0e-13

#ifdef STANDALONE
#define PROTONMASS 1.6726e-24
#define BOLTZMANN 1.3806e-16
#endif

void eos_init( char* datafile, char* speciesfile );
void eos_deinit();
int eos_calc_egiven( double rho, double *xnuc, double e, double *temp, double *p, double *dpdr );
int eos_calc_tgiven( double rho, double *xnuc, double temp, double *e, double *dedt );
int eos_calc_egiven_v( double rho, double *xnuc, double *dxnuc, double fac, double e, double *temp, double *p, double *dpdr );
int eos_calc_tgiven_v( double rho, double *xnuc, double *dxnuc, double fac, double temp, double *e, double *dedt );
void eos_trilinear_e( double temp, double rho, double ye, double *e, double *dedt );
void eos_trilinear( double temp, double rho, double ye, double *e, double *dedt, double *p, double *dpdr, double *dpdt );

double eos_SwapDouble( double Val );
int eos_SwapInt( int Val );
double eos_calcYe( double *xnuc );
void eos_checkswap( char* fname, int *swap );


struct eos_table {
  double *nuclearmasses;
  double *nuclearcharges;
  int ntemp, nrho, nye;
  double tempMin, tempMax, ltempMin, ltempMax, ltempDelta, ltempDeltaI;
  double rhoMin, rhoMax, lrhoMin, lrhoMax, lrhoDelta, lrhoDeltaI;
  double yeMin, yeMax, yeDelta, yeDeltaI;
  double *ltemp, *lrho, *ye; 
   
  double *p;           /* pressure */
  double *dpdt;        /* derivative of pressure with respect to temperature */
  double *dpdr;        /* derivative of pressure with respect to density */
  double *e;           /* energy per mass */
  double *dedt;        /* derivative of energy with respect to temperature */
} eos_table;

#endif /* EOS_DEGENERATE */

#endif /* EOS_H */
