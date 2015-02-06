#ifdef EOS_DEGENERATE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eos.h"

#ifndef STANDALONE
#include "allvars.h"
#include "proto.h"
#endif

#ifndef STANDALONE
#define EXIT endrun(1337)
#else
#define EXIT exit(1337)
#endif

void eos_init( char* datafile, char* speciesfile ) {
  int swap, skip;
  FILE *fd;
  int i, entries, bytes, species;
  char cdummy[5];
  int mass, charge;

  fd = 0;
#ifndef STANDALONE
  if (ThisTask == 0) {
#endif
    eos_checkswap( datafile, &swap );

    if (swap == 2) {
      printf( "eos table `%s' is corrupt.\n", datafile );
      EXIT;
    }

    if (!(fd = fopen(datafile, "r"))) {
      printf( "can't open file `%s' for reading eos table.\n", datafile );
      EXIT;
    }

    fread( &skip, sizeof(int), 1, fd );
    fread( &eos_table.ntemp, sizeof(int), 1, fd );
    fread( &eos_table.nrho, sizeof(int), 1, fd );
    fread( &eos_table.nye, sizeof(int), 1, fd );
    fread( &skip, sizeof(int), 1, fd );

    fread( &skip, sizeof(int), 1, fd );
    fread( &eos_table.ltempMin, sizeof(double), 1, fd );
    fread( &eos_table.ltempMax, sizeof(double), 1, fd );
    fread( &eos_table.lrhoMin, sizeof(double), 1, fd );
    fread( &eos_table.lrhoMax, sizeof(double), 1, fd );
    fread( &eos_table.yeMin, sizeof(double), 1, fd );
    fread( &eos_table.yeMax, sizeof(double), 1, fd );
    fread( &skip, sizeof(int), 1, fd );

    if (swap) {
      eos_table.ntemp = eos_SwapInt( eos_table.ntemp );
      eos_table.nrho = eos_SwapInt( eos_table.nrho );
      eos_table.nye = eos_SwapInt( eos_table.nye );

      eos_table.ltempMin = eos_SwapDouble( eos_table.ltempMin );
      eos_table.ltempMax = eos_SwapDouble( eos_table.ltempMax );
      eos_table.lrhoMin = eos_SwapDouble( eos_table.lrhoMin );
      eos_table.lrhoMax = eos_SwapDouble( eos_table.lrhoMax );
      eos_table.yeMin = eos_SwapDouble( eos_table.yeMin );
      eos_table.yeMax = eos_SwapDouble( eos_table.yeMax );
    }
    
#ifndef STANDALONE
  }
#endif

#ifndef STANDALONE
  MPI_Bcast( &eos_table.ntemp, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &eos_table.nrho , 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &eos_table.nye  , 1, MPI_INT, 0, MPI_COMM_WORLD );

  MPI_Bcast( &eos_table.ltempMin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( &eos_table.ltempMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( &eos_table.lrhoMin , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( &eos_table.lrhoMax , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( &eos_table.yeMin  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( &eos_table.yeMax  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif

  eos_table.ltempDelta = (eos_table.ltempMax - eos_table.ltempMin) / (eos_table.ntemp-1);
  eos_table.lrhoDelta = (eos_table.lrhoMax - eos_table.lrhoMin) / (eos_table.nrho-1);
  eos_table.yeDelta = (eos_table.yeMax - eos_table.yeMin) / (eos_table.nye-1);

  eos_table.tempMin = pow( 10.0, eos_table.ltempMin );
  eos_table.tempMax = pow( 10.0, eos_table.ltempMax );
  eos_table.rhoMin = pow( 10.0, eos_table.lrhoMin );
  eos_table.rhoMax = pow( 10.0, eos_table.lrhoMax );  

#ifndef STANDALONE
  if (ThisTask == 0) {
#endif
    printf( "EOS table spans rho [%e,%e], temp [%e,%e], ye [%e,%e].\n", eos_table.rhoMin, eos_table.rhoMax, eos_table.tempMin, eos_table.tempMax, eos_table.yeMin, eos_table.yeMax );
    printf( "Resolution is rho: %d, temp: %d, ye: %d points.\n", eos_table.nrho, eos_table.ntemp, eos_table.nye );
#ifndef STANDALONE
  }
#endif

  eos_table.ltemp = (double *)malloc( eos_table.ntemp * sizeof(double) );
  eos_table.lrho = (double *)malloc( eos_table.nrho * sizeof(double) );
  eos_table.ye = (double *)malloc( eos_table.nye * sizeof(double) );

  if (!(eos_table.ltemp && eos_table.lrho && eos_table.ye)) {
    printf( "not enough memory to allocate eos arrays.\n" );
    EXIT;
  }

  eos_table.ltemp[0] = eos_table.ltempMin;
  for ( i=1; i<eos_table.ntemp; i++ ) eos_table.ltemp[i] = eos_table.ltemp[i-1] + eos_table.ltempDelta;
  eos_table.lrho[0] = eos_table.lrhoMin;
  for ( i=1; i<eos_table.nrho; i++ ) eos_table.lrho[i] = eos_table.lrho[i-1] + eos_table.lrhoDelta;
  eos_table.ye[0] = eos_table.yeMin;
  for ( i=1; i<eos_table.nye; i++ ) eos_table.ye[i] = eos_table.ye[i-1] + eos_table.yeDelta;

  eos_table.ltempDeltaI = 1. / eos_table.ltempDelta;
  eos_table.lrhoDeltaI = 1. / eos_table.lrhoDelta;
  eos_table.yeDeltaI = 1. / eos_table.yeDelta;

  entries = eos_table.ntemp * eos_table.nrho * eos_table.nye;
  bytes = entries * sizeof(double);
  eos_table.p = (double *)malloc( bytes );
  eos_table.dpdt = (double *)malloc( bytes );
  eos_table.dpdr = (double *)malloc( bytes );
  eos_table.e = (double *)malloc( bytes );
  eos_table.dedt = (double *)malloc( bytes );

  if (!(eos_table.p && eos_table.e && eos_table.dedt)) {
    printf( "not enough memory to allocate eos arrays.\n" );
    EXIT;
  }

#ifndef STANDALONE
  if (ThisTask == 0) {
#endif
    printf( "Reading grid containing %d points for each quantity.\n", entries );

    fread( &skip, sizeof(int), 1, fd );
    fread( eos_table.p, sizeof(double), entries, fd );
    fread( eos_table.dpdt, sizeof(double), entries, fd );
    fread( eos_table.dpdr, sizeof(double), entries, fd );
    fread( eos_table.e, sizeof(double), entries, fd );
    fread( eos_table.dedt, sizeof(double), entries, fd );
    fread( &skip, sizeof(int), 1, fd );

    if (swap) {
      for ( i=0; i<entries; i++ ) {
        eos_table.p[i] = eos_SwapDouble( eos_table.p[i] );
        eos_table.dpdt[i] = eos_SwapDouble( eos_table.dpdt[i] );
        eos_table.dpdr[i] = eos_SwapDouble( eos_table.dpdr[i] );
        eos_table.e[i] = eos_SwapDouble( eos_table.e[i] );
        eos_table.dedt[i] = eos_SwapDouble( eos_table.dedt[i] );
      }
    }

    fclose( fd );
#ifndef STANDALONE
  }
#endif

#ifndef STANDALONE
  MPI_Bcast( eos_table.p, entries, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( eos_table.dpdt, entries, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( eos_table.dpdr, entries, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( eos_table.e, entries, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( eos_table.dedt, entries, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif

  eos_table.nuclearmasses = (double *)malloc( EOS_NSPECIES * sizeof(double) );
  eos_table.nuclearcharges = (double *)malloc( EOS_NSPECIES * sizeof(double) );

#ifndef STANDALONE
  if (ThisTask == 0) {
#endif
    if (!(fd = fopen(speciesfile, "r"))) {
      printf( "can't open file `%s' for reading species information.\n", speciesfile );
      EXIT;
    }
    
    fscanf( fd, "%d", &species );
    if (species != EOS_NSPECIES) {
      printf( "code compiled for wrong number of species (compiled for %d species, speciesfile contains %d species).\n", EOS_NSPECIES, species );
      EXIT;
	}
    
    for ( i=0; i<EOS_NSPECIES; i++ ) {
      fscanf( fd, "%5s%d%d", cdummy, &mass, &charge );
      eos_table.nuclearmasses[i] = (double)mass;
      eos_table.nuclearcharges[i] = (double)charge;
	}
	    
    fclose( fd );
#ifndef STANDALONE
  }
#endif

#ifndef STANDALONE
  MPI_Bcast( eos_table.nuclearmasses, EOS_NSPECIES, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( eos_table.nuclearcharges, EOS_NSPECIES, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif

#ifndef STANDALONE
  if (ThisTask == 0) {
    printf( "EOS init done.\n" );
  }
#else
    printf( "EOS init done.\n" );
#endif
}

void eos_deinit() {
  free( eos_table.p );
  free( eos_table.dpdt );
  free( eos_table.dpdr );
  free( eos_table.e );
  free( eos_table.dedt );
  free( eos_table.nuclearmasses );
}

int eos_calc_egiven_v( double rho, double *xnuc, double *dxnuc, double fac, double e, double *temp, double *p, double *dpdr ) {
	double x[EOS_NSPECIES];
	int i;
	
	for (i=0; i<EOS_NSPECIES; i++) {
		x[i] = xnuc[i] + dxnuc[i]*fac;
	}
	
	return eos_calc_egiven( rho, x, e, temp, p, dpdr );
}

int eos_calc_tgiven_v( double rho, double *xnuc, double *dxnuc, double fac, double temp, double *e, double *dedt ) {
	double x[EOS_NSPECIES];
	int i;
	
	for (i=0; i<EOS_NSPECIES; i++) {
		x[i] = xnuc[i] + dxnuc[i]*fac;
	}
	
	return eos_calc_tgiven( rho, x, temp, e, dedt );
}

int eos_calc_egiven( double rho, double *xnuc, double e, double *temp, double *p, double *dpdr ) {
  double *n, ni, nn, ne;        /* number densities [ particles per volume ] */
  double ye;                    /* electron fraction ( electrons per baryon ) */
  double _temp, _tempold;       /* temperatures [ K ] */
  double _ee, _deedt;           /* energy per mass of electrons + positrons + radiation */
  double _ie, _diedt;           /* energy per mass of ions */
  double _e, _ep;               /* total energy per mass, pressure of electrons */
  double _dpedr, _dpedt;		/* partial deriverative of pressure with respect to density and temperature of electrons */
  double dpdt, dedt;
  int iter, i;

  if (*temp == 0) return 0;

  if (rho > eos_table.rhoMax) {
    printf( "Density exceeds allowed maximum, rho: %g, maximum: %g\n", rho, eos_table.rhoMax );
    return -1;
  }

  n = (double *)malloc( EOS_NSPECIES * sizeof( double ) );
  ni = 0.0;
  ne = 0.0;
  nn = 0.0;
  for ( i=0; i<EOS_NSPECIES; i++ ) {
    n[i] = xnuc[i] * rho / PROTONMASS / eos_table.nuclearmasses[i];
    ni += n[i];
    nn += n[i] * eos_table.nuclearmasses[i];
    ne += n[i] * eos_table.nuclearcharges[i];
  }

  if (rho <= eos_table.rhoMin) {
    /* assume ideal fully ionised gas */
    *temp = 2.0/3.0 * e * rho / (ni+ne) / BOLTZMANN;
    *p = (ni + ne) * BOLTZMANN * (*temp);
    return 0;
  }

  ye = ne / nn;

  if (ye < eos_table.yeMin || ye > eos_table.yeMax) {
    printf( "Electron fraction out of table, ye: %g, table: [%g,%g]\n", ye, eos_table.yeMin, eos_table.yeMax );
    return -1;
  }

  _diedt = 1.5 * ni * BOLTZMANN / rho;

  if (*temp == 1.0) {
    /* lets make a guess for the temperatures, assuming an ideal gas */
    _temp = 2.0/3.0 * e * rho / (ni+ne) / BOLTZMANN;
  } else {
    _temp = *temp;
  }

  _tempold = 0.0;
  iter = 0;
  while (iter < EOS_MAXITER) {
    eos_trilinear_e( _temp, rho, ye, &_ee, &_deedt );

    _ie = _diedt * _temp;
    _e = _ee + _ie;
    if ( fabs(_e - e) <= (EOS_EPS * e) ) {
      break;
    }

    _tempold = _temp;
    _temp = _temp - (_e - e) / (_deedt + _diedt);
    if (_temp <= eos_table.tempMin && _tempold == eos_table.tempMin) {
      _temp = eos_table.tempMin;
      break;
    }

    if (_temp >= eos_table.tempMax && _tempold >= eos_table.tempMax) {
      _temp = eos_table.tempMax;
      break;
    }

    _temp = _temp < eos_table.tempMax ?  _temp : eos_table.tempMax;
    _temp = _temp > eos_table.tempMin ?  _temp : eos_table.tempMin;
    iter++;
  }

  eos_trilinear( _temp, rho, ye, &_ee, &_deedt, &_ep, &_dpedr, &_dpedt );

  *temp = _temp;
  *p = ni * BOLTZMANN * (*temp) + _ep;
  dpdt = ni * BOLTZMANN + _dpedt;
  dedt = 1.5 * ni * BOLTZMANN / rho + _deedt;
  *dpdr = ni * BOLTZMANN * (*temp) / rho + _dpedr + (*temp) * dpdt * dpdt / dedt / rho / rho;

  free( n );
  return 0;
}

int eos_calc_tgiven( double rho, double *xnuc, double temp, double *e, double *dedt ) {
  double *n, ni, ne, nn;        /* number densities [ particles per volume ] */
  double ye;                    /* electron fraction ( electrons per baryon ) */
  double _ee;                   /* energy per mass of electrons + positrons + radiation */
  double _deedt;                /* derivative of energy per mass of electrons + positrons + radiation with respect to temperature */
  int i;

  if (temp == 0) return -1;

  if (rho > eos_table.rhoMax) {
    printf( "Density exceeds allowed maximum, rho: %g, maximum: %g\n", rho, eos_table.rhoMax );
    return -1;
  }

  n = (double *)malloc( EOS_NSPECIES * sizeof( double ) );
  ni = 0.0;
  ne = 0.0;
  nn = 0.0;
  for ( i=0; i<EOS_NSPECIES; i++ ) {
    n[i] = xnuc[i] * rho / PROTONMASS / eos_table.nuclearmasses[i];
    ni += n[i];
    nn += n[i] * eos_table.nuclearmasses[i];
    ne += n[i] * eos_table.nuclearcharges[i];
  }

  if (rho <= eos_table.rhoMin) {
    /* assume ideal fully ionised gas */
    *e = 1.5 * temp / rho * (ni+ne) * BOLTZMANN;
    return 0;
  }

  ye = ne / nn;

  if (ye < eos_table.yeMin || ye > eos_table.yeMax) {
    printf( "Electron fraction out of table, ye: %g, table: [%g,%g]\n", ye, eos_table.yeMin, eos_table.yeMax );
    return -1;
  }

  eos_trilinear_e( temp, rho, ye, &_ee, &_deedt );
  *e = 1.5 * ni * BOLTZMANN * temp / rho + _ee;
  *dedt = 1.5 * ni * BOLTZMANN / rho + _deedt;

  free( n );
  return 0;
}

void eos_trilinear_e( double temp, double rho, double ye, double *e, double *dedt ) {
  double logtemp, logrho;
  int itemp, irho, iye;
  unsigned int entry111, entry211, entry121, entry221, entry112, entry212, entry122, entry222;
  double w111, w211, w121, w221, w112, w212, w122, w222;
  double dx1, dx2, dy1, dy2, dz1, dz2;

  logtemp = log10( temp );
  logrho = log10( rho );

  itemp = (logtemp-eos_table.ltempMin) / eos_table.ltempDelta;
  irho = (logrho-eos_table.lrhoMin) / eos_table.lrhoDelta;
  iye = (ye-eos_table.yeMin) / eos_table.yeDelta;

  dx1 = (logtemp - eos_table.ltemp[ itemp ]) * eos_table.ltempDeltaI;
  dy1 = (logrho - eos_table.lrho[ irho ]) * eos_table.lrhoDeltaI;
  dz1 = (ye - eos_table.ye[ iye ]) * eos_table.yeDeltaI;
  dx2 = 1. - dx1;
  dy2 = 1. - dy1;
  dz2 = 1. - dz1;

  w111 = dx2 * dy2 * dz2;
  w112 = dx2 * dy2 * dz1;
  w121 = dx2 * dy1 * dz2;
  w122 = dx2 * dy1 * dz1;
  w211 = dx1 * dy2 * dz2;
  w212 = dx1 * dy2 * dz1;
  w221 = dx1 * dy1 * dz2;
  w222 = dx1 * dy1 * dz1;

  entry111 = ( iye    *eos_table.nrho + irho   )*eos_table.ntemp + itemp;
  entry112 = ( (iye+1)*eos_table.nrho + irho   )*eos_table.ntemp + itemp;
  entry121 = ( iye    *eos_table.nrho + irho+1 )*eos_table.ntemp + itemp;
  entry122 = ( (iye+1)*eos_table.nrho + irho+1 )*eos_table.ntemp + itemp;
  entry211 = ( iye    *eos_table.nrho + irho   )*eos_table.ntemp + itemp+1;
  entry212 = ( (iye+1)*eos_table.nrho + irho   )*eos_table.ntemp + itemp+1;
  entry221 = ( iye    *eos_table.nrho + irho+1 )*eos_table.ntemp + itemp+1;
  entry222 = ( (iye+1)*eos_table.nrho + irho+1 )*eos_table.ntemp + itemp+1;

  *e = eos_table.e[ entry111 ]*w111 + eos_table.e[ entry112 ]*w112 +
       eos_table.e[ entry121 ]*w121 + eos_table.e[ entry122 ]*w122 +
       eos_table.e[ entry211 ]*w211 + eos_table.e[ entry212 ]*w212 +
       eos_table.e[ entry221 ]*w221 + eos_table.e[ entry222 ]*w222;

  *dedt = eos_table.dedt[ entry111 ]*w111 + eos_table.dedt[ entry112 ]*w112 +
          eos_table.dedt[ entry121 ]*w121 + eos_table.dedt[ entry122 ]*w122 +
          eos_table.dedt[ entry211 ]*w211 + eos_table.dedt[ entry212 ]*w212 +
          eos_table.dedt[ entry221 ]*w221 + eos_table.dedt[ entry222 ]*w222;
}

void eos_trilinear( double temp, double rho, double ye, double *e, double *dedt, double *p, double *dpdr, double *dpdt ) {
  double logtemp, logrho;
  int itemp, irho, iye;
  unsigned int entry111, entry211, entry121, entry221, entry112, entry212, entry122, entry222;
  double w111, w211, w121, w221, w112, w212, w122, w222;
  double dx1, dx2, dy1, dy2, dz1, dz2;

  logtemp = log10( temp );
  logrho = log10( rho );

  itemp = (logtemp-eos_table.ltempMin) / eos_table.ltempDelta;
  irho = (logrho-eos_table.lrhoMin) / eos_table.lrhoDelta;
  iye = (ye-eos_table.yeMin) / eos_table.yeDelta;
  
  dx1 = (logtemp - eos_table.ltemp[ itemp ]) * eos_table.ltempDeltaI;
  dy1 = (logrho - eos_table.lrho[ irho ]) * eos_table.lrhoDeltaI;
  dz1 = (ye - eos_table.ye[ iye ]) * eos_table.yeDeltaI;
  dx2 = 1. - dx1;
  dy2 = 1. - dy1;
  dz2 = 1. - dz1;

  w111 = dx2 * dy2 * dz2;
  w112 = dx2 * dy2 * dz1;
  w121 = dx2 * dy1 * dz2;
  w122 = dx2 * dy1 * dz1;
  w211 = dx1 * dy2 * dz2;
  w212 = dx1 * dy2 * dz1;
  w221 = dx1 * dy1 * dz2;
  w222 = dx1 * dy1 * dz1;

  entry111 = ( iye    *eos_table.nrho + irho   )*eos_table.ntemp + itemp;
  entry112 = ( (iye+1)*eos_table.nrho + irho   )*eos_table.ntemp + itemp;
  entry121 = ( iye    *eos_table.nrho + irho+1 )*eos_table.ntemp + itemp;
  entry122 = ( (iye+1)*eos_table.nrho + irho+1 )*eos_table.ntemp + itemp;
  entry211 = ( iye    *eos_table.nrho + irho   )*eos_table.ntemp + itemp+1;
  entry212 = ( (iye+1)*eos_table.nrho + irho   )*eos_table.ntemp + itemp+1;
  entry221 = ( iye    *eos_table.nrho + irho+1 )*eos_table.ntemp + itemp+1;
  entry222 = ( (iye+1)*eos_table.nrho + irho+1 )*eos_table.ntemp + itemp+1;

  *e = eos_table.e[ entry111 ]*w111 + eos_table.e[ entry112 ]*w112 +
       eos_table.e[ entry121 ]*w121 + eos_table.e[ entry122 ]*w122 +
       eos_table.e[ entry211 ]*w211 + eos_table.e[ entry212 ]*w212 +
       eos_table.e[ entry221 ]*w221 + eos_table.e[ entry222 ]*w222;

  *dedt = eos_table.dedt[ entry111 ]*w111 + eos_table.dedt[ entry112 ]*w112 +
          eos_table.dedt[ entry121 ]*w121 + eos_table.dedt[ entry122 ]*w122 +
          eos_table.dedt[ entry211 ]*w211 + eos_table.dedt[ entry212 ]*w212 +
          eos_table.dedt[ entry221 ]*w221 + eos_table.dedt[ entry222 ]*w222;

  *p = eos_table.p[ entry111 ]*w111 + eos_table.p[ entry112 ]*w112 +
       eos_table.p[ entry121 ]*w121 + eos_table.p[ entry122 ]*w122 +
       eos_table.p[ entry211 ]*w211 + eos_table.p[ entry212 ]*w212 +
       eos_table.p[ entry221 ]*w221 + eos_table.p[ entry222 ]*w222;

  *dpdr = eos_table.dpdr[ entry111 ]*w111 + eos_table.dpdr[ entry112 ]*w112 +
          eos_table.dpdr[ entry121 ]*w121 + eos_table.dpdr[ entry122 ]*w122 +
          eos_table.dpdr[ entry211 ]*w211 + eos_table.dpdr[ entry212 ]*w212 +
          eos_table.dpdr[ entry221 ]*w221 + eos_table.dpdr[ entry222 ]*w222;

  *dpdt = eos_table.dpdt[ entry111 ]*w111 + eos_table.dpdt[ entry112 ]*w112 +
          eos_table.dpdt[ entry121 ]*w121 + eos_table.dpdt[ entry122 ]*w122 +
          eos_table.dpdt[ entry211 ]*w211 + eos_table.dpdt[ entry212 ]*w212 +
          eos_table.dpdt[ entry221 ]*w221 + eos_table.dpdt[ entry222 ]*w222;
}

double eos_SwapDouble( double Val ) {
  double nVal;
  int i;
  const char* readFrom = ( const char* ) &Val;
  char * writeTo = ( ( char* ) &nVal ) + sizeof( nVal );
  for (i=0; i<sizeof( Val ); ++i) {
    *( --writeTo ) = *( readFrom++ );
  }
  return nVal;
}

int eos_SwapInt( int Val ) {
  int nVal;
  int i;
  const char* readFrom = ( const char* ) &Val;
  char * writeTo = ( ( char* ) &nVal ) + sizeof( nVal );
  for (i=0; i<sizeof( Val ); ++i) {
    *( --writeTo ) = *( readFrom++ );
  }
  return nVal;
}

double eos_calcYe( double *xnuc ) {
  double ye, atot, ztot;
  int i;

  atot = 0.0;
  ztot = 0.0;
  for (i=0; i<EOS_NSPECIES; i++) {
    atot += xnuc[i] * eos_table.nuclearmasses[i];
    ztot += xnuc[i] * eos_table.nuclearcharges[i];
  }
  ye = ztot / atot;  

  return ye;
}

void eos_checkswap( char* fname, int *swap ) {
  FILE *fd;
  size_t fsize, fpos;
  int blocksize, blockend;

  if (!(fd = fopen(fname, "r"))) {
    printf( "can't open file `%s' for reading eos table.\n", fname );
    EXIT;
  }

  fseek( fd, 0, SEEK_END );
  fsize = ftell( fd );

  *swap = 0;
  fpos = 0;
  fseek( fd, 0, SEEK_SET );
  fread( &blocksize, sizeof(int), 1, fd );
  while (!feof(fd)) {
    if (fpos + blocksize + 4 > fsize) {
      *swap += 1;
      break;
    }
    fpos += 4 + blocksize;
    fseek( fd, fpos, SEEK_SET );
    fread( &blockend, sizeof(int), 1, fd );
    if (blocksize != blockend) {
      *swap += 1;
      break;
    }
    fpos += 4;
    fread( &blocksize, sizeof(int), 1, fd );
  }

  if (*swap == 0) { fclose( fd ); return; }

  fpos = 0;
  fseek( fd, 0, SEEK_SET );
  fread( &blocksize, sizeof(int), 1, fd );
  while (!feof(fd)) {
    blocksize = eos_SwapInt( blocksize );
    if (fpos + blocksize + 4 > fsize) {
      *swap += 1;
      break;
    }
    fpos += 4 + blocksize;
    fseek( fd, fpos, SEEK_SET ); 
    fread( &blockend, sizeof(int), 1, fd );
    blockend = eos_SwapInt( blockend );
    if (blocksize != blockend) {
      *swap += 1;
      break;
    }
    fpos += 4;
    fread( &blocksize, sizeof(int), 1, fd );
  }

  fclose( fd );
}

#endif
