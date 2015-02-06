#ifdef NUCLEAR_NETWORK

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "eos.h"
#include "network_thiel.h"
#include "network_solver.h"

#ifndef STANDALONE
#include "allvars.h"
#include "proto.h"
#endif

#ifndef STANDALONE
#define EXIT endrun(1338)
#else
#define EXIT exit(1338)
#endif

/*
  routines to read the thielemann library
 */

void network_init() {
	int swap, skip;
	FILE *f_data, *f_rates, *f_nucdata;
	int nm;
	char *nuc_names;
	int nrates[3];	/* amount of one, two and three body reaction rates */
	int i, j, idummy;
	char cdummy[100];
	float fdummy, na, sp, me;
	
#ifndef STANDALONE
  if (ThisTask == 0) {
#endif

#ifndef STANDALONE
		sprintf( cdummy, "%s/nets4", All.NetworkData );
		network_checkswap( cdummy, &swap );
		f_data  = fopen( cdummy, "rb" );
		sprintf( cdummy, "%s/nets3", All.NetworkData );
		f_rates = fopen( cdummy, "rb" );
		sprintf( cdummy, "%s/netwinv", All.NetworkData );
		f_nucdata = fopen( cdummy, "rb" );
#else
		network_checkswap( "../netiso13/nets4", &swap );
		f_data  = fopen( "../netiso13/nets4", "rb" );
		f_rates = fopen( "../netiso13/nets3", "rb" );
		f_nucdata = fopen( "../netiso13/netwinv", "rb" );
#endif
		
		fseek( f_data, 0, SEEK_SET );
		fseek( f_rates, 0, SEEK_SET );
		fseek( f_nucdata, 0, SEEK_SET );
		
		fread( &skip, sizeof(int), 1, f_data );
		fread( &network_data.nuc_count, sizeof(int), 1, f_data );
		fread( &skip, sizeof(int), 1, f_data );
		if (swap) network_data.nuc_count = network_SwapInt( network_data.nuc_count );	
		
		printf( "Network: Loading %d nuclear species for network.\n", network_data.nuc_count );
		
		nuc_names = (char*)malloc( 5 * network_data.nuc_count * sizeof(char) );
		fread( &skip, sizeof(int), 1, f_data );
		fread( nuc_names, sizeof(char), 5*network_data.nuc_count, f_data );
		fread( &skip, sizeof(int), 1, f_data );
		
		fread( &skip, sizeof(int), 1, f_data );
		fread( &nm, sizeof(int), 1, f_data );
		fread( &skip, sizeof(int), 1, f_data );
		if (swap) nm = network_SwapInt( nm );
		
		fread( &skip, sizeof(int), 1, f_data );
		fread( &nrates, sizeof(int), 3, f_data );
		fread( &skip, sizeof(int), 1, f_data );
		if (swap) for (i=0; i<3; i++) nrates[i] = network_SwapInt( nrates[i] );
		
		printf( "Network: Loading %d|%d|%d 1b|2b|3b reaction rates.\n", nrates[0], nrates[1], nrates[2] );		
#ifndef STANDALONE
	}
#endif

#ifndef STANDALONE
	MPI_Bcast( &network_data.nuc_count, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( nrates, 3, MPI_INT, 0, MPI_COMM_WORLD );
#endif

	network_rates_1b = (struct network_ratedata_1b*)malloc( nrates[0] * sizeof( struct network_ratedata_1b ) );
	network_rates_2b = (struct network_ratedata_2b*)malloc( nrates[1] * sizeof( struct network_ratedata_2b ) );
	network_rates_3b = (struct network_ratedata_3b*)malloc( nrates[2] * sizeof( struct network_ratedata_3b ) );
	
#ifndef STANDALONE
	if (ThisTask == 0) {
#endif		

		for (i=0; i<nrates[0]; i++) {
			fread( &skip, sizeof(int), 1, f_rates );
			fread( &network_rates_1b[i], sizeof(struct network_ratedata_1b), 1, f_rates );
			fread( &skip, sizeof(int), 1, f_rates );
			if (swap) {
				network_rates_1b[i].n = network_SwapInt( network_rates_1b[i].n );
				for (j=0; j<3; j++) network_rates_1b[i].ni[j] = network_SwapInt( network_rates_1b[i].ni[j] );
				network_rates_1b[i].ic = network_SwapInt( network_rates_1b[i].ic );
				network_rates_1b[i].ir = network_SwapInt( network_rates_1b[i].ir );
				network_rates_1b[i].iv = network_SwapInt( network_rates_1b[i].iv );
				for (j=0; j<7; j++) network_rates_1b[i].p[j] = network_SwapDouble( network_rates_1b[i].p[j] );
				network_rates_1b[i].q = network_SwapDouble( network_rates_1b[i].q );
			}
			for (j=0; j<3; j++) network_rates_1b[i].ni[j] -= 1; /* fortran arrays start with index 1 */
		}
		
		for (i=0; i<nrates[1]; i++) {
			fread( &skip, sizeof(int), 1, f_rates );
			fread( &network_rates_2b[i], sizeof(struct network_ratedata_2b), 1, f_rates );
			fread( &skip, sizeof(int), 1, f_rates );
			if (swap) {
				network_rates_2b[i].n = network_SwapInt( network_rates_2b[i].n );
				for (j=0; j<4; j++) network_rates_2b[i].ni[j] = network_SwapInt( network_rates_2b[i].ni[j] );
				network_rates_2b[i].ic = network_SwapInt( network_rates_2b[i].ic );
				network_rates_2b[i].ir = network_SwapInt( network_rates_2b[i].ir );
				network_rates_2b[i].iv = network_SwapInt( network_rates_2b[i].iv );
				for (j=0; j<7; j++) network_rates_2b[i].p[j] = network_SwapDouble( network_rates_2b[i].p[j] );
				network_rates_2b[i].dummy = network_SwapDouble( network_rates_2b[i].dummy );
			}
			for (j=0; j<4; j++) network_rates_2b[i].ni[j] -= 1; /* fortran arrays start with index 1 */
		}
		
		for (i=0; i<nrates[2]; i++) {
			fread( &skip, sizeof(int), 1, f_rates );
			fread( &network_rates_3b[i], sizeof(struct network_ratedata_3b), 1, f_rates );
			fread( &skip, sizeof(int), 1, f_rates );
			if (swap) {
				network_rates_3b[i].n = network_SwapInt( network_rates_3b[i].n );
				for (j=0; j<3; j++) network_rates_3b[i].ni[j] = network_SwapInt( network_rates_3b[i].ni[j] );
				network_rates_3b[i].ic = network_SwapInt( network_rates_3b[i].ic );
				network_rates_3b[i].ir = network_SwapInt( network_rates_3b[i].ir );
				network_rates_3b[i].iv = network_SwapInt( network_rates_3b[i].iv );
				for (j=0; j<7; j++) network_rates_3b[i].p[j] = network_SwapDouble( network_rates_3b[i].p[j] );
				network_rates_3b[i].dummy = network_SwapDouble( network_rates_3b[i].dummy );
			}
			for (j=0; j<3; j++) network_rates_3b[i].ni[j] -= 1; /* fortran arrays start with index 1 */
		}
#ifndef STANDALONE
	}
#endif

#ifndef STANDALONE
	MPI_Bcast( network_rates_1b, nrates[0]*sizeof(struct network_ratedata_1b), MPI_BYTE, 0, MPI_COMM_WORLD );
	MPI_Bcast( network_rates_2b, nrates[1]*sizeof(struct network_ratedata_2b), MPI_BYTE, 0, MPI_COMM_WORLD );
	MPI_Bcast( network_rates_3b, nrates[2]*sizeof(struct network_ratedata_3b), MPI_BYTE, 0, MPI_COMM_WORLD );
#endif

	network_nucdata = (struct network_nucdata*)malloc( network_data.nuc_count * sizeof( struct network_nucdata ) );
	
	for (i=0; i<network_data.nuc_count; i++) {
		network_nucdata[i].la = (int*)malloc( 3 * sizeof( int ) );
		network_nucdata[i].le = (int*)malloc( 3 * sizeof( int ) );
#ifndef STANDALONE
		if (ThisTask == 0) {
#endif
			fread( &skip, sizeof(int), 1, f_data );
			fread( &idummy, sizeof(int), 1, f_data );
			for (j=0; j<3; j++) {
				fread( &network_nucdata[i].la[j], sizeof(int), 1, f_data );
				fread( &network_nucdata[i].le[j], sizeof(int), 1, f_data );
				if (swap) {
					network_nucdata[i].la[j] = network_SwapInt( network_nucdata[i].la[j] );
					network_nucdata[i].le[j] = network_SwapInt( network_nucdata[i].le[j] );
				}
				network_nucdata[i].la[j] -= 1; /* fortran arrays start with index 1 */
				network_nucdata[i].le[j] -= 1; /* fortran arrays start with index 1 */
			}
			fread( &skip, sizeof(int), 1, f_data );
#ifndef STANDALONE
		}
		MPI_Bcast( network_nucdata[i].la, 3, MPI_INT, 0, MPI_COMM_WORLD );
		MPI_Bcast( network_nucdata[i].le, 3, MPI_INT, 0, MPI_COMM_WORLD );
#endif
	}

	for (i=0; i<3; i++) {
		network_data.a[i] = (double*)malloc( (network_nucdata[network_data.nuc_count-1].le[i]+1) * sizeof( double ) );
		network_data.mu[i] = (int*)malloc( (network_nucdata[network_data.nuc_count-1].le[i]+1) * sizeof( int ) );
#ifndef STANDALONE
		if (ThisTask == 0) {
#endif
			network_data.mumax[i] = 0;
			for (j=0; j<network_nucdata[network_data.nuc_count-1].le[i]+1; j++) {
				fread( &skip, sizeof(int), 1, f_rates );
				fread( &network_data.a[i][j], sizeof(double), 1, f_rates );
				fread( &network_data.mu[i][j], sizeof(int), 1, f_rates );
				fread( &skip, sizeof(int), 1, f_rates );
				if (swap) {
					network_data.a[i][j] = network_SwapDouble( network_data.a[i][j] );
					network_data.mu[i][j] = network_SwapInt( network_data.mu[i][j] );
				}
				network_data.mu[i][j] -= 1; /* fortran arrays start with index 1 */
				
				if (network_data.mu[i][j] > network_data.mumax[i]) {
					network_data.mumax[i] = network_data.mu[i][j];
				}
			}
#ifndef STANDALONE
		}
		MPI_Bcast( network_data.a[i],  network_nucdata[network_data.nuc_count-1].le[i]+1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		MPI_Bcast( network_data.mu[i], network_nucdata[network_data.nuc_count-1].le[i]+1, MPI_INT, 0, MPI_COMM_WORLD );
		MPI_Bcast( &network_data.mumax[i], 1, MPI_INT, 0, MPI_COMM_WORLD );
#endif
	}
		
#ifndef STANDALONE
	if (ThisTask == 0) {
#endif
		fscanf( f_nucdata, "%5c", (char*)&cdummy );	/* skip line */
		/* temperature table */
		for (i=0; i<24; i++) {
			fscanf( f_nucdata, "%3d", &idummy );
			network_data.t9[i] = (double)idummy * 0.01;
		}
		network_data.t9[23] = network_data.t9[23] * 10.0;
		
		for (i=0; i<network_data.nuc_count; i++) {
			fscanf( f_nucdata, "%s", network_nucdata[i].name );
			network_nucdata[i].name[5] = 0;
		}
		
		for (i=0; i<network_data.nuc_count; i++) {
			fscanf( f_nucdata, "%5s%f%d%d%f%f", cdummy, &na, &network_nucdata[i].np, &network_nucdata[i].nn, &sp, &me );
					
			network_nucdata[i].na = (double)na;	/* atomic weight */
			network_nucdata[i].bi = 8.07144 * network_nucdata[i].nn + 7.28899 * network_nucdata[i].np - me; /* constants are mass excess of neutron and proton in MeV; me = binding energy ? */
			network_nucdata[i].exm = me;
			network_nucdata[i].angm = 2.0 * sp + 1.0; /* spin ? */
			/* skip dummy, its again the name of the nucleus */

			for (j=0; j<24; j++) {
				fscanf( f_nucdata, "%f", &fdummy );
				network_nucdata[i].dlg[j] = log( (double)fdummy );
			}
		}
		
		fclose( f_data );
		fclose( f_rates );
		fclose( f_nucdata );
#ifndef STANDALONE
	}

	MPI_Bcast( network_data.t9, 24, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	for (i=0; i<network_data.nuc_count; i++) {
		MPI_Bcast( &network_nucdata[i].na, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		MPI_Bcast( &network_nucdata[i].nn, 1, MPI_INT, 0, MPI_COMM_WORLD );
		MPI_Bcast( &network_nucdata[i].np, 1, MPI_INT, 0, MPI_COMM_WORLD );
		MPI_Bcast( &network_nucdata[i].bi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		MPI_Bcast( &network_nucdata[i].exm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		MPI_Bcast( &network_nucdata[i].angm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		MPI_Bcast( network_nucdata[i].dlg, 24, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		MPI_Bcast( network_nucdata[i].name, 5, MPI_BYTE, 0, MPI_COMM_WORLD );
	}
#endif

	/* allocate some memory that is needed to calculate rates */
	network_data.gg = (double*)malloc( network_data.nuc_count * sizeof(double) );
	network_data.sig1b = (double*)malloc( (network_data.mumax[0]+1) * sizeof(double) );
	network_data.sig2b = (double*)malloc( (network_data.mumax[1]+1) * sizeof(double) );
	network_data.sig3b = (double*)malloc( (network_data.mumax[2]+1) * sizeof(double) );
	network_data.x = (double*)malloc( network_data.nuc_count * sizeof(double) );
	network_data.y = (double*)malloc( (network_data.nuc_count+1) * sizeof(double) );
	network_data.conv = 1.602177e-12 * 1.0e6 * 6.0221367e23; /* eV2erg * 1.0e6 * avogadro */
	network_data.na = (double*)malloc( network_data.nuc_count * sizeof(double) );
	for (i=0; i<network_data.nuc_count; i++) { network_data.na[i] = network_nucdata[i].na; }

#ifndef FIXED_TEMPERATURE
	network_solver_init( &network_getrhs, &network_getjacob, 1e-6, network_data.nuc_count+1, network_data.nuc_count, network_data.na );
#else
	network_solver_init( &network_getrhs, &network_getjacob, 1e-6, network_data.nuc_count, network_data.nuc_count, network_data.na );
#endif
}

/* do not call this routine if you want to run test cases with constant temperature */
void network_integrate( double temp, double rho, double *x, double *dx, double dt, double *dedt ) {
	double *y;
	int i;

	if (dt == 0) {
	  for (i=0; i<network_data.nuc_count; i++) dx[i] = 0;
	  *dedt = 0;
	  return;
	}
	
	/* calculate number densities */
	y = network_data.y;
	for (i=0; i<network_data.nuc_count; i++) {
		y[i] = x[i] / network_nucdata[i].na;
	}
	y[network_data.nuc_count] = temp;
	
	/* run network */
	network_solver_integrate( temp, rho, y, dt );
	
	/* calculate change of mass fractions and energy release */
	*dedt = 0;
	for (i=0; i<network_data.nuc_count; i++) {
	  dx[i] = ( y[i] * network_nucdata[i].na - x[i] ) / dt;
	  *dedt -= dx[i] / network_nucdata[i].na * network_nucdata[i].exm;
	}
	*dedt *= network_data.conv;
}

void network_part( double temp ) {
	/* interpolates partition functions, given the temperature */
	int index, i;
	double tempLeft, tempRight;
	double dlgLeft, dlgRight;
	double grad;
	
	index = 0;
	temp = min( max( temp/1e9, network_data.t9[0] ), network_data.t9[23] );
	
	while (temp > network_data.t9[index]) {
		index++;
	}
	index--;
	
	tempLeft  = network_data.t9[index];
	tempRight = network_data.t9[index+1];
    
    for (i=0; i<network_data.nuc_count; i++) {
    	dlgLeft = network_nucdata[i].dlg[index];
    	dlgRight = network_nucdata[i].dlg[index+1];
    	
    	grad = (dlgRight-dlgLeft) / (tempRight-tempLeft);
    	network_data.gg[i] = exp( dlgLeft + (temp - tempLeft)*grad );
    }
}

void network_getrates( double temp, double rho ) {
	double t9, t9l, temp9[9];
	int i, j;
	double sig;
	
	network_part( temp );
	
	t9 = temp * 1.0e-9;
	t9l = log( t9 );
	
	temp9[0] = 1.0;
	temp9[1] = 1.0 / t9;
	temp9[2] = exp( -t9l / 3.0 );
	temp9[3] = 1.0 / temp9[2];
	temp9[4] = t9;
	temp9[5] = exp( t9l * 5.0 / 3.0 );
	temp9[6] = t9l;
	
	/* 1 body reactions */
	for (i=0; i<=network_data.mumax[0]; i++) {
		sig = 0;
		for (j=0; j<7; j++) {
			sig += temp9[j] * network_rates_1b[i].p[j];
		}
		network_data.sig1b[i] = exp( sig );
	}
	
	for (i=0; i<=network_data.mumax[0]; i++) {
		/* account for inverse rate */
		if (network_rates_1b[i].iv == 1) {
			network_data.sig1b[i] *= network_data.gg[ network_rates_1b[i].ni[1] ] * network_data.gg[ network_rates_1b[i].ni[2] ] / network_data.gg[ network_rates_1b[i].ni[0] ];
		}
		
		/* we do not use weak rates at the moment */
		if (network_rates_1b[i].ic > 0) {
			network_data.sig1b[i] = 0;
		}
	}

	/* 2 body reactions */
	for (i=0; i<=network_data.mumax[1]; i++) {
		sig = 0;
		for (j=0; j<7; j++) {
			sig += temp9[j] * network_rates_2b[i].p[j];
		}
		network_data.sig2b[i] = exp( sig ) * rho;
	}
	
	for (i=0; i<=network_data.mumax[1]; i++) {
		/* account for inverse rate */
		if (network_rates_2b[i].iv == 1) {
			network_data.sig2b[i] *= network_data.gg[ network_rates_2b[i].ni[2] ] * network_data.gg[ network_rates_2b[i].ni[3] ] / (network_data.gg[ network_rates_2b[i].ni[0] ] * network_data.gg[ network_rates_2b[i].ni[1] ] );
		}
		
		/* we do not use weak rates at the moment */
		if (network_rates_2b[i].ic > 0) {
			network_data.sig2b[i] = 0;
		}
	}
	
	/* 3 body reactions */
	for (i=0; i<=network_data.mumax[2]; i++) {
		sig = 0;
		for (j=0; j<7; j++) {
			sig += temp9[j] * network_rates_3b[i].p[j];
		}
		network_data.sig3b[i] = exp( sig ) * rho * rho;
	}
	
	/* inverse reactions do not occur... */
	
	for (i=0; i<=network_data.mumax[2]; i++) {
		/* we do not use weak rates at the moment */
		if (network_rates_3b[i].ic > 0) {
			network_data.sig3b[i] = 0;
		}
	}
}

/* calculates the right hand dy/dt of the set of linear equations */
void network_getrhs( double temp, double rho, double *y, double *rhs ) {
	int iTemp;
	int i, j;
	double deriv;
	double e, dedT, dy;
	double newTemp, p, dpdr;

	for (i=0; i<network_data.nuc_count; i++) {
		if (y[i] < 1e-30) y[i] = 1e-30;
		if (y[i] > 1.0) y[i] = 1.0;
	}

#ifndef FIXED_TEMPERATURE	
	/* i = network_solver_data.nelements */
    y[network_data.nuc_count] = max( 1e7, min( y[network_data.nuc_count], 1e10 ) );
    
    temp = y[network_data.nuc_count];
#endif
	
	network_getrates( temp, rho );

	/* dy_i/dt */
	for (i=0; i<network_data.nuc_count; i++) {
		deriv = 0;
		
		/* 1 body reactions */
		for (j=network_nucdata[i].la[0]; j<= network_nucdata[i].le[0]; j++) {
			deriv += network_data.a[0][j] * network_data.sig1b[ network_data.mu[0][j] ] * y[ network_rates_1b[ network_data.mu[0][j] ].ni[ 0 ] ];
		}
		
		/* 2 body reactions */
		for (j=network_nucdata[i].la[1]; j<= network_nucdata[i].le[1]; j++) {
			deriv += network_data.a[1][j] * network_data.sig2b[ network_data.mu[1][j] ] * y[ network_rates_2b[ network_data.mu[1][j] ].ni[ 0 ] ]
			                                                                            * y[ network_rates_2b[ network_data.mu[1][j] ].ni[ 1 ] ];
		}
		
		/* 3 body reactions */
		for (j=network_nucdata[i].la[2]; j<= network_nucdata[i].le[2]; j++) {
			deriv += network_data.a[2][j] * network_data.sig3b[ network_data.mu[2][j] ] * y[ network_rates_3b[ network_data.mu[2][j] ].ni[ 0 ] ]
						                                                                * y[ network_rates_3b[ network_data.mu[2][j] ].ni[ 1 ] ]
			                                                                            * y[ network_rates_3b[ network_data.mu[2][j] ].ni[ 2 ] ];
		}
		
		rhs[i] = deriv;
	}

#ifndef FIXED_TEMPERATURE
	iTemp = network_data.nuc_count;
	
	/* dT/dt = dT/dE * dE/dt + sum_i ( dT/dy_i * dy_i/dt )
	 * dE/dt = sum_i ( ebind_i * dy_i/dt ) */
	for (i=0; i<network_data.nuc_count; i++) { 
		network_data.x[i] = y[i] * network_nucdata[i].na;
	}
	eos_calc_tgiven( rho, network_data.x, temp, &e, &dedT );
	rhs[iTemp] = 0;
	for (i=0; i<network_data.nuc_count; i++) {
		rhs[iTemp] += rhs[i] * network_nucdata[i].exm;
	}
	rhs[iTemp] *= -network_data.conv / dedT;

#ifndef NEGLECT_DTDY_TERMS
	/* sum_i ( dT/dy_i * dy_i/dt ) */
	for (i=0; i<network_data.nuc_count; i++) {
		dy = max( NETWORK_DIFFVAR, y[i]*NETWORK_DIFFVAR );
		network_data.x[i] = (y[i]+dy) * network_nucdata[i].na;
		newTemp = temp;
		eos_calc_egiven( rho, network_data.x, e, &newTemp, &p, &dpdr );
		rhs[iTemp] += (newTemp-temp) / dy * rhs[i];
		network_data.x[i] = y[i] * network_nucdata[i].na;
	}	
#endif
#endif
}

void network_getjacob( double temp, double rho, double *y, double *rhs, double *jacob ) {
	int iTemp, nMatrix;
	int i, j;
	double dTemp, yold, dy;
	double *deriv;

#ifndef FIXED_TEMPERATURE
	temp = y[network_data.nuc_count];
#endif
	
	network_getrates( temp, rho );
	
#ifndef FIXED_TEMPERATURE
	nMatrix = network_data.nuc_count + 1;
#else
	nMatrix = network_data.nuc_count;
#endif
	deriv = (double*)malloc( nMatrix * sizeof(double) );
	
	/* dy_i/dy_j */
	for (i=0; i<network_data.nuc_count; i++) {
		/* do row by row */
		for (j=0; j<network_data.nuc_count; j++) {
			deriv[j] = 0;
		}
		
		/* 1 body reactions */
		for (j=network_nucdata[i].la[0]; j<= network_nucdata[i].le[0]; j++) {
			deriv[ network_rates_1b[ network_data.mu[0][j] ].ni[0] ] += network_data.a[0][j] * network_data.sig1b[ network_data.mu[0][j] ];
		}
		
		/* 2 body reactions */
		for (j=network_nucdata[i].la[1]; j<= network_nucdata[i].le[1]; j++) {
			deriv[ network_rates_2b[ network_data.mu[1][j] ].ni[0] ] += network_data.a[1][j] * network_data.sig2b[ network_data.mu[1][j] ] * y[ network_rates_2b[ network_data.mu[1][j] ].ni[ 1 ] ];
			deriv[ network_rates_2b[ network_data.mu[1][j] ].ni[1] ] += network_data.a[1][j] * network_data.sig2b[ network_data.mu[1][j] ] * y[ network_rates_2b[ network_data.mu[1][j] ].ni[ 0 ] ];
		}
		
		/* 3 body reactions */
		for (j=network_nucdata[i].la[2]; j<= network_nucdata[i].le[2]; j++) {
			deriv[ network_rates_3b[ network_data.mu[2][j] ].ni[0] ] += network_data.a[2][j] * network_data.sig3b[ network_data.mu[2][j] ] * y[ network_rates_3b[ network_data.mu[2][j] ].ni[ 1 ] ] * y[ network_rates_3b[ network_data.mu[2][j] ].ni[ 2 ] ];
			deriv[ network_rates_3b[ network_data.mu[2][j] ].ni[0] ] += network_data.a[2][j] * network_data.sig3b[ network_data.mu[2][j] ] * y[ network_rates_3b[ network_data.mu[2][j] ].ni[ 0 ] ] * y[ network_rates_3b[ network_data.mu[2][j] ].ni[ 2 ] ];
			deriv[ network_rates_3b[ network_data.mu[2][j] ].ni[0] ] += network_data.a[2][j] * network_data.sig3b[ network_data.mu[2][j] ] * y[ network_rates_3b[ network_data.mu[2][j] ].ni[ 0 ] ] * y[ network_rates_3b[ network_data.mu[2][j] ].ni[ 1 ] ];
		}
		
		/* move entries into the matrix */
		for (j=0; j<network_data.nuc_count; j++) {
			jacob[ i*nMatrix + j ] = deriv[j];
		}
	}

#ifndef FIXED_TEMPERATURE
	iTemp = network_data.nuc_count;

	/* dy_i/dT & dT/dT */
	dTemp = max( fabs( temp ) * NETWORK_DIFFVAR, NETWORK_DIFFVAR );
	y[iTemp] = temp + dTemp;
	network_getrhs( temp + dTemp, rho, y, deriv );
	y[iTemp] = temp;
	for (i=0; i<nMatrix; i++) {
		jacob[ i*nMatrix + iTemp ] = ( deriv[i] - rhs[i] ) / dTemp;
	}

	/* dT/dy_i */
	for (i=0; i<network_data.nuc_count; i++) {
		yold = y[i];
		dy = max( fabs( yold ) * NETWORK_DIFFVAR, NETWORK_DIFFVAR );
		y[i] = yold + dy;
		network_getrhs( temp, rho, y, deriv );
		jacob[ iTemp*nMatrix + i ] = ( deriv[iTemp] - rhs[iTemp] ) / dy;
		y[i] = yold;
	}
#endif

	free( deriv );
}

double network_SwapDouble( double Val ) {
  double nVal;
  int i;
  const char* readFrom = ( const char* ) &Val;
  char * writeTo = ( ( char* ) &nVal ) + sizeof( nVal );
  for (i=0; i<sizeof( Val ); ++i) {
    *( --writeTo ) = *( readFrom++ );
  }
  return nVal;
}

int network_SwapInt( int Val ) {
  int nVal;
  int i;
  const char* readFrom = ( const char* ) &Val;
  char * writeTo = ( ( char* ) &nVal ) + sizeof( nVal );
  for (i=0; i<sizeof( Val ); ++i) {
    *( --writeTo ) = *( readFrom++ );
  }
  return nVal;
}

void network_checkswap( char* fname, int *swap ) {
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
    blocksize = network_SwapInt( blocksize );
    if (fpos + blocksize + 4 > fsize) {
      *swap += 1;
      break;
    }
    fpos += 4 + blocksize;
    fseek( fd, fpos, SEEK_SET ); 
    fread( &blockend, sizeof(int), 1, fd );
    blockend = network_SwapInt( blockend );
    if (blocksize != blockend) {
      *swap += 1;
      break;
    }
    fpos += 4;
    fread( &blocksize, sizeof(int), 1, fd );
  }

  fclose( fd );
}

#endif /* NUCLEAR_NETWORK */
