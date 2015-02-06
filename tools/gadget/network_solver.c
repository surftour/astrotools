#ifdef NUCLEAR_NETWORK

#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

#include "network_solver.h"
#include "network_thiel.h"
#include "eos.h"

#include "proto.h"
#include "allvars.h"

void network_solver_init( func_network_getrhs getrhs, func_network_getjacob getjacob, double tolerance, int matrixsize, int nelements, double *aion ) {
	int i, k;
	
	network_solver_data.getrhs = getrhs;
	network_solver_data.getjacob = getjacob;
	network_solver_data.tolerance = tolerance;
	network_solver_data.matrixsize = matrixsize;
	network_solver_data.matrixsize2 = matrixsize * matrixsize;
	network_solver_data.nelements = nelements;
	network_solver_data.aion = aion;
	
	network_solver_data.nsteps = 8;
	network_solver_data.maxstep = network_solver_data.nsteps-1;
	network_solver_data.steps = (int*)malloc( network_solver_data.nsteps * sizeof( int ) );
	network_solver_data.steps[0] =  2;
	network_solver_data.steps[1] =  6;
	network_solver_data.steps[2] = 10;
	network_solver_data.steps[3] = 14;
	network_solver_data.steps[4] = 22;
	network_solver_data.steps[5] = 34;
	network_solver_data.steps[6] = 50;
	network_solver_data.steps[7] = 70;
	
	network_solver_data.dy = (double*)malloc( network_solver_data.matrixsize * sizeof( double ) );
	network_solver_data.ynew = (double*)malloc( network_solver_data.matrixsize * sizeof( double ) );
	network_solver_data.yscale = (double*)malloc( network_solver_data.matrixsize * sizeof( double ) );
	network_solver_data.rhs = (double*)malloc( network_solver_data.matrixsize * sizeof( double ) );
	network_solver_data.jacob = (double*)malloc( network_solver_data.matrixsize2 * sizeof( double ) );
	network_solver_data.x = (double*)malloc( network_solver_data.maxstep * sizeof( double ) );
	network_solver_data.err = (double*)malloc( network_solver_data.maxstep * sizeof( double ) );
	network_solver_data.qcol = (double*)malloc( network_solver_data.matrixsize * network_solver_data.nsteps * sizeof( double ) );
	network_solver_data.d = (double*)malloc( network_solver_data.matrixsize * sizeof( double ) );
	network_solver_data.a = (double*)malloc( network_solver_data.nsteps * sizeof( double ) );
	network_solver_data.alf = (double*)malloc( (network_solver_data.nsteps-1) * (network_solver_data.nsteps-1) * sizeof( double ) );	
	
	/* performance checker */
	network_solver_data.a[0] = network_solver_data.steps[0] + 1.0;
	for (i=0; i<network_solver_data.maxstep; i++) network_solver_data.a[i+1] = network_solver_data.a[i] + network_solver_data.steps[i+1];

	for (i=1; i<network_solver_data.maxstep; i++) {
    	for (k=0; k<i; k++) {
      		network_solver_data.alf[i*network_solver_data.maxstep+k] = pow( 1.0e-6 * 0.25, (network_solver_data.a[k+1] - network_solver_data.a[i+1]) /
      			((network_solver_data.a[i+1]-network_solver_data.a[0] + 1.0) * (2.*(k+1.0) + 1.0)) );
    	}
	}

	network_solver_data.a[0] += network_solver_data.matrixsize;
	for (i=0; i<network_solver_data.maxstep; i++) network_solver_data.a[i+1] = network_solver_data.a[i] + network_solver_data.steps[i+1];

	for (i=1; i<network_solver_data.maxstep-1; i++) {
		if (network_solver_data.a[i+1] > network_solver_data.a[i]*network_solver_data.alf[i*network_solver_data.maxstep+i-1]) break;
	}
	network_solver_data.maxiter = i+1;

#ifndef STANDALONE
	if (ThisTask == 0)
#endif  
		printf( "Network solver initialization done\n" );
}

void network_solver_integrate( double temp, double rho, double *y, double dt ) {
  double *dy, *ynew, *yscale;
  double *x, *err, *qcol;
  double *a, *alf;
  double t, dttry, dtnext;
  double errmax, dum, red;
  double workmin, work, fac, scale;
  double sum;
  int i, m, iter, lastiter, iterm, iteropt;
  int first, reduce, step;
#ifdef NETWORK_OUTPUT_TIMEEVOLUTION
  FILE *fp;
  
  fp = fopen( "network_output.dat", "w" );
  fwrite( &network_solver_data.nelements, sizeof(int), 1, fp );
  fwrite( &network_solver_data.matrixsize, sizeof(int), 1, fp );
  for (i=0; i<network_solver_data.nelements; i++) fwrite( &network_solver_data.aion[i], sizeof(double), 1, fp );
#endif

  dy = network_solver_data.dy;
  ynew = network_solver_data.ynew;
  yscale = network_solver_data.yscale;
  x = network_solver_data.x;
  err = network_solver_data.err;
  qcol = network_solver_data.qcol;
  a = network_solver_data.a;
  alf = network_solver_data.alf;

  t = 0;
  step = 0;
  first = 1;
  dttry = dt;
  iteropt = network_solver_data.maxiter;

  while (t < dt) {
#ifdef NETWORK_OUTPUT_TIMEEVOLUTION
    printf( "t(%03d): %11.5e, dt (%02d): %11.5e, x:", step, t, m, dttry ); for (i=0; i<network_solver_data.nelements; i++) { printf( " %8.1e", y[i]*network_solver_data.aion[i] ); } if (network_solver_data.nelements<network_solver_data.matrixsize) { printf( " %8.1e", y[network_solver_data.nelements] ); } printf( "\n" );
  	fwrite( &t, sizeof(double), 1, fp );
  	fwrite( y, sizeof(double), network_solver_data.matrixsize, fp );
#endif
  	
    for (i=0; i<network_solver_data.nelements; i++) {
      if (y[i] > 1.0) y[i] = 1.0;
      if (y[i] < 1e-30) y[i] = 1e-30;
    }
    
	/* check temperature if present */
    if (network_solver_data.nelements < network_solver_data.matrixsize) {
		/* i = network_solver_data.nelements */
      	y[i] = max( 1e7, min( y[i], 1e10 ) );
	}

	sum = 0;
    for (i=0; i<network_solver_data.nelements; i++) { sum += y[i]*network_solver_data.aion[i]; }
    for (i=0; i<network_solver_data.nelements; i++) y[i] /= sum;

    for (i=0; i<network_solver_data.matrixsize; i++) {
      yscale[i] = max( fabs( y[i] ), network_solver_data.tolerance );
    }
 
    reduce = 0;
    for (iter=0; iter<network_solver_data.maxiter; iter++) {
      if (dttry == 0) {
        printf( "dt is zero.\n" );
        return;
      }

      if (reduce == 1) reduce = 2;

      network_solver_calc_dy( temp, rho, dttry, network_solver_data.steps[iter], y, ynew );
      x[iter] = ( dttry / network_solver_data.steps[iter] )*( dttry / network_solver_data.steps[iter] );
      network_solver_extrapolate( iter, x, ynew, dy, qcol );

      /* compute normalized error estimate */
      if (iter > 0) {
        errmax = 1.0e-30;
        for (i=0; i<network_solver_data.matrixsize; i++) {
          dum = fabs( dy[i] / yscale[i] );
          if (dum > errmax) {errmax = dum; m = i;}
        }

        errmax /= network_solver_data.tolerance;
        iterm = iter-1;
        err[iterm] = pow( errmax / 0.25, 1.0 / (2.*(iterm+1.)+1.) );
      }

      if (iter > 0 && ( iter >= iteropt-1 || first ) ) {
        /* if converged, leave */
        if (errmax < 1.0) break;

        if (iter == network_solver_data.maxiter-1 || iter == iteropt+1) {
          red = 0.7 / err[iterm];
          reduce = 1;
        } else if ( iter == iteropt ) {
          if ( alf[iteropt*network_solver_data.maxstep + iteropt-1] < err[iterm] ) {
            red = 1.0 / err[iterm];
            reduce = 1;
          }
        } else if ( iteropt == network_solver_data.maxiter-1 ) {
          if ( alf[(network_solver_data.maxiter-2)*network_solver_data.maxstep + iterm] < err[iterm] ) {
            red = alf[(network_solver_data.maxiter-2)*network_solver_data.maxstep + iterm] * 0.7 / err[iterm];
            reduce = 1;
          }
        } else if ( alf[iteropt*network_solver_data.maxstep + iterm] < err[iterm] ) {
          red = alf[(iteropt-1)*network_solver_data.maxstep + iterm] / err[iterm];
          reduce = 1;
        }
        
         /* reduce step size */
        if (reduce == 1) {
          /* step not successful */
          red = red > 1e-5 ? red : 1e-5;
          red = red < 0.7 ? red : 0.7;
          dttry = dttry * red;
          iter = -1;          
        }
      }
    }

    lastiter = iter;
    /* step successful */
    for (i=0; i<network_solver_data.matrixsize; i++) {
      y[i] = ynew[i];
    }

    first = 0;

    scale = 1.0;
    workmin = 1.0e35;
    for (iter = 0; iter < iterm+1; iter++) {
      fac = err[iter] > 0.1 ? err[iter] : 0.1;
      work = fac * a[iter+1];
      if (work < workmin) {
        scale = fac;
        workmin = work;
        iteropt = iter + 1;
      }
    }

    dtnext = dttry / scale;
    if (iteropt >= lastiter && iteropt != network_solver_data.maxiter-1 && !reduce) {
      fac = scale / alf[iteropt*network_solver_data.maxstep + iteropt-1];
      fac = fac > 0.1 ? fac : 0.1;
      if (a[iteropt+1]*fac <= workmin) {
        dtnext = dttry / fac;
        iteropt += 1;
      }
    } 

    t += dttry;

    if (t + dtnext > dt) {
      first = 1;
      iteropt = network_solver_data.maxiter-1;
      dttry = dt - t;
    } else {
      dttry = dtnext;
    }
	step++;
  }
  
#ifdef NETWORK_OUTPUT_TIMEEVOLUTION
  printf( "t(%03d): %11.5e, dt (%02d): %11.5e, x:", step, t, m, dttry ); for (i=0; i<network_solver_data.nelements; i++) { printf( " %8.1e", y[i]*network_solver_data.aion[i] ); } if (network_solver_data.nelements<network_solver_data.matrixsize) { printf( " %8.1e", y[network_solver_data.nelements] ); } printf( "\n" );
  fwrite( &t, sizeof(double), 1, fp );
  fwrite( y, sizeof(double), network_solver_data.matrixsize, fp );
  fclose( fp );
#endif
}

void network_solver_calc_dy( double temp, double rho, double dt, int steps, double *y, double *yn ) {
  double *dy, *rhs, *jacob;
  int matrixsize, matrixsize2;
  int i, j, s;
  double h;
  gsl_matrix_view A;
  gsl_vector_view b;
  gsl_vector *x;
  gsl_permutation *p;

  /* take variables from global struct */
  dy = network_solver_data.dy;
  rhs = network_solver_data.rhs;
  jacob = network_solver_data.jacob;
  matrixsize = network_solver_data.matrixsize;
  matrixsize2 = network_solver_data.matrixsize2;
  
  h = dt / steps; /* subtimestep */
  (*network_solver_data.getrhs)( temp, rho, y, rhs );
  (*network_solver_data.getjacob)( temp, rho, y, rhs, jacob );

  for ( i=0; i<matrixsize2; i++ ) {
    jacob[i] = - h * jacob[i];
    if (i % (matrixsize+1) == 0) jacob[i] = jacob[i] + 1.0;
  }

  A = gsl_matrix_view_array( jacob, matrixsize, matrixsize );
  x = gsl_vector_alloc( matrixsize );
  p = gsl_permutation_alloc( matrixsize );
  gsl_linalg_LU_decomp( &A.matrix, p, &s );

  for ( i=0; i<matrixsize; i++ ) { rhs[i] = h*rhs[i]; }
  b = gsl_vector_view_array( rhs, matrixsize );
  gsl_linalg_LU_solve( &A.matrix, p, &b.vector, x );
  for ( i=0; i<matrixsize; i++ ) {
    dy[i] = gsl_vector_get(x,i);
    yn[i] = y[i] + dy[i];
  }

  for (j=1; j<steps; j++) {
    (*network_solver_data.getrhs)( temp, rho, yn, rhs );
    for ( i=0; i<matrixsize; i++ ) { rhs[i] = h*rhs[i] - dy[i]; }
    b = gsl_vector_view_array( rhs, matrixsize );
    gsl_linalg_LU_solve( &A.matrix, p, &b.vector, x );
    for ( i=0; i<matrixsize; i++ ) {
      dy[i] = dy[i] + 2.0 * gsl_vector_get(x,i);
      yn[i] = yn[i] + dy[i];
    }
  }

  (*network_solver_data.getrhs)( temp, rho, yn, rhs );
  for ( i=0; i<matrixsize; i++ ) { rhs[i] = h*rhs[i] - dy[i]; }
  b = gsl_vector_view_array( rhs, matrixsize );
  gsl_linalg_LU_solve( &A.matrix, p, &b.vector, x );
  for ( i=0; i<matrixsize; i++ ) {
    yn[i] = yn[i] + gsl_vector_get(x,i);
  }

  gsl_vector_free( x );
  gsl_permutation_free( p );
}

void network_solver_extrapolate( int iter, double *x, double *y, double *dy, double *qcol ) {
  int matrixsize;
  double *d;
  double delta, f1, f2, q;
  int i, k;
  
  matrixsize = network_solver_data.matrixsize;
  d = network_solver_data.d;
  
  for (i=0; i<matrixsize; i++) {
    dy[i] = y[i];
  }

  if (iter == 0) {
    for (i=0; i<matrixsize; i++) qcol[i] = y[i];
  } else {
    for (i=0; i<matrixsize; i++) d[i] = y[i];

    for (k=0; k<iter; k++) {
      delta = 1.0 / ( x[iter-k-1] - x[iter] );
      f1 = x[iter] * delta;
      f2 = x[iter-k-1] * delta;

      for (i=0; i<matrixsize; i++) {
        q = qcol[ k*matrixsize + i ];
        qcol[ k*matrixsize + i ] = dy[i];
        delta = d[i] - q;
        dy[i] = f1 * delta;
        d[i] = f2 * delta;
        y[i] = y[i] + dy[i];
      }
    }

    for (i=0; i<matrixsize; i++) qcol[ iter*matrixsize + i ] = dy[i];
  }
}

#endif /* NUCLEAR_NETWORK */
