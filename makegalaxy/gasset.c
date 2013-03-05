#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"
#include "gsl/gsl_qrng.h"

#include "prototypes.h"
#include "globvars.h"


static double R, z;



void set_gas_velocities(void)
{
  int i;
  //long dum;
  int iz, ir;
  double ur, uz;
  double vstream_phi;
  double vr, vphi;
  double vx, vy, vz;



  if(N_GAS == 0)
    return;

  //dum = drand48() * 1e8;

  printf("set gas velocities...\t");
  fflush(stdout);

  for(i = 1; i <= N_GAS; i++)
    {
      R = sqrt(xp_gas[i] * xp_gas[i] + yp_gas[i] * yp_gas[i]);
      z = fabs(zp_gas[i]);

      if(R < Baselen)
	ir = 0;
      else
	ir = find_idx(R, list_R, RSIZE);
      /* ir = (log(R) - log(Baselen)) / (log(LL) - log(Baselen)) * (RSIZE - 1) + 1; */

      ur = (R - list_R[ir]) / (list_R[ir + 1] - list_R[ir]);

      if(z < Baselen)
	iz = 0;
      else
	iz = find_idx(z, list_z, ZSIZE);
      /* iz = (log(z) - log(Baselen)) / (log(LL) - log(Baselen)) * (ZSIZE - 1) + 1; */

      uz = (z - list_z[iz]) / (list_z[iz + 1] - list_z[iz]);

      if(ir < 0 || ir >= RSIZE)
	{
	  printf("ir=%d out of range\n", ir);
	}

      if(iz < 0 || iz >= ZSIZE)
	{
	  printf("iz=%d out of range\n", iz);
	}


      vstream_phi = VelStreamGas[ir][iz] * (1 - ur) * (1 - uz)
	+ VelStreamGas[ir + 1][iz] * (ur) * (1 - uz)
	+ VelStreamGas[ir][iz + 1] * (1 - ur) * (uz) + VelStreamGas[ir + 1][iz + 1] * (ur) * (uz);



      vr = 0;
      vz = 0;

      vphi = vstream_phi;

      vx = vr * xp_gas[i] / R - vphi * yp_gas[i] / R;
      vy = vr * yp_gas[i] / R + vphi * xp_gas[i] / R;

      vxp_gas[i] = vx;
      vyp_gas[i] = vy;
      vzp_gas[i] = vz;

      u_gas[i] = U4;

    }

  printf("done.\n");
  fflush(stdout);
}



void set_gas_positions(void)
{
  int i;
  double q[3];
  double pos[3];
  int n;

  if(N_GAS == 0)
    return;

  printf("set gas positions...\t");
  fflush(stdout);

  // we want a 3d qrng here because we draw three numbers per particle
  gsl_qrng* qrng;
  if(UseQrng) {
    printf("Using quasi-random numbers for gas positions\n");
    qrng = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 3);
  }

  n = N_GAS - (WriteDensity ? (int)(0.005*N_GAS)*2 : 0);  

  for(i = 1; i <= n; i++)
    {
      if(UseQrng)
	gsl_qrng_get(qrng, q);
      else {
        q[0] = drand48();
        q[1] = drand48();
        q[2] = drand48();
      };
      q[0]*=.999;
      q[1]= q[1]*2-1;

      rho_gas[i] = sample_gas_position(q, pos);
  
      xp_gas[i]=pos[0];
      yp_gas[i]=pos[1];
      zp_gas[i]=pos[2];

      assert(!isnan(xp_gas[i]));
      assert(!isnan(yp_gas[i]));
      assert(!isnan(zp_gas[i]));

      mp_gas[i] = M_GAS / n;
    }


  // if we are using the rho option, we need a "shell" of particles
  // that prevents extrapolation of disk densities into the
  // background. We do this by allocating 1% of the particles to a
  // series with very large z.
  if(WriteDensity) {
    if(UseQrng) {
      qrng = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 2);
    }

    for(i=n; i <= N_GAS; i+=2)
      {
	if(UseQrng) {
	  // only fills first 2 components of q now
	  gsl_qrng_get(qrng, q);
	  q[2]=q[1];
	}
	else {
	  q[0] = drand48();
	  q[2] = drand48();
	};
	q[0]*=.999;
	q[1]= .9999;

	rho_gas[i] = sample_gas_position(q, pos);

	xp_gas[i]=pos[0];
	yp_gas[i]=pos[1];
	zp_gas[i]=pos[2];

	mp_gas[i] = 0;

	// duplicate particle to opposite side of disk
	rho_gas[i+1]=rho_gas[i];
	xp_gas[i+1]=pos[0];
	yp_gas[i+1]=pos[1];
	zp_gas[i+1]=-pos[2];
	mp_gas[i+1] = mp_gas[i];
      
      }
  }

  if(UseQrng)
    gsl_qrng_free(qrng);

  printf("done.\n");
}


/** This function samples the position of the gas particles based on
    the 3 random numbers in q. Position returned in pos. The volume
    density at the position is returned. */
double sample_gas_position(double* q, double *pos)
{
  double R, phi, ur, uz; 
  int rbin, zbin;
  double CumMax, CumCur, CumNext, target;

  R = gas_q_to_R(q[0]);
  
  rbin = find_idx(R, list_R, RSIZE);
  assert(rbin<RSIZE-1);
  
  // ur is fractional bin size between rbin and our point
  ur = (R - list_R[rbin]) / (list_R[rbin + 1] - list_R[rbin]);
  

  // if this is zero we have no density here and consequently can't
  // have a particle there either. This would indicate that gas_q_to_R
  // has done something crazy.
  // printf("q= %g  R= %g    rbin= %d   CumulMassGas[rbin][ZSIZE-1]= %g\n", q[0], R, rbin, CumulMassGas[rbin][ZSIZE-1]);
  assert(CumulMassGas[rbin][ZSIZE-1]>0);

  // find a first guess to what bin we want
  zbin = find_idx(fabs(q[1])*CumulMassGas[rbin][ZSIZE-1],
		  CumulMassGas[rbin], ZSIZE);
  
  assert(zbin>=0);
  assert(zbin<ZSIZE);

  // must use ZSIZE-2 to ensure we don't overrun
  CumMax = interpolate_table_2d(ur, 1, rbin, ZSIZE-2, CumulMassGas); 
  target = fabs(q[1])*CumMax;

  while(1) {
    CumCur = interpolate_table_2d(ur, 0, rbin, zbin, CumulMassGas);
    CumNext = interpolate_table_2d(ur, 0, rbin, zbin+1, CumulMassGas);
    // We need to include the equalities in both comparisons so that we find a bin for q==1, where it will be exactly 
    if( target >= CumCur ) {
      // target is potentially in this bin
      if (target < CumNext)
	// simple, found bin
	break;

      // special test for if the target is exactly on a degenerate bin
      // edge (notably happens when q==+-1)
      if ( (target == CumNext) &&
	   // this is *potentially* a degenerate find, but only if zbin-1 has a different value
	   CumCur > interpolate_table_2d(ur, 0, rbin, zbin-1, 
					 CumulMassGas) ) {
	// we now know that the bin below is the correct bin (with uz=1)
	--zbin;
	CumCur = interpolate_table_2d(ur, 0, rbin, zbin, CumulMassGas);
	CumNext = interpolate_table_2d(ur, 0, rbin, zbin+1, CumulMassGas);
	break;
      }
    }

    if (target <= CumCur)
      // the '=' applies only if we've tested above and concluded
      // it's not in this bin
      --zbin;
    else
      ++zbin;
    assert(zbin>=0);
    assert(zbin<ZSIZE);
  }
  assert(zbin>=0);
  assert(zbin<ZSIZE);

  uz = (CumMax*fabs(q[1])-CumCur)/(CumNext-CumCur);

  z = interpolate_table(uz, zbin, list_z);

  if(q[1]<0)
    z = -z;


  phi = q[2] * PI * 2;

  pos[0] = R * cos(phi);
  pos[1] = R * sin(phi);
  pos[2] = z;

  double rho = interpolate_table_2d(ur, uz, rbin, zbin, RhoGas);
  return rho;
}

double interpolate_table(double f, int bin, double* array)
{
  assert(bin>=0);
  assert(f>=0);
  assert(f<=1);

  return (1-f)*array[bin]+f*array[bin+1];
}

double interpolate_table_2d(double f1, double f2, int bin1, int bin2,
			    double** array)
{
  assert(bin1>=0);
  assert(bin2>=0);
  assert(f1>=0);
  assert(f1<=1);
  assert(f2>=0);
  assert(f2<=1);
  
  return
    (1-f1)* ((1-f2)*array[bin1][bin2]   + f2*array[bin1][bin2+1]) +
    f1*     ((1-f2)*array[bin1+1][bin2] + f2*array[bin1+1][bin2+1]);
}

/* Find a value in a sorted array of size SIZE using bisection
   search. If value is exactly on a bin value, that bin is returned,
   otherwise the return value points to the lower side bin. */
int find_bin(double value, double* array, int size)
{
  int min=0, max=size-1, bin=size/2;
  assert(value>=array[0]);
  assert(value<array[size-1]);
  
  while(value >= array[bin+1] || value < array[bin]) {
    if(value < array[bin])
      max=bin;
    else
      min=bin;
    bin=(min+max)/2;
  }

  assert(value>=array[bin]);
  assert(value<array[bin+1]);
  return bin;
}
