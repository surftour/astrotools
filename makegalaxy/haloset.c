#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"



static double R, z, r;



double set_halo_velocities(void)
{
  int i, rangecount=0;
  long dum;
  int iz, ir;
  double ur, uz;
  double vx, vy, vz;
  double vr, vphi;

#ifdef MAXWELLIAN
  double vdisp_rz, vdisp_phi, vstream_phi;
#else
  int iter;
  double frand=1.0, f0=0.0, fmax;
  double iphi, E, L;
  double vmax, v2, vtheta;
#endif


  if(N_HALO == 0)
    return 0;

  //srand48(85);
  dum = (long) -1.0 * drand48() * 1e8;

  printf("set halo velocities...\t");
  fflush(stdout);

  for(i = 1; i <= N_HALO; i++)
    {
      R = sqrt(xp_halo[i] * xp_halo[i] + yp_halo[i] * yp_halo[i]);
      z = fabs(zp_halo[i]);
      r = sqrt(xp_halo[i] * xp_halo[i] + yp_halo[i] * yp_halo[i] + zp_halo[i] * zp_halo[i]);

      //if(R < Baselen)
      //ir = 0;
      //else
	ir = find_idx(R, list_R, RSIZE);

      ur = (R - list_R[ir]) / (list_R[ir + 1] - list_R[ir]);

      //if(z < Baselen)
      //iz = 0;
      //else
	iz = find_idx(z, list_z, ZSIZE);

      uz = (z - list_z[iz]) / (list_z[iz + 1] - list_z[iz]);

      if(ir < 0 || ir >= RSIZE || ur < 0 || ur > 1)
	{
	  rangecount++;
	  /* 
	     printf("ir=%d out of range R=%g\n", ir, R);
	  */
	  vxp_halo[i] = 0;
	  vyp_halo[i] = 0;
	  vzp_halo[i] = 0;
	  continue;
	}

      if(iz < 0 || iz >= ZSIZE || uz < 0 || uz > 1)
	{
	  rangecount++;
	  /*
	  printf("iz=%d out of range Z=%g\n", iz, z);
	  */
	  vxp_halo[i] = 0;
	  vyp_halo[i] = 0;
	  vzp_halo[i] = 0;
	  continue;
	}

#ifdef MAXWELLIAN

      vdisp_rz = VelDispRz_halo[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispRz_halo[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispRz_halo[ir][iz + 1] * (1 - ur) * (uz) + VelDispRz_halo[ir + 1][iz + 1] * (ur) * (uz);

      vdisp_phi = VelDispPhi_halo[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispPhi_halo[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispPhi_halo[ir][iz + 1] * (1 - ur) * (uz) + VelDispPhi_halo[ir + 1][iz + 1] * (ur) * (uz);
/*
if((i==104317) || (i==104320))
{
	printf("ir= %d  iz= %d  ur= %g  uz= %g\n",ir,iz,ur,uz);
	printf("VelDispPhi_halo[ir][iz]= %g\n",VelDispPhi_halo[ir][iz]);
	printf("VelDispPhi_halo[ir+1][iz]= %g\n",VelDispPhi_halo[ir+1][iz]);
	printf("VelDispPhi_halo[ir][iz+1]= %g\n",VelDispPhi_halo[ir][iz+1]);
	printf("VelDispPhi_halo[ir+1][iz+1]= %g\n",VelDispPhi_halo[ir+1][iz+1]);
	printf("---\n");
}
*/

      vstream_phi = VelStreamPhi_halo[ir][iz] * (1 - ur) * (1 - uz)
	+ VelStreamPhi_halo[ir + 1][iz] * (ur) * (1 - uz)
	+ VelStreamPhi_halo[ir][iz + 1] * (1 - ur) * (uz) + VelStreamPhi_halo[ir + 1][iz + 1] * (ur) * (uz);
/*
if((i==104317) || (i==104320))
{
        printf("VelStreamPhi_halo[ir][iz]= %g\n",VelStreamPhi_halo[ir][iz]);
        printf("VelStreamPhi_halo[ir+1][iz]= %g\n",VelStreamPhi_halo[ir+1][iz]);
        printf("VelStreamPhi_halo[ir][iz+1]= %g\n",VelStreamPhi_halo[ir][iz+1]);
        printf("VelStreamPhi_halo[ir+1][iz+1]= %g\n",VelStreamPhi_halo[ir+1][iz+1]);
        fflush(stdout);
}
*/

      if(vdisp_rz < 0)
	{
	  printf("in halo: vdisp_rz:%g   %g %g %d %d \n", vdisp_rz, ur, uz, ir, iz);
	  vdisp_rz = -vdisp_rz;
	}
      if(vdisp_phi < 0)
	{
	  printf("in halo: vdisp_phi:%g  %g %g %d %d\n", vdisp_phi, ur, uz, ir, iz);

	  vdisp_phi = -vdisp_phi;
	}

      vr = gasdev(&dum) * sqrt(vdisp_rz);
      vz = gasdev(&dum) * sqrt(vdisp_rz);

      vphi = vstream_phi + gasdev(&dum) * sqrt(vdisp_phi);
/*
if(i==104317)
{
	printf("vphi= %g    vstream_phi= %g  vdisp_phi= %g  gasdev(&dum)= %g\n", vphi, vstream_phi, vdisp_phi, gasdev(&dum));
	fflush(stdout);
}
*/


      vx = vr * xp_halo[i] / R - vphi * yp_halo[i] / R;
      vy = vr * yp_halo[i] / R + vphi * xp_halo[i] / R;

#else
      iter= 0;

      /* what's the potential at this point? */
      iphi= comp_phi(R,z);

      vmax= sqrt(2. * fabs(iphi));
      fmax= comp_DF_halo(iphi, 0.0);

      do {

	iter++;

	/* select velocities at random (but make sure it's bound) */
	do {
	  /* easy to show that if we define beta = 1 - sig_theta^2/sig_r^2 and
	     v^2 = sig_r^2 + sig_theta^2 + sig_phi^2 then

		sig_r^2 = v^2 / (3-2*beta)
		sig_theta^2 = sig_phi^2 = v^2 * (1-beta) / (3-2*beta)

	     and then we get the 3 (or sqrt(3), since we're selecting velocities)
	     factor requiring this to be normalized the same as a purely isotropic
	     distribution.
	  */
	  vr = vmax * (drand48()*2.-1.) * sqrt(3.) / sqrt(3.-2.*compute_ani_beta(r));
	  vtheta = vmax * (drand48()*2.-1.) * sqrt(3.) * sqrt(1.-compute_ani_beta(r)) / sqrt(3.-2.*compute_ani_beta(r));
	  vphi = vmax * (drand48()*2.-1.) * sqrt(3.) * sqrt(1.-compute_ani_beta(r)) / sqrt(3.-2.*compute_ani_beta(r));

	  vx = vr * (R/r) * (xp_halo[i]/R) + vtheta * (zp_halo[i]/r) * (xp_halo[i]/R) - vphi * (yp_halo[i]/R);
	  vy = vr * (R/r) * (yp_halo[i]/R) + vtheta * (zp_halo[i]/r) * (yp_halo[i]/R) + vphi * (xp_halo[i]/R);
	  vz = vr * (zp_halo[i]/r) - vtheta * (R/r);

	  v2= vx*vx + vy*vy + vz*vz;
	} while (v2 > (vmax*vmax));


	E= 0.5*v2 + iphi;

	L= r * sqrt(vtheta*vtheta + vphi*vphi);

	if(iter>5e5)
	  {
	    printf("PROBLEM: run out of iterations!\n");
	    printf("i=%d   ****   R,z= %g|%g   ", i, R, z);
	    printf("v= %g vmax=%g   iphi= %g   E=%g   ", sqrt(v2), vmax, iphi, E);
	    printf("f0= %g  fmax= %g  frand= %g   iter= %d\n", f0, fmax, frand, iter);
	    fflush(stdout);
	    break;
	  }

	/* this should never occur (provided we check v2 above) */
	if(E>=0)
	  {
	    frand= 1.0;
	    f0= 0.0;
	    continue;
	  }

	/* Now, we're ready to do the acceptance or rejection */
	f0= comp_DF_halo(E,L);
	frand= fmax * drand48();
      
      } while (frand > f0);


/*
if(r>3.0 && r<3.5)
  printf("inner i=%d   ****   R|z= %g|%g  v= %g (%g|%g|%g) (%g|%g|%g) vmax=%g   iphi= %g   E=%g  beta=%g  iter= %d\n", i, R, z, sqrt(v2), vx, vy, vz, vr, vtheta, vphi, vmax, iphi, E, compute_ani_beta(r),iter); fflush(stdout);
if(r>1000.0)
  printf("outer i=%d   ****   R|z= %g|%g  v= %g (%g|%g|%g) (%g|%g|%g) vmax=%g   iphi= %g   E=%g  beta=%g  iter= %d\n", i, R, z, sqrt(v2), vx, vy, vz, vr, vtheta, vphi, vmax, iphi, E, compute_ani_beta(r),iter); fflush(stdout);
*/

#endif

      vxp_halo[i] = vx;
      vyp_halo[i] = vy;
      vzp_halo[i] = vz;
/*
if((i > 104312) && (i < 104322))
{
	printf("i= %d   R= %g vr= %g  vz= %g  vphi= %g  xyz= %g|%g|%g  v_xyz= %g|%g|%g\n",i,R,vr,vz,vphi,xp_halo[i],yp_halo[i],zp_halo[i],vx,vy,vz); fflush(stdout);
}
*/
    }

  /* printf("rangecount=%d\n", rangecount); */
  printf("done.\n");
  fflush(stdout);

  return 0;
}






double set_halo_positions(void)
{
  int i, countr, countz;
  double q, R, phi, theta;
  double halo_q_to_r(double q);

  if(N_HALO == 0)
    return 0;

  //srand48(22);

  printf("set halo positions...\t");

  for(i = 1, countr = countz = 0; i <= N_HALO;)
    {

      do
	{
	  q = drand48();

	  R = halo_q_to_r(q);
	}
      while((R > 5 * R200) || (R < list_R[2]));

      phi = drand48() * PI * 2;
      theta = acos(drand48() * 2 - 1);

      xp_halo[i] = R * sin(theta) * cos(phi);
      yp_halo[i] = R * sin(theta) * sin(phi);
      zp_halo[i] = R * cos(theta);

      i++;
    }

  for(i = 1; i <= N_HALO; i++)
    mp_halo[i] = M_HALO / N_HALO;

  printf("done.\n");

  return 0;
}


