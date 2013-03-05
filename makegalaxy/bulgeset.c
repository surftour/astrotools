#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"



static double R, z, r;


double set_bulge_velocities(void)
{
  int i;
  long dum;
  int iz, ir;
  double ur, uz;
  double vr, vphi;
  double vx, vy, vz;

#ifdef MAXWELLIAN
  double vdisp_rz, vdisp_phi, vstream_phi;
#else
  int iter;
  double frand=1.0, f0=0.0, fmax;
  double iphi, E, L;
  double vmax, v2, vtheta;
#endif


  if(N_BULGE == 0)
    return 0;

  //srand48(7285);
  dum = drand48() * 1e8;

  printf("set bulge velocities...\t");
  fflush(stdout);

  for(i = 1; i <= N_BULGE; i++)
    {
      R = sqrt(xp_bulge[i] * xp_bulge[i] + yp_bulge[i] * yp_bulge[i]);
      z = fabs(zp_bulge[i]);
      r = sqrt(xp_bulge[i] * xp_bulge[i] + yp_bulge[i] * yp_bulge[i] + zp_bulge[i] * zp_bulge[i]);

      //if(R < Baselen)
      //ir = 0;
      //else
	ir = find_idx(R, list_R, RSIZE);
	/* ir = (log(R) - log(Baselen)) / (log(LL) - log(Baselen)) * (RSIZE - 1) + 1; */

      ur = (R - list_R[ir]) / (list_R[ir + 1] - list_R[ir]);

      //if(z < Baselen)
      //iz = 0;
      //else
	iz = find_idx(z, list_z, ZSIZE);
	/* iz = (log(z) - log(Baselen)) / (log(LL) - log(Baselen)) * (ZSIZE - 1) + 1; */

      uz = (z - list_z[iz]) / (list_z[iz + 1] - list_z[iz]);

      if(ir < 0 || ir >= RSIZE || ur < 0 || ur > 1)
	{
	  printf("\nR= %g ( in bounds? list_R[0]= %g  list_R[RSIZE]= %g)\nin: set_bulge_velocities\nproblem: ir=%d out of range ur=%g\n", R, list_R[0], list_R[RSIZE], ir, ur);
	  exit(0);
	}

      if(iz < 0 || iz >= ZSIZE || uz < 0 || uz > 1)
	{
	  printf("iz=%d out of range uz=%g\n", iz, uz);
	  exit(0);
	}

#ifdef MAXWELLIAN


      vdisp_rz = VelDispRz_bulge[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispRz_bulge[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispRz_bulge[ir][iz + 1] * (1 - ur) * (uz) + VelDispRz_bulge[ir + 1][iz + 1] * (ur) * (uz);

      vdisp_phi = VelDispPhi_bulge[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispPhi_bulge[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispPhi_bulge[ir][iz + 1] * (1 - ur) * (uz) + VelDispPhi_bulge[ir + 1][iz + 1] * (ur) * (uz);

      vstream_phi = VelStreamPhi_bulge[ir][iz] * (1 - ur) * (1 - uz)
	+ VelStreamPhi_bulge[ir + 1][iz] * (ur) * (1 - uz)
	+ VelStreamPhi_bulge[ir][iz + 1] * (1 - ur) * (uz) + VelStreamPhi_bulge[ir + 1][iz + 1] * (ur) * (uz);


      if(vdisp_rz < 0)
	{
	  printf("in bulge: vdisp_rz:%g   %g %g %d %d \n", vdisp_rz, ur, uz, ir, iz);
	  vdisp_rz = -vdisp_rz;
	}
      if(vdisp_phi < 0)
	{
	  printf("in bulge: vdisp_phi:%g  %g %g %d %d\n", vdisp_phi, ur, uz, ir, iz);

	  vdisp_phi = -vdisp_phi;
	}

      vr = gasdev(&dum) * sqrt(vdisp_rz);
      vz = gasdev(&dum) * sqrt(vdisp_rz);

      vphi = vstream_phi + gasdev(&dum) * sqrt(vdisp_phi);


      vx = vr * xp_bulge[i] / R - vphi * yp_bulge[i] / R;
      vy = vr * yp_bulge[i] / R + vphi * xp_bulge[i] / R;

#else

      iter= 0;

      /* what's the potential at this point? */
      iphi= comp_phi(R,z);

      vmax= sqrt(2. * fabs(iphi));
      fmax= comp_DF_bulge(iphi, 0.0);

      do {

        iter++;

        do {
          /* see haloset.c for more details about what's going on here */
          vr = vmax * (drand48()*2.-1.) * sqrt(3.) / sqrt(3.-2.*compute_ani_beta(r));
          vtheta = vmax * (drand48()*2.-1.) * sqrt(3.) * sqrt(1.-compute_ani_beta(r)) / sqrt(3.-2.*compute_ani_beta(r));
          vphi = vmax * (drand48()*2.-1.) * sqrt(3.) * sqrt(1.-compute_ani_beta(r)) / sqrt(3.-2.*compute_ani_beta(r));

          vx = vr * (R/r) * (xp_bulge[i]/R) + vtheta * (zp_bulge[i]/r) * (xp_bulge[i]/R) - vphi * (yp_bulge[i]/R);
          vy = vr * (R/r) * (yp_bulge[i]/R) + vtheta * (zp_bulge[i]/r) * (yp_bulge[i]/R) + vphi * (xp_bulge[i]/R);
          vz = vr * (zp_bulge[i]/r) - vtheta * (R/r);

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
        f0= comp_DF_bulge(E,L);
        frand= fmax * drand48();

      } while (frand > f0);

#endif

      vxp_bulge[i] = vx;
      vyp_bulge[i] = vy;
      vzp_bulge[i] = vz;
    }

  printf("done.\n");
  fflush(stdout);

  return 0;
}





double set_bulge_positions(void)
{
  int i, countr, countz;
  double q, r, phi, theta;

  if(N_BULGE == 0)
    return 0;

  //srand48(666);

  printf("set bulge positions...\t");

  for(i = 1, countr = countz = 0; i <= N_BULGE;)
    {
      q = drand48();
      r = bulge_q_to_r(q);

      phi = drand48() * PI * 2;
      theta = acos(drand48() * 2 - 1);

      xp_bulge[i] = r * sin(theta) * cos(phi);
      yp_bulge[i] = r * sin(theta) * sin(phi);
      zp_bulge[i] = r * cos(theta);

      /* now stretch */

      xp_bulge[i] *= A;
      yp_bulge[i] *= A;
      zp_bulge[i] *= A;


      mp_bulge[i] = M_BULGE / N_BULGE;

      if((r * A) > LL)
	countr++;
      else
	i++;
    }

  /*
     printf("bulge discarded:  %d  \n",countr);
   */

  printf("done.\n");

  return 0;
}
