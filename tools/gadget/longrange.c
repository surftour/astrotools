#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


/*! \file longrange.c
 *  \brief driver routines for computation of long-range gravitational PM force
 */

#ifdef PMGRID

/*! Driver routine to call initializiation of periodic or/and non-periodic FFT
 *  routines.
 */
void long_range_init(void)
{
#ifdef PERIODIC
  pm_init_periodic();
#ifdef PLACEHIGHRESREGION
  pm_init_nonperiodic();
#endif
#else
  pm_init_nonperiodic();
#endif
}


void long_range_init_regionsize(void)
{
#ifdef PERIODIC
#ifdef PLACEHIGHRESREGION
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif
#else
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif
}


/*! This function computes the long-range PM force for all particles.
 */
void long_range_force(void)
{
  int i;

#ifndef PERIODIC
  int j;
  double fac;
#endif


  for(i = 0; i < NumPart; i++)
    {
      P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;
#ifdef EVALPOTENTIAL
      P[i].PM_Potential = 0;
#endif
#ifdef DISTORTIONTENSORPS
      P[i].tidal_tensorpsPM[0][0] = 0;
      P[i].tidal_tensorpsPM[0][1] = 0;
      P[i].tidal_tensorpsPM[0][2] = 0;
      P[i].tidal_tensorpsPM[1][0] = 0;
      P[i].tidal_tensorpsPM[1][1] = 0;
      P[i].tidal_tensorpsPM[1][2] = 0;
      P[i].tidal_tensorpsPM[2][0] = 0;
      P[i].tidal_tensorpsPM[2][1] = 0;
      P[i].tidal_tensorpsPM[2][2] = 0;
#endif
    }

#ifdef NOGRAVITY
  return;
#endif


#ifdef PERIODIC
  pmforce_periodic();
#ifdef DISTORTIONTENSORPS
  pmtidaltensor_periodic(0);
  pmtidaltensor_periodic(1);
  pmtidaltensor_periodic(2);
  pmtidaltensor_periodic(3);
  pmtidaltensor_periodic(4);
  pmtidaltensor_periodic(5);
  /* check tidal tensor of given particle ID */
  /* check_tidaltensor(500); */
#endif
#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);
#ifdef DISTORTIONTENSORPS
  pmtidaltensor_nonperiodic(1, 0);
  pmtidaltensor_nonperiodic(1, 1);
  pmtidaltensor_nonperiodic(1, 2);
  pmtidaltensor_nonperiodic(1, 3);
  pmtidaltensor_nonperiodic(1, 4);
  pmtidaltensor_nonperiodic(1, 5);
#endif
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(1);	/* try again */
#ifdef DISTORTIONTENSORPS
      pmtidaltensor_nonperiodic(1, 0);
      pmtidaltensor_nonperiodic(1, 1);
      pmtidaltensor_nonperiodic(1, 2);
      pmtidaltensor_nonperiodic(1, 3);
      pmtidaltensor_nonperiodic(1, 4);
      pmtidaltensor_nonperiodic(1, 5);
#endif

    }
  if(i == 1)
    endrun(68686);
#endif
#else
  i = pmforce_nonperiodic(0);
#ifdef DISTORTIONTENSORPS
  pmtidaltensor_nonperiodic(0, 0);
  pmtidaltensor_nonperiodic(0, 1);
  pmtidaltensor_nonperiodic(0, 2);
  pmtidaltensor_nonperiodic(0, 3);
  pmtidaltensor_nonperiodic(0, 4);
  pmtidaltensor_nonperiodic(0, 5);
#endif

  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(0);	/* try again */
#ifdef DISTORTIONTENSORPS
      pmtidaltensor_nonperiodic(0, 0);
      pmtidaltensor_nonperiodic(0, 1);
      pmtidaltensor_nonperiodic(0, 2);
      pmtidaltensor_nonperiodic(0, 3);
      pmtidaltensor_nonperiodic(0, 4);
      pmtidaltensor_nonperiodic(0, 5);
#endif
    }
  if(i == 1)
    endrun(68687);
#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);
#ifdef DISTORTIONTENSORPS
  pmtidaltensor_nonperiodic(1, 0);
  pmtidaltensor_nonperiodic(1, 1);
  pmtidaltensor_nonperiodic(1, 2);
  pmtidaltensor_nonperiodic(1, 3);
  pmtidaltensor_nonperiodic(1, 4);
  pmtidaltensor_nonperiodic(1, 5);
#endif
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();

      /* try again */

      for(i = 0; i < NumPart; i++)
	P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;

#ifdef DISTORTIONTENSORPS
      P[i].tidal_tensorpsPM[0][0] = 0;
      P[i].tidal_tensorpsPM[0][1] = 0;
      P[i].tidal_tensorpsPM[0][2] = 0;
      P[i].tidal_tensorpsPM[1][0] = 0;
      P[i].tidal_tensorpsPM[1][1] = 0;
      P[i].tidal_tensorpsPM[1][2] = 0;
      P[i].tidal_tensorpsPM[2][0] = 0;
      P[i].tidal_tensorpsPM[2][1] = 0;
      P[i].tidal_tensorpsPM[2][2] = 0;
#endif
      i = pmforce_nonperiodic(0) + pmforce_nonperiodic(1);

#ifdef DISTORTIONTENSORPS
      pmtidaltensor_nonperiodic(0, 0);
      pmtidaltensor_nonperiodic(0, 1);
      pmtidaltensor_nonperiodic(0, 2);
      pmtidaltensor_nonperiodic(0, 3);
      pmtidaltensor_nonperiodic(0, 4);
      pmtidaltensor_nonperiodic(0, 5);

      pmtidaltensor_nonperiodic(1, 0);
      pmtidaltensor_nonperiodic(1, 1);
      pmtidaltensor_nonperiodic(1, 2);
      pmtidaltensor_nonperiodic(1, 3);
      pmtidaltensor_nonperiodic(1, 4);
      pmtidaltensor_nonperiodic(1, 5);
#endif
    }
  if(i != 0)
    endrun(68688);
#endif
#endif


#ifndef PERIODIC
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0;

      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].GravPM[j] += fac * P[i].Pos[j];
    }


  /* Finally, the following factor allows a computation of cosmological simulation 
     with vacuum energy in physical coordinates */

  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].GravPM[j] += fac * P[i].Pos[j];
    }
#endif

}


#endif
