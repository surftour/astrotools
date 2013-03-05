#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "allvars.h"
#include "proto.h"


/* This routine computes the accelerations for all active particles. 
 * First, the gravitational forces are computed (this also reconstructs
 * the tree, if needed). Also note that the gas-particle tree will
 * in any case be updated in its geometrical properties.

 * Then the density for active SPH particles is computed,
 * and finally, the hydrodynamical forces are added.
 */
void compute_accelerations(int mode) 
{
  double tstart, tend;
  
  printf("Start force computation...\n");


  tstart=second();   /* measure the time for the full force computation */

#ifdef GRAPE
  gravity_grape();     /* computes gravity accel. & potential */
#else
  gravity_tree();     /* computes gravity accel. & potential  */
#endif

  tend=second();
  All.CPU_Gravity+= timediff(tstart,tend);

#ifdef VELDISP
  tstart=second();
  veldisp();
  tend=second();
  All.CPU_Hydro+= timediff(tstart, tend);
#endif

  if(N_gas>0) 
    {
      tstart=second();

      density();       /* computes density, and pressure */
      
      hydro_force();   /* adds hydrodynamical accelerations 
			  and computes du/dt  */
      tend=second();
      All.CPU_Hydro+= timediff(tstart,tend);
    }

  printf("done force.\n"); fflush(stdout);
}


