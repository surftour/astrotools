#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


/*
 *  init() reads in the initial conditions,
 *  and allocates storage for the tree(s).
 *  An intial force computation is done, 
 *  then the first particle timesteps
 *  are determined. The simulation is set up for
 *  the timestep iteration in run().
 */
void init_2d(void)
{
  int i,j;

  for(i=1; i<=NumPart_2d; i++) /*  start-up initialization */
    {
      for(j=0;j<3;j++)
	{
	  Pn[i].PosPred[j]=Pn[i].Pos[j];
	  Pn[i].VelPred[j]=Pn[i].Vel[j];
	  Pn[i].Accel[j]=0;
	}
      Pn[i].OldAcc=0;
      Pn[i].Potential=0;
    }


  
  for(i=1; i<=N_gas_2d; i++)     /* initialize sph properties */
    {
      SphPn[i].EgySpecPred  = SphPn[i].EgySpec;
      SphPn[i].DtEgySpec=0;
      SphPn[i].DtDensity = SphPn[i].Density=0;
      SphPn[i].DtHsml = SphPn[i].Hsml=0;
#ifdef COOLING
      SphPn[i].Ne= 1.0;
#endif
    }
  
  //o.k.  
  force_treeallocate_2d(All.TreeAllocFactor*All.MaxPart, All.MaxPart);

  //o.k.
  ngb_treeallocate_2d(MAX_NGB);

  for(i=1; i<=NumPart_2d; i++) 
    {
      Pn[i].CurrentTime=All.TimeBegin;
      Pn[i].ForceFlag=i+1;
    }
  Pn[NumPart_2d].ForceFlag=1; IndFirstUpdate=1; NumForceUpdate=NumPart_2d; NumSphUpdate=N_gas_2d;


  //o.k.
  ngb_treebuild_2d();

  All.NumForcesSinceLastTreeConstruction= All.TreeUpdateFrequency*All.TotNumPart_2d ; /* ensures that new tree will be constructed */

}
