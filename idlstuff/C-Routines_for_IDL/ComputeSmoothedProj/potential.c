#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


/* This function computes the gravitational potential for ALL the particles.
 * It expects that the particles are predicted to the current time.
 * The routine constructs a new force-tree. 
 */
void compute_potential() 
{
#ifdef GRAPE
  compute_potential_grape();
  return;
#else


  double  t0,t1;
  int     i,j;
  double  r2,fac;



  t0=second();

  if(All.ComovingIntegrationOn)
    set_softenings();

  printf("Start potential computation..."); fflush(stdout);


  force_treebuild();
  All.NumForcesSinceLastTreeConstruction = All.TreeUpdateFrequency*All.TotNumPart ; /* ensures that new tree will be constructed next time*/


  for(i=1;i<=NumPart;i++)
    {
      force_treeevaluate_potential(i-1);    
      P[i].Potential += P[i].Mass/All.SofteningTable[P[i].Type];  /* removes self energy */
    }


  if(All.ComovingIntegrationOn)
    {
      fac=0.5*All.Omega0*All.Hubble*All.Hubble;

      for(i=1;i<=NumPart;i++)
	{
#ifdef PERIODIC
	  P[i].Potential = All.G*P[i].Potential;
#else
	  for(j=0, r2=0; j<3; j++)
	    r2 += P[i].PosPred[j]*P[i].PosPred[j];

	  P[i].Potential = All.G*P[i].Potential - fac*r2;
#endif
	}
    }
  else
    {
      fac= -0.5*All.OmegaLambda*All.Hubble*All.Hubble;

      for(i=1;i<=NumPart;i++)
	{
	  P[i].Potential *= All.G;

	  if(fac!=0)
	     {
	       for(j=0,r2=0;j<3;j++)
		 r2 += P[i].PosPred[j]*P[i].PosPred[j];
	       
	       P[i].Potential += fac*r2;
	     }
	}
    }

  printf("done.\n"); fflush(stdout);

  t1=second();

  All.CPU_Potential+= timediff(t0,t1);

#endif
}






















