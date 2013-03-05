#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "allvars.h"
#include "proto.h"



/*  This function computes the gravitational forces for all active particles.
 *  A new tree is constructed, if the number of force computations since
 *  it's last construction exceeds some fraction of the total
 *  particle number, otherwise tree nodes are dynamically updated if needed.
 */
void gravity_tree(void)
{
  int     i,j;
  double  tstart,tend,timebuild,timewalk;
  static int numnodes;
  int     count;
  double  dt,dt_h0;
  double  s_a, s_a_inverse,  a, a2, fac1, fac2, fac3;
#ifdef DIAG
  int  costtotal, costtotal_quadru;
#endif


  if(All.ComovingIntegrationOn)
    {
      set_softenings();   /* set new softening lengths */

      s_a_inverse=1/(All.Hubble*sqrt(All.Omega0 + All.Time*(1-All.Omega0-All.OmegaLambda)+ All.Time*All.Time*All.Time*All.OmegaLambda));
    }
  else
    s_a_inverse=1;
  All.NumForcesSinceLastTreeConstruction+= NumForceUpdate;


  tstart=second();
  if(All.NumForcesSinceLastTreeConstruction >= All.TreeUpdateFrequency*NumPart)
    {
#ifdef PERIODIC
      do_box_wrapping();    /* map the particles back onto the box */
#endif	  
      predict_collisionless_only(All.Time);  /* gas particles have already been predicted */
      numnodes= force_treebuild();
      All.NumForcesSinceLastTreeConstruction=0;
    }
  else
    {
#ifdef VELDISP
      predict_collisionless_only(All.Time); 
#endif
      ngb_update_nodes();  /* need to make sure that we have outer boundaries for all nodes */
    }

  tend=second();
  timebuild=timediff(tstart,tend);
  All.CPU_TreeConstruction += timebuild;


  force_resetcost();  /* resets counters for statistics of number of particle-node interactions */
 
  tstart=second();
  for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
    {
      /* need to predict the particle, since it may not have been done yet */
      /* NOTE: the predicted velocity is needed in the comoving integration below */
      dt = (All.Time - P[i].CurrentTime);  
      dt_h0 = dt*s_a_inverse; 

      for(j=0;j<3;j++)  
	{
	  P[i].PosPred[j] = P[i].Pos[j] + P[i].Vel[j]*dt_h0;
	  P[i].VelPred[j] = P[i].Vel[j] + P[i].Accel[j]*dt;
	}

      force_treeevaluate(i-1, s_a_inverse);
    }
  tend=second();
  timewalk=timediff(tstart,tend);
  All.CPU_TreeWalk+=timewalk;


 
  if(All.ComovingIntegrationOn)
    {
      if(All.TypeOfOpeningCriterion==1)
	{
	  fac3= 0.5*All.Hubble*All.Hubble*All.Omega0 / All.G;
	  for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag,count++)   
	    {
	      for(j=0,a2=0;j<3;j++)
		{
#ifdef PERIODIC
		  a = P[i].Accel[j];
#else
		  a = P[i].Accel[j]+ fac3*P[i].PosPred[j];
#endif
		  a2 += a*a;
		}
	      P[i].OldAcc= sqrt(a2);
	    }
	}


      s_a=sqrt(All.Omega0 + All.Time*(1-All.Omega0-All.OmegaLambda)+ All.Time*All.Time*All.Time*All.OmegaLambda);

      fac1= All.G/ ( All.Hubble * All.Time * All.Time* s_a );

      fac2= -1.5/All.Time;

      fac3= 0.5*All.Hubble * All.Omega0 / ( All.Time * All.Time * s_a );
	
      for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag,count++)   
	{
	  for(j=0;j<3;j++)
#ifdef PERIODIC
	    P[i].Accel[j] = fac1*P[i].Accel[j] 
       	                  + fac2*P[i].VelPred[j];
#else
	    P[i].Accel[j] = fac1*P[i].Accel[j] 
	                  + fac2*P[i].VelPred[j] 
	                  + fac3*P[i].PosPred[j];  
#endif
	}
    }
  else
    {
      if(All.TypeOfOpeningCriterion==1)
	for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag,count++)   
	  {
	    P[i].OldAcc= sqrt(P[i].Accel[0]*P[i].Accel[0] + 
				  P[i].Accel[1]*P[i].Accel[1] +
				  P[i].Accel[2]*P[i].Accel[2] );
	  }

      fac1= All.OmegaLambda*All.Hubble*All.Hubble;
             /* this factor allows a computation of cosmological simulation 
                 with vacuum energy in physical coordinates */

      for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag,count++)   
	{
	  for(j=0;j<3;j++)
	    P[i].Accel[j] = All.G*P[i].Accel[j]  + fac1*P[i].PosPred[j];
	}
    }

#ifdef DIAG
  /*  gather some diagnostic information */
  costtotal= force_getcost_single();
  costtotal_quadru= force_getcost_quadru();

  All.TotNumOfForces+= NumForceUpdate;

  fprintf(FdTimings,"Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
  fprintf(FdTimings,"Nf= %d  total-Nf= %d\n", NumForceUpdate, All.TotNumOfForces);
  fprintf(FdTimings,"nodeupdates: %d (%d part)\n", Num_nodeupdates, Num_nodeupdate_particles);
  fprintf(FdTimings,"nodes: %d, filled: %g\n", numnodes, numnodes/(All.TreeAllocFactor*All.MaxPart));
  fprintf(FdTimings,"part/sec=%g  ia/part=%g (QP-frac %g)\n", 
	  NumForceUpdate/(timewalk+1.0e-20), 
	  ((double)(costtotal+costtotal_quadru))/NumForceUpdate, 
	  costtotal_quadru/((double)(costtotal+costtotal_quadru)));
  fprintf(FdTimings,"\n");
  fflush(FdTimings);
#endif
 
}




void set_softenings(void)
{
  if(All.SofteningGas*All.Time > All.SofteningGasMaxPhys)
    All.SofteningTable[0] = All.SofteningGasMaxPhys/All.Time;
  else
    All.SofteningTable[0] = All.SofteningGas;
  
  if(All.SofteningHalo*All.Time > All.SofteningHaloMaxPhys)
    All.SofteningTable[1] = All.SofteningHaloMaxPhys/All.Time;
  else
    All.SofteningTable[1] = All.SofteningHalo;
  
  if(All.SofteningDisk*All.Time > All.SofteningDiskMaxPhys)
    All.SofteningTable[2] = All.SofteningDiskMaxPhys/All.Time;
  else
    All.SofteningTable[2] = All.SofteningDisk;

  if(All.SofteningBulge*All.Time > All.SofteningBulgeMaxPhys)
    All.SofteningTable[3] = All.SofteningBulgeMaxPhys/All.Time;
  else
    All.SofteningTable[3] = All.SofteningBulge;

  if(All.SofteningStars*All.Time > All.SofteningStarsMaxPhys)
    All.SofteningTable[4] = All.SofteningStarsMaxPhys/All.Time;
  else
    All.SofteningTable[4] = All.SofteningStars;
}











