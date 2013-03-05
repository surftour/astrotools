#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


/* This function linearly predicts all particles to the time passed as 
 * an argument (will usually be the current time, or a snapshot time).
 * 
 * NOTE: For periodic boundary conditions we do the mapping of coordinates
 * onto the interval [0, All.BoxSize] only before each tree construction.
 * This simplifies dynamic tree updates.
 * In the force computation and neighbour search the correct nearest image
 * will nevertheless always be selected. 
 */
void predict(double time)
{
  int i,j;
  double dt,dt_h0;
  double s_a,s_a_inverse;
  double t0, t1;
  
  t0=second();

  if(All.ComovingIntegrationOn)
    {
      s_a=All.Hubble*sqrt(All.Omega0 + time*(1-All.Omega0-All.OmegaLambda)+pow(time,3)*All.OmegaLambda);
      s_a_inverse=1/s_a;

      for(i=1;i<=NumPart;i++)   
	{
	  dt = time - P[i].CurrentTime; 
	  dt_h0 = dt*s_a_inverse;  

	  for(j=0;j<3;j++)  
	    {
	      P[i].PosPred[j] = P[i].Pos[j] + P[i].Vel[j]*dt_h0;
	      P[i].VelPred[j] = P[i].Vel[j] + P[i].Accel[j]*dt;
	    }

	  if(P[i].Type==0)
	    {
	      SphP[i].DensityPred =  dmax(0, SphP[i].Density + SphP[i].DtDensity*dt);
	      SphP[i].EgySpecPred =  dmax(0, SphP[i].EgySpec + SphP[i].DtEgySpec*dt);
	      SphP[i].Pressure = GAMMA_MINUS1*SphP[i].EgySpecPred*SphP[i].DensityPred;
	    }
	}
    }
  else
    {
      for(i=1;i<=NumPart;i++)   
	{
	  dt = time - P[i].CurrentTime; 
	  
	  for(j=0;j<3;j++)  
	    {
	      P[i].PosPred[j] = P[i].Pos[j] + P[i].Vel[j]*dt;
	      P[i].VelPred[j] = P[i].Vel[j] + P[i].Accel[j]*dt;
	    }
	  
	  if(P[i].Type==0)
	    {
	      SphP[i].DensityPred =  dmax(0, SphP[i].Density + SphP[i].DtDensity*dt);
	      SphP[i].EgySpecPred =  dmax(0, SphP[i].EgySpec + SphP[i].DtEgySpec*dt);
	      SphP[i].Pressure = GAMMA_MINUS1*SphP[i].EgySpecPred*SphP[i].DensityPred;
	    }
	}
    }

  t1=second();
  
  All.CPU_Predict+= timediff(t0,t1);
}


/* This function restricts the prediction to the collisionless particles.
 * Since SPH particles are predicted every timestep, this routine is 
 * called before the tree-construction in order not to save some time
 * by not redoing the SPH prediction a second time.
 */
void predict_collisionless_only(double time)
{
  int i,j;
  double dt,dt_h0;
  double s_a,s_a_inverse;
  double t0, t1;

  t0=second();

  if(All.ComovingIntegrationOn)
    {
      s_a=All.Hubble*sqrt(All.Omega0 + time*(1-All.Omega0-All.OmegaLambda)+pow(time,3)*All.OmegaLambda);
      s_a_inverse=1/s_a;

      for(i=1+N_gas;i<=NumPart;i++)   
	{
	  dt = time - P[i].CurrentTime; 
	  dt_h0 = dt*s_a_inverse;  

	  for(j=0;j<3;j++)  
	    {
	      P[i].PosPred[j] = P[i].Pos[j] + P[i].Vel[j]*dt_h0;
	      P[i].VelPred[j] = P[i].Vel[j] + P[i].Accel[j]*dt;
	    }
	}
    }
  else
    {
      for(i=1+N_gas;i<=NumPart;i++)   
	{
	  dt = time - P[i].CurrentTime; 
	  
	  for(j=0;j<3;j++)  
	    {
	      P[i].PosPred[j] = P[i].Pos[j] + P[i].Vel[j]*dt;
	      P[i].VelPred[j] = P[i].Vel[j] + P[i].Accel[j]*dt;
	    }
	}
    }

  t1=second();
  
  All.CPU_Predict+= timediff(t0,t1);
}



/* This function predicts only the SPH particles.
 */
void predict_sph_particles(double time)
{
  int i,j;
  double dt,dt_h0;
  double s_a,s_a_inverse;
  double t0, t1;

  t0=second();
 
  
  if(All.ComovingIntegrationOn)
    {
      s_a=All.Hubble*sqrt(All.Omega0 + time*(1-All.Omega0-All.OmegaLambda)+pow(time,3)*All.OmegaLambda);
      s_a_inverse=1/s_a;
      
      for(i=1;i<=N_gas;i++)   
	{
	  dt = time - P[i].CurrentTime; 
	  dt_h0 = dt*s_a_inverse;  
	  
	  for(j=0;j<3;j++)  
	    {
	      P[i].PosPred[j] = P[i].Pos[j] + P[i].Vel[j]*dt_h0;
	      P[i].VelPred[j] = P[i].Vel[j] + P[i].Accel[j]*dt;
	    }

	  SphP[i].DensityPred =  dmax(0, SphP[i].Density + SphP[i].DtDensity*dt);
	  SphP[i].EgySpecPred =  dmax(0, SphP[i].EgySpec + SphP[i].DtEgySpec*dt);
	  SphP[i].Pressure = GAMMA_MINUS1*SphP[i].EgySpecPred*SphP[i].DensityPred;
	  
	}
    }
  else
    {
      for(i=1;i<=N_gas;i++)   
	{
	  dt = time - P[i].CurrentTime; 
	  
	  for(j=0;j<3;j++)  
	    {
	      P[i].PosPred[j] = P[i].Pos[j] + P[i].Vel[j]*dt;
	      P[i].VelPred[j] = P[i].Vel[j] + P[i].Accel[j]*dt;
	    }
	  
	  SphP[i].DensityPred =  dmax(0, SphP[i].Density + SphP[i].DtDensity*dt);
	  SphP[i].EgySpecPred =  dmax(0, SphP[i].EgySpec + SphP[i].DtEgySpec*dt);
	  SphP[i].Pressure = GAMMA_MINUS1*SphP[i].EgySpecPred*SphP[i].DensityPred;
	}
    }

  t1=second();
  
  All.CPU_Predict+= timediff(t0,t1);
}


/* This routine advances the active particles for their timestep.
 * Note that no periodic wrapping onto the interval [0, BoxSize]
 * is done here (see comments above).
 */
void advance(void)
{
  int i,j,count;
  double dt,dt_h0,s_a;
  double t0, t1;

  t0=second();

  if(All.ComovingIntegrationOn)
    {
      s_a=All.Hubble*sqrt(All.Omega0 + All.Time*(1-All.Omega0-All.OmegaLambda)+pow(All.Time,3)*All.OmegaLambda);

      for(i=IndFirstUpdate,count=0;count<NumForceUpdate; i=P[i].ForceFlag,count++)   
	{
	  dt = 2*(All.Time - P[i].CurrentTime);  /*  the actual time-step */
	  dt_h0 = dt/s_a; 
      
	  for(j=0;j<3;j++)  
	    {
	      /* update velocities */
	      P[i].Pos[j] = P[i].Pos[j] + 0.5*P[i].Vel[j]*dt_h0;
	      P[i].Vel[j] = P[i].Vel[j] + P[i].Accel[j]*dt;
	      P[i].Pos[j] = P[i].Pos[j] + 0.5*P[i].Vel[j]*dt_h0;
	    }

	  if(P[i].Type==0)
	    {
	      SphP[i].EgySpec   = dmax(All.MinEgySpec, SphP[i].EgySpec   + SphP[i].DtEgySpec*dt);

	      
	      /* Note: density has been updated on the middle of the timestep.
		 Now, predict them at the end of the timestep */
	      
	      SphP[i].Density = SphP[i].Density + SphP[i].DtDensity*0.5*dt;

	      SphP[i].Hsml    = SphP[i].Hsml    + SphP[i].DtHsml*dt; /*estimated new smoothing length */
	      
	      if(SphP[i].Hsml<All.MinGasHsml)
		SphP[i].Hsml= All.MinGasHsml;
	    }
	      
	  P[i].CurrentTime = All.Time + 0.5*dt;
	}
    }
  else
    {
  
      for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag,count++)   
	{
	  dt = 2*(All.Time - P[i].CurrentTime);  /*  the actual time-step */
	  
	  for(j=0;j<3;j++)  
	    {
	      P[i].Pos[j] = P[i].Pos[j] + 0.5*P[i].Vel[j]*dt;
	      P[i].Vel[j] = P[i].Vel[j] + P[i].Accel[j]*dt; 	      /* advance velocities */
	      P[i].Pos[j] = P[i].Pos[j] + 0.5*P[i].Vel[j]*dt;
	    }
	  
	  if(P[i].Type==0)
	    {
	      SphP[i].EgySpec   = dmax(All.MinEgySpec, SphP[i].EgySpec   + SphP[i].DtEgySpec*dt);
	      
	      /* Note: density and h have been updated on the middle of the timestep.
		 Now, predict them at the end of the timestep */
	      
	      SphP[i].Density = SphP[i].Density + SphP[i].DtDensity*0.5*dt;
	      SphP[i].Hsml    = SphP[i].Hsml    + SphP[i].DtHsml*dt;
	      
	      if(SphP[i].Hsml<All.MinGasHsml)
		SphP[i].Hsml= All.MinGasHsml;
	    }

	  P[i].CurrentTime = All.Time + 0.5*dt;
	}
    }

  t1=second();
  
  All.CPU_Predict+= timediff(t0,t1);
}


/*  This function makes sure that all particles coordinates (Pos) are
 *  mapped onto the interval [0, BoxSize]. 
 *  After this function has been called, a new domain decomposition
 *  should be done, or at the very least a new force-tree needs to
 *  be constructed.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i,j;

  for(i=1; i<=NumPart; i++)   
    for(j=0; j<3; j++)  
      {
	while(P[i].Pos[j] < 0)
	  {
	    P[i].Pos[j] += All.BoxSize;
	    P[i].PosPred[j] += All.BoxSize;
	  }
	
	while(P[i].Pos[j] > All.BoxSize)
	  {
	    P[i].Pos[j] -= All.BoxSize;
	    P[i].PosPred[j] -= All.BoxSize;
	  }
      }
}    
#endif

