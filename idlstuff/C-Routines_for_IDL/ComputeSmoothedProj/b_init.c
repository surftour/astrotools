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
void init(void)
{
  int i,j;

  All.Time = All.TimeBegin;

  All.SofteningTable[0] = All.SofteningGas;
  All.SofteningTable[1] = All.SofteningHalo;
  All.SofteningTable[2] = All.SofteningDisk;
  All.SofteningTable[3] = All.SofteningBulge;
  All.SofteningTable[4] = All.SofteningStars;
  All.MinGasHsml= All.MinGasHsmlFractional*All.SofteningTable[0]; 


  All.NumCurrentTiStep= 0;   /* ... setup counters */
  All.SnapshotFileCount=0;
  if(RestartFlag==2)
    All.SnapshotFileCount= atoi(All.InitCondFile+strlen(All.InitCondFile)-3)+1;
 
  All.TotNumOfForces=0;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn==1)
      check_omega();

 if(RestartFlag==2)
    {
      if(All.ComovingIntegrationOn)
	All.TimeOfFirstSnapshot = All.Time*All.TimeBetSnapshot;
      else
	All.TimeOfFirstSnapshot = All.Time+All.TimeBetSnapshot;
    }


  All.TimeLastStatistics= All.TimeBegin-All.TimeBetStatistics;


  for(i=1; i<=NumPart; i++) /*  start-up initialization */
    {
      for(j=0;j<3;j++)
	{
	  P[i].PosPred[j]=P[i].Pos[j];
	  P[i].VelPred[j]=P[i].Vel[j];
	  P[i].Accel[j]=0;
	}
      P[i].OldAcc=0;
      P[i].Potential=0;
    }


  
  for(i=1; i<=(N_gas); i++)     /* initialize sph properties */
    {
      SphP[i].EgySpecPred  = SphP[i].EgySpec;
      SphP[i].DtEgySpec=0;
      SphP[i].DtDensity = SphP[i].Density=0;
      SphP[i].DtHsml = SphP[i].Hsml=0;
#ifdef COOLING
      SphP[i].Ne= 1.0;
#endif
    }
  
  if(All.ComovingIntegrationOn==0)
    compute_global_quantities_of_system(); 
   
  force_treeallocate(All.TreeAllocFactor*All.MaxPart, All.MaxPart);
  ngb_treeallocate(MAX_NGB);

  for(i=1; i<=NumPart; i++) 
    {
      P[i].CurrentTime=All.TimeBegin;
      P[i].ForceFlag=i+1;
    }
  P[NumPart].ForceFlag=1; IndFirstUpdate=1; NumForceUpdate=NumPart; NumSphUpdate=N_gas;

  setup_smoothinglengths(All.DesNumNgb);

#ifdef VELDISP
  setup_smoothinglengths_veldisp(All.DesNumNgb);
  veldisp();
#endif


  All.NumForcesSinceLastTreeConstruction= All.TreeUpdateFrequency*All.TotNumPart ; /* ensures that new tree will be constructed */

  //compute_accelerations(1);        /* ... compute accelerations */
  

  //find_timesteps(2);               /* ... set particles to initial timesteps   */
                                   /* ordered timeline will be constructed there */

  //compute_global_quantities_of_system();
  SysStateAtStart = SysState;     /* ... remember global initial state */
}


/* This routine computes the mass content of the box and
 * compares it to the specified value of Omega.
 * If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass=0, omega;
  int    i;

  for(i=1; i<=NumPart; i++)
    mass+= P[i].Mass;

  omega= mass/(All.BoxSize*All.BoxSize*All.BoxSize)/ (3*All.Hubble*All.Hubble/(8*PI*All.G));

  if(fabs(omega-All.Omega0) > 1.0e-3)
    {
      printf("\n\nI've found something odd!\n");
      printf("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n", 
	     omega, All.Omega0);
      printf("\nI better stop.\n");
      endrun(0);
    }
}


/*
 *  This function is used to find an initial smoothing length 
 *  for each SPH particle. 
 */
void setup_smoothinglengths(int desired_ngb)
{
  int    i;
  float  *r2list;
  int    *ngblist;

  ngb_treebuild();

  for(i=1; i<=N_gas; i++) 
    SphP[i].Hsml= sqrt(ngb_treefind( P[i].PosPred, desired_ngb, 0, 0, &ngblist, &r2list));  
}



#ifdef VELDISP 
/*
 *  This function is used to find an initial smoothing length for each 
 *  dark matter particle. 
 */
void setup_smoothinglengths_veldisp(int desired_ngb)
{
  int    i;
  float  *r2list;
  int    *ngblist;

  for(i=1+N_gas; i<=NumPart; i++) 
    P[i].HsmlVelDisp= sqrt(ngb_treefind( P[i].PosPred, desired_ngb, 0, P[i].Type, &ngblist, &r2list));  
}

#endif  







