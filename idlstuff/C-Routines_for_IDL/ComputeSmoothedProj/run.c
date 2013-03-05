#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "allvars.h"
#include "proto.h"



/* This routine contains the main simulation loop that iterates over
 * the single timesteps. The loop terminates when the cpu-time
 * limit is reached, when a `stop' file is found in the output 
 * directory, or when the simulation ends because we arrived
 * at TimeMax
 */
void run(void) 
{
  FILE   *fd;
  int    stopflag=0;
  double savetime;
  char   buf[200],stopfname[200];
  double t0,t1;


  sprintf(stopfname,"%sstop",All.OutputDir);


  do       /* main loop */
    {
      t0=second();

      find_next_time();   /* increase Time to smallest prediction time and
                             determine which particles are grouped together
                             for force evaluation */


      if(All.Time > All.TimeMax)
	All.Time = All.TimeMax;

#ifdef COOLING
      IonizeParams();
#endif

      every_timestep_stuff();
 
      if((All.Time-All.TimeLastStatistics)>=All.TimeBetStatistics)
	{
	  savetime=All.Time;
	  predict(All.Time=All.TimeLastStatistics+All.TimeBetStatistics);
	  compute_potential();  /* needed for potential energy */
 	  energy_statistics();
	  All.TimeLastStatistics += All.TimeBetStatistics;
	  All.Time=savetime;
	}

      if((All.Time-All.TimeOfFirstSnapshot)>=0) /* check whether it's time for a snapshot file */
	{ 
	  savetime=All.Time;
	  predict(All.Time=All.TimeOfFirstSnapshot);
	  savepositions(All.SnapshotFileCount++);   /* write snapshot file */
	  if(All.OutputListOn)
	    All.TimeOfFirstSnapshot= find_next_outputtime(savetime);
	  else
	    if(All.ComovingIntegrationOn)
	      All.TimeOfFirstSnapshot *= All.TimeBetSnapshot;
	    else
	      All.TimeOfFirstSnapshot += All.TimeBetSnapshot;
	  All.Time=savetime;
	}


      predict_sph_particles(All.Time);  /* SPH particles are alwyas predicted, while
                                           this is done for collisionless particles
                                           either before the tree construction,
                                           or on the fly while the tree is walked */

   
      compute_accelerations(0);     /* ... compute accelerations for 
                                           the particles that are to be advanced 
                                           Note: the particle positions and tree nodes are 
                                           predicted for the current time during the force 
					   computation */


      advance();                    /* advance the active particles */



      find_timesteps(0);               /* compute new timesteps for these particles, 
                                          and set their maximum prediction times, 
                                          and reinsert them into the timeline as needed */


      All.NumCurrentTiStep++;


      if((fd=fopen(stopfname,"r")))  /* test for stop-file */
	{
	  fclose(fd);
	  stopflag=1;
	  sprintf(buf,"rm -f %s",stopfname); /* remove it */
	  system(buf);
	}

      if(CPUThisRun > 0.90*All.TimeLimitCPU)  /* produce a restart file if CPU-time expires soon */
	{
	  printf("reaching time-limit. stopping.\n");
	  stopflag=2;
	}

      if(stopflag)
	{
	  restart(0);  /* write restart file */

	  if(stopflag==2 && All.ResubmitOn)
	    {
	      close_outputfiles();
	      sprintf(buf,"%s", All.ResubmitCommand);
	      system(buf); 
	    }
	  return;
	}
      
      if((CPUThisRun - All.TimeLastRestartFile)>= All.CpuTimeBetRestartFile) 
	{
	  All.TimeLastRestartFile= CPUThisRun;
	  restart(0);
	}

      t1=second();

      All.CPU_Total+= timediff(t0,t1);
      CPUThisRun+= timediff(t0,t1);
    }
  while(All.Time < All.TimeMax);

  restart(0);  /* write a restart file at end of simulation (it can be continued later on) */


 /* write a last snapshot file at final time (will be overwritten if 
     All.TimeMax is increased and the run is continued !)*/

  savetime=All.Time;
  predict(All.Time=All.TimeMax);
  compute_potential();
  energy_statistics();
  savepositions(All.SnapshotFileCount++);   /* write snapshot file */
  All.Time=savetime;
}



/* This routine writes one line for every timestep to two log-files.
 * In FdInfo, we just list the timesteps that have been done,
 * while in FdCPU the cumulative cpu-time consumption in various parts
 * of the code is stored.
 */
void every_timestep_stuff(void)
{
  double z;


  if(All.ComovingIntegrationOn)
    {
      z=1/All.Time-1;
      fprintf(FdInfo,"\nBegin Timestep %d, Time: %g, Redshift: %g, Timestep: %g\n",All.NumCurrentTiStep,All.Time,z,All.TimeStep);
      printf(        "\nBegin Timestep %d, Time: %g, Redshift: %g, Timestep: %g\n",All.NumCurrentTiStep,All.Time,z,All.TimeStep);
      fflush(FdInfo);
    }
  else
    {
      fprintf(FdInfo,"\nBegin Timestep %d, Time: %g, Timestep: %g\n",All.NumCurrentTiStep,All.Time,All.TimeStep);
      printf(        "\nBegin Timestep %d, Time: %g, Timestep: %g\n",All.NumCurrentTiStep,All.Time,All.TimeStep);
      fflush(FdInfo);
    }


  fprintf(FdCPU,"Timestep %d, Time: %g\n",All.NumCurrentTiStep,All.Time);
  fprintf(FdCPU,"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
	  All.CPU_Total,
	  All.CPU_Gravity,
	  All.CPU_Hydro,
	  All.CPU_Potential,
	  All.CPU_Predict,
	  All.CPU_TimeLine,
	  All.CPU_Snapshot,
	  All.CPU_TreeWalk,
	  All.CPU_TreeConstruction);

  fflush(FdCPU);
}



/* This routine first calls a computation of various global quantities
 * of the particle distribution, and then writes some statistics
 * about the energies in the various particle components to the 
 * file FdEnergy.
 */
void energy_statistics(void)
{
  compute_global_quantities_of_system();
  
  fprintf(FdEnergy,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	  All.Time,
	  SysState.EnergyInt,
	  SysState.EnergyPot,
	  SysState.EnergyKin,
	  SysState.EnergyIntComp[0],
	  SysState.EnergyPotComp[0],
	  SysState.EnergyKinComp[0],
	  SysState.EnergyIntComp[1],
	  SysState.EnergyPotComp[1],
	  SysState.EnergyKinComp[1],
	  SysState.EnergyIntComp[2],
	  SysState.EnergyPotComp[2],
	  SysState.EnergyKinComp[2],
	  SysState.EnergyIntComp[3],
	  SysState.EnergyPotComp[3],
	  SysState.EnergyKinComp[3],
	  SysState.EnergyIntComp[4],
	  SysState.EnergyPotComp[4],
	  SysState.EnergyKinComp[4],
	  SysState.MassComp[0],
	  SysState.MassComp[1],
	  SysState.MassComp[2],
	  SysState.MassComp[3],
	  SysState.MassComp[4]
	  ); 

  fflush(FdEnergy);
}

