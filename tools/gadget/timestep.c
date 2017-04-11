#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif

/*! \file timestep.c 
 *  \brief routines for 'kicking' particles in
 *  momentum space and assigning new timesteps
 */

static double fac1, fac2, fac3, hubble_a, atime, a3inv;
static double dt_displacement = 0;

#ifdef PMGRID
static double dt_gravkickA, dt_gravkickB;
#endif


/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */
void advance_and_find_timesteps(void)
{
  int i, ti_step, ti_step_old, ti_min, tend, tstart, bin, binold, prev, next;
  double aphys;

#ifdef PMGRID
  int j, dt_step;
  double dt_gravkick, dt_hydrokick;

#ifdef DISTORTIONTENSORPS
  /* for the distortion 'velocity part', so only the lower two 3x3 submatrices will be != 0 */
  MyDouble dv_distortion_tensorps[6][6];
  int j1, j2;
#endif
#endif
#ifdef MAKEGLASS
  double disp, dispmax, globmax, dmean, fac, disp2sum, globdisp2sum;
#endif
#if defined(TIME_DEP_ART_VISC) || defined(TIME_DEP_MAGN_DISP)
  double dmin1, dmin2;
#endif
#ifdef CHEMISTRY
  int ifunc, mode;
  double a_start, a_end;
#endif
#ifdef WAKEUP
  int n, k, dt_bin, dt_step, ti_next_for_bin, ti_next_kick, ti_next_kick_global, max_time_bin_active;
  double dt_entr;
  long long ntot;
#endif

  CPU_Step[CPU_MISC] += measure_time();

  if(All.ComovingIntegrationOn)
    {
      fac1 = 1 / (All.Time * All.Time);
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
      fac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      hubble_a = hubble_function(All.Time);
      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    fac1 = fac2 = fac3 = hubble_a = a3inv = atime = 1;

  if(Flag_FullStep || dt_displacement == 0)
    find_dt_displacement_constraint(hubble_a * atime * atime);

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkickB = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkickB = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
#endif


#ifdef MAKEGLASS
  for(i = 0, dispmax = 0, disp2sum = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].g.GravAccel[j] *= -1;
#ifdef PMGRID
	  P[i].GravPM[j] *= -1;
	  P[i].g.GravAccel[j] += P[i].GravPM[j];
	  P[i].GravPM[j] = 0;
#endif
	}

      disp = sqrt(P[i].g.GravAccel[0] * P[i].g.GravAccel[0] +
		  P[i].g.GravAccel[1] * P[i].g.GravAccel[1] + P[i].g.GravAccel[2] * P[i].g.GravAccel[2]);

      disp *= 2.0 / (3 * All.Hubble * All.Hubble);

      disp2sum += disp * disp;

      if(disp > dispmax)
	dispmax = disp;
    }

  MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  dmean = pow(P[0].Mass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  if(globmax > dmean)
    fac = dmean / globmax;
  else
    fac = 1.0;

  if(ThisTask == 0)
    {
      printf("\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n",
	     dmean, globmax, sqrt(globdisp2sum / All.TotNumPart));
      fflush(stdout);
    }

  for(i = 0, dispmax = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Vel[j] = 0;
	  P[i].Pos[j] += fac * P[i].g.GravAccel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
	  P[i].g.GravAccel[j] = 0;
	}
    }
#endif




  All.DoDynamicUpdate = ShouldWeDoDynamicUpdate();

  /* Now assign new timesteps and kick */

  if(All.DoDynamicUpdate)
    {
      GlobFlag++;
      DomainNumChanged = 0;
      DomainList = (int *) mymalloc(NTopleaves * sizeof(int));
      if(ThisTask == 0)
	printf("kicks will prepare for dynamic update of tree\n");
    }

#ifdef FORCE_EQUAL_TIMESTEPS
  for(i = FirstActiveParticle, ti_min = TIMEBASE; i >= 0; i = NextActiveParticle[i])
    {
      ti_step = get_timestep(i, &aphys, 0);

      if(ti_step < ti_min)
	ti_min = ti_step;
    }

  int ti_min_glob;

  MPI_Allreduce(&ti_min, &ti_min_glob, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef FORCE_EQUAL_TIMESTEPS
      ti_step = ti_min_glob;
#else
      ti_step = get_timestep(i, &aphys, 0);
#endif

      /* make it a power 2 subdivision */
      ti_min = TIMEBASE;
      while(ti_min > ti_step)
	ti_min >>= 1;
      ti_step = ti_min;

      bin = get_timestep_bin(ti_step);
      binold = P[i].TimeBin;

      if(bin > binold)		/* timestep wants to increase */
	{
	  if(TimeBinActive[bin] == 0)	/* leave at old step if not synchronized */
	    {
	      bin = binold;
	      ti_step = bin ? (1 << bin) : 0;
	    }
	}

      if(All.Ti_Current >= TIMEBASE)	/* we here finish the last timestep. */
	{
	  ti_step = 0;
	  bin = 0;
	}

      if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
	{
	  fprintf(stderr, "\n @ /* should not happen */ \n");
	  endrun(888);		/* should not happen */
	  fprintf(stderr, "\n @ /* should not happen */ \n");
	  ti_step = TIMEBASE - All.Ti_Current;
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;
	}

      if(bin != binold)
	{
	  TimeBinCount[binold]--;
	  if(P[i].Type == 0)
	    {
	      TimeBinCountSph[binold]--;
#ifdef SFR
	      TimeBinSfr[binold] -= SphP[i].Sfr;
	      TimeBinSfr[bin] += SphP[i].Sfr;
#endif
	    }

#ifdef BLACK_HOLES
	  if(P[i].Type == 5)
	    {
	      TimeBin_BH_mass[binold] -= P[i].BH_Mass;
	      TimeBin_BH_dynamicalmass[binold] -= P[i].Mass;
	      TimeBin_BH_Mdot[binold] -= P[i].BH_Mdot;
	      if(P[i].BH_Mass > 0)
		TimeBin_BH_Medd[binold] -= P[i].BH_Mdot / P[i].BH_Mass;
	      TimeBin_BH_mass[bin] += P[i].BH_Mass;
	      TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
	      TimeBin_BH_Mdot[bin] += P[i].BH_Mdot;
	      if(P[i].BH_Mass > 0)
		TimeBin_BH_Medd[bin] += P[i].BH_Mdot / P[i].BH_Mass;
	    }
#endif

	  prev = PrevInTimeBin[i];
	  next = NextInTimeBin[i];

	  if(FirstInTimeBin[binold] == i)
	    FirstInTimeBin[binold] = next;
	  if(LastInTimeBin[binold] == i)
	    LastInTimeBin[binold] = prev;
	  if(prev >= 0)
	    NextInTimeBin[prev] = next;
	  if(next >= 0)
	    PrevInTimeBin[next] = prev;

	  if(TimeBinCount[bin] > 0)
	    {
	      PrevInTimeBin[i] = LastInTimeBin[bin];
	      NextInTimeBin[LastInTimeBin[bin]] = i;
	      NextInTimeBin[i] = -1;
	      LastInTimeBin[bin] = i;
	    }
	  else
	    {
	      FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
	      PrevInTimeBin[i] = NextInTimeBin[i] = -1;
	    }
	  TimeBinCount[bin]++;
	  if(P[i].Type == 0)
	    TimeBinCountSph[bin]++;

	  P[i].TimeBin = bin;
	}

#ifndef WAKEUP
      ti_step_old = binold ? (1 << binold) : 0;
#else
      ti_step_old = P[i].dt_step;
#endif

      tstart = P[i].Ti_begstep + ti_step_old / 2;	/* midpoint of old step */
      tend = P[i].Ti_begstep + ti_step_old + ti_step / 2;	/* midpoint of new step */

      P[i].Ti_begstep += ti_step_old;
#ifdef WAKEUP
      P[i].dt_step = ti_step;
#endif

      do_the_kick(i, tstart, tend, P[i].Ti_begstep);
    }

  if(All.DoDynamicUpdate)
    {
      force_finish_kick_nodes();
      myfree(DomainList);
    }


#ifdef CONDUCTION
  if(All.Conduction_Ti_endstep == All.Ti_Current)
    {
      ti_step = TIMEBASE;
      while(ti_step > (All.MaxSizeConductionStep / All.Timebase_interval))
	ti_step >>= 1;
      while(ti_step > (All.MaxSizeTimestep / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.Conduction_Ti_endstep - All.Conduction_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.Conduction_Ti_endstep) % ti_step) > 0)
	    ti_step = All.Conduction_Ti_endstep - All.Conduction_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      All.Conduction_Ti_begstep = All.Conduction_Ti_endstep;
      All.Conduction_Ti_endstep = All.Conduction_Ti_begstep + ti_step;
    }
#endif

#ifdef CR_DIFFUSION
  if(All.CR_Diffusion_Ti_endstep == All.Ti_Current)
    {
      if(All.CR_Diffusion_Ti_endstep < All.Ti_Current)
	endrun(1231);

      ti_step = TIMEBASE;
      while(ti_step > (All.CR_DiffusionMaxSizeTimestep / All.Timebase_interval))
	ti_step >>= 1;
      while(ti_step > (All.MaxSizeTimestep / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.CR_Diffusion_Ti_endstep - All.CR_Diffusion_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.CR_Diffusion_Ti_endstep) % ti_step) > 0)
	    ti_step = All.CR_Diffusion_Ti_endstep - All.CR_Diffusion_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      All.CR_Diffusion_Ti_begstep = All.CR_Diffusion_Ti_endstep;
      All.CR_Diffusion_Ti_endstep = All.CR_Diffusion_Ti_begstep + ti_step;
    }
#endif


#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      All.DoDynamicUpdate = 0;

      ti_step = TIMEBASE;
      while(ti_step > (dt_displacement / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.PM_Ti_endstep) % ti_step) > 0)
	    ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      tstart = (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2;
      tend = All.PM_Ti_endstep + ti_step / 2;

      if(All.ComovingIntegrationOn)
	dt_gravkick = get_gravkick_factor(tstart, tend);
      else
	dt_gravkick = (tend - tstart) * All.Timebase_interval;

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;

      if(All.ComovingIntegrationOn)
	dt_gravkickB = -get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
      else
	dt_gravkickB =
	  -((All.PM_Ti_begstep + All.PM_Ti_endstep) / 2 - All.PM_Ti_begstep) * All.Timebase_interval;

      for(i = 0; i < NumPart; i++)
	{
	  for(j = 0; j < 3; j++)	/* do the kick */
	    P[i].Vel[j] += P[i].GravPM[j] * dt_gravkick;

#ifdef DISTORTIONTENSORPS
	  /* add long range tidal forces calculated on mesh */
	  /* now we do the distortiontensor kick */
	  for(j1 = 0; j1 < 3; j1++)
	    for(j2 = 0; j2 < 3; j2++)
	      {
		dv_distortion_tensorps[j1 + 3][j2] = 0.0;
		dv_distortion_tensorps[j1 + 3][j2 + 3] = 0.0;

		/* the 'acceleration' is given by the product of tidaltensor and distortiontensor */
		for(j = 0; j < 3; j++)
		  {
		    dv_distortion_tensorps[j1 + 3][j2] +=
		      P[i].tidal_tensorpsPM[j1][j] * P[i].distortion_tensorps[j][j2];
		    dv_distortion_tensorps[j1 + 3][j2 + 3] +=
		      P[i].tidal_tensorpsPM[j1][j] * P[i].distortion_tensorps[j][j2 + 3];
		  }
		dv_distortion_tensorps[j1 + 3][j2] *= dt_gravkick;
		dv_distortion_tensorps[j1 + 3][j2 + 3] *= dt_gravkick;
		/* add it to the distortiontensor 'velocities' */
		P[i].distortion_tensorps[j1 + 3][j2] += dv_distortion_tensorps[j1 + 3][j2];
		P[i].distortion_tensorps[j1 + 3][j2 + 3] += dv_distortion_tensorps[j1 + 3][j2 + 3];
	      }
#endif



	  if(P[i].Type == 0)
	    {
#ifndef WAKEUP
	      dt_step = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);
#else
	      dt_step = P[i].dt_step;
#endif

	      if(All.ComovingIntegrationOn)
		{
		  dt_gravkickA = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_gravkick_factor(P[i].Ti_begstep, P[i].Ti_begstep + dt_step / 2);
		  dt_hydrokick = get_hydrokick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_hydrokick_factor(P[i].Ti_begstep, P[i].Ti_begstep + dt_step / 2);
		}
	      else
		dt_gravkickA = dt_hydrokick =
		  (All.Ti_Current - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;

	      for(j = 0; j < 3; j++)
		SphP[i].VelPred[j] = P[i].Vel[j]
		  + P[i].g.GravAccel[j] * dt_gravkickA
		  + SphP[i].a.HydroAccel[j] * dt_hydrokick + P[i].GravPM[j] * dt_gravkickB;
	    }
	}
    }
#endif

#ifdef WAKEUP
  /* find the next kick time */
  for(n = 0, ti_next_kick = TIMEBASE; n < TIMEBINS; n++)
    {
      if(TimeBinCount[n])
	{
	  if(n > 0)
	    {
	      dt_bin = (1 << n);
	      ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;	/* next kick time for this timebin */
	    }
	  else
	    {
	      dt_bin = 0;
	      ti_next_for_bin = All.Ti_Current;
	    }

	  if(ti_next_for_bin < ti_next_kick)
	    ti_next_kick = ti_next_for_bin;
	}
    }

  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("predicting next timestep: %g\n", (ti_next_kick_global - All.Ti_Current) * All.Timebase_interval);

  max_time_bin_active = 0;
  /* get the highest bin, that is active next time */
  for(n = 0; n < TIMEBINS; n++)
    {
      dt_bin = (1 << n);

      if((ti_next_kick_global % dt_bin) == 0)
	max_time_bin_active = n;
    }

  /* move the particle on the highest bin, that is active in the next timestep and that is lower than its last timebin */
  bin = 0;
  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinCount[n] > 0)
	{
	  bin = n;
	  break;
	}
    }
  n = 0;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type != 0)
	continue;

      if(!SphP[i].wakeup)
	continue;

      binold = P[i].TimeBin;
      if(TimeBinActive[binold])
	continue;

      bin = max_time_bin_active < binold ? max_time_bin_active : binold;

      if(bin != binold)
	{
	  TimeBinCount[binold]--;
	  if(P[i].Type == 0)
	    TimeBinCountSph[binold]--;

	  prev = PrevInTimeBin[i];
	  next = NextInTimeBin[i];

	  if(FirstInTimeBin[binold] == i)
	    FirstInTimeBin[binold] = next;
	  if(LastInTimeBin[binold] == i)
	    LastInTimeBin[binold] = prev;
	  if(prev >= 0)
	    NextInTimeBin[prev] = next;
	  if(next >= 0)
	    PrevInTimeBin[next] = prev;

	  if(TimeBinCount[bin] > 0)
	    {
	      PrevInTimeBin[i] = LastInTimeBin[bin];
	      NextInTimeBin[LastInTimeBin[bin]] = i;
	      NextInTimeBin[i] = -1;
	      LastInTimeBin[bin] = i;
	    }
	  else
	    {
	      FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
	      PrevInTimeBin[i] = NextInTimeBin[i] = -1;
	    }
	  TimeBinCount[bin]++;
	  if(P[i].Type == 0)
	    TimeBinCountSph[bin]++;

	  P[i].TimeBin = bin;

	  /* correct quantities predicted for a longer timestep */
	  ti_step_old = P[i].dt_step;
	  dt_step = ti_next_kick_global - P[i].Ti_begstep;
	  P[i].dt_step = dt_step;
	  dt_entr = (-ti_step_old / 2 + dt_step / 2) * All.Timebase_interval;

	  /* works only in non-comoving run at the moment */
	  /* WARNING: this velocity correction is inconsistent, as the position of the particle was calculated with a "wrong" velocity before  */
	  for(k = 0; k < 3; k++)
	    {
	      P[i].Vel[k] += P[i].g.GravAccel[k] * dt_entr;
	    }

	  for(k = 0; k < 3; k++)
	    {
	      P[i].Vel[k] += SphP[i].a.HydroAccel[k] * dt_entr;
	    }

#if !defined(EOS_DEGENERATE)
	  SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr;
#else
	  SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr * All.UnitTime_in_s;
#endif

#ifdef NUCLEAR_NETWORK
	  for(k = 0; k < EOS_NSPECIES; k++)
	    {
	      SphP[i].xnuc[k] += SphP[i].dxnuc[k] * dt_entr * All.UnitTime_in_s;
	    }
#endif

	  n++;
	}
    }

  sumup_large_ints(1, &n, &ntot);
  if(ThisTask == 0)
    printf("%d%09d particles woken up.\n", (int) (ntot / 1000000000), (int) (ntot % 1000000000));
#endif

  CPU_Step[CPU_TIMELINE] += measure_time();
}



void do_the_kick(int i, int tstart, int tend, int tcurrent)
{
  int j;
  MyFloat dv[3];
  double minentropy;
  double dt_entr, dt_gravkick, dt_hydrokick, dt_gravkick2, dt_hydrokick2, dt_entr2;

#if defined(TIME_DEP_ART_VISC) || defined(TIME_DEP_MAGN_DISP)
  double dmin1, dmin2;
#endif

#ifdef DISTORTIONTENSORPS
  /* for the distortion 'velocity part', so only the lower two 3x3 submatrices will be != 0 */
  MyDouble dv_distortion_tensorps[6][6];
  int j1, j2;
#endif


  if(All.ComovingIntegrationOn)
    {
      dt_entr = (tend - tstart) * All.Timebase_interval;
      dt_entr2 = (tend - tcurrent) * All.Timebase_interval;
      dt_gravkick = get_gravkick_factor(tstart, tend);
      dt_hydrokick = get_hydrokick_factor(tstart, tend);
      dt_gravkick2 = get_gravkick_factor(tcurrent, tend);
      dt_hydrokick2 = get_hydrokick_factor(tcurrent, tend);
    }
  else
    {
      dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
      dt_gravkick2 = dt_hydrokick2 = dt_entr2 = (tend - tcurrent) * All.Timebase_interval;
    }


  /* do the kick */

  for(j = 0; j < 3; j++)
    {
      dv[j] = P[i].g.GravAccel[j] * dt_gravkick;
      P[i].Vel[j] += dv[j];
    }


#ifdef DISTORTIONTENSORPS
  /* now we do the distortiontensor kick */
  for(j1 = 0; j1 < 3; j1++)
    for(j2 = 0; j2 < 3; j2++)
      {
	dv_distortion_tensorps[j1 + 3][j2] = 0.0;
	dv_distortion_tensorps[j1 + 3][j2 + 3] = 0.0;

	/* the 'acceleration' is given by the product of tidaltensor and distortiontensor */
	for(j = 0; j < 3; j++)
	  {
	    dv_distortion_tensorps[j1 + 3][j2] +=
	      P[i].tidal_tensorps[j1][j] * P[i].distortion_tensorps[j][j2];
	    dv_distortion_tensorps[j1 + 3][j2 + 3] +=
	      P[i].tidal_tensorps[j1][j] * P[i].distortion_tensorps[j][j2 + 3];
	  }
	dv_distortion_tensorps[j1 + 3][j2] *= dt_gravkick;
	dv_distortion_tensorps[j1 + 3][j2 + 3] *= dt_gravkick;

	/* add it to the distortiontensor 'velocities' */
	P[i].distortion_tensorps[j1 + 3][j2] += dv_distortion_tensorps[j1 + 3][j2];
	P[i].distortion_tensorps[j1 + 3][j2 + 3] += dv_distortion_tensorps[j1 + 3][j2 + 3];
      }
#endif

  if(P[i].Type == 0)		/* SPH stuff */
    {
      for(j = 0; j < 3; j++)
	{
//        SphP[i].a.HydroAccel[0] = 0.0;
	  dv[j] += SphP[i].a.HydroAccel[j] * dt_hydrokick;
	  P[i].Vel[j] += SphP[i].a.HydroAccel[j] * dt_hydrokick;

	  SphP[i].VelPred[j] =
	    P[i].Vel[j] - dt_gravkick2 * P[i].g.GravAccel[j] - dt_hydrokick2 * SphP[i].a.HydroAccel[j];
#ifdef PMGRID
	  SphP[i].VelPred[j] += P[i].GravPM[j] * dt_gravkickB;
#endif
#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
	  SphP[i].B[j] += SphP[i].DtB[j] * dt_entr;
	  SphP[i].BPred[j] = SphP[i].B[j] - SphP[i].DtB[j] * dt_entr2;
#endif
	}
#if defined(MAGNETIC) && defined(DIVBCLEANING_DEDNER)
      SphP[i].Phi += SphP[i].DtPhi * dt_entr;
      SphP[i].PhiPred = SphP[i].Phi - SphP[i].DtPhi * dt_entr2;
#endif
#ifdef TIME_DEP_ART_VISC
      SphP[i].alpha += SphP[i].Dtalpha * dt_entr;
      SphP[i].alpha = DMIN(SphP[i].alpha, All.ArtBulkViscConst);
      if(SphP[i].alpha < All.AlphaMin)
	SphP[i].alpha = All.AlphaMin;
#endif
#ifdef VORONOI_TIME_DEP_ART_VISC
      SphP[i].alpha += SphP[i].Dtalpha * dt_entr;
#ifdef VORONOI_RELAX_VIA_VISC
      if(SphP[i].alpha < All.ArtBulkViscConst / 128.0 / 128.0)
	SphP[i].alpha = All.ArtBulkViscConst / 128.0 / 128.0;
#else
      if(SphP[i].alpha < All.ArtBulkViscConst / 128.0)
	SphP[i].alpha = All.ArtBulkViscConst / 128.0;
#endif
//       if(SphP[i].alpha > 16.0*All.ArtBulkViscConst)          SphP[i].alpha = 16.0*All.ArtBulkViscConst;
//       fprintf(stderr,"\n All.ArtBulkViscConst %f SphP[i].alpha %f SphP[i].Dtalpha * dt_entr %f dt_entr %f   \n", All.ArtBulkViscConst, SphP[i].alpha, SphP[i].Dtalpha * dt_entr, dt_entr);}
//       SphP[i].alpha = DMIN(SphP[i].alpha, All.ArtBulkViscConst);
#endif
#ifdef TIME_DEP_MAGN_DISP
      SphP[i].Balpha += SphP[i].DtBalpha * dt_entr;
      SphP[i].Balpha = DMIN(SphP[i].Balpha, All.ArtMagDispConst);
      if(SphP[i].Balpha < All.ArtMagDispMin)
	SphP[i].Balpha = All.ArtMagDispMin;
#endif
      /* In case of cooling, we prevent that the entropy (and
         hence temperature decreases by more than a factor 0.5 */

      if(SphP[i].e.DtEntropy * dt_entr > -0.5 * SphP[i].Entropy)
#if !defined(EOS_DEGENERATE)
	SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr;
#else
	SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr * All.UnitTime_in_s;
#endif
      else
	SphP[i].Entropy *= 0.5;

#ifdef NUCLEAR_NETWORK
      for(j = 0; j < EOS_NSPECIES; j++)
	{
	  SphP[i].xnuc[j] += SphP[i].dxnuc[j] * dt_entr * All.UnitTime_in_s;
	}
#endif

#ifdef CHEMISTRY
      /* update the chemical abundances for the new density and temperature */
      a_start = All.TimeBegin * exp(tstart * All.Timebase_interval);
      a_end = All.TimeBegin * exp(tend * All.Timebase_interval);

      /* time in cosmic expansion parameter */
      ifunc = compute_abundances(mode = 1, i, a_start, a_end);
#endif

      if(All.MinEgySpec)
	{
#ifndef TRADITIONAL_SPH_FORMULATION
	  minentropy = All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
#else
	  minentropy = All.MinEgySpec;
#endif
	  if(SphP[i].Entropy < minentropy)
	    {
	      SphP[i].Entropy = minentropy;
	      SphP[i].e.DtEntropy = 0;
	    }
	}

      /* In case the timestep increases in the new step, we
         make sure that we do not 'overcool'. */
#ifndef WAKEUP
      dt_entr = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) / 2 * All.Timebase_interval;
#else
      dt_entr = P[i].dt_step / 2 * All.Timebase_interval;
#endif
      if(SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr < 0.5 * SphP[i].Entropy)
	SphP[i].e.DtEntropy = -0.5 * SphP[i].Entropy / dt_entr;
    }

  if(All.DoDynamicUpdate)
    force_kick_node(i, dv);
}



/*! This function normally (for flag==0) returns the maximum allowed timestep of a particle, expressed in
 *  terms of the integer mapping that is used to represent the total simulated timespan. The physical
 *  acceleration is returned in aphys. The latter is used in conjunction with the PSEUDOSYMMETRIC integration
 *  option, which also makes of the second function of get_timestep. When it is called with a finite timestep
 *  for flag, it returns the physical acceleration that would lead to this timestep, assuming timestep
 *  criterion 0.
 */
int get_timestep(int p,		/*!< particle index */
		 double *aphys,	/*!< acceleration (physical units) */
		 int flag	/*!< either 0 for normal operation, or finite timestep to get corresponding
				   aphys */ )
{
  double ax, ay, az, ac;
  double csnd = 0, dt = 0, dt_courant = 0;
  int ti_step;

#ifdef BLACK_HOLES
  double dt_accr;

#ifdef UNIFIED_FEEDBACK
  double meddington = 0;
#endif
#endif

#ifdef NS_TIMESTEP
  double dt_NS = 0;
#endif

#ifdef NONEQUILIBRIUM
  double dt_cool, dt_elec;
#endif

#ifdef COSMIC_RAYS
  int CRpop;
#endif

  if(flag == 0)
    {
      ax = fac1 * P[p].g.GravAccel[0];
      ay = fac1 * P[p].g.GravAccel[1];
      az = fac1 * P[p].g.GravAccel[2];

#ifdef PMGRID
      ax += fac1 * P[p].GravPM[0];
      ay += fac1 * P[p].GravPM[1];
      az += fac1 * P[p].GravPM[2];
#endif

      if(P[p].Type == 0)
	{
	  ax += fac2 * SphP[p].a.HydroAccel[0];
	  ay += fac2 * SphP[p].a.HydroAccel[1];
	  az += fac2 * SphP[p].a.HydroAccel[2];
	}

      ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
      *aphys = ac;
    }
  else
    ac = *aphys;

  if(ac == 0)
    ac = 1.0e-30;


  switch (All.TypeOfTimestepCriterion)
    {
    case 0:
      if(flag > 0)
	{
	  dt = flag * All.Timebase_interval;

	  dt /= hubble_a;	/* convert dloga to physical timestep  */

	  ac = 2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / (dt * dt);
	  *aphys = ac;
	  return flag;
	}
      dt = sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac);
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
      if(P[p].Type == 0)
	dt = sqrt(2 * All.ErrTolIntAccuracy * atime * PPP[p].Hsml / 2.8 / ac);
#else
      if(P[p].Type == 0)
	dt =
	  sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] *
	       pow(P[p].Mass / All.ReferenceGasMass, 1.0 / 3) / ac);
#endif
#endif
      break;

    default:
      fprintf(stderr, "\n !!!2@@@!!! \n");
      endrun(888);
      fprintf(stderr, "\n !!!2@@@!!! \n");
      break;
    }


  if(P[p].Type == 0)
    {
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].d.Density);

#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
      double dmax1, dmax2;

      if(All.ComovingIntegrationOn)
	dt_courant = All.CourantFac * All.Time * DMAX(PPP[p].Hsml, All.SofteningTable[0]) / (fac3 * csnd);
      else
	dt_courant = All.CourantFac * DMAX(PPP[p].Hsml, All.SofteningTable[0]) / csnd;

      if(dt_courant > 2.0 * All.CourantFac * SphP[p].MinViscousDt)
	dt_courant = 2.0 * All.CourantFac * SphP[p].MinViscousDt;
#else
      if(All.ComovingIntegrationOn)
      {
	dt_courant = 2.0 * All.CourantFac * All.Time * PPP[p].Hsml / (fac3 * SphP[p].MaxSignalVel);
      }else{
	if(SphP[p].MaxSignalVel<0.1)
	  SphP[p].MaxSignalVel= sqrt(GAMMA*SphP[p].Pressure/SphP[p].d.Density);
	dt_courant = 2.0 * All.CourantFac * PPP[p].Hsml / SphP[p].MaxSignalVel;
      }
#endif

#ifdef VORONOI
#ifdef TWODIMS
      dt_courant = All.CourantFac * sqrt(SphP[p].Volume / M_PI) / SphP[p].MaxSignalVel;
#else
      dt_courant = All.CourantFac * pow(SphP[p].Volume / (4.0 / 3 * M_PI), 1.0 / 3) / SphP[p].MaxSignalVel;
#endif
#endif

  if(P[p].ID == 1341940)
{
printf("xyz=(%g|%g|%g)\n",P[p].Pos[0], P[p].Pos[1], P[p].Pos[2]);
printf("vxvyvz=(%g|%g|%g)\n",P[p].Vel[0],P[p].Vel[1], P[p].Vel[2]);
printf("ent= %g   rho=  %g\n",SphP[p].Entropy, SphP[p].d.Density);
printf("GasAge = %g\n",SphP[p].GasAge);

} 

/*  if(P[p].GasAge > 0.0)
  {
    printf("\nP[p].Type       = %d\n",P[p].Type);
    printf("P[p].GasAge       = %g\n",P[p].GasAge);
    printf("SphP[p].GasAge    = %g\n",SphP[p].GasAge);
    printf("SphP[p].Pressure  = %g\n",SphP[p].Pressure);
    printf("SphP[p].d.Density = %g\n",SphP[p].d.Density);
    printf("SphP[p].Entropy   = %g\n\n",SphP[p].Entropy);
  }
*/

  if(dt_courant < All.MinSizeTimestep)
    {
      printf("\n\nWe might be in trouble...\n");
      printf("SphP[p].MaxSignalVel = %g\n",SphP[p].MaxSignalVel);
      printf("GAMMA             = %g\n",GAMMA);
      printf("SphP[p].Pressure  = %g\n",SphP[p].Pressure);
      printf("SphP[p].d.Density = %g\n",SphP[p].d.Density);
      printf("SphP[p].Entropy   = %g\n",SphP[p].Entropy);
      printf("PPP[p].Hsml       = %g\n",PPP[p].Hsml);
      printf("SphP[p].GasAge    = %g\n",SphP[p].GasAge);
      printf("P[p].StellarAge   = %g\n",P[p].StellarAge);
      printf("P[p].ID           = %d\n",P[p].ID);
      printf("P[p].Vel[0]       = %g\n",P[p].Vel[0]);
      printf("P[p].Vel[1]       = %g\n",P[p].Vel[1]);
      printf("P[p].Vel[2]       = %g\n",P[p].Vel[2]);
      printf("Courant Time      = %g\n\n",2.0*All.CourantFac*PPP[p].Hsml/SphP[p].MaxSignalVel);

      printf("All.AGBGasTemp    = %g\n",All.AGBGasTemp);
      printf("All.AGBGasDensity = %g\n",All.AGBGasDensity);
      printf("All.AGBGasEnergy  = %g\n",All.AGBGasEnergy);
      printf("All.AGBGasEntropy = %g\n",All.AGBGasEntropy);

/*      printf("SphP[p].MaxSignalVel = %g\n",1.0*SphP[p].MaxSignalVel);
      printf("GAMMA = %g\n",1.0*GAMMA);
      printf("SphP[p].Pressure = %g\n",1.0*SphP[p].Pressure);
      printf("SphP[p].d.Density = %g\n",1.0*SphP[p].d.Density);
      printf("PPP[p].Hsml = %g\n",1.0*PPP[p].Hsml);
      printf("Courant Time = %g\n\n",2.0*All.CourantFac*1.0*PPP[p].Hsml/(1.0*SphP[p].MaxSignalVel));
*/
     } 
     if(dt_courant < dt)
	dt = dt_courant;

#ifdef MYFALSE
      dt_viscous = All.CourantFac * SphP[p].MaxViscStep / hubble_a;	/* to convert dloga to physical dt */

      if(dt_viscous < dt)
	dt = dt_viscous;
#endif

#ifdef NS_TIMESTEP
      if(fabs(SphP[p].ViscEntropyChange))
	{
	  dt_NS = VISC_TIMESTEP_PARAMETER * SphP[p].Entropy / SphP[p].ViscEntropyChange / hubble_a;

	  if(dt_NS < dt)
	    dt = dt_NS;
	}
#endif

    }

#ifdef BLACK_HOLES
  if(P[p].Type == 5)
    {
      if(P[p].BH_Mdot > 0 && P[p].BH_Mass > 0)
	{
	  dt_accr = 0.25 * P[p].BH_Mass / P[p].BH_Mdot;
	  if(dt_accr < dt)
	    dt = dt_accr;
	}
    }
#endif

#ifdef BH_BUBBLES
  if(P[p].Type == 5)
    {
      if(P[p].BH_Mdot > 0 && P[p].BH_Mass > 0)
	{
#ifdef UNIFIED_FEEDBACK
	  meddington = (4 * M_PI * GRAVITY * C * PROTONMASS /
			(0.1 * C * C * THOMPSON)) * P[p].BH_Mass * All.UnitTime_in_s;
	  if(P[p].BH_Mdot < All.RadioThreshold * meddington)
#endif
	    dt_accr = (All.BlackHoleRadioTriggeringFactor - 1) * P[p].BH_Mass / P[p].BH_Mdot;
	  if(dt_accr < dt)
	    dt = dt_accr;
	}
    }
#endif

#ifdef NONEQUILIBRIUM
  /* another criterion given by the local cooling time */

  if(P[p].Type == 0)
    {
      dt_cool = fabs(SphP[p].t_cool);	/* still in yrs */
      dt_cool *= SEC_PER_YEAR;	/* in seconds */
      dt_cool /= All.UnitTime_in_s;
      dt_cool *= All.HubbleParam;	/* internal units */

      dt_cool = All.Epsilon * dt_cool;


      if(dt_cool > 0 && dt_cool < dt)
	dt = dt_cool;


      /* yet another criterion given by the electron number density change */

      dt_elec = fabs(SphP[p].t_elec);	/* still in yrs */
      dt_elec *= SEC_PER_YEAR;	/* in seconds */
      dt_elec /= All.UnitTime_in_s;
      dt_elec *= All.HubbleParam;	/* internal units */

      dt_elec = All.Epsilon * dt_elec;

      if(dt_elec > 0 && dt_elec < dt)
	dt = dt_elec;
    }
#endif



  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     hubble_a=1.
   */
  dt *= hubble_a;

#ifdef ONLY_PM
  dt = All.MaxSizeTimestep;
#endif



  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;


  if(dt >= dt_displacement)
    dt = dt_displacement;


#ifdef CONDUCTION
  if(P[p].Type == 0)
    if(dt >= All.MaxSizeConductionStep)
      dt = All.MaxSizeConductionStep;
#endif
#ifdef CR_DIFFUSION
  if(P[p].Type == 0)
    if(dt >= All.CR_DiffusionMaxSizeTimestep)
      dt = All.CR_DiffusionMaxSizeTimestep;
#endif

  if(dt < All.MinSizeTimestep)
    {
#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");

      if(P[p].Type == 0)
	{
	  printf
	    ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g\n",
	     (int) P[p].ID, dt, dt_courant * hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], PPP[p].Hsml,
	     csnd,
	     sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac) * hubble_a,
	     All.SofteningTable[P[p].Type]);

         printf
	    ("Part-Dens=%g \t Part-Pressure=%g \t Part-Ent=%g\n",
		SphP[p].d.Density,SphP[p].Pressure,SphP[p].Entropy);

#ifdef NS_TIMESTEP
	  printf
	    ("Part-ID=%d  dt_NS=%g  A=%g  rho=%g  dotAvisc=%g  dtold=%g, meanpath=%g \n",
	     (int) P[p].ID, dt_NS * hubble_a, SphP[p].Entropy, SphP[p].d.Density,
	     SphP[p].ViscEntropyChange, (P[p].TimeBin ? (1 << P[p].TimeBin) : 0) * All.Timebase_interval,
	     All.IonMeanFreePath *
	     pow((SphP[p].Entropy * pow(SphP[p].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1),
		 2.0) / SphP[p].d.Density);

	  printf("Stressd=(%g|%g|%g) \n", SphP[p].u.s.StressDiag[0], SphP[p].u.s.StressDiag[1],
		 SphP[p].u.s.StressDiag[2]);
	  printf("Stressoffd=(%g|%g|%g) \n", SphP[p].u.s.StressOffDiag[0], SphP[p].u.s.StressOffDiag[1],
		 SphP[p].u.s.StressOffDiag[2]);
#endif


	}
      else
	{
	  printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
		 P[p].Pos[2]);
	}
      fflush(stdout);
      fprintf(stderr, "\n @ fflush \n");
      endrun(888);
#endif
      dt = All.MinSizeTimestep;
    }

  ti_step = (int) (dt / All.Timebase_interval);

#ifdef CHEMISTRY
  if(ti_step == 0)
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%d dt=%g dt_elec=%g dt_cool=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g)\n\n",
	     ThisTask, P[p].ID, dt, SphP[p].t_elec, SphP[p].t_cool, All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2]);
      fflush(stdout);
      endrun(818);
    }
#endif


  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%d dt=%g dtc=%g dtv=%g dtdis=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) tree=(%g|%g%g)\n\n",
	     ThisTask, (int) P[p].ID, dt, dt_courant, dt, dt_displacement,
	     All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].g.GravAccel[0], P[p].g.GravAccel[1],
	     P[p].g.GravAccel[2]);
#ifdef PMGRID
      printf("pm_force=(%g|%g|%g)\n", P[p].GravPM[0], P[p].GravPM[1], P[p].GravPM[2]);
#endif
      if(P[p].Type == 0)
	printf("hydro-frc=(%g|%g|%g) dens=%g hsml=%g\n", SphP[p].a.HydroAccel[0], SphP[p].a.HydroAccel[1],
	       SphP[p].a.HydroAccel[2], SphP[p].d.Density, PPP[p].Hsml);

#ifdef COSMIC_RAYS
      if(P[p].Type == 0)
	for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	  printf("Cosmic Ray Properties: C0: %g -- q0  : %g -- P  : %g\n"
		 "                       Rho: %g\n\n",
		 SphP[p].CR_C0[CRpop], SphP[p].CR_q0[CRpop], CR_Particle_Pressure(SphP + p, CRpop),
		 SphP[p].d.Density);
#endif

      fflush(stdout);
      endrun(818);
    }

  return ti_step;
}


/*! This function computes an upper limit ('dt_displacement') to the global timestep of the system based on
 *  the rms velocities of particles. For cosmological simulations, the criterion used is that the rms
 *  displacement should be at most a fraction MaxRMSDisplacementFac of the mean particle separation. Note that
 *  the latter is estimated using the assigned particle masses, separately for each particle type. If comoving
 *  integration is not used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
  int i, type;
  int count[6];
  long long count_sum[6];
  double v[6], v_sum[6], mim[6], min_mass[6];
  double dt, dmean, asmth = 0;

  dt_displacement = All.MaxSizeTimestep;

  if(All.ComovingIntegrationOn)
    {
      for(type = 0; type < 6; type++)
	{
	  count[type] = 0;
	  v[type] = 0;
	  mim[type] = 1.0e30;
	}

      for(i = 0; i < NumPart; i++)
	{
	  v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
	  if(P[i].Mass > 0)
	    {
	      if(mim[P[i].Type] > P[i].Mass)
		mim[P[i].Type] = P[i].Mass;
	    }
	  count[P[i].Type]++;
	}

      MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      sumup_large_ints(6, count, count_sum);

#ifdef SFR
      /* add star and gas particles together to treat them on equal footing, using the original gas particle
         spacing. */
      v_sum[0] += v_sum[4];
      count_sum[0] += count_sum[4];
      v_sum[4] = v_sum[0];
      count_sum[4] = count_sum[0];
#ifdef BLACK_HOLES
      v_sum[0] += v_sum[5];
      count_sum[0] += count_sum[5];
      v_sum[5] = v_sum[0];
      count_sum[5] = count_sum[0];
      min_mass[5] = min_mass[0];
#endif
#endif

      for(type = 0; type < 6; type++)
	{
	  if(count_sum[type] > 0)
	    {
	      if(type == 0 || (type == 4 && All.StarformationOn))
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
	      else
		dmean =
		  pow(min_mass[type] /
		      ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);

#ifdef BLACK_HOLES
	      if(type == 5)
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
#endif
	      dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

#ifdef PMGRID
	      asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
	      if(((1 << type) & (PLACEHIGHRESREGION)))
		asmth = All.Asmth[1];
#endif
	      if(asmth < dmean)
		dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
#endif

	      if(ThisTask == 0)
		printf("type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
		       type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);


#ifdef NEUTRINOS
	      if(type != 2)	/* don't constrain the step to the neutrinos */
#endif
		if(dt < dt_displacement)
		  dt_displacement = dt;
	    }
	}

      if(ThisTask == 0)
	printf("displacement time constraint: %g  (%g)\n", dt_displacement, All.MaxSizeTimestep);
    }
}

int get_timestep_bin(int ti_step)
{
  int bin = -1;

  if(ti_step == 0)
    return 0;

  if(ti_step == 1)
    {
      printf("time-step of integer size 1 not allowed\n");
      endrun(112313);
    }

  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  return bin;
}
