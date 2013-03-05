#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef CS_MODEL
#include "cs_metals.h"
#endif


/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */


/*! This routine computes the accelerations for all active particles.  First, the gravitational forces are
 * computed. This also reconstructs the tree, if needed, otherwise the drift/kick operations have updated the
 * tree to make it fullu usable at the current time.
 *
 * If gas particles are presented, the `interior' of the local domain is determined. This region is guaranteed
 * to contain only particles local to the processor. This information will be used to reduce communication in
 * the hydro part.  The density for active SPH particles is computed next. If the number of neighbours should
 * be outside the allowed bounds, it will be readjusted by the function ensure_neighbours(), and for those
 * particle, the densities are recomputed accordingly. Finally, the hydrodynamical forces are added.
 */
void compute_accelerations(int mode)
{
#ifdef RADTRANSFER
  double timeeach = 0, timeall = 0, tstart = 0, tend = 0;
#endif
#if defined(BUBBLES) || defined(MULTI_BUBBLES)
  double hubble_a;
#endif

  if(ThisTask == 0)
    {
      printf("Start force computation...\n");
      fflush(stdout);
    }

#ifdef REIONIZATION
  heating();
#endif

  CPU_Step[CPU_MISC] += measure_time();

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      long_range_force();

      CPU_Step[CPU_MESH] += measure_time();

    }
#endif



#ifndef ONLY_PM

  gravity_tree();		/* computes gravity accel. */

  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
    gravity_tree();		/* For the first timestep, we redo it
				 * to allow usage of relative opening
				 * criterion for consistent accuracy.
				 */
#endif


  if(All.Ti_Current == 0 && RestartFlag == 0 && header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
    second_order_ics();		/* produces the actual ICs from the special second order IC file */


#ifdef FORCETEST
  gravity_forcetest();
#endif


#ifdef GASRETURN
  stochastic_gas_return();
//  MPI_Barrier(MPI_COMM_WORLD);
  rearrange_particle_sequence();
//  MPI_Barrier(MPI_COMM_WORLD);
#endif


  if(All.TotN_gas > 0)
    {
      /***** density *****/
/*      MPI_Barrier(MPI_COMM_WORLD);  */
      if(ThisTask == 0)
	{
	  printf("Start density computation...\n");
	  fflush(stdout);
	}

#ifdef CS_MODEL
      CPU_Step[CPU_MISC] += measure_time();

#if defined(CS_SNI) || defined(CS_SNII)
      cs_flag_SN_starparticles();	/* mark SNI star particles */
#endif
      cs_copy_densities();
      cs_find_low_density_tail();

      CPU_Step[CPU_CSMISC] += measure_time();
#endif

#ifndef VORONOI
      density();		/* computes density, and pressure */
#else
      voronoi_mesh();
      voronoi_density();
#endif

#if (defined(SMOOTH_PHI) || defined(SMOOTH_ROTB) || defined(BSMOOTH))
      smoothed_values();
#endif


#if defined(SNIA_HEATING)
      snIa_heating();
#endif


#if defined(CS_MODEL) && defined(CS_FEEDBACK)

      CPU_Step[CPU_CSMISC] += measure_time();

      cs_find_hot_neighbours();

      cs_promotion();
      cs_copy_densities();
      CPU_Step[CPU_CSMISC] += measure_time();
      density();		/* recalculate densities again */
      CPU_Step[CPU_CSMISC] += measure_time();
#endif


      /***** update smoothing lengths in tree *****/
#ifndef VORONOI
      force_update_hmax();
#endif

      /***** hydro forces *****/
      if(ThisTask == 0)
	{
	  printf("Start hydro-force computation...\n");
	  fflush(stdout);
	}

#ifndef VORONOI
      hydro_force();		/* adds hydrodynamical accelerations  and computes du/dt  */
#else
      voronoi_hydro_force();
#endif

#ifdef CONDUCTION
      if(All.Conduction_Ti_endstep == All.Ti_Current)
	conduction();
#endif

#ifdef CR_DIFFUSION
      if(All.CR_Diffusion_Ti_endstep == All.Ti_Current)
	cosmic_ray_diffusion();
#endif


#ifdef RADTRANSFER
      if(Flag_FullStep)		/* only do it for full timesteps */
	{
	  All.Radiation_Ti_endstep = All.Ti_Current;
	      

	  if(ThisTask == 0)
	    {
	      printf("Start Eddington tensor computation...\n");
	      fflush(stdout);
	    }
	  
	  eddington();
	  
	  if(ThisTask == 0)
	    {
	      printf("done Eddington tensor! \n");
	      fflush(stdout);
	    }
	  
	  star_density();
	  

	  /***** set simple initial conditions *****/
	  if(All.Time == All.TimeBegin)
	    {
	      if(ThisTask == 0)
		{
		  printf("Setting simple inits...\n");
		  fflush(stdout);
		}

	      radtransfer_set_simple_inits();

	      if(ThisTask == 0)
		{
		  printf("done with simple inits! \n");
		  fflush(stdout);
		}
	    }

	  /***** evolve the transport of radiation *****/
	  if(ThisTask == 0)
	    {
	      printf("start radtransfer...\n");
	      fflush(stdout);
	    }

	  tstart = second();

	  radtransfer();

	  tend = second();
	  timeeach = timediff(tstart, tend);
	  MPI_Allreduce(&timeeach, &timeall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	  if(ThisTask == 0)
	    {
	      printf("time consumed is %g \n", timeall);
	      printf("done with radtransfer! \n");
	      fflush(stdout);
	    }

	  All.Radiation_Ti_begstep = All.Radiation_Ti_endstep;
	}
#endif

#ifdef MHM
      /***** kinetic feedback *****/
      kinetic_feedback_mhm();
#endif


#ifdef BLACK_HOLES
      /***** black hole accretion and feedback *****/
      blackhole_accretion();
#ifdef FOF
      /* this will find new black hole seed halos */
      if(All.Time >= All.TimeNextBlackHoleCheck)
	{
	  fof_fof(-1);

	  if(All.ComovingIntegrationOn)
	    All.TimeNextBlackHoleCheck *= All.TimeBetBlackHoleSearch;
	  else
	    All.TimeNextBlackHoleCheck += All.TimeBetBlackHoleSearch;
	}
#endif
#endif


#ifdef COOLING	/**** radiative cooling and star formation *****/

#ifdef CS_MODEL
      cs_cooling_and_starformation();
#else

#ifdef SFR
      cooling_and_starformation();
#else
      cooling_only();
#endif

#endif
      CPU_Step[CPU_COOLINGSFR] += measure_time();
#endif /*ends COOLING */


#if defined(CS_MODEL) && defined(CS_ENRICH)
#ifndef CS_FEEDBACK
      Flag_phase = 0;		/* no destinction between phases */

      cs_update_weights();
      CPU_Step[CPU_WEIGHTS_HOT] += measure_time();
      cs_enrichment();
      CPU_Step[CPU_ENRICH_HOT] += measure_time();
#else

      Flag_phase = 1;		/* COLD phase Flag_phase  = 1 */

      cs_update_weights();
      CPU_Step[CPU_WEIGHTS_HOT] += measure_time();
      cs_enrichment();
      CPU_Step[CPU_ENRICH_HOT] += measure_time();

      Flag_phase = 2;		/* HOT phase Flag_phase = 2 */

      cs_update_weights();
      CPU_Step[CPU_WEIGHTS_COLD] += measure_time();
      cs_enrichment();
      CPU_Step[CPU_ENRICH_COLD] += measure_time();

      Flag_phase = 0;
#endif
#endif

#ifdef CS_TESTS
      cs_energy_test();
#endif



#ifndef BH_BUBBLES
#ifdef BUBBLES
      /**** bubble feedback *****/
      if(All.Time >= All.TimeOfNextBubble)
	{
#ifdef FOFs
	  fof_fof(-1);
	  bubble();
#else
	  bubble();
#endif
	  if(All.ComovingIntegrationOn)
	    {
	      hubble_a = hubble_function(All.Time);
	      All.TimeOfNextBubble *= (1.0 + All.BubbleTimeInterval * hubble_a);
	    }
	  else
	    All.TimeOfNextBubble += All.BubbleTimeInterval / All.UnitTime_in_Megayears;

	  if(ThisTask == 0)
	    printf("Time of the bubble generation: %g\n", 1. / All.TimeOfNextBubble - 1.);
	}
#endif
#endif

#if defined(MULTI_BUBBLES) && defined(FOF)
      if(All.Time >= All.TimeOfNextBubble)
	{
	  fof_fof(-1);

	  if(All.ComovingIntegrationOn)
	    {
	      hubble_a = hubble_func(All.Time);
	      All.TimeOfNextBubble *= (1.0 + All.BubbleTimeInterval * hubble_a);
	    }
	  else
	    All.TimeOfNextBubble += All.BubbleTimeInterval / All.UnitTime_in_Megayears;

	  if(ThisTask == 0)
	    printf("Time of the bubble generation: %g\n", 1. / All.TimeOfNextBubble - 1.);
	}
#endif

    }



  if(ThisTask == 0)
    {
      printf("force computation done.\n");
      fflush(stdout);
    }
}
