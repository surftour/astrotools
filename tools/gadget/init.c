#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif
#ifdef CS_MODEL
#include "cs_metals.h"
#endif

#include "voronoi.h"

#ifdef MACHNUM
#ifdef COSMIC_RAYS
#define h  All.HubbleParam
#define cm (h/All.UnitLength_in_cm)
#define s  (h/All.UnitTime_in_s)
#define LightSpeed (2.9979e10*cm/s)
#define c2   ( LightSpeed * LightSpeed )
#endif
#endif



/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the inial SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3, atime;

#ifdef COSMIC_RAYS
  int CRpop;
#endif

#if defined(COSMIC_RAYS) && defined(MACHNUM)
  double Pth1, PCR1[NUMCRPOP], rBeta[NUMCRPOP], C_phys[NUMCRPOP], q_phys[NUMCRPOP];
#endif
#ifdef CR_INITPRESSURE
  double cr_pressure, q_phys, C_phys[NUMCRPOP];
#endif
#ifdef BLACK_HOLES
  int count_holes = 0;
#endif
#ifdef CHEMISTRY
  int ifunc;
  double min_t_cool, max_t_cool;
  double min_t_elec, max_t_elec;
  double a_start, a_end;
#endif

#ifdef EULERPOTENTIALS
  double a0, a1, a2;
  double b0, b1, b2;
#endif

#ifdef START_WITH_EXTRA_NGBDEV
  double MaxNumNgbDeviationMerk;
#endif

#ifdef DISTORTIONTENSORPS
  int i1, i2;
#endif

  All.Time = All.TimeBegin;

  if(RestartFlag == 3 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf("Need to give the snapshot number if FOF/SUBFIND is selected for output\n");
      endrun(0);
    }

  if(RestartFlag == 4 && RestartSnapNum < 0)
    {
      if(ThisTask == 0)
	printf("Need to give the snapshot number if snapshot should be converted\n");
      endrun(0);
    }

  switch (All.ICFormat)
    {
    case 1:
    case 2:
    case 3:
    case 4:
      if(RestartFlag >= 2 && RestartSnapNum >= 0)
	{
	  char fname[1000];

	  if(All.NumFilesPerSnapshot > 1)
	    sprintf(fname, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase,
		    RestartSnapNum);
	  else
	    sprintf(fname, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);
	  read_ic(fname);
	  if(RestartFlag == 4)
	    {
	      sprintf(All.SnapshotFileBase, "%s_converted", All.SnapshotFileBase);
	      if(ThisTask == 0)
		printf("Start writing file %s\n", All.SnapshotFileBase);
	      printf("RestartSnapNum %d\n", RestartSnapNum);
	      savepositions(RestartSnapNum);
	      endrun(0);
	    }

	}
      else
	{
	  if(All.ICFormat == 4)
	    read_ic_cluster(All.InitCondFile);
	  else
	    read_ic(All.InitCondFile);
	}
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;

#ifdef COOLING
#ifdef CS_MODEL
  if(RestartFlag == 0)
    XH = HYDROGEN_MASSFRAC;
#endif
  IonizeParams();
#endif

#ifdef CHEMISTRY
  InitChem();
#endif

#ifdef GASRETURN
  ReadMassReturnParams("MASSRETURN");
#endif

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = All.Time * All.Time * All.Time;
      atime = All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = 1;
      atime = 1;
    }

#ifdef RADTRANSFER
  All.Radiation_Ti_begstep = 0;
#endif


  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    {
      if(RestartSnapNum < 0)
	All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;
      else
	All.SnapshotFileCount = RestartSnapNum + 1;
    }

#ifdef OUTPUTLINEOFSIGHT
  All.Ti_nextlineofsight = (int) (log(All.TimeFirstLineOfSight / All.TimeBegin) / All.Timebase_interval);
  if(RestartFlag == 2)
    endrun(78787);
#endif

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;
#if defined(MAGNETIC) && defined(BSMOOTH)
#ifdef SETMAINTIMESTEPCOUNT
  All.MainTimestepCounts = All.MainTimestepCountIni;
#else
  All.MainTimestepCounts = 0;
#endif
#endif

  All.TopNodeAllocFactor = 0.008;
  All.TreeAllocFactor = 0.7;

  All.Cadj_Cost = 1.0e-30;
  All.Cadj_Cpu = 1.0e-3;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#ifdef BLACK_HOLES
  All.TimeNextBlackHoleCheck = All.TimeBegin;
#endif



#ifdef BUBBLES
  if(All.ComovingIntegrationOn)
    All.TimeOfNextBubble = 1. / (1. + All.FirstBubbleRedshift);
  else
    All.TimeOfNextBubble = All.TimeBegin + All.BubbleTimeInterval / All.UnitTime_in_Megayears;
  if(ThisTask == 0)
    printf("Initial time: %g and first bubble time %g \n", All.TimeBegin, All.TimeOfNextBubble);

  if(RestartFlag == 2 && All.TimeBegin > All.TimeOfNextBubble)
    {
      printf("Restarting from the snapshot file with the wrong FirstBubbleRedshift! \n");
      endrun(0);
    }
#endif

#ifdef MULTI_BUBBLES
  if(All.ComovingIntegrationOn)
    All.TimeOfNextBubble = 1. / (1. + All.FirstBubbleRedshift);
  else
    All.TimeOfNextBubble = All.TimeBegin + All.BubbleTimeInterval / All.UnitTime_in_Megayears;
  if(ThisTask == 0)
    printf("Initial time: %g and time of the first bubbles %g \n", All.TimeBegin, All.TimeOfNextBubble);
  if(RestartFlag == 2 && All.TimeBegin > All.TimeOfNextBubble)
    {
      printf("Restarting from the snapshot file with the wrong FirstBubbleRedshift! \n");
      endrun(0);
    }
#endif

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }

#ifdef REINIT_AT_TURNAROUND
  double turnaround_radius_local = 0.0, turnaround_radius_global = 0.0, v_part, r_part;

  /* get local turnaroundradius */
  for(i = 0; i < NumPart; i++)
    {
      r_part = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);
      v_part = (P[i].Pos[0] * P[i].Vel[0] + P[i].Pos[1] * P[i].Vel[1] + P[i].Pos[2] * P[i].Vel[2]);
      if((v_part < 0.0) && (r_part > turnaround_radius_local))
	turnaround_radius_local = r_part;
    }

  /* find global turnaround radius by taking maximum of all CPUs */
  MPI_Allreduce(&turnaround_radius_local, &turnaround_radius_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  All.CurrentTurnaroundRadius = turnaround_radius_global;
  All.InitialTurnaroundRadius = turnaround_radius_global;

  if(ThisTask == 0)
    {
      printf("REINIT_AT_TURNAROUND: current (initial) turnaround radius = %g\n", All.CurrentTurnaroundRadius);
      fflush(stdout);
    }
#endif


  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].g.GravAccel[j] = 0;

/* here comes all the DISTORTIONTENSORPS setup */
#ifdef DISTORTIONTENSORPS
      /*init tidal tensor for first output (this is not used for any calculation) */
      for(i1 = 0; i1 < 3; i1++)
	for(i2 = 0; i2 < 3; i2++)
	  P[i].tidal_tensorps[i1][i2] = 0.0;

      /* find caustics on the fly by sign analysis of configuration space distortion */
      P[i].last_stream_determinant = 1.0;

      /* distortion tensor was read in from ICs, so get correct last_stream_determinant */
#ifdef DISTORTION_READALL
      MyDouble product_matrix[3][3];

      product_matrix[0][0] = P[i].distortion_tensorps[0][0] +
	P[i].distortion_tensorps[0][3] * P[i].V_matrix[0][0] +
	P[i].distortion_tensorps[0][4] * P[i].V_matrix[1][0] +
	P[i].distortion_tensorps[0][5] * P[i].V_matrix[2][0];
      product_matrix[0][1] = P[i].distortion_tensorps[0][1] +
	P[i].distortion_tensorps[0][3] * P[i].V_matrix[0][1] +
	P[i].distortion_tensorps[0][4] * P[i].V_matrix[1][1] +
	P[i].distortion_tensorps[0][5] * P[i].V_matrix[2][1];
      product_matrix[0][2] = P[i].distortion_tensorps[0][2] +
	P[i].distortion_tensorps[0][3] * P[i].V_matrix[0][2] +
	P[i].distortion_tensorps[0][4] * P[i].V_matrix[1][2] +
	P[i].distortion_tensorps[0][5] * P[i].V_matrix[2][2];
      product_matrix[1][0] = P[i].distortion_tensorps[1][0] +
	P[i].distortion_tensorps[1][3] * P[i].V_matrix[0][0] +
	P[i].distortion_tensorps[1][4] * P[i].V_matrix[1][0] +
	P[i].distortion_tensorps[1][5] * P[i].V_matrix[2][0];
      product_matrix[1][1] = P[i].distortion_tensorps[1][1] +
	P[i].distortion_tensorps[1][3] * P[i].V_matrix[0][1] +
	P[i].distortion_tensorps[1][4] * P[i].V_matrix[1][1] +
	P[i].distortion_tensorps[1][5] * P[i].V_matrix[2][1];
      product_matrix[1][2] = P[i].distortion_tensorps[1][2] +
	P[i].distortion_tensorps[1][3] * P[i].V_matrix[0][2] +
	P[i].distortion_tensorps[1][4] * P[i].V_matrix[1][2] +
	P[i].distortion_tensorps[1][5] * P[i].V_matrix[2][2];
      product_matrix[2][0] = P[i].distortion_tensorps[2][0] +
	P[i].distortion_tensorps[2][3] * P[i].V_matrix[0][0] +
	P[i].distortion_tensorps[2][4] * P[i].V_matrix[1][0] +
	P[i].distortion_tensorps[2][5] * P[i].V_matrix[2][0];
      product_matrix[2][1] = P[i].distortion_tensorps[2][1] +
	P[i].distortion_tensorps[2][3] * P[i].V_matrix[0][1] +
	P[i].distortion_tensorps[2][4] * P[i].V_matrix[1][1] +
	P[i].distortion_tensorps[2][5] * P[i].V_matrix[2][1];
      product_matrix[2][2] = P[i].distortion_tensorps[2][2] +
	P[i].distortion_tensorps[2][3] * P[i].V_matrix[0][2] +
	P[i].distortion_tensorps[2][4] * P[i].V_matrix[1][2] +
	P[i].distortion_tensorps[2][5] * P[i].V_matrix[2][2];

      /* this determinant will change sign when we pass through a caustic -> criterion for caustics */
      P[i].last_stream_determinant =
	((product_matrix[0][0]) * (product_matrix[1][1]) * (product_matrix[2][2]) +
	 (product_matrix[0][1]) * (product_matrix[1][2]) * (product_matrix[2][0]) +
	 (product_matrix[0][2]) * (product_matrix[1][0]) * (product_matrix[2][1]) -
	 (product_matrix[0][2]) * (product_matrix[1][1]) * (product_matrix[2][0]) -
	 (product_matrix[0][0]) * (product_matrix[1][2]) * (product_matrix[2][1]) -
	 (product_matrix[0][1]) * (product_matrix[1][0]) * (product_matrix[2][2]));


#endif

#ifdef REINIT_AT_TURNAROUND
      All.SIM_epsilon = REINIT_AT_TURNAROUND;
      /* no phase-space analysis for particles that are initially within the turnaround radius */
      if(sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]) <
	 All.CurrentTurnaroundRadius)
#ifdef DISTORTION_READALL
	P[i].turnaround_flag = 1;
#else
	P[i].turnaround_flag = -1;
#endif
      else
	P[i].turnaround_flag = 0;
#endif

#ifdef OUTPUT_LAST_CAUSTIC
      /* all entries zero -> no caustic yet */
      P[i].lc_Time = 0.0;
      P[i].lc_Pos[0] = 0.0;
      P[i].lc_Pos[1] = 0.0;
      P[i].lc_Pos[2] = 0.0;
      P[i].lc_Vel[0] = 0.0;
      P[i].lc_Vel[1] = 0.0;
      P[i].lc_Vel[2] = 0.0;
      P[i].lc_rho_normed_cutoff = 0.0;

      P[i].lc_Dir_x[0] = 0.0;
      P[i].lc_Dir_x[1] = 0.0;
      P[i].lc_Dir_x[2] = 0.0;
      P[i].lc_Dir_y[0] = 0.0;
      P[i].lc_Dir_y[1] = 0.0;
      P[i].lc_Dir_y[2] = 0.0;
      P[i].lc_Dir_z[0] = 0.0;
      P[i].lc_Dir_z[1] = 0.0;
      P[i].lc_Dir_z[2] = 0.0;

      P[i].lc_smear_x = 0.0;
      P[i].lc_smear_y = 0.0;
      P[i].lc_smear_z = 0.0;
#endif


#ifdef PMGRID
      /* long range tidal field init */
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

#ifndef DISTORTION_READALL
      /* 
         6D phase space distortion tensor in the beginning.
         Both in physical and comoving coordinates this equals unity.
       */
      for(i1 = 0; i1 < 6; i1++)
	for(i2 = 0; i2 < 6; i2++)
	  {
	    if((i1 == i2))
	      {
		P[i].distortion_tensorps[i1][i2] = 1.0;
	      }
	    else
	      {
		P[i].distortion_tensorps[i1][i2] = 0.0;
	      }
	  }
#endif

#ifdef COSMIC_DISTORTION
      /* for cosmological simulations we do init here, not read from ICs */

      /* no caustic passages in the beginning */
      P[i].caustic_counter = 0.0;

      /* 
         Approximation: perfect Hubble Flow -> peculiar sheet orientation is exactly zero. 
         the stream density is then given by the upper left 3x3 submatrix
         determinant of the phase-space distortiontensor.
       */
      for(i1 = 0; i1 < 3; i1++)
	for(i2 = 0; i2 < 3; i2++)
	  P[i].V_matrix[i1][i2] = 0.0;

      /* 
         Approximation: initial sream density equals background density. 
         this is comoving density. the scalefactor is taken into
         account when handling the normed stream density, so that
         initial stream density * normed stream density gives the
         correct physical stream density at every time.
       */
      P[i].init_density = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
#endif

      /* annihilation stuff */

      P[i].s_1_last = 1.0;
      P[i].s_2_last = 1.0;
      P[i].s_3_last = 1.0;
      P[i].second_deriv_last = 0.0;
      P[i].rho_normed_cutoff_last = 1.0;

      P[i].s_1_current = 1.0;
      P[i].s_2_current = 1.0;
      P[i].s_3_current = 1.0;
      P[i].second_deriv_current = 0.0;
      P[i].rho_normed_cutoff_current = 1.0;

      P[i].annihilation = 0.0;
#ifdef COSMIC_DISTORTION
      P[i].stream_density = P[i].init_density / (All.TimeBegin * All.TimeBegin * All.TimeBegin);
#else
      P[i].stream_density = P[i].init_density;
#endif
      P[i].analytic_caustics = 0.0;
#endif /* DISTORTIONTENSORPS */

#ifdef KEEP_DM_HSML_AS_GUESS
      if(RestartFlag != 1)
	P[i].DM_Hsml = -1;
#endif

#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_begstep = 0;
      P[i].Ti_current = 0;
      P[i].TimeBin = 0;

      if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
	P[i].OldAcc = 0;	/* Do not zero in 2lpt case as masses are stored here */

      P[i].GravCost = 1;
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)
      P[i].p.Potential = 0;
#endif
#ifdef STELLARAGE
      if(RestartFlag == 0 && (P[i].Type < 2 || P[i].Type > 3))
	P[i].StellarAge = 0;
#endif

#ifdef GASRETURN
      //if(RestartFlag == 0 && (P[i].Type < 2 || P[i].Type > 3)) - TJ says it maybe should be this instead
      if(RestartFlag == 0)
	{
//	  P[i].StellarAge = - 5.0*get_random_number(P[i].ID + 3);
	  P[i].StellarAge = - 5.0; 
	  P[i].GasAge 	  = - 5.0;
	}
#endif 	

#ifdef METALS
#ifndef CS_MODEL
      if(RestartFlag == 0 && (P[i].Type == 1 || P[i].Type > 3))
	P[i].Metallicity = 0;
#else
      if(RestartFlag == 0)
	{
	  for(j = 0; j < 12; j++)
	    {
	      P[i].Zm[j] = 0;
	      P[i].ZmReservoir[j] = 0;
	    }

	  P[i].Zm[6] = HYDROGEN_MASSFRAC * (P[i].Mass);
	  P[i].Zm[0] = (1 - HYDROGEN_MASSFRAC) * (P[i].Mass);
	  /*     Order of chemical elements:   He, Carbon,Mg,O,Fe,Si,H,N,Ne,S,Ca,Zn */

	  PPP[i].Hsml = 0;
#ifdef CS_FEEDBACK
	  PPP[i].EnergySN = 0;
	  PPP[i].EnergySNCold = 0;
#endif
	}
#endif

#endif

#ifdef BLACK_HOLES
      if(P[i].Type == 5)
	{
	  count_holes++;

	  if(RestartFlag == 0)
	    P[i].BH_Mass = All.SeedBlackHoleMass;
#ifdef BH_BUBBLES
	  if(RestartFlag == 0)
	    {
	      P[i].BH_Mass_bubbles = All.SeedBlackHoleMass;
	      P[i].BH_Mass_ini = All.SeedBlackHoleMass;
#ifdef UNIFIED_FEEDBACK
	      P[i].BH_Mass_radio = All.SeedBlackHoleMass;
#endif
	    }
#endif
	}
#endif
    }

#ifdef BLACK_HOLES
  MPI_Allreduce(&count_holes, &All.TotBHs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  for(i = 0; i < TIMEBINS; i++)
    TimeBinActive[i] = 1;

  reconstruct_timebins();

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

#ifdef CONDUCTION
  All.Conduction_Ti_endstep = All.Conduction_Ti_begstep = 0;
#endif
#ifdef CR_DIFFUSION
  All.CR_Diffusion_Ti_endstep = All.CR_Diffusion_Ti_begstep = 0;
#endif

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].a.HydroAccel[j] = 0;
	}

      SphP[i].e.DtEntropy = 0;

#ifdef CHEMISTRY
      SphP[i].Gamma = GAMMA;	/* set universal value */
      SphP[i].t_cool = 0;
      SphP[i].t_elec = 0;
#endif

#ifdef GASRETURN
      SphP[i].GasAge = -1.0;
#endif

#ifdef CS_MODEL
      SphP[i].DensityOld = 0;
#endif
#ifdef CS_FEEDBACK
      SphP[i].da.DensityAvg = 0;
      SphP[i].ea.EntropyAvg = 0;
      SphP[i].DensPromotion = 0;
      SphP[i].TempPromotion = 0;
      SphP[i].HotHsml = 0;
#endif

      if(RestartFlag == 0)
	{
#ifndef READ_HSML
	  PPP[i].Hsml = 0;
#endif
	  SphP[i].d.Density = -1;
#ifdef VOLUME_CORRECTION
	  SphP[i].DensityOld = 1;
#endif
#ifdef COOLING
	  SphP[i].Ne = 1.0;
#endif
	  SphP[i].v.DivVel = 0;
	}
#ifdef WINDS
      SphP[i].DelayTime = 0;
#endif
#ifdef SFR
      SphP[i].Sfr = 0;
#endif
#ifdef MAGNETIC
#ifdef BINISET
      SphP[i].BPred[0] = All.BiniX;
      SphP[i].BPred[1] = All.BiniY;
      SphP[i].BPred[2] = All.BiniZ;
#ifdef EULERPOTENTIALS
      if(All.BiniY == 0 && All.BiniZ == 0)
	{
	  a0 = 0;
	  a1 = All.BiniX;
	  a2 = 0;
	  b0 = 0;
	  b1 = 1;
	  b2 = 1;
	}
      else
	{
	  if(All.BiniX != 0 && (All.BiniY != 0 || All.BiniZ != 0))
	    {
	      b0 = -(All.BiniZ + All.BiniY) / All.BiniX;
	      b1 = 1;
	      b2 = 1;
	      a0 = 0;
	      a1 = -All.BiniZ / b0;
	      a2 = All.BiniY / b0;
	    }
	  else
	    {
	      a0 = a1 = a2 = b0 = b1 = b2 = 0;
	      if(ThisTask == 0)
		{
		  printf("Can not reconstruct Euler potentials from Bini values !\n");
		  endrun(6723);
		}
	    }
	}
      SphP[i].EulerA = (a0 * P[i].Pos[0] + a1 * P[i].Pos[1] + a2 * P[i].Pos[2]) * atime * All.HubbleParam;
      SphP[i].EulerB = (b0 * P[i].Pos[0] + b1 * P[i].Pos[1] + b2 * P[i].Pos[2]) * atime * All.HubbleParam;
#endif
#endif
#ifndef EULERPOTENTIALS
      for(j = 0; j < 3; j++)
	{
	  SphP[i].DtB[j] = 0;
	  SphP[i].B[j] = SphP[i].BPred[j];
	}
#endif
#ifdef TIME_DEP_MAGN_DISP
#ifdef HIGH_MAGN_DISP_START
      SphP[i].Balpha = All.ArtMagDispConst;
#else
      SphP[i].Balpha = All.ArtMagDispMin;
#endif
      SphP[i].DtBalpha = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
      SphP[i].Phi = SphP[i].PhiPred = SphP[i].DtPhi = 0;
#endif
#endif

#ifdef TIME_DEP_ART_VISC
#ifdef HIGH_ART_VISC_START
      if(HIGH_ART_VISC_START == 0)
	SphP[i].alpha = All.ArtBulkViscConst;
      if(HIGH_ART_VISC_START > 0)
	if(P[i].Pos[0] > HIGH_ART_VISC_START)
	  SphP[i].alpha = All.ArtBulkViscConst;
	else
	  SphP[i].alpha = All.AlphaMin;
      if(HIGH_ART_VISC_START < 0)
	if(P[i].Pos[0] < -HIGH_ART_VISC_START)
	  SphP[i].alpha = All.ArtBulkViscConst;
	else
	  SphP[i].alpha = All.AlphaMin;
#else
      SphP[i].alpha = All.AlphaMin;
#endif
      SphP[i].Dtalpha = 0.0;
#endif

#ifdef VORONOI_TIME_DEP_ART_VISC
      SphP[i].Dtalpha = 0.0;
      SphP[i].alpha = All.ArtBulkViscConst / 128.0;
#ifdef VORONOI_RELAX_VIA_VISC
      SphP[i].alpha = All.ArtBulkViscConst;
#endif
#endif
#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
      SphP[i].i.Injected_BH_Energy = 0;
#endif
    }

#ifdef TWODIMS
  for(i = 0; i < NumPart; i++)
    {
      P[i].Pos[2] = 0;
      P[i].Vel[2] = 0;

      P[i].g.GravAccel[2] = 0;

      if(P[i].Type == 0)
	{
	  SphP[i].VelPred[2] = 0;
	  SphP[i].a.HydroAccel[2] = 0;
	}
    }
#endif

#ifdef ASSIGN_NEW_IDS
  assign_unique_ids();
#endif

#ifndef NOTEST_FOR_IDUNIQUENESS
  test_id_uniqueness();
#endif


  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  TreeReconstructFlag = 1;

  domain_Decomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  set_softenings();

  /* will build tree */
  ngb_treebuild();

  All.Ti_Current = 0;

#ifdef START_WITH_EXTRA_NGBDEV
  MaxNumNgbDeviationMerk = All.MaxNumNgbDeviation;
  All.MaxNumNgbDeviation = All.MaxNumNgbDeviationStart;
#endif

  if(RestartFlag != 3)
    setup_smoothinglengths();

#ifdef VORONOI
  for(i = 0; i < N_gas; i++)
    SphP[i].MaxDelaunayRadius = SphP[i].Hsml;

  voronoi_mesh();

#ifdef MESHRELAX_DENSITY_IN_INPUT
  if(RestartFlag == 0)
    for(i = 0; i < N_gas; i++)
      P[i].Mass *= SphP[i].Volume;
#endif

#ifdef VORONOI_MESHRELAX
  double mass, masstot, egy, egytot;

  for(i = 0, mass = 0, egy = 0; i < N_gas; i++)
    {
      mass += P[i].Mass;
      egy += P[i].Mass * SphP[i].Entropy;
      SphP[i].Pressure = GAMMA_MINUS1 * SphP[i].Entropy * P[i].Mass / SphP[i].Volume;
    }
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&egy, &egytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  All.MeanMass = masstot / All.TotN_gas;
#ifdef TWODIMS
  All.MeanPressure = GAMMA_MINUS1 * egytot / (boxSize_X * boxSize_Y);
#else
  All.MeanPressure = GAMMA_MINUS1 * egytot / (boxSize_X * boxSize_Y * boxSize_Z);
#endif
#endif

  voronoi_density();

  myfree(List_P);
  myfree(ListExports);

  myfree(DT);
  myfree(DP - 5);
  myfree(VF);
#endif

#ifdef START_WITH_EXTRA_NGBDEV
  All.MaxNumNgbDeviation = MaxNumNgbDeviationMerk;
#endif

#ifdef HEALPIX
  //this should be readed in the parameterfile
  All.Nside = 8;
  //
  if(ThisTask == 0)
    printf(" First calculation of Healpix %i with %i \n", All.Nside, NSIDE2NPIX(All.Nside));
  // initialize the healpix array (just in case)
  All.healpixmap = (float *) malloc(NSIDE2NPIX(All.Nside) * sizeof(float));
  for(i = 0; i < NSIDE2NPIX(All.Nside); i++)
    All.healpixmap[i] = 0;

  healpix_halo(All.healpixmap);

#endif
  /* at this point, the entropy variable actually contains the 
   * internal energy, read in from the initial conditions file. 
   * Once the density has been computed, we can convert to entropy.
   */

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      if(header.flag_entropy_instead_u == 0)
	{
#ifndef EOS_DEGENERATE

#if !defined(VORONOI_MESHRELAX) && !defined(TRADITIONAL_SPH_FORMULATION)

	  if(ThisTask == 0 && i == 0)
	    printf("Converting u -> entropy !\n");

	  SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].d.Density / a3, GAMMA_MINUS1);
#endif

#else
	  for(j = 0; j < EOS_NSPECIES; j++)
	    {
	      SphP[i].dxnuc[j] = 0;
	    }

	  SphP[i].Entropy *= All.UnitEnergy_in_cgs;
	  /* call eos with physical units, energy and entropy are always stored in physical units */
	  SphP[i].temp = 1.0;
	  eos_calc_egiven_v(SphP[i].d.Density * All.UnitDensity_in_cgs, SphP[i].xnuc, SphP[i].dxnuc,
			    0, SphP[i].Entropy, &SphP[i].temp, &SphP[i].Pressure, &SphP[i].dpdr);
	  SphP[i].Pressure /= All.UnitPressure_in_cgs;
#endif
	}

      SphP[i].e.DtEntropy = 0;

      SphP[i].v.DivVel = 0;

#ifdef MACHNUM
      SphP[i].Shock_MachNumber = 1.0;
#ifdef COSMIC_RAYS
      Pth1 = SphP[i].Entropy * pow(SphP[i].d.Density / a3, GAMMA);

#ifdef CR_IC_PHYSICAL
      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  C_phys[CRpop] = SphP[i].CR_C0[CRpop];
	  q_phys[CRpop] = SphP[i].CR_q0[CRpop];
	}
#else
      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  C_phys[CRpop] = SphP[i].CR_C0[CRpop] * pow(SphP[i].d.Density, (All.CR_Alpha[CRpop] - 1.0) / 3.0);
	  q_phys[CRpop] = SphP[i].CR_q0[CRpop] * pow(SphP[i].d.Density, 1.0 / 3.0);
	}
#endif
      SphP[i].PreShock_XCR = 0.0;
      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  rBeta[CRpop] = gsl_sf_beta((All.CR_Alpha[CRpop] - 2.0) * 0.5, (3.0 - All.CR_Alpha[CRpop]) * 0.5) *
	    gsl_sf_beta_inc((All.CR_Alpha[CRpop] - 2.0) * 0.5, (3.0 - All.CR_Alpha[CRpop]) * 0.5,
			    1.0 / (1.0 + q_phys[CRpop] * q_phys[CRpop]));

	  PCR1[CRpop] = C_phys[CRpop] * c2 * SphP[i].d.Density * rBeta[CRpop] / 6.0;
	  PCR1[CRpop] *= pow(atime, -3.0 * GAMMA);
	  SphP[i].PreShock_XCR += PCR1[CRpop] / Pth1;
	}

      SphP[i].PreShock_PhysicalDensity = SphP[i].d.Density / a3;
      SphP[i].PreShock_PhysicalEnergy =
	SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].d.Density / a3, GAMMA_MINUS1);

      SphP[i].Shock_DensityJump = 1.0001;
      SphP[i].Shock_EnergyJump = 1.0;
#endif /* COSMIC_RAYS */
#ifdef OUTPUT_PRESHOCK_CSND
      Pth1 = SphP[i].Entropy * pow(SphP[i].d.Density / a3, GAMMA);
      SphP[i].PreShock_PhysicalSoundSpeed =
	sqrt(GAMMA * Pth1 / SphP[i].d.Density) * pow(atime, -3. / 2. * GAMMA_MINUS1);
      SphP[i].PreShock_PhysicalDensity = SphP[i].d.Density / a3;
#endif /* OUTPUT_PRESHOCK_CSND */
#endif /* MACHNUM */

#ifdef REIONIZATION
      All.not_yet_reionized = 1;
#endif

#ifdef CR_IC_PHYSICAL
      /* Scale CR variables so that values from IC file are now the
       * physical values, not the adiabatic invariants
       */
      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  SphP[i].CR_C0[CRpop] *= pow(SphP[i].d.Density, (1.0 - All.CR_Alpha[CRpop]) / 3.0);
	  SphP[i].CR_q0[CRpop] *= pow(SphP[i].d.Density, -1.0 / 3.0);
	}
#endif

#ifdef CR_INITPRESSURE

      cr_pressure = CR_INITPRESSURE * SphP[i].Entropy * pow(SphP[i].d.Density / a3, GAMMA);
      SphP[i].Entropy *= (1 - CR_INITPRESSURE);
      q_phys = 1.685;

      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  C_phys[CRpop] =
	    cr_pressure / (SphP[i].d.Density / a3 * CR_Tab_Beta(q_phys, CRpop) *
			   (C / All.UnitVelocity_in_cm_per_s) * (C / All.UnitVelocity_in_cm_per_s) / 6.0);

	  SphP[i].CR_C0[CRpop] = C_phys[CRpop] * pow(SphP[i].d.Density, (1.0 - All.CR_Alpha[CRpop]) / 3.0);
	  SphP[i].CR_q0[CRpop] = q_phys * pow(SphP[i].d.Density, -1.0 / 3.0);
	}
#endif
    }

  if(RestartFlag == 3)
    {
#ifdef FOF
      fof_fof(RestartSnapNum);
#endif
      endrun(0);
    }

#ifdef CHEMISTRY
  if(ThisTask == 0)
    {
      printf("Initial abundances: \n");
      printf("HI=%g, HII=%g, HeI=%g, HeII=%g, HeIII=%g \n",
	     SphP[1].HI, SphP[1].HII, SphP[1].HeI, SphP[1].HeII, SphP[1].HeIII);

      printf("HM=%g, H2I=%g, H2II=%g, elec=%g, %d\n",
	     SphP[1].HM, SphP[1].H2I, SphP[1].H2II, SphP[1].elec, P[1].ID);

      printf("x=%g, y=%g, z=%g, vx=%g, vy=%g, vz=%g, density=%g, entropy=%g\n",
	     P[N_gas - 1].Pos[0], P[N_gas - 1].Pos[1], P[N_gas - 1].Pos[2], P[N_gas - 1].Vel[0],
	     P[N_gas - 1].Vel[1], P[N_gas - 1].Vel[2], SphP[N_gas - 1].Density, SphP[N_gas - 1].Entropy);
    }

  /* need predict the cooling time and elec_dot here */
  min_t_cool = min_t_elec = 1.0e30;
  max_t_cool = max_t_elec = -1.0e30;

  for(i = 0; i < N_gas; i++)
    {
      a_start = All.Time;
      a_end = All.Time + 0.001;	/* 0.001 as an arbitrary value */

      ifunc = compute_abundances(0, i, a_start, a_end);


      if(fabs(SphP[i].t_cool) < min_t_cool)
	min_t_cool = fabs(SphP[i].t_cool);
      if(fabs(SphP[i].t_cool) > max_t_cool)
	max_t_cool = fabs(SphP[i].t_cool);

      if(fabs(SphP[i].t_elec) < min_t_elec)
	min_t_elec = fabs(SphP[i].t_elec);
      if(fabs(SphP[i].t_elec) > max_t_elec)
	max_t_elec = fabs(SphP[i].t_elec);

    }

  fprintf(stdout, "PE %d t_cool min= %g, max= %g in yrs \n", ThisTask, min_t_cool, max_t_cool);
  fflush(stdout);
  fprintf(stdout, "PE %d t_elec min= %g, max= %g in yrs \n", ThisTask, min_t_elec, max_t_elec);
  fflush(stdout);

#endif



}


/*! This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));
#ifdef TIMEDEPGRAV
  omega *= All.Gini / All.G;
#endif
  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;

  if(RestartFlag == 0 || RestartFlag == 2)
    {
#ifdef GASRETURN
      for(i = 0; i < NumPart; i++)
#else
      for(i = 0; i < N_gas; i++)
#endif
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

#ifndef READ_HSML
#ifndef TWODIMS
	  PPP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  PPP[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
         if(PPP[i].Hsml > 200.0 * All.SofteningTable[0])
           PPP[i].Hsml = 200.0 * All.SofteningTable[0];

         if(P[i].Mass <= 0.0 || Nodes[no].u.d.mass <= 0.0 || Nodes[no].len == 0)
           PPP[i].Hsml = All.SofteningTable[0];

         if(PPP[i].Hsml < 0.025 * All.SofteningTable[0])
           PPP[i].Hsml = 0.025 * All.SofteningTable[0];
#endif
	}
    }

#ifdef BLACK_HOLES
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 5)
	  PPP[i].Hsml = All.SofteningTable[5];
    }
#endif

#ifdef RADTRANSFER
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 4)
	  PPP[i].Hsml = All.SofteningTable[4];
    }
#endif

  density();

#if defined(MAGNETIC) && defined(BFROMROTA)
  if(RestartFlag == 0)
    {
      if(ThisTask == 0)
	printf("Converting: Vector Potential -> Bfield\n");
      rot_a();
    }
#endif

}


void assign_unique_ids(void)
{
  int i, *numpartlist;
  MyIDType idfirst;

  numpartlist = mymalloc(NTask * sizeof(int));

  MPI_Allgather(&NumPart, 1, MPI_INT, numpartlist, 1, MPI_INT, MPI_COMM_WORLD);

  idfirst = 1;

  for(i = 0; i < ThisTask; i++)
    idfirst += numpartlist[i];

  for(i = 0; i < NumPart; i++)
    {
      P[i].ID = idfirst;
      idfirst++;
    }

  myfree(numpartlist);
}



void test_id_uniqueness(void)
{
  int i;
  double t0, t1;
  MyIDType *ids, *ids_first;

  if(ThisTask == 0)
    {
      printf("Testing ID uniqueness...\n");
      fflush(stdout);
    }

  if(NumPart == 0)
    {
      printf("need at least one particle per cpu\n");
      endrun(8);
    }

  t0 = second();

#ifndef SPH_BND_PARTICLES
  ids = (MyIDType *) mymalloc(NumPart * sizeof(MyIDType));
  ids_first = (MyIDType *) mymalloc(NTask * sizeof(MyIDType));

  for(i = 0; i < NumPart; i++)
    ids[i] = P[i].ID;

  parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);

  for(i = 1; i < NumPart; i++)
    if(ids[i] == ids[i - 1])
      {
#ifdef LONGIDS
	printf("non-unique ID=%d%09d found on task=%d (i=%d NumPart=%d)\n",
	       (int) (ids[i] / 1000000000), (int) (ids[i] % 1000000000), ThisTask, i, NumPart);

#else
	printf("non-unique ID=%d found on task=%d   (i=%d NumPart=%d)\n", (int) ids[i], ThisTask, i, NumPart);
#endif
	endrun(12);
      }

  MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

  if(ThisTask < NTask - 1)
    if(ids[NumPart - 1] == ids_first[ThisTask + 1])
      {
	printf("non-unique ID=%d found on task=%d\n", (int) ids[NumPart - 1], ThisTask);
	endrun(13);
      }

  myfree(ids_first);
  myfree(ids);
#endif

  t1 = second();

  if(ThisTask == 0)
    {
      printf("success.  took=%g sec\n", timediff(t0, t1));
      fflush(stdout);
    }
}

int compare_IDs(const void *a, const void *b)
{
  if(*((MyIDType *) a) < *((MyIDType *) b))
    return -1;

  if(*((MyIDType *) a) > *((MyIDType *) b))
    return +1;

  return 0;
}
