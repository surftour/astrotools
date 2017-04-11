#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file blackhole.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
 */

#ifdef BLACK_HOLES


static struct blackholedata_in
{
  MyDouble Pos[3];
  MyFloat Density;
  MyFloat Mdot;
  MyFloat Dt;
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat BH_Mass;
  MyFloat Vel[3];
  MyFloat Csnd;
  MyIDType ID;
  int Index;
  int NodeList[NODELISTLENGTH];
}
 *BlackholeDataIn, *BlackholeDataGet;

static struct blackholedata_out
{
  MyLongDouble Mass;
  MyLongDouble BH_Mass;
  MyLongDouble AccretedMomentum[3];
#ifdef REPOSITION_ON_POTMIN
  MyFloat BH_MinPotPos[3];
  MyFloat BH_MinPot;
#endif
#ifdef BH_BUBBLES
  MyLongDouble BH_Mass_bubbles;
#ifdef UNIFIED_FEEDBACK
  MyLongDouble BH_Mass_radio;
#endif
#endif
#ifdef BH_COUNTPROGS
  int BH_CountProgs;
#endif
}
 *BlackholeDataResult, *BlackholeDataOut;



static double hubble_a, ascale;

static int N_gas_swallowed, N_BH_swallowed;

void blackhole_accretion(void)
{
  int i, j, k, n, bin;
  int ndone_flag, ndone;
  int ngrp, sendTask, recvTask, place, nexport, nimport, dummy;
  int Ntot_gas_swallowed, Ntot_BH_swallowed;
  double mdot, rho, bhvel, soundspeed, meddington, dt, mdot_in_msun_per_year;
  double mass_real, total_mass_real, medd, total_mdoteddington;
  double mass_holes, total_mass_holes, total_mdot;

#ifdef BH_BUBBLES
  MyFloat bh_center[3];
  double *bh_dmass, *tot_bh_dmass;
  float *bh_posx, *bh_posy, *bh_posz;
  float *tot_bh_posx, *tot_bh_posy, *tot_bh_posz;
  int l, num_activebh = 0, total_num_activebh = 0;
  int *common_num_activebh, *disp;
  int *bh_id, *tot_bh_id;
#endif
  MPI_Status status;

  if(ThisTask == 0)
    {
      printf("Beginning black-hole accretion\n");
      fflush(stdout);
    }

  CPU_Step[CPU_MISC] += measure_time();

  if(All.ComovingIntegrationOn)
    {
      ascale = All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    hubble_a = ascale = 1;


  /* Let's first compute the Mdot values */

  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5)
      {
	mdot = 0;		/* if no accretion model is enabled, we have mdot=0 */

	rho = P[n].b1.BH_Density;

	bhvel = sqrt(pow(P[n].Vel[0] - P[n].b3.BH_SurroundingGasVel[0], 2) +
		     pow(P[n].Vel[1] - P[n].b3.BH_SurroundingGasVel[1], 2) +
		     pow(P[n].Vel[2] - P[n].b3.BH_SurroundingGasVel[2], 2));

	if(All.ComovingIntegrationOn)
	  {
	    bhvel /= All.Time;
	    rho /= pow(All.Time, 3);
	  }

	soundspeed = sqrt(GAMMA * P[n].b2.BH_Entropy * pow(rho, GAMMA_MINUS1));

	/* Note: we take here a radiative efficiency of 0.1 for Eddington accretion */

	meddington = (4 * M_PI * GRAVITY * C * PROTONMASS / (0.1 * C * C * THOMPSON)) * P[n].BH_Mass
	  * All.UnitTime_in_s / All.HubbleParam;

#ifdef BONDI
	mdot = 4. * M_PI * All.BlackHoleAccretionFactor * All.G * All.G *
	  P[n].BH_Mass * P[n].BH_Mass * rho / pow((pow(soundspeed, 2) + pow(bhvel, 2)), 1.5);
#endif


#ifdef ENFORCE_EDDINGTON_LIMIT
	if(mdot > All.BlackHoleEddingtonFactor * meddington)
	  mdot = All.BlackHoleEddingtonFactor * meddington;
#endif
	P[n].BH_Mdot = mdot;

	if(P[n].BH_Mass > 0)
	  {
#ifndef LONGIDS
	    fprintf(FdBlackHolesDetails, "BH=%u %g %g %g %g %g\n",
		    P[n].ID, All.Time, P[n].BH_Mass, mdot, rho, soundspeed);
#else
	    fprintf(FdBlackHolesDetails, "BH=%llu %g %g %g %g %g\n",
		    P[n].ID, All.Time, P[n].BH_Mass, mdot, rho, soundspeed);
#endif
	  }

	dt = (P[n].TimeBin ? (1 << P[n].TimeBin) : 0) * All.Timebase_interval / hubble_a;

#ifdef BH_DRAG
	/* add a drag force for the black-holes,
	   accounting for the accretion */
	double fac;

	if(P[n].BH_Mass > 0)
	  {
	    /*
	       fac = P[n].BH_Mdot * dt / P[n].BH_Mass;
	     */
	    fac = meddington * dt / P[n].BH_Mass;

	    if(fac > 1)
	      fac = 1;

	    if(dt > 0)
	      for(k = 0; k < 3; k++)
		P[n].g.GravAccel[k] +=
		  -ascale * ascale * fac / dt * (P[n].Vel[k] - P[n].b3.BH_SurroundingGasVel[k]) / ascale;
	  }
#endif

	P[n].BH_Mass += P[n].BH_Mdot * dt;

#ifdef BH_BUBBLES
	P[n].BH_Mass_bubbles += P[n].BH_Mdot * dt;
#ifdef UNIFIED_FEEDBACK
	if(P[n].BH_Mdot < All.RadioThreshold * meddington)
	  P[n].BH_Mass_radio += P[n].BH_Mdot * dt;
#endif
#endif
      }


  /* Now let's invoke the functions that stochasticall swallow gas
   * and deal with black hole mergers.
   */

  if(ThisTask == 0)
    {
      printf("Start swallowing of gas particles and black holes\n");
      fflush(stdout);
    }


  N_gas_swallowed = N_BH_swallowed = 0;


  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct blackholedata_in) +
					     sizeof(struct blackholedata_out) +
					     sizemax(sizeof(struct blackholedata_in),
						     sizeof(struct blackholedata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  /** Let's first spread the feedback energy, and determine which particles may be swalled by whom */

  i = FirstActiveParticle;	/* first particle for this task */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type == 5)
	  if(blackhole_evaluate(i, 0, &nexport, Send_count) < 0)
	    break;

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      BlackholeDataGet = (struct blackholedata_in *) mymalloc(nimport * sizeof(struct blackholedata_in));
      BlackholeDataIn = (struct blackholedata_in *) mymalloc(nexport * sizeof(struct blackholedata_in));

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
	      BlackholeDataIn[j].Pos[k] = P[place].Pos[k];
	      BlackholeDataIn[j].Vel[k] = P[place].Vel[k];
	    }

	  BlackholeDataIn[j].Hsml = PPP[place].Hsml;
	  BlackholeDataIn[j].Mass = P[place].Mass;
	  BlackholeDataIn[j].BH_Mass = P[place].BH_Mass;
	  BlackholeDataIn[j].Density = P[place].b1.BH_Density;
	  BlackholeDataIn[j].Mdot = P[place].BH_Mdot;
	  BlackholeDataIn[j].Csnd =
	    sqrt(GAMMA * P[place].b2.BH_Entropy *
		 pow(P[place].b1.BH_Density / (ascale * ascale * ascale), GAMMA_MINUS1));
	  BlackholeDataIn[j].Dt =
	    (P[place].TimeBin ? (1 << P[place].TimeBin) : 0) * All.Timebase_interval / hubble_a;
	  BlackholeDataIn[j].ID = P[place].ID;

	  memcpy(BlackholeDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}


      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &BlackholeDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}


      myfree(BlackholeDataIn);
      BlackholeDataResult = (struct blackholedata_out *) mymalloc(nimport * sizeof(struct blackholedata_out));
      BlackholeDataOut = (struct blackholedata_out *) mymalloc(nexport * sizeof(struct blackholedata_out));


      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
	blackhole_evaluate(j, 1, &dummy, &dummy);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &BlackholeDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the result to the particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

#ifdef REPOSITION_ON_POTMIN
	  if(P[place].BH_MinPot > BlackholeDataOut[j].BH_MinPot)
	    {
	      P[place].BH_MinPot = BlackholeDataOut[j].BH_MinPot;
	      for(k = 0; k < 3; k++)
		P[place].BH_MinPotPos[k] = BlackholeDataOut[j].BH_MinPotPos[k];
	    }
#endif
	}

      myfree(BlackholeDataOut);
      myfree(BlackholeDataResult);
      myfree(BlackholeDataGet);
    }
  while(ndone < NTask);





  /* Now do the swallowing of particles */

  i = FirstActiveParticle;	/* first particle for this task */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type == 5)
	  if(P[i].SwallowID == 0)
	    if(blackhole_evaluate_swallow(i, 0, &nexport, Send_count) < 0)
	      break;


      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

      MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      BlackholeDataGet = (struct blackholedata_in *) mymalloc(nimport * sizeof(struct blackholedata_in));
      BlackholeDataIn = (struct blackholedata_in *) mymalloc(nexport * sizeof(struct blackholedata_in));

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    BlackholeDataIn[j].Pos[k] = P[place].Pos[k];

	  BlackholeDataIn[j].Hsml = PPP[place].Hsml;
	  BlackholeDataIn[j].BH_Mass = P[place].BH_Mass;
	  BlackholeDataIn[j].ID = P[place].ID;

	  memcpy(BlackholeDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}


      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &BlackholeDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
		}
	    }
	}


      myfree(BlackholeDataIn);
      BlackholeDataResult = (struct blackholedata_out *) mymalloc(nimport * sizeof(struct blackholedata_out));
      BlackholeDataOut = (struct blackholedata_out *) mymalloc(nexport * sizeof(struct blackholedata_out));


      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
	blackhole_evaluate_swallow(j, 1, &dummy, &dummy);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &BlackholeDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct blackholedata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the result to the particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  P[place].b4.dBH_accreted_Mass += BlackholeDataOut[j].Mass;
	  P[place].b5.dBH_accreted_BHMass += BlackholeDataOut[j].BH_Mass;
#ifdef BH_BUBBLES
	  P[place].b7.dBH_accreted_BHMass_bubbles += BlackholeDataOut[j].BH_Mass_bubbles;
#ifdef UNIFIED_FEEDBACK
	  P[place].b8.dBH_accreted_BHMass_radio += BlackholeDataOut[j].BH_Mass_radio;
#endif
#endif
	  for(k = 0; k < 3; k++)
	    P[place].b6.dBH_accreted_momentum[k] += BlackholeDataOut[j].AccretedMomentum[k];

#ifdef BH_COUNTPROGS
	  P[place].BH_CountProgs += BlackholeDataOut[j].BH_CountProgs;
#endif
	}

      myfree(BlackholeDataOut);
      myfree(BlackholeDataResult);
      myfree(BlackholeDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);


  MPI_Reduce(&N_gas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Accretion done: %d gas particles swallowed, %d BH particles swallowed\n",
	     Ntot_gas_swallowed, Ntot_BH_swallowed);
      fflush(stdout);
    }



#ifdef REPOSITION_ON_POTMIN
  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5)
      for(k = 0; k < 3; k++)
	P[n].Pos[k] = P[n].BH_MinPotPos[k];
#endif


  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	{
	  TimeBin_BH_mass[n] = 0;
	  TimeBin_BH_dynamicalmass[n] = 0;
	  TimeBin_BH_Mdot[n] = 0;
	  TimeBin_BH_Medd[n] = 0;
	}
    }

  for(n = FirstActiveParticle; n >= 0; n = NextActiveParticle[n])
    if(P[n].Type == 5)
      {
#ifdef FLTROUNDOFFREDUCTION
	P[n].b4.BH_accreted_Mass = FLT(P[n].b4.dBH_accreted_Mass);
	P[n].b5.BH_accreted_BHMass = FLT(P[n].b5.dBH_accreted_BHMass);
#ifdef BH_BUBBLES
	P[n].b7.BH_accreted_BHMass_bubbles = FLT(P[n].b7.dBH_accreted_BHMass_bubbles);
#ifdef UNIFIED_FEEDBACK
	P[n].b8.BH_accreted_BHMass_radio = FLT(P[n].b8.dBH_accreted_BHMass_radio);
#endif
#endif
	for(k = 0; k < 3; k++)
	  P[n].b6.BH_accreted_momentum[k] = FLT(P[n].b6.dBH_accreted_momentum[k]);
#endif
	if(P[n].b4.BH_accreted_Mass > 0)
	  {
	    for(k = 0; k < 3; k++)
	      P[n].Vel[k] =
		(P[n].Vel[k] * P[n].Mass + P[n].b6.BH_accreted_momentum[k]) /
		(P[n].Mass + P[n].b4.BH_accreted_Mass);

	    P[n].Mass += P[n].b4.BH_accreted_Mass;
	    P[n].BH_Mass += P[n].b5.BH_accreted_BHMass;
#ifdef BH_BUBBLES
	    P[n].BH_Mass_bubbles += P[n].b7.BH_accreted_BHMass_bubbles;
#ifdef UNIFIED_FEEDBACK
	    P[n].BH_Mass_radio += P[n].b8.BH_accreted_BHMass_radio;
#endif
#endif
	    P[n].b4.BH_accreted_Mass = 0;
	  }
	bin = P[n].TimeBin;
	TimeBin_BH_mass[bin] += P[n].BH_Mass;
	TimeBin_BH_dynamicalmass[bin] += P[n].Mass;
	TimeBin_BH_Mdot[bin] += P[n].BH_Mdot;
	if(P[n].BH_Mass > 0)
	  TimeBin_BH_Medd[bin] += P[n].BH_Mdot / P[n].BH_Mass;

#ifdef BH_BUBBLES
	if(P[n].BH_Mass_bubbles > 0
	   && P[n].BH_Mass_bubbles > All.BlackHoleRadioTriggeringFactor * P[n].BH_Mass_ini)
	  num_activebh++;
#endif
      }

#ifdef BH_BUBBLES
  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  MPI_Allreduce(&num_activebh, &total_num_activebh, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("The total number of active BHs is: %d\n", total_num_activebh);
      fflush(stdout);
    }

  if(total_num_activebh > 0)
    {
      bh_dmass = mymalloc(num_activebh * sizeof(double));
      tot_bh_dmass = mymalloc(total_num_activebh * sizeof(double));
      bh_posx = mymalloc(num_activebh * sizeof(float));
      bh_posy = mymalloc(num_activebh * sizeof(float));
      bh_posz = mymalloc(num_activebh * sizeof(float));
      tot_bh_posx = mymalloc(total_num_activebh * sizeof(float));
      tot_bh_posy = mymalloc(total_num_activebh * sizeof(float));
      tot_bh_posz = mymalloc(total_num_activebh * sizeof(float));
      bh_id = mymalloc(num_activebh * sizeof(int));
      tot_bh_id = mymalloc(total_num_activebh * sizeof(int));

      for(n = 0; n < num_activebh; n++)
	{
	  bh_dmass[n] = 0.0;
	  bh_posx[n] = 0.0;
	  bh_posy[n] = 0.0;
	  bh_posz[n] = 0.0;
	  bh_id[n] = 0;
	}

      for(n = 0; n < total_num_activebh; n++)
	{
	  tot_bh_dmass[n] = 0.0;
	  tot_bh_posx[n] = 0.0;
	  tot_bh_posy[n] = 0.0;
	  tot_bh_posz[n] = 0.0;
	  tot_bh_id[n] = 0;
	}

      for(n = FirstActiveParticle, l = 0; n >= 0; n = NextActiveParticle[n])
	if(P[n].Type == 5)
	  {
	    if(P[n].BH_Mass_bubbles > 0
	       && P[n].BH_Mass_bubbles > All.BlackHoleRadioTriggeringFactor * P[n].BH_Mass_ini)
	      {
#ifndef UNIFIED_FEEDBACK
		bh_dmass[l] = P[n].BH_Mass_bubbles - P[n].BH_Mass_ini;
#else
		bh_dmass[l] = P[n].BH_Mass_radio - P[n].BH_Mass_ini;
		P[n].BH_Mass_radio = P[n].BH_Mass;
#endif
		P[n].BH_Mass_ini = P[n].BH_Mass;
		P[n].BH_Mass_bubbles = P[n].BH_Mass;

		bh_posx[l] = P[n].Pos[0];
		bh_posy[l] = P[n].Pos[1];
		bh_posz[l] = P[n].Pos[2];
		bh_id[l] = P[n].ID;

		l++;
	      }
	  }
      common_num_activebh = mymalloc(NTask * sizeof(int));
      disp = mymalloc(NTask * sizeof(int));

      MPI_Allgather(&num_activebh, 1, MPI_INT, common_num_activebh, 1, MPI_INT, MPI_COMM_WORLD);

      for(k = 1, disp[0] = 0; k < NTask; k++)
	disp[k] = disp[k - 1] + common_num_activebh[k - 1];


      MPI_Allgatherv(bh_dmass, num_activebh, MPI_DOUBLE, tot_bh_dmass, common_num_activebh, disp, MPI_DOUBLE,
		     MPI_COMM_WORLD);
      MPI_Allgatherv(bh_posx, num_activebh, MPI_FLOAT, tot_bh_posx, common_num_activebh, disp, MPI_FLOAT,
		     MPI_COMM_WORLD);
      MPI_Allgatherv(bh_posy, num_activebh, MPI_FLOAT, tot_bh_posy, common_num_activebh, disp, MPI_FLOAT,
		     MPI_COMM_WORLD);
      MPI_Allgatherv(bh_posz, num_activebh, MPI_FLOAT, tot_bh_posz, common_num_activebh, disp, MPI_FLOAT,
		     MPI_COMM_WORLD);

      MPI_Allgatherv(bh_id, num_activebh, MPI_INT, tot_bh_id, common_num_activebh, disp, MPI_INT, MPI_COMM_WORLD);

      for(l = 0; l < total_num_activebh; l++)
	{
	  bh_center[0] = tot_bh_posx[l];
	  bh_center[1] = tot_bh_posy[l];
	  bh_center[2] = tot_bh_posz[l];

	  if(tot_bh_dmass[l] > 0)
	    bh_bubble(tot_bh_dmass[l], bh_center, tot_bh_id[l]);	  
	  
	}

      myfree(disp);
      myfree(common_num_activebh);
      myfree(tot_bh_id);
      myfree(bh_id);
      myfree(tot_bh_posz);
      myfree(tot_bh_posy);
      myfree(tot_bh_posx);
      myfree(bh_posz);
      myfree(bh_posy);
      myfree(bh_posx);
      myfree(tot_bh_dmass);
      myfree(bh_dmass);
    }
  myfree(Ngblist);
#endif

  mdot = 0;
  mass_holes = 0;
  mass_real = 0;
  medd = 0;

  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      {
	mass_holes += TimeBin_BH_mass[bin];
	mass_real += TimeBin_BH_dynamicalmass[bin];
	mdot += TimeBin_BH_Mdot[bin];
	medd += TimeBin_BH_Medd[bin];
      }

  MPI_Reduce(&mass_holes, &total_mass_holes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mass_real, &total_mass_real, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mdot, &total_mdot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&medd, &total_mdoteddington, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      /* convert to solar masses per yr */
      mdot_in_msun_per_year =
	total_mdot * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      total_mdoteddington *= 1.0 / ((4 * M_PI * GRAVITY * C * PROTONMASS /
				     (0.1 * C * C * THOMPSON)) * All.UnitTime_in_s);

      fprintf(FdBlackHoles, "%g %d %g %g %g %g %g\n",
	      All.Time, All.TotBHs, total_mass_holes, total_mdot, mdot_in_msun_per_year,
	      total_mass_real, total_mdoteddington);
      fflush(FdBlackHoles);
    }


  fflush(FdBlackHolesDetails);

  CPU_Step[CPU_BLACKHOLES] += measure_time();
}






int blackhole_evaluate(int target, int mode, int *nexport, int *nSend_local)
{
  int startnode, numngb, j, k, n, index, id, listindex = 0;
  MyFloat *pos, *velocity, h_i, dt, mdot, rho, mass, bh_mass, csnd;
  double dx, dy, dz, h_i2, r2, r, u, hinv, hinv3, wk, vrel;

#ifdef UNIFIED_FEEDBACK
  double meddington;
#endif

#ifdef BH_KINETICFEEDBACK
  /*  double deltavel; */
  double activetime, activeenergy;
#endif
#ifdef BH_THERMALFEEDBACK
  double energy;
#endif
#ifdef REPOSITION_ON_POTMIN
  MyFloat minpotpos[3] = { 0, 0, 0 }, minpot = 1.0e30;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      rho = P[target].b1.BH_Density;
      mdot = P[target].BH_Mdot;
      dt = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0) * All.Timebase_interval / hubble_a;
      h_i = PPP[target].Hsml;
      mass = P[target].Mass;
      bh_mass = P[target].BH_Mass;
      velocity = P[target].Vel;
      csnd =
	sqrt(GAMMA * P[target].b2.BH_Entropy *
	     pow(P[target].b1.BH_Density / (ascale * ascale * ascale), GAMMA_MINUS1));
      index = target;
      id = P[target].ID;
#ifdef BH_KINETICFEEDBACK
      activetime = P[target].ActiveTime;
      activeenergy = P[target].ActiveEnergy;
#endif

    }
  else
    {
      pos = BlackholeDataGet[target].Pos;
      rho = BlackholeDataGet[target].Density;
      mdot = BlackholeDataGet[target].Mdot;
      dt = BlackholeDataGet[target].Dt;
      h_i = BlackholeDataGet[target].Hsml;
      mass = BlackholeDataGet[target].Mass;
      bh_mass = BlackholeDataGet[target].BH_Mass;
      velocity = BlackholeDataGet[target].Vel;
      csnd = BlackholeDataGet[target].Csnd;
      index = BlackholeDataGet[target].Index;
      id = BlackholeDataGet[target].ID;
#ifdef BH_KINETICFEEDBACK
      activetime = BlackholeDataGet[target].ActiveTime;
      activeenergy = BlackholeDataGet[target].ActiveEnergy;
#endif
    }

  /* initialize variables before SPH loop is started */
  h_i2 = h_i * h_i;

  /* Now start the actual SPH computation for this particle */
  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = BlackholeDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_blackhole(pos, h_i, target, &startnode, mode, nexport, nSend_local);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

	      if(P[j].Mass > 0)
		{
		  if(mass > 0)
		    {
		      dx = pos[0] - P[j].Pos[0];
		      dy = pos[1] - P[j].Pos[1];
		      dz = pos[2] - P[j].Pos[2];
#ifdef PERIODIC			/*  now find the closest image in the given box size  */
		      if(dx > boxHalf_X)
			dx -= boxSize_X;
		      if(dx < -boxHalf_X)
			dx += boxSize_X;
		      if(dy > boxHalf_Y)
			dy -= boxSize_Y;
		      if(dy < -boxHalf_Y)
			dy += boxSize_Y;
		      if(dz > boxHalf_Z)
			dz -= boxSize_Z;
		      if(dz < -boxHalf_Z)
			dz += boxSize_Z;
#endif
		      r2 = dx * dx + dy * dy + dz * dz;

		      if(r2 < h_i2)
			{
#ifdef REPOSITION_ON_POTMIN
			  /* if this option is switched on, we may also encounter dark matter particles or stars */
			  if(P[j].p.Potential < minpot)
			    {
			      if(P[j].Type == 0 || P[j].Type == 1 || P[j].Type == 4 || P[j].Type == 5)
				{
				  /* compute relative velocities */

				  for(k = 0, vrel = 0; k < 3; k++)
				    vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);

				  vrel = sqrt(vrel) / ascale;

				  if(vrel <= 0.25 * csnd)
				    {
				      minpot = P[j].p.Potential;
				      for(k = 0; k < 3; k++)
					minpotpos[k] = P[j].Pos[k];
				    }
				}
			    }
#endif
			  if(P[j].Type == 5)	/* we have a black hole merger */
			    {
			      if(id != P[j].ID)
				{
				  /* compute relative velocity of BHs */

				  for(k = 0, vrel = 0; k < 3; k++)
				    vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);

				  vrel = sqrt(vrel) / ascale;

				  if(vrel > 0.5 * csnd)
				    {
#ifndef LONGIDS
				      fprintf(FdBlackHolesDetails,
					      "ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g csnd=%g\n",
					      ThisTask, All.Time, id, P[j].ID, vrel, csnd);
#else
				      fprintf(FdBlackHolesDetails,
					      "ThisTask=%d, time=%g: id=%llu would like to swallow %llu, but vrel=%g csnd=%g\n",
					      ThisTask, All.Time, id, P[j].ID, vrel, csnd);
#endif
				    }
				  else
				    {
				      if(P[j].SwallowID < id)
					P[j].SwallowID = id;
				    }
				}
			    }
			  if(P[j].Type == 0)
			    {
			      /* here we have a gas particle */

			      r = sqrt(r2);
			      hinv = 1 / h_i;
#ifndef  TWODIMS
			      hinv3 = hinv * hinv * hinv;
#else
			      hinv3 = hinv * hinv / boxSize_Z;
#endif

			      u = r * hinv;

			      if(u < 0.5)
				wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			      else
				wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

#ifdef SWALLOWGAS
			      /* compute accretion probability */
			      double p, w;

			      if((bh_mass - mass) > 0)
				p = (bh_mass - mass) * wk / rho;
			      else
				p = 0;

			      /* compute random number, uniform in [0,1] */
			      w = get_random_number(P[j].ID);
			      if(w < p)
				{
				  if(P[j].SwallowID < id)
				    P[j].SwallowID = id;
				}
#endif

			      if(P[j].Mass > 0)
				{
#ifdef BH_THERMALFEEDBACK
#ifndef UNIFIED_FEEDBACK
				  energy = All.BlackHoleFeedbackFactor * 0.1 * mdot * dt *
				    pow(C / All.UnitVelocity_in_cm_per_s, 2);

				  SphP[j].i.dInjected_BH_Energy += FLT(energy * P[j].Mass * wk / rho);
#else
				  meddington = (4 * M_PI * GRAVITY * C *
						PROTONMASS / (0.1 * C * C * THOMPSON)) * bh_mass *
				    All.UnitTime_in_s;

				  if(mdot > All.RadioThreshold * meddington)
				    {
				      energy =
					All.BlackHoleFeedbackFactor * 0.1 * mdot * dt * pow(C /
											    All.
											    UnitVelocity_in_cm_per_s,
											    2);

				      SphP[j].i.dInjected_BH_Energy += FLT(energy * P[j].Mass * wk / rho);
				    }
#endif
#endif
				}

			    }
			}
		    }
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = BlackholeDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }



  /* Now collect the result at the right place */
  if(mode == 0)
    {
#ifdef REPOSITION_ON_POTMIN
      for(k = 0; k < 3; k++)
	P[target].BH_MinPotPos[k] = minpotpos[k];
      P[target].BH_MinPot = minpot;
#endif
    }
  else
    {
#ifdef REPOSITION_ON_POTMIN
      for(k = 0; k < 3; k++)
	BlackholeDataResult[target].BH_MinPotPos[k] = minpotpos[k];
      BlackholeDataResult[target].BH_MinPot = minpot;
#endif
    }

  return 0;
}


int blackhole_evaluate_swallow(int target, int mode, int *nexport, int *nSend_local)
{
  int startnode, numngb, j, k, n, id, bin, listindex = 0;
  MyLongDouble accreted_mass, accreted_BH_mass, accreted_momentum[3];
  MyFloat *pos, h_i, bh_mass;

#ifdef BH_BUBBLES
  MyLongDouble accreted_BH_mass_bubbles = 0;
  MyLongDouble accreted_BH_mass_radio = 0;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      id = P[target].ID;
      bh_mass = P[target].BH_Mass;
    }
  else
    {
      pos = BlackholeDataGet[target].Pos;
      h_i = BlackholeDataGet[target].Hsml;
      id = BlackholeDataGet[target].ID;
      bh_mass = BlackholeDataGet[target].BH_Mass;
    }

  accreted_mass = 0;
  accreted_BH_mass = 0;
  accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;
#ifdef BH_COUNTPROGS
  int accreted_BH_progs = 0;
#endif


  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = BlackholeDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_blackhole(pos, h_i, target, &startnode, mode, nexport, nSend_local);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

	      if(P[j].SwallowID == id)
		{
		  if(P[j].Type == 5)	/* we have a black hole merger */
		    {
#ifndef LONGIDS
		      fprintf(FdBlackHolesDetails,
			      "ThisTask=%d, time=%g: id=%u swallows %u (%g %g)\n",
			      ThisTask, All.Time, id, P[j].ID, bh_mass, P[j].BH_Mass);
#else
		      fprintf(FdBlackHolesDetails,
			      "ThisTask=%d, time=%g: id=%llu swallows %llu (%g %g)\n",
			      ThisTask, All.Time, id, P[j].ID, bh_mass, P[j].BH_Mass);
#endif

		      accreted_mass += FLT(P[j].Mass);
		      accreted_BH_mass += FLT(P[j].BH_Mass);
#ifdef BH_BUBBLES
		      accreted_BH_mass_bubbles += FLT(P[j].BH_Mass_bubbles - P[j].BH_Mass_ini);
#ifdef UNIFIED_FEEDBACK
		      accreted_BH_mass_radio += FLT(P[j].BH_Mass_radio - P[j].BH_Mass_ini);
#endif
#endif
		      for(k = 0; k < 3; k++)
			accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);

#ifdef BH_COUNTPROGS
		      accreted_BH_progs += P[j].BH_CountProgs;
#endif

                      bin = P[j].TimeBin;

                      TimeBin_BH_mass[bin] -= P[j].BH_Mass;
                      TimeBin_BH_dynamicalmass[bin] -= P[j].Mass;
                      TimeBin_BH_Mdot[bin] -= P[j].BH_Mdot;
                      if(P[j].BH_Mass > 0)
                        TimeBin_BH_Medd[bin] -= P[j].BH_Mdot / P[j].BH_Mass;

		      P[j].Mass = 0;
		      P[j].BH_Mass = 0;
                      P[j].BH_Mdot = 0;

#ifdef BH_BUBBLES
		      P[j].BH_Mass_bubbles = 0;
		      P[j].BH_Mass_ini = 0;
#ifdef UNIFIED_FEEDBACK
		      P[j].BH_Mass_radio = 0;
#endif
#endif
		      N_BH_swallowed++;
		    }
		}

	      if(P[j].Type == 0)
		{
		  if(P[j].SwallowID == id)
		    {
		      accreted_mass += FLT(P[j].Mass);

		      for(k = 0; k < 3; k++)
			accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);

		      P[j].Mass = 0;

		      N_gas_swallowed++;
		    }
		}
	    }
	}
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = BlackholeDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  /* Now collect the result at the right place */
  if(mode == 0)
    {
      P[target].b4.dBH_accreted_Mass = accreted_mass;
      P[target].b5.dBH_accreted_BHMass = accreted_BH_mass;
      for(k = 0; k < 3; k++)
	P[target].b6.dBH_accreted_momentum[k] = accreted_momentum[k];
#ifdef BH_BUBBLES
      P[target].b7.dBH_accreted_BHMass_bubbles = accreted_BH_mass_bubbles;
#ifdef UNIFIED_FEEDBACK
      P[target].b8.dBH_accreted_BHMass_radio = accreted_BH_mass_radio;
#endif
#endif
#ifdef BH_COUNTPROGS
      P[target].BH_CountProgs += accreted_BH_progs;
#endif
    }
  else
    {
      BlackholeDataResult[target].Mass = accreted_mass;
      BlackholeDataResult[target].BH_Mass = accreted_BH_mass;
      for(k = 0; k < 3; k++)
	BlackholeDataResult[target].AccretedMomentum[k] = accreted_momentum[k];
#ifdef BH_BUBBLES
      BlackholeDataResult[target].BH_Mass_bubbles = accreted_BH_mass_bubbles;
#ifdef UNIFIED_FEEDBACK
      BlackholeDataResult[target].BH_Mass_radio = accreted_BH_mass_radio;
#endif
#endif
#ifdef BH_COUNTPROGS
      BlackholeDataResult[target].BH_CountProgs = accreted_BH_progs;
#endif
    }

  return 0;
}




int ngb_treefind_blackhole(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			   int *nexport, int *nsend_local)
{
  int numngb, no, p, task, nexport_save;
  struct NODE *current;
  MyDouble dx, dy, dz, dist;

#ifdef PERIODIC
  MyDouble xtmp;
#endif
  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

#ifndef REPOSITION_ON_POTMIN
	  if(P[p].Type != 0 && P[p].Type != 5)
	    continue;
#endif
	  dist = hsml;
	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
		{
		  Exportflag[task] = target;
		  Exportnodecount[task] = NODELISTLENGTH;
		}

	      if(Exportnodecount[task] == NODELISTLENGTH)
		{
		  if(*nexport >= All.BunchSize)
		    {
		      *nexport = nexport_save;
		      for(task = 0; task < NTask; task++)
			nsend_local[task] = 0;
		      for(no = 0; no < nexport_save; no++)
			nsend_local[DataIndexTable[no].Task]++;
		      return -1;
		    }
		  Exportnodecount[task] = 0;
		  Exportindex[task] = *nexport;
		  DataIndexTable[*nexport].Task = task;
		  DataIndexTable[*nexport].Index = target;
		  DataIndexTable[*nexport].IndexGet = *nexport;
		  *nexport = *nexport + 1;
		  nsend_local[task]++;
		}

	      DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

	      if(Exportnodecount[task] < NODELISTLENGTH)
		DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
		  return numngb;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;;
	  dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}

#ifdef BH_BUBBLES
void bh_bubble(double bh_dmass, MyFloat center[3], int BH_id)
{
  double phi, theta;
  double dx, dy, dz, rr, r2, dE;
  double E_bubble, totE_bubble;
  double BubbleDistance = 0.0, BubbleRadius = 0.0, BubbleEnergy = 0.0;
  double ICMDensity;
  double Mass_bubble, totMass_bubble;
  double u_to_temp_fac;
  MyDouble pos[3];
  int numngb, tot_numngb, startnode, numngb_inbox;
  int n, i, j, dummy;

#ifdef CR_BUBBLES
  double tinj = 0.0, instant_reheat = 0.0;
  double sum_instant_reheat = 0.0, tot_instant_reheat = 0.0;
#endif

  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS /
    BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  if(All.ComovingIntegrationOn)
    {

      BubbleDistance = All.BubbleDistance;
      BubbleRadius = All.BubbleRadius;

      /*switch to comoving if it is assumed that Rbub should be constant with redshift */

      /* BubbleDistance = All.BubbleDistance / All.Time; 
         BubbleRadius = All.BubbleRadius / All.Time; */
    }
  else
    {
      BubbleDistance = All.BubbleDistance;
      BubbleRadius = All.BubbleRadius;
    }

  BubbleEnergy = All.RadioFeedbackFactor * 0.1 * bh_dmass * All.UnitMass_in_g / All.HubbleParam * pow(C, 2);	/*in cgs units */

  phi = 2 * M_PI * get_random_number(BH_id);
  theta = acos(2 * get_random_number(BH_id + 1) - 1);
  rr = pow(get_random_number(BH_id + 2), 1. / 3.) * BubbleDistance;

  pos[0] = sin(theta) * cos(phi);
  pos[1] = sin(theta) * sin(phi);
  pos[2] = cos(theta);

  for(i = 0; i < 3; i++)
    pos[i] *= rr;

  for(i = 0; i < 3; i++)
    pos[i] += center[i];


  /* First, let's see how many particles are in the bubble of the default radius */

  numngb = 0;
  E_bubble = 0.;
  Mass_bubble = 0.;

  startnode = All.MaxPart;
  do
    {
      numngb_inbox = ngb_treefind_variable(pos, BubbleRadius, -1, &startnode, 0, &dummy, &dummy);

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];
	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < BubbleRadius * BubbleRadius)
	    {
	      if(P[j].Type == 0)
		{
		  numngb++;

		  if(All.ComovingIntegrationOn)
		    E_bubble +=
		      SphP[j].Entropy * P[j].Mass * pow(SphP[j].d.Density / pow(All.Time, 3),
							GAMMA_MINUS1) / GAMMA_MINUS1;
		  else
		    E_bubble +=
		      SphP[j].Entropy * P[j].Mass * pow(SphP[j].d.Density, GAMMA_MINUS1) / GAMMA_MINUS1;

		  Mass_bubble += P[j].Mass;
		}
	    }
	}
    }
  while(startnode >= 0);


  MPI_Allreduce(&numngb, &tot_numngb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&E_bubble, &totE_bubble, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Mass_bubble, &totMass_bubble, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  totE_bubble *= All.UnitEnergy_in_cgs;

  if(totMass_bubble > 0)
    {
      if(ThisTask == 0)
	{
	  printf("found %d particles in bubble with energy %g and total mass %g \n",
		 tot_numngb, totE_bubble, totMass_bubble);
	  fflush(stdout);
	}


      /*calculate comoving density of ICM inside the bubble */

      ICMDensity = totMass_bubble / (4.0 * M_PI / 3.0 * pow(BubbleRadius, 3));

      if(All.ComovingIntegrationOn)
	ICMDensity = ICMDensity / (pow(All.Time, 3));	/*now physical */

      /*Rbub=R0*[(Ejet/Ejet,0)/(rho_ICM/rho_ICM,0)]^(1./5.) - physical */

      rr = rr / BubbleDistance;

      BubbleRadius =
	All.BubbleRadius * pow((BubbleEnergy * All.DefaultICMDensity / (All.BubbleEnergy * ICMDensity)),
			       1. / 5.);

      BubbleDistance =
	All.BubbleDistance * pow((BubbleEnergy * All.DefaultICMDensity / (All.BubbleEnergy * ICMDensity)),
				 1. / 5.);

      if(All.ComovingIntegrationOn)
	{
	  /*switch to comoving if it is assumed that Rbub should be constant with redshift */
	  /* BubbleRadius = BubbleRadius / All.Time;
	     BubbleDistance = BubbleDistance / All.Time; */
	}

      /*recalculate pos */
      rr = rr * BubbleDistance;

      pos[0] = sin(theta) * cos(phi);
      pos[1] = sin(theta) * sin(phi);
      pos[2] = cos(theta);

      for(i = 0; i < 3; i++)
	pos[i] *= rr;

      for(i = 0; i < 3; i++)
	pos[i] += center[i];

      /* now find particles in Bubble again, 
         and recalculate number, mass and energy */

      numngb = 0;
      E_bubble = 0.;
      Mass_bubble = 0.;
      tot_numngb = 0;
      totE_bubble = 0.;
      totMass_bubble = 0.;

      startnode = All.MaxPart;

      do
	{
	  numngb_inbox = ngb_treefind_variable(pos, BubbleRadius, -1, &startnode, 0, &dummy, &dummy);

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];
	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      if(dx > boxHalf_X)
		dx -= boxSize_X;
	      if(dx < -boxHalf_X)
		dx += boxSize_X;
	      if(dy > boxHalf_Y)
		dy -= boxSize_Y;
	      if(dy < -boxHalf_Y)
		dy += boxSize_Y;
	      if(dz > boxHalf_Z)
		dz -= boxSize_Z;
	      if(dz < -boxHalf_Z)
		dz += boxSize_Z;
#endif
	      r2 = dx * dx + dy * dy + dz * dz;

	      if(r2 < BubbleRadius * BubbleRadius)
		{
		  if(P[j].Type == 0 && P[j].Mass > 0)
		    {
		      numngb++;

		      if(All.ComovingIntegrationOn)
			E_bubble +=
			  SphP[j].Entropy * P[j].Mass * pow(SphP[j].d.Density / pow(All.Time, 3),
							    GAMMA_MINUS1) / GAMMA_MINUS1;
		      else
			E_bubble +=
			  SphP[j].Entropy * P[j].Mass * pow(SphP[j].d.Density, GAMMA_MINUS1) / GAMMA_MINUS1;

		      Mass_bubble += P[j].Mass;
		    }
		}
	    }
	}
      while(startnode >= 0);


      MPI_Allreduce(&numngb, &tot_numngb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&E_bubble, &totE_bubble, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&Mass_bubble, &totMass_bubble, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      totE_bubble *= All.UnitEnergy_in_cgs;

      if(totMass_bubble > 0)
	{
	  if(ThisTask == 0)
	    {
	      printf("found %d particles in bubble of rescaled radius with energy %g and total mass %g \n",
		     tot_numngb, totE_bubble, totMass_bubble);
	      printf("energy shall be increased by: (Eini+Einj)/Eini = %g \n",
		     (BubbleEnergy + totE_bubble) / totE_bubble);
	      fflush(stdout);
	    }
	}

      /* now find particles in Bubble again, and inject energy */

#ifdef CR_BUBBLES
      sum_instant_reheat = 0.0;
      tot_instant_reheat = 0.0;
#endif

      startnode = All.MaxPart;

      do
	{
	  numngb_inbox = ngb_treefind_variable(pos, BubbleRadius, -1, &startnode, 0, &dummy, &dummy);

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];
	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      if(dx > boxHalf_X)
		dx -= boxSize_X;
	      if(dx < -boxHalf_X)
		dx += boxSize_X;
	      if(dy > boxHalf_Y)
		dy -= boxSize_Y;
	      if(dy < -boxHalf_Y)
		dy += boxSize_Y;
	      if(dz > boxHalf_Z)
		dz -= boxSize_Z;
	      if(dz < -boxHalf_Z)
		dz += boxSize_Z;
#endif
	      r2 = dx * dx + dy * dy + dz * dz;

	      if(r2 < BubbleRadius * BubbleRadius)
		{
		  if(P[j].Type == 0 && P[j].Mass > 0)
		    {
		      /* energy we want to inject in this particle */

		      if(All.StarformationOn)
			dE = ((BubbleEnergy / All.UnitEnergy_in_cgs) / totMass_bubble) * P[j].Mass;
		      else
			dE = (BubbleEnergy / All.UnitEnergy_in_cgs) / tot_numngb;

		      if(u_to_temp_fac * dE / P[j].Mass > 5.0e9)
			dE = 5.0e9 * P[j].Mass / u_to_temp_fac;

#ifndef CR_BUBBLES
		      if(All.ComovingIntegrationOn)
			SphP[j].Entropy +=
			  GAMMA_MINUS1 * dE / P[j].Mass / pow(SphP[j].d.Density / pow(All.Time, 3),
							      GAMMA_MINUS1);
		      else
			SphP[j].Entropy +=
			  GAMMA_MINUS1 * dE / P[j].Mass / pow(SphP[j].d.Density, GAMMA_MINUS1);
#else

		      if(All.ComovingIntegrationOn)
			tinj = 10.0 * All.HubbleParam * hubble_a / All.UnitTime_in_Megayears;
		      else
			tinj = 10.0 * All.HubbleParam / All.UnitTime_in_Megayears;

		      instant_reheat =
			CR_Particle_SupernovaFeedback(&SphP[j], dE / P[j].Mass * All.CR_AGNEff, tinj);

		      if(instant_reheat > 0)
			{
			  if(All.ComovingIntegrationOn)
			    SphP[j].Entropy +=
			      instant_reheat * GAMMA_MINUS1 / pow(SphP[j].d.Density / pow(All.Time, 3),
								  GAMMA_MINUS1);
			  else
			    SphP[j].Entropy +=
			      instant_reheat * GAMMA_MINUS1 / pow(SphP[j].d.Density, GAMMA_MINUS1);
			}

		      if(All.CR_AGNEff < 1)
			{
			  if(All.ComovingIntegrationOn)
			    SphP[j].Entropy +=
			      (1 -
			       All.CR_AGNEff) * dE * GAMMA_MINUS1 / P[j].Mass / pow(SphP[j].d.Density /
										    pow(All.Time, 3),
										    GAMMA_MINUS1);
			  else
			    SphP[j].Entropy +=
			      (1 - All.CR_AGNEff) * dE * GAMMA_MINUS1 / P[j].Mass / pow(SphP[j].d.Density,
											GAMMA_MINUS1);
			}


		      sum_instant_reheat += instant_reheat * P[j].Mass;
#endif

		    }
		}
	    }
	}
      while(startnode >= 0);

#ifdef CR_BUBBLES
      MPI_Allreduce(&sum_instant_reheat, &tot_instant_reheat, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if(ThisTask == 0)
	{
	  printf("Total BubbleEnergy %g Thermalized Energy %g \n", BubbleEnergy,
		 tot_instant_reheat * All.UnitEnergy_in_cgs);
	  fflush(stdout);

	}
#endif
    }
  else
    {
      if(ThisTask == 0)
	{
	  printf("No particles in bubble found! \n");
	  fflush(stdout);
	}

    }

}
#endif /* end of BH_BUBBLE */

#endif
