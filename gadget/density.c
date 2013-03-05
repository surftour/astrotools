#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif
#ifdef CS_MODEL
#include "cs_metals.h"
#endif


/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct densdata_in
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
#ifdef VOLUME_CORRECTION
  MyFloat DensityOld;
#endif
#ifdef WINDS
  MyFloat DelayTime;
#endif
  int NodeList[NODELISTLENGTH];
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS) || defined(TRACEDIVB)
  MyFloat BPred[3];
#endif
#ifdef EULERPOTENTIALS
  MyFloat EulerA, EulerB;
#endif
#ifdef CS_MODEL
  MyFloat DensityOld, Entropy;
#endif
}
 *DensDataIn, *DensDataGet;


static struct densdata_out
{
  MyLongDouble Rho;
  MyLongDouble DhsmlDensity;
  MyLongDouble Ngb;
#ifndef NAVIERSTOKES
  MyLongDouble Div, Rot[3];
#else
  MyFloat DV[3][3];
#endif
#ifdef MAGNETIC
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
  MyFloat RotB[3];
#endif
#ifdef TRACEDIVB
  MyFloat divB;
#endif
#ifdef VECT_PRO_CLEAN
  MyFloat BPredVec[3];
#endif
#endif

#ifdef RADTRANSFER_FLUXLIMITER
  MyFloat Grad_ngamma[3];
#endif

#if defined(BLACK_HOLES)
  MyLongDouble SmoothedEntr;
#endif
#ifdef CONDUCTION_SATURATION
  MyFloat GradEntr[3];
#endif

#ifdef BLACK_HOLES
  MyLongDouble GasVel[3];
#endif

#ifdef HYDRO_COST_FACTOR
  int Ninteractions;
#endif

#ifdef VOLUME_CORRECTION
  MyFloat DensityStd;
#endif

#ifdef EULERPOTENTIALS
  MyFloat dEulerA[3], dEulerB[3];
#endif

}
 *DensDataResult, *DensDataOut;


/*! \file density.c 
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  corrects the smoothing length if needed.
 */


/*! This function computes the local density for each active SPH particle, the
 * number of neighbours in the current smoothing radius, and the divergence
 * and rotation of the velocity field.  The pressure is updated as well.  If a
 * particle with its smoothing region is fully inside the local domain, it is
 * not exported to the other processors. The function also detects particles
 * that have a number of neighbours outside the allowed tolerance range. For
 * these particles, the smoothing length is adjusted accordingly, and the
 * density() computation is called again.  Note that the smoothing length is
 * not allowed to fall below the lower bound set by MinGasHsml (this may mean
 * that one has to deal with substantially more than normal number of
 * neighbours.)
 */
void density(void)
{
  MyFloat *Left, *Right;
  int i, j, ndone, ndone_flag, npleft, dt_step, dummy, iter = 0;
  int ngrp, sendTask, recvTask, place, nexport, nimport;
  long long ntot;
  double dmax1, dmax2, fac;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double dt_entr, tstart, tend, t0, t1;
  double desnumngb;

#ifdef COSMIC_RAYS
  int CRpop;
#endif

#if defined(NAVIERSTOKES)
  int k;
#endif
#ifdef NAVIERSTOKES
  double dvel[3][3];
  double rotx, roty, rotz;
#endif

#if defined(SOFTEREQS)
  double a3inv;

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;
#endif

#if defined(EULERPOTENTIALS) || defined(VECT_PRO_CLEAN)
  double efak;

  if(All.ComovingIntegrationOn)
    efak = 1. / All.Time / All.HubbleParam;
  else
    efak = 1;
#endif

  CPU_Step[CPU_DENSMISC] += measure_time();

  Left = (MyFloat *) mymalloc(NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc(NumPart * sizeof(MyFloat));

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(density_isactive(i))
	{
	  Left[i] = Right[i] = 0;

#ifdef BLACK_HOLES
	  P[i].SwallowID = 0;
#endif
#if defined(BLACK_HOLES) && defined(FLTROUNDOFFREDUCTION)
	  if(P[i].Type == 0)
	    SphP[i].i.dInjected_BH_Energy = SphP[i].i.Injected_BH_Energy;
#endif
	}
    }

  /* allocate buffers to arrange communication */


  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct densdata_in) + sizeof(struct densdata_out) +
					     sizemax(sizeof(struct densdata_in),
						     sizeof(struct densdata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));

  t0 = second();

  desnumngb = All.DesNumNgb;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      i = FirstActiveParticle;	/* begin with this index */

      do
	{
	  for(j = 0; j < NTask; j++)
	    {
	      Send_count[j] = 0;
	      Exportflag[j] = -1;
	    }

	  /* do local particles and prepare export list */
	  tstart = second();
	  for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	    {
	      if(density_isactive(i))
		{
		  if(density_evaluate(i, 0, &nexport, Send_count) < 0)
		    break;
		}
	    }
	  tend = second();
	  timecomp1 += timediff(tstart, tend);

#ifdef MYSORT
	  mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
	  qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

	  tstart = second();

	  MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);

	  tend = second();
	  timewait1 += timediff(tstart, tend);

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

	  DensDataGet = (struct densdata_in *) mymalloc(nimport * sizeof(struct densdata_in));
	  DensDataIn = (struct densdata_in *) mymalloc(nexport * sizeof(struct densdata_in));

	  /* prepare particle data for export */
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      DensDataIn[j].Pos[0] = P[place].Pos[0];
	      DensDataIn[j].Pos[1] = P[place].Pos[1];
	      DensDataIn[j].Pos[2] = P[place].Pos[2];
	      DensDataIn[j].Hsml = PPP[place].Hsml;

#ifdef CS_MODEL
	      DensDataIn[j].DensityOld = SphP[place].DensityOld;
	      DensDataIn[j].Entropy = SphP[place].Entropy;
#endif
	      memcpy(DensDataIn[j].NodeList,
		     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));

#if defined(BLACK_HOLES)
	      if(P[place].Type != 0)
		{
		  DensDataIn[j].Vel[0] = 0;
		  DensDataIn[j].Vel[1] = 0;
		  DensDataIn[j].Vel[2] = 0;
		}
	      else
#endif
		{
		  DensDataIn[j].Vel[0] = SphP[place].VelPred[0];
		  DensDataIn[j].Vel[1] = SphP[place].VelPred[1];
		  DensDataIn[j].Vel[2] = SphP[place].VelPred[2];
		}
#ifdef VOLUME_CORRECTION
	      DensDataIn[j].DensityOld = SphP[place].DensityOld;
#endif

#ifdef EULERPOTENTIALS
	      DensDataIn[j].EulerA = SphP[place].EulerA;
	      DensDataIn[j].EulerB = SphP[place].EulerB;
#endif


#ifdef WINDS
	      DensDataIn[j].DelayTime = SphP[place].DelayTime;
#endif

#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS) || defined(TRACEDIVB)
	      DensDataIn[j].BPred[0] = SphP[place].BPred[0];
	      DensDataIn[j].BPred[1] = SphP[place].BPred[1];
	      DensDataIn[j].BPred[2] = SphP[place].BPred[2];
#endif
	    }
	  /* exchange particle data */
	  tstart = second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&DensDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &DensDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }
	  tend = second();
	  timecommsumm1 += timediff(tstart, tend);

	  myfree(DensDataIn);
	  DensDataResult = (struct densdata_out *) mymalloc(nimport * sizeof(struct densdata_out));
	  DensDataOut = (struct densdata_out *) mymalloc(nexport * sizeof(struct densdata_out));


	  /* now do the particles that were sent to us */

	  tstart = second();
	  for(j = 0; j < nimport; j++)
	    density_evaluate(j, 1, &dummy, &dummy);
	  tend = second();
	  timecomp2 += timediff(tstart, tend);

	  if(i < 0)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

	  tstart = second();
	  MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  tend = second();
	  timewait2 += timediff(tstart, tend);


	  /* get the result */
	  tstart = second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&DensDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &DensDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}

	    }
	  tend = second();
	  timecommsumm2 += timediff(tstart, tend);


	  /* add the result to the local particles */
	  tstart = second();
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      PPP[place].n.dNumNgb += DensDataOut[j].Ngb;
#ifdef HYDRO_COST_FACTOR
	      P[place].GravCost += HYDRO_COST_FACTOR * DensDataOut[j].Ninteractions;
#endif

	      if(P[place].Type == 0)
		{
		  SphP[place].d.dDensity += DensDataOut[j].Rho;
		  SphP[place].h.dDhsmlDensityFactor += DensDataOut[j].DhsmlDensity;

#ifndef NAVIERSTOKES
		  SphP[place].v.dDivVel += DensDataOut[j].Div;
		  SphP[place].r.dRot[0] += DensDataOut[j].Rot[0];
		  SphP[place].r.dRot[1] += DensDataOut[j].Rot[1];
		  SphP[place].r.dRot[2] += DensDataOut[j].Rot[2];
#else
		  for(k = 0; k < 3; k++)
		    {
		      SphP[place].u.DV[k][0] += DensDataOut[j].DV[k][0];
		      SphP[place].u.DV[k][1] += DensDataOut[j].DV[k][1];
		      SphP[place].u.DV[k][2] += DensDataOut[j].DV[k][2];
		    }
#endif

#ifdef VOLUME_CORRECTION
		  SphP[place].DensityStd += DensDataOut[j].DensityStd;
#endif

#ifdef CONDUCTION_SATURATION
		  SphP[place].GradEntr[0] += DensDataOut[j].GradEntr[0];
		  SphP[place].GradEntr[1] += DensDataOut[j].GradEntr[1];
		  SphP[place].GradEntr[2] += DensDataOut[j].GradEntr[2];
#endif

#ifdef RADTRANSFER_FLUXLIMITER
		  SphP[place].Grad_ngamma[0] += DensDataOut[j].Grad_ngamma[0];
		  SphP[place].Grad_ngamma[1] += DensDataOut[j].Grad_ngamma[1];
		  SphP[place].Grad_ngamma[2] += DensDataOut[j].Grad_ngamma[2];
#endif


#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
		  SphP[place].RotB[0] += DensDataOut[j].RotB[0];
		  SphP[place].RotB[1] += DensDataOut[j].RotB[1];
		  SphP[place].RotB[2] += DensDataOut[j].RotB[2];
#endif

#ifdef TRACEDIVB
		  SphP[place].divB += DensDataOut[j].divB;
#endif

#ifdef VECT_PRO_CLEAN
		  SphP[place].BPredVec[0] += DensDataOut[j].BPredVec[0];
		  SphP[place].BPredVec[1] += DensDataOut[j].BPredVec[1];
		  SphP[place].BPredVec[2] += DensDataOut[j].BPredVec[2];
		  //DO here the change or COrrection
#endif
#ifdef EULERPOTENTIALS
		  SphP[place].dEulerA[0] += DensDataOut[j].dEulerA[0];
		  SphP[place].dEulerA[1] += DensDataOut[j].dEulerA[1];
		  SphP[place].dEulerA[2] += DensDataOut[j].dEulerA[2];
		  SphP[place].dEulerB[0] += DensDataOut[j].dEulerB[0];
		  SphP[place].dEulerB[1] += DensDataOut[j].dEulerB[1];
		  SphP[place].dEulerB[2] += DensDataOut[j].dEulerB[2];
#endif
		}

#if defined(RADTRANSFER) || defined(SNIA_HEATING)
	      if(P[place].Type == 4)
		P[place].DensAroundStar += DensDataOut[j].Rho;
#endif


#ifdef BLACK_HOLES
	      if(P[place].Type == 5)
		{
		  P[place].b1.dBH_Density += DensDataOut[j].Rho;
		  P[place].b2.dBH_Entropy += DensDataOut[j].SmoothedEntr;
		  P[place].b3.dBH_SurroundingGasVel[0] += DensDataOut[j].GasVel[0];
		  P[place].b3.dBH_SurroundingGasVel[1] += DensDataOut[j].GasVel[1];
		  P[place].b3.dBH_SurroundingGasVel[2] += DensDataOut[j].GasVel[2];
		}
#endif
	    }
	  tend = second();
	  timecomp1 += timediff(tstart, tend);


	  myfree(DensDataOut);
	  myfree(DensDataResult);
	  myfree(DensDataGet);
	}
      while(ndone < NTask);

#ifdef FLTROUNDOFFREDUCTION
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	if(density_isactive(i))
	  {
	    PPP[i].n.NumNgb = FLT(PPP[i].n.dNumNgb);

	    if(P[i].Type == 0)
	      {
		SphP[i].d.Density = FLT(SphP[i].d.dDensity);
		SphP[i].h.DhsmlDensityFactor = FLT(SphP[i].h.dDhsmlDensityFactor);
		SphP[i].v.DivVel = FLT(SphP[i].v.dDivVel);
		for(j = 0; j < 3; j++)
		  SphP[i].r.Rot[j] = FLT(SphP[i].r.dRot[j]);
	      }

#ifdef BLACK_HOLES
	    if(P[i].Type == 5)
	      {
		P[i].b1.BH_Density = FLT(P[i].b1.dBH_Density);
		P[i].b2.BH_Entropy = FLT(P[i].b2.dBH_Entropy);
		for(j = 0; j < 3; j++)
		  P[i].b3.BH_SurroundingGasVel[j] = FLT(P[i].b3.dBH_SurroundingGasVel[j]);
	      }
#endif
	  }
#endif


      /* do final operations on results */
      tstart = second();
      for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
	{
	  if(density_isactive(i))
	    {
	      if(P[i].Type == 0)
		{
		  if(SphP[i].d.Density > 0)
		    {
#ifdef VOLUME_CORRECTION
		      SphP[i].DensityOld = SphP[i].DensityStd;
#endif
		      SphP[i].h.DhsmlDensityFactor *= PPP[i].Hsml / (NUMDIMS * SphP[i].d.Density);
		      if(SphP[i].h.DhsmlDensityFactor > -0.9)	/* note: this would be -1 if only a single particle at zero lag is found */
			SphP[i].h.DhsmlDensityFactor = 1 / (1 + SphP[i].h.DhsmlDensityFactor);
		      else
			SphP[i].h.DhsmlDensityFactor = 1;
#ifndef NAVIERSTOKES
		      SphP[i].r.CurlVel = sqrt(SphP[i].r.Rot[0] * SphP[i].r.Rot[0] +
					       SphP[i].r.Rot[1] * SphP[i].r.Rot[1] +
					       SphP[i].r.Rot[2] * SphP[i].r.Rot[2]) / SphP[i].d.Density;

		      SphP[i].v.DivVel /= SphP[i].d.Density;
#else
		      for(k = 0; k < 3; k++)
			{
			  dvel[k][0] = SphP[i].u.DV[k][0] / SphP[i].d.Density;
			  dvel[k][1] = SphP[i].u.DV[k][1] / SphP[i].d.Density;
			  dvel[k][2] = SphP[i].u.DV[k][2] / SphP[i].d.Density;
			}
		      SphP[i].u.s.DivVel = dvel[0][0] + dvel[1][1] + dvel[2][2];

		      SphP[i].u.s.StressDiag[0] = 2 * dvel[0][0] - 2.0 / 3 * SphP[i].u.s.DivVel;
		      SphP[i].u.s.StressDiag[1] = 2 * dvel[1][1] - 2.0 / 3 * SphP[i].u.s.DivVel;
		      SphP[i].u.s.StressDiag[2] = 2 * dvel[2][2] - 2.0 / 3 * SphP[i].u.s.DivVel;

		      SphP[i].u.s.StressOffDiag[0] = dvel[0][1] + dvel[1][0];	/* xy */
		      SphP[i].u.s.StressOffDiag[1] = dvel[0][2] + dvel[2][0];	/* xz */
		      SphP[i].u.s.StressOffDiag[2] = dvel[1][2] + dvel[2][1];	/* yz */

#ifdef NAVIERSTOKES_BULK
		      SphP[i].u.s.StressBulk = All.NavierStokes_BulkViscosity * SphP[i].u.s.DivVel;
#endif
		      rotx = dvel[1][2] - dvel[2][1];
		      roty = dvel[2][0] - dvel[0][2];
		      rotz = dvel[0][1] - dvel[1][0];
		      SphP[i].u.s.CurlVel = sqrt(rotx * rotx + roty * roty + rotz * rotz);
#endif


#ifdef CONDUCTION_SATURATION
		      SphP[i].GradEntr[0] /= SphP[i].d.Density;
		      SphP[i].GradEntr[1] /= SphP[i].d.Density;
		      SphP[i].GradEntr[2] /= SphP[i].d.Density;
#endif


#ifdef RADTRANSFER_FLUXLIMITER
		      SphP[i].Grad_ngamma[0] /= SphP[i].d.Density;
		      SphP[i].Grad_ngamma[1] /= SphP[i].d.Density;
		      SphP[i].Grad_ngamma[2] /= SphP[i].d.Density;
#endif


#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
		      SphP[i].RotB[0] /= SphP[i].d.Density;
		      SphP[i].RotB[1] /= SphP[i].d.Density;
		      SphP[i].RotB[2] /= SphP[i].d.Density;
#endif

#ifdef TRACEDIVB
		      SphP[i].divB /= SphP[i].d.Density;
#endif






#ifdef VECT_PRO_CLEAN
		      // dooo the changes to Bpred
		      SphP[i].BPred[0] += efak * SphP[i].BPredVec[0];	//SphP[i].d.Density; 
		      SphP[i].BPred[1] += efak * SphP[i].BPredVec[1];	//SphP[i].d.Density; 
		      SphP[i].BPred[2] += efak * SphP[i].BPredVec[2];	//SphP[i].d.Density;
#endif
#ifdef EULERPOTENTIALS
		      SphP[i].dEulerA[0] *= efak / SphP[i].d.Density;
		      SphP[i].dEulerA[1] *= efak / SphP[i].d.Density;
		      SphP[i].dEulerA[2] *= efak / SphP[i].d.Density;
		      SphP[i].dEulerB[0] *= efak / SphP[i].d.Density;
		      SphP[i].dEulerB[1] *= efak / SphP[i].d.Density;
		      SphP[i].dEulerB[2] *= efak / SphP[i].d.Density;
#ifdef	PERIODIC
		      SphP[i].dEulerB[2] = 1.0;
#endif
		      SphP[i].BPred[0] =
			SphP[i].dEulerA[1] * SphP[i].dEulerB[2] - SphP[i].dEulerA[2] * SphP[i].dEulerB[1];
		      SphP[i].BPred[1] =
			SphP[i].dEulerA[2] * SphP[i].dEulerB[0] - SphP[i].dEulerA[0] * SphP[i].dEulerB[2];
		      SphP[i].BPred[2] =
			SphP[i].dEulerA[0] * SphP[i].dEulerB[1] - SphP[i].dEulerA[1] * SphP[i].dEulerB[0];
#endif
		    }

#ifndef WAKEUP
		  dt_step = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);
#else
		  dt_step = P[i].dt_step;
#endif
		  dt_entr = (All.Ti_Current - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;

#ifndef EOS_DEGENERATE
#ifndef MHM
#ifndef SOFTEREQS
#ifndef TRADITIONAL_SPH_FORMULATION
		  SphP[i].Pressure =
		    (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr) * pow(SphP[i].d.Density, GAMMA);
#else
		  SphP[i].Pressure =
		    GAMMA_MINUS1 * (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr) * SphP[i].d.Density;
#endif

#else
		  /* use an intermediate EQS, between isothermal and the full multiphase model */
		  if(SphP[i].d.Density * a3inv >= All.PhysDensThresh)
		    SphP[i].Pressure = All.FactorForSofterEQS *
		      (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr) * pow(SphP[i].d.Density, GAMMA) +
		      (1 - All.FactorForSofterEQS) * GAMMA_MINUS1 * SphP[i].d.Density * All.InitGasU;
		  else
		    SphP[i].Pressure =
		      (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr) * pow(SphP[i].d.Density, GAMMA);
#endif
#else
		  /* Here we use an isothermal equation of state */
		  SphP[i].Pressure = GAMMA_MINUS1 * SphP[i].d.Density * All.InitGasU;
		  SphP[i].Entropy = SphP[i].Pressure / pow(SphP[i].d.Density, GAMMA);
#endif
#else
		  /* call tabulated eos with physical units */
		  eos_calc_egiven_v(SphP[i].d.Density * All.UnitDensity_in_cgs, SphP[i].xnuc, SphP[i].dxnuc,
				    dt_entr * All.UnitTime_in_s,
				    SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr * All.UnitTime_in_s,
				    &SphP[i].temp, &SphP[i].Pressure, &SphP[i].dpdr);
		  SphP[i].Pressure /= All.UnitPressure_in_cgs;
#endif

#ifdef COSMIC_RAYS
		  for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
		    {
		      CR_Particle_Update(SphP + i, CRpop);
#ifndef CR_NOPRESSURE
		      SphP[i].Pressure += CR_Comoving_Pressure(SphP + i, CRpop);
#endif
		    }
#endif
		}

#ifdef BLACK_HOLES
	      if(P[i].Type == 5)
		{
		  if(P[i].b1.BH_Density > 0)
		    {
		      P[i].b2.BH_Entropy /= P[i].b1.BH_Density;
		      P[i].b3.BH_SurroundingGasVel[0] /= P[i].b1.BH_Density;
		      P[i].b3.BH_SurroundingGasVel[1] /= P[i].b1.BH_Density;
		      P[i].b3.BH_SurroundingGasVel[2] /= P[i].b1.BH_Density;
		    }
		}
#endif

	      /* now check whether we had enough neighbours */

	      desnumngb = All.DesNumNgb;

#ifdef BLACK_HOLES
	      if(P[i].Type == 5)
		desnumngb = All.DesNumNgb * All.BlackHoleNgbFactor;
#endif

#ifdef RADTRANSFER
	      if(P[i].Type == 4)
		desnumngb = 64;	//NORM_COEFF * KERNEL_COEFF_1;   /* will assign the stellar luminosity to very few (one actually) gas particles */
#endif

	      if(PPP[i].n.NumNgb < (desnumngb - All.MaxNumNgbDeviation) ||
		 (PPP[i].n.NumNgb > (desnumngb + All.MaxNumNgbDeviation)
		  && PPP[i].Hsml > (1.01 * All.MinGasHsml)))
		{
		  /* need to redo this particle */
		  npleft++;

		  if(Left[i] > 0 && Right[i] > 0)
		    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
		      {
			/* this one should be ok */
			npleft--;
			P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
			continue;
		      }

		  if(PPP[i].n.NumNgb < (desnumngb - All.MaxNumNgbDeviation))
		    Left[i] = DMAX(PPP[i].Hsml, Left[i]);
		  else
		    {
		      if(Right[i] != 0)
			{
			  if(PPP[i].Hsml < Right[i])
			    Right[i] = PPP[i].Hsml;
			}
		      else
			Right[i] = PPP[i].Hsml;
		    }

		  if(iter >= MAXITER - 10)
		    {
		      printf
			("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (int) P[i].ID, PPP[i].Hsml, Left[i], Right[i],
			 (float) PPP[i].n.NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
		      fflush(stdout);
		    }

		  if(Right[i] > 0 && Left[i] > 0)
		    PPP[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
		  else
		    {
		      if(Right[i] == 0 && Left[i] == 0)
			endrun(8188);	/* can't occur */

		      if(Right[i] == 0 && Left[i] > 0)
			{
			  if(P[i].Type == 0 && fabs(PPP[i].n.NumNgb - desnumngb) < 0.5 * desnumngb)
			    {
			      fac = 1 - (PPP[i].n.NumNgb -
					 desnumngb) / (NUMDIMS * PPP[i].n.NumNgb) *
				SphP[i].h.DhsmlDensityFactor;

			      if(fac < 1.26)
				PPP[i].Hsml *= fac;
			      else
				PPP[i].Hsml *= 1.26;
			    }
			  else
			    PPP[i].Hsml *= 1.26;
			}

		      if(Right[i] > 0 && Left[i] == 0)
			{
			  if(P[i].Type == 0 && fabs(PPP[i].n.NumNgb - desnumngb) < 0.5 * desnumngb)
			    {
			      fac = 1 - (PPP[i].n.NumNgb -
					 desnumngb) / (NUMDIMS * PPP[i].n.NumNgb) *
				SphP[i].h.DhsmlDensityFactor;

			      if(fac > 1 / 1.26)
				PPP[i].Hsml *= fac;
			      else
				PPP[i].Hsml /= 1.26;
			    }
			  else
			    PPP[i].Hsml /= 1.26;
			}
		    }

		  if(PPP[i].Hsml < All.MinGasHsml)
		    PPP[i].Hsml = All.MinGasHsml;

#ifdef BLACK_HOLES
		  if(P[i].Type == 5)
		    if(Left[i] > All.BlackHoleMaxAccretionRadius)
		      {
			/* this will stop the search for a new BH smoothing length in the next iteration */
			PPP[i].Hsml = Left[i] = Right[i] = All.BlackHoleMaxAccretionRadius;
		      }
#endif
		}
	      else
		P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
	    }
	}
      tend = second();
      timecomp1 += timediff(tstart, tend);

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
	{
	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
    }
  while(ntot > 0);


  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
  myfree(Right);
  myfree(Left);



#if defined(CS_MODEL) && defined(CS_FEEDBACK)
  double xhyd, yhel, ne, mu, energy, temp, a3inv;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0 && P[i].EnergySN < 0)	/* only gas particles should enter here */
	{

	  printf("prom 2i=%d type=%d energySN%g\n", i, P[i].Type, P[i].EnergySN);
	  fflush(0);

	  P[i].EnergySN = 0;

	}

      if(All.ComovingIntegrationOn)
	a3inv = 1 / (All.Time * All.Time * All.Time);
      else
	a3inv = 1;

      if(P[i].Type == 0 && (SphP[i].TempPromotion > 0 || SphP[i].DensPromotion > 0))
	{
	  xhyd = P[i].Zm[6] / P[i].Mass;
	  yhel = (1 - xhyd) / (4. * xhyd);
	  ne = SphP[i].Ne;
	  mu = (1 + 4 * yhel) / (1 + yhel + ne);
	  energy = SphP[i].Entropy * P[i].Mass / GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);	/* Total Energys */
	  temp = GAMMA_MINUS1 / BOLTZMANN * energy / P[i].Mass * PROTONMASS * mu;
	  temp *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;	/* Temperature in Kelvin */

	  fprintf(FdPromotion, "%g %d %g %g %g %g %g %g\n", All.Time, P[i].ID, SphP[i].DensPromotion,
		  SphP[i].TempPromotion, SphP[i].d.Density, temp, SphP[i].da.DensityAvg,
		  SphP[i].ea.EntropyAvg);
	  fflush(FdPromotion);
	}
    }
#endif

  /* mark as active again */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].TimeBin < 0)
      P[i].TimeBin = -P[i].TimeBin - 1;

  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_DENSCOMPUTE] += timecomp;
  CPU_Step[CPU_DENSWAIT] += timewait;
  CPU_Step[CPU_DENSCOMM] += timecomm;
  CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm);
}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int density_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;
  int startnode, numngb, numngb_inbox, listindex = 0;
  double h, h2, fac, hinv, hinv3, hinv4;
  MyLongDouble rho;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  double dvx, dvy, dvz;
  MyLongDouble weighted_numngb;
  MyLongDouble dhsmlrho;

#ifdef HYDRO_COST_FACTOR
  int ninteractions = 0;
#endif

#ifdef CS_MODEL
  double densityold, entropy = 0;
#endif
#ifdef BLACK_HOLES
  MyLongDouble gasvel[3];
#endif
#ifndef NAVIERSTOKES
  MyLongDouble divv, rotv[3];
#else
  int k;
  double dvel[3][3];
#endif

  MyDouble *pos;
  MyFloat *vel;
  static MyFloat veldummy[3] = { 0, 0, 0 };

#ifdef EULERPOTENTIALS
  MyDouble deulera[3], deulerb[3], eulera, eulerb, dea, deb;

  deulera[0] = deulera[1] = deulera[2] = deulerb[0] = deulerb[1] = deulerb[2] = 0;
#endif

#if defined(BLACK_HOLES)
  MyLongDouble smoothentr;

  smoothentr = 0;
#endif

#ifdef WINDS
  double delaytime;
#endif

#ifdef VOLUME_CORRECTION
  double densityold, densitystd = 0;
#endif

#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS) || defined(TRACEDIVB)
  MyFloat *bflt;
  MyDouble dbx, dby, dbz;
#endif

#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
  MyDouble rotb[3];
#endif

#ifdef TRACEDIVB
  double divB = 0;
#endif
#ifdef VECT_PRO_CLEAN
  double BVec[3];

  BVec[0] = BVec[1] = BVec[2] = 0;
#endif

#ifdef CONDUCTION_SATURATION
  double gradentr[3];

  gradentr[0] = gradentr[1] = gradentr[2] = 0;
#endif

#ifdef RADTRANSFER_FLUXLIMITER
  double grad_ngamma[3];

  grad_ngamma[0] = grad_ngamma[1] = grad_ngamma[2] = 0;
#endif


#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
  rotb[0] = rotb[1] = rotb[2] = 0;
#endif

#ifndef NAVIERSTOKES
  divv = rotv[0] = rotv[1] = rotv[2] = 0;
#else
  for(k = 0; k < 3; k++)
    dvel[k][0] = dvel[k][1] = dvel[k][2] = 0;
#endif
#ifdef BLACK_HOLES
  gasvel[0] = gasvel[1] = gasvel[2] = 0;
#endif
  rho = weighted_numngb = dhsmlrho = 0;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
#ifdef VOLUME_CORRECTION
      densityold = SphP[target].DensityOld;
#endif
      if(P[target].Type == 0)
	{
	  vel = SphP[target].VelPred;
#ifdef WINDS
	  delaytime = SphP[target].DelayTime;
#endif
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS) || defined(TRACEDIVB)
	  bflt = SphP[target].BPred;
#endif
#ifdef EULERPOTENTIALS
	  eulera = SphP[target].EulerA;
	  eulerb = SphP[target].EulerB;
#endif
	}
      else
	vel = veldummy;
#ifdef CS_MODEL
      if(P[target].Type == 0)
	{
	  densityold = SphP[target].DensityOld;
	  entropy = SphP[target].Entropy;
	}
      else
	densityold = 0;
#endif
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
#ifdef WINDS
      delaytime = DensDataGet[target].DelayTime;
#endif
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS) || defined(TRACEDIVB)
      bflt = DensDataGet[target].BPred;
#endif
#ifdef VOLUME_CORRECTION
      densityold = DensDataGet[target].DensityOld;
#endif
#ifdef EULERPOTENTIALS
      eulera = DensDataGet[target].EulerA;
      eulerb = DensDataGet[target].EulerB;
#endif
#ifdef CS_MODEL
      densityold = DensDataGet[target].DensityOld;
      entropy = DensDataGet[target].Entropy;
#endif
    }


  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;


  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = DensDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  numngb = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
#ifdef CS_MODEL
	  numngb_inbox = cs_ngb_treefind_variable_decoupling(&pos[0], h, target, &startnode,
							     densityold, entropy, &vel[0], mode, nexport,
							     nsend_local);
#else
	  numngb_inbox = ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, nsend_local);
#endif

	  if(numngb_inbox < 0)
	    return -1;

	  for(n = 0; n < numngb_inbox; n++)
	    {
#ifdef HYDRO_COST_FACTOR
	      ninteractions++;
#endif
	      j = Ngblist[n];
#ifdef WINDS
	      if(SphP[j].DelayTime > 0)	/* partner is a wind particle */
		if(!(delaytime > 0))	/* if I'm not wind, then ignore the wind particle */
		  continue;
#endif
#ifdef BLACK_HOLES
	      if(P[j].Mass == 0)
		continue;
#endif
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

	      if(r2 < h2)
		{
		  numngb++;

		  r = sqrt(r2);

		  u = r * hinv;

		  if(u < 0.5)
		    {
		      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		      dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		    }
		  else
		    {
		      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		      dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		    }

		  mass_j = P[j].Mass;

#ifdef VOLUME_CORRECTION
		  rho += FLT(mass_j * wk * pow(densityold / SphP[j].DensityOld, VOLUME_CORRECTION));
		  densitystd += FLT(mass_j * wk);
#else
		  rho += FLT(mass_j * wk);
#endif

		  weighted_numngb += FLT(NORM_COEFF * wk / hinv3);	/* 4.0/3 * PI = 4.188790204786 */

		  dhsmlrho += FLT(-mass_j * (NUMDIMS * hinv * wk + u * dwk));


#ifdef BLACK_HOLES
		  smoothentr += FLT(mass_j * wk * SphP[j].Entropy);
		  gasvel[0] += FLT(mass_j * wk * SphP[j].VelPred[0]);
		  gasvel[1] += FLT(mass_j * wk * SphP[j].VelPred[1]);
		  gasvel[2] += FLT(mass_j * wk * SphP[j].VelPred[2]);
#endif

#ifdef CONDUCTION_SATURATION
		  if(r > 0)
		    {
		      gradentr[0] += mass_j * dwk * dx / r * SphP[j].Entropy;
		      gradentr[1] += mass_j * dwk * dy / r * SphP[j].Entropy;
		      gradentr[2] += mass_j * dwk * dz / r * SphP[j].Entropy;
		    }
#endif


#ifdef RADTRANSFER_FLUXLIMITER
		  if(r > 0)
		    {
		      grad_ngamma[0] += mass_j * dwk * dx / r * SphP[j].n_gamma;
		      grad_ngamma[1] += mass_j * dwk * dy / r * SphP[j].n_gamma;
		      grad_ngamma[2] += mass_j * dwk * dz / r * SphP[j].n_gamma;
		    }
#endif


		  if(r > 0)
		    {
		      fac = mass_j * dwk / r;

		      dvx = vel[0] - SphP[j].VelPred[0];
		      dvy = vel[1] - SphP[j].VelPred[1];
		      dvz = vel[2] - SphP[j].VelPred[2];

#ifndef NAVIERSTOKES
		      divv += FLT(-fac * (dx * dvx + dy * dvy + dz * dvz));

		      rotv[0] += FLT(fac * (dz * dvy - dy * dvz));
		      rotv[1] += FLT(fac * (dx * dvz - dz * dvx));
		      rotv[2] += FLT(fac * (dy * dvx - dx * dvy));
#else
		      dvel[0][0] -= fac * dx * dvx;
		      dvel[0][1] -= fac * dx * dvy;
		      dvel[0][2] -= fac * dx * dvz;
		      dvel[1][0] -= fac * dy * dvx;
		      dvel[1][1] -= fac * dy * dvy;
		      dvel[1][2] -= fac * dy * dvz;
		      dvel[2][0] -= fac * dz * dvx;
		      dvel[2][1] -= fac * dz * dvy;
		      dvel[2][2] -= fac * dz * dvz;
#endif
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS) || defined(TRACEDIVB)
		      dbx = bflt[0] - SphP[j].BPred[0];
		      dby = bflt[1] - SphP[j].BPred[1];
		      dbz = bflt[2] - SphP[j].BPred[2];
#endif
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
		      rotb[0] += FLT(fac * (dz * dby - dy * dbz));
		      rotb[1] += FLT(fac * (dx * dbz - dz * dbx));
		      rotb[2] += FLT(fac * (dy * dbx - dx * dby));
#endif
#ifdef TRACEDIVB
		      divB += FLT(-fac * (dbx * dx + dby * dy + dbz * dz));
#endif
#ifdef VECT_PRO_CLEAN
		      BVec[0] +=
			FLT(fac * r2 * (SphP[j].RotB[1] * dz - SphP[j].RotB[2] * dy) / SphP[j].d.Density);
		      BVec[1] +=
			FLT(fac * r2 * (SphP[j].RotB[2] * dx - SphP[j].RotB[0] * dz) / SphP[j].d.Density);
		      BVec[2] +=
			FLT(fac * r2 * (SphP[j].RotB[0] * dy - SphP[j].RotB[1] * dx) / SphP[j].d.Density);
#endif
#ifdef EULERPOTENTIALS
		      dea = eulera - SphP[j].EulerA;
		      deb = eulerb - SphP[j].EulerB;
		      deulera[0] -= fac * dx * dea;
		      deulera[1] -= fac * dy * dea;
		      deulera[2] -= fac * dz * dea;
		      deulerb[0] -= fac * dx * deb;
		      deulerb[1] -= fac * dy * deb;
		      deulerb[2] -= fac * dz * deb;
#endif
		    }
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = DensDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    {
      PPP[target].n.dNumNgb = weighted_numngb;
#ifdef HYDRO_COST_FACTOR
      P[target].GravCost += HYDRO_COST_FACTOR * ninteractions;
#endif
      if(P[target].Type == 0)
	{
	  SphP[target].d.dDensity = rho;
#ifdef VOLUME_CORRECTION
	  SphP[target].DensityStd = densitystd;
#endif
	  SphP[target].h.dDhsmlDensityFactor = dhsmlrho;
#ifndef NAVIERSTOKES
	  SphP[target].v.dDivVel = divv;
	  SphP[target].r.dRot[0] = rotv[0];
	  SphP[target].r.dRot[1] = rotv[1];
	  SphP[target].r.dRot[2] = rotv[2];
#else
	  for(k = 0; k < 3; k++)
	    {
	      SphP[target].u.DV[k][0] = dvel[k][0];
	      SphP[target].u.DV[k][1] = dvel[k][1];
	      SphP[target].u.DV[k][2] = dvel[k][2];
	    }
#endif

#ifdef CONDUCTION_SATURATION
	  SphP[target].GradEntr[0] = gradentr[0];
	  SphP[target].GradEntr[1] = gradentr[1];
	  SphP[target].GradEntr[2] = gradentr[2];
#endif

#ifdef RADTRANSFER_FLUXLIMITER
	  SphP[target].Grad_ngamma[0] = grad_ngamma[0];
	  SphP[target].Grad_ngamma[1] = grad_ngamma[1];
	  SphP[target].Grad_ngamma[2] = grad_ngamma[2];
#endif

#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
	  SphP[target].RotB[0] = rotb[0];
	  SphP[target].RotB[1] = rotb[1];
	  SphP[target].RotB[2] = rotb[2];
#endif
#ifdef TRACEDIVB
	  SphP[target].divB = divB;
#endif
#ifdef EULERPOTENTIALS
	  SphP[target].dEulerA[0] = deulera[0];
	  SphP[target].dEulerA[1] = deulera[1];
	  SphP[target].dEulerA[2] = deulera[2];
	  SphP[target].dEulerB[0] = deulerb[0];
	  SphP[target].dEulerB[1] = deulerb[1];
	  SphP[target].dEulerB[2] = deulerb[2];
#endif
#ifdef VECT_PRO_CLEAN
	  SphP[target].BPredVec[0] = BVec[0];
	  SphP[target].BPredVec[1] = BVec[1];
	  SphP[target].BPredVec[2] = BVec[2];
#endif
	}
#ifdef BLACK_HOLES
      P[target].b1.dBH_Density = rho;
      P[target].b2.dBH_Entropy = smoothentr;
      P[target].b3.dBH_SurroundingGasVel[0] = gasvel[0];
      P[target].b3.dBH_SurroundingGasVel[1] = gasvel[1];
      P[target].b3.dBH_SurroundingGasVel[2] = gasvel[2];
#endif
#if defined(RADTRANSFER) || defined(SNIA_HEATING)
      if(P[target].Type == 4)
	P[target].DensAroundStar = rho;
#endif
    }
  else
    {
#ifdef HYDRO_COST_FACTOR
      DensDataResult[target].Ninteractions = ninteractions;
#endif
      DensDataResult[target].Rho = rho;
#ifdef VOLUME_CORRECTION
      DensDataResult[target].DensityStd = densitystd;
#endif
      DensDataResult[target].Ngb = weighted_numngb;
      DensDataResult[target].DhsmlDensity = dhsmlrho;
#ifndef NAVIERSTOKES
      DensDataResult[target].Div = divv;
      DensDataResult[target].Rot[0] = rotv[0];
      DensDataResult[target].Rot[1] = rotv[1];
      DensDataResult[target].Rot[2] = rotv[2];
#else
      for(k = 0; k < 3; k++)
	{
	  DensDataResult[target].DV[k][0] = dvel[k][0];
	  DensDataResult[target].DV[k][1] = dvel[k][1];
	  DensDataResult[target].DV[k][2] = dvel[k][2];
	}
#endif

#if defined(BLACK_HOLES)
      DensDataResult[target].SmoothedEntr = smoothentr;
#endif
#ifdef CONDUCTION_SATURATION
      DensDataResult[target].GradEntr[0] = gradentr[0];
      DensDataResult[target].GradEntr[1] = gradentr[1];
      DensDataResult[target].GradEntr[2] = gradentr[2];
#endif

#ifdef RADTRANSFER_FLUXLIMITER
      DensDataResult[target].Grad_ngamma[0] = grad_ngamma[0];
      DensDataResult[target].Grad_ngamma[1] = grad_ngamma[1];
      DensDataResult[target].Grad_ngamma[2] = grad_ngamma[2];
#endif

#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
      DensDataResult[target].RotB[0] = rotb[0];
      DensDataResult[target].RotB[1] = rotb[1];
      DensDataResult[target].RotB[2] = rotb[2];
#endif
#ifdef TRACEDIVB
      DensDataResult[target].divB = divB;
#endif
#ifdef BLACK_HOLES
      DensDataResult[target].GasVel[0] = gasvel[0];
      DensDataResult[target].GasVel[1] = gasvel[1];
      DensDataResult[target].GasVel[2] = gasvel[2];
#endif

#ifdef EULERPOTENTIALS
      DensDataResult[target].dEulerA[0] = deulera[0];
      DensDataResult[target].dEulerA[1] = deulera[1];
      DensDataResult[target].dEulerA[2] = deulera[2];
      DensDataResult[target].dEulerB[0] = deulerb[0];
      DensDataResult[target].dEulerB[1] = deulerb[1];
      DensDataResult[target].dEulerB[2] = deulerb[2];
#endif
#ifdef VECT_PRO_CLEAN
      DensDataResult[target].BPredVec[0] = BVec[0];
      DensDataResult[target].BPredVec[1] = BVec[1];
      DensDataResult[target].BPredVec[2] = BVec[2];
#endif
    }

  return 0;
}





int density_isactive(int n)
{
  if(P[n].TimeBin < 0)
    return 0;

#if defined(RADTRANSFER) || defined(SNIA_HEATING)
  if(P[n].Type == 4)
    return 1;
#endif

#ifdef BLACK_HOLES
  if(P[n].Type == 5)
    return 1;
#endif

  if(P[n].Type == 0)
    return 1;

  return 0;
}


#ifdef NAVIERSTOKES
double get_shear_viscosity(int i)
{
  return All.NavierStokes_ShearViscosity;
}
#endif
