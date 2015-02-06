#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef RADTRANSFER
#ifdef CG

#define NSTEP 1
#define MAX_ITER 200
#define ACCURACY 1.0e-5
#define EPSILON 1.0e-5

/*structures for radtransfer*/
struct radtransferdata_in
{
  int NodeList[NODELISTLENGTH];
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat ET[6];
  double Kappa, Lambda;
  MyFloat Mass, Density;
}
 *RadTransferDataIn, *RadTransferDataGet;

struct radtransferdata_out
{
  double Out, Sum;
}
 *RadTransferDataResult, *RadTransferDataOut;

static int *NextActiveParticleRT;
static double *XVec;
static double *QVec, *DVec, *Residue, *Zvec;
static double *Kappa, *Lambda, *Diag, *Diag2;
static double c_light, dt, a3inv, hubble_a;
static double sigma;

void radtransfer(void)
{
  int i, j, iter;
  double alpha_cg, beta, delta_old, delta_new, sum, old_sum, min_diag, glob_min_diag, max_diag, glob_max_diag;
  double timestep, nH, je=0;
  double rel, res, maxrel, glob_maxrel;
  double DQ;

  NextActiveParticleRT = (int *) mymalloc(N_gas * sizeof(int));
  XVec = (double *) mymalloc(N_gas * sizeof(double));
  QVec = (double *) mymalloc(N_gas * sizeof(double));
  DVec = (double *) mymalloc(N_gas * sizeof(double));
  Residue = (double *) mymalloc(N_gas * sizeof(double));
  Kappa = (double *) mymalloc(N_gas * sizeof(double));
  Lambda = (double *) mymalloc(N_gas * sizeof(double));
  Diag = (double *) mymalloc(N_gas * sizeof(double));
  Zvec = (double *) mymalloc(N_gas * sizeof(double));
  Diag2 = (double *) mymalloc(N_gas * sizeof(double));

  c_light = C / All.UnitVelocity_in_cm_per_s;
  sigma = 6.3e-18 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;;


  /*  the actual time-step we need to do */
  timestep = (All.Radiation_Ti_endstep - All.Radiation_Ti_begstep) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      /* in comoving case, timestep is dloga at this point. Convert to dt */
      timestep /= hubble_a;
    }
  else
    {
      a3inv = hubble_a = 1.0;
    }

  //radtransfer_get_active();

  dt = timestep / NSTEP;
  
  for(i = 0; i < NSTEP; i++)
    {
      if(ThisTask == 0)
	{
	  printf("%s %i\n", "the step is ", i);
	  printf("%s %g\n", "c_light is ", c_light);
	  printf("%s %g\n", "dt is ", dt);
	  fflush(stdout);
	}
      
      /* initialization for the CG method */
      
      for(j = 0; j < N_gas; j++)
	if(P[j].Type == 0)
#ifdef RT_VAR_DT
	  if(SphP[j].rt_flag == 1)
#endif
	    {
	      XVec[j] = SphP[j].n_gamma;
	      
	      nH = (SphP[j].d.Density * a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam); //physical
	      
	      Kappa[j] = (SphP[j].nHI + 1.0e-8) * nH * sigma; //physical
	      
	      /* physical kappa = 1/distance => comoving is 1/(distance/a) = a/distance = a * kappa */
	      if(All.ComovingIntegrationOn)
		Kappa[j] = (Kappa[j] + 3.0 * hubble_a / c_light) * All.Time; //comoving
	      
#ifdef RADTRANSFER_FLUXLIMITER
	      /* now calculate flux limiter */
	      
	      if(SphP[j].n_gamma > 0)
		{
		  double R =  sqrt(SphP[j].Grad_ngamma[0] * SphP[j].Grad_ngamma[0] +
				   SphP[j].Grad_ngamma[1] * SphP[j].Grad_ngamma[1] +
				   SphP[j].Grad_ngamma[2] * SphP[j].Grad_ngamma[2]) / (SphP[j].n_gamma * Kappa[j]);
		  
		  if(All.ComovingIntegrationOn)
		    R /= All.Time;
		  
		  Lambda[j] = (1+R)/(1+R+R*R);
		  if(Lambda[j] <  1e-100)
		    Lambda[j] = 0;
		}
	      else
		Lambda[j] = 1.0;
#endif 	    
	    }
      
      /* add the source term */
      
      for(j = 0; j < N_gas; j++)
	if(P[j].Type == 0)
#ifdef RT_VAR_DT
	  if(SphP[j].rt_flag == 1)
#endif
	    {
	      SphP[j].n_gamma += dt * SphP[j].Je * P[j].Mass;
	      if(SphP[j].Je)
		je+=SphP[j].Je;
	  }
      // printf("glowing particles lum %g \n", je);
      
      radtransfer_matrix_multiply(XVec, Residue, Diag);
      
      for(j = 0; j < N_gas; j++)
	if(P[j].Type == 0)
#ifdef RT_VAR_DT
	  if(SphP[j].rt_flag == 1)
#endif
	    Residue[j] = SphP[j].n_gamma - Residue[j];
      
      /* Let's take the diagonal matrix elements as Jacobi preconditioner */
      
      for(j = 0, min_diag = MAX_REAL_NUMBER, max_diag = -MAX_REAL_NUMBER; j < N_gas; j++)
	if(P[j].Type == 0)
#ifdef RT_VAR_DT
	  if(SphP[j].rt_flag == 1)
#endif
	    {
	      /* note: in principle we would have to substract the w_ii term, but this is always zero */
	      
	      if(Diag[j] < min_diag)
		min_diag = Diag[j];
	      if(Diag[j] > max_diag)
		max_diag = Diag[j];
	      
	      Zvec[j] = Residue[j] / Diag[j];
	      DVec[j] = Zvec[j];
	    }
      
      MPI_Allreduce(&min_diag, &glob_min_diag, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&max_diag, &glob_max_diag, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      
      delta_new = radtransfer_vector_multiply(Zvec, Residue);
      delta_old = delta_new;
	
      old_sum = radtransfer_vector_sum(XVec);
      
      if(ThisTask == 0)
	printf("\nBegin CG iteration\nold |x|=%g, min-diagonal=%g, max-diagonal=%g\n", old_sum,
	       glob_min_diag, glob_max_diag);
      
      
      /* begin the CG method iteration */
      iter = 0;
      
      do
	{
	  radtransfer_matrix_multiply(DVec, QVec, Diag2);
	  
	  DQ = radtransfer_vector_multiply(DVec, QVec);
	  if(DQ == 0)
	    alpha_cg = 0;
	  else
	    alpha_cg = delta_new / DQ;
	  
	  
	  for(j = 0, maxrel = 0; j < N_gas; j++)
	    {
	      XVec[j] += alpha_cg * DVec[j];
	      Residue[j] -= alpha_cg * QVec[j];
	      
	      Zvec[j] = Residue[j] / Diag[j];
	      
	      rel = fabs(alpha_cg * DVec[j]) / (XVec[j] + 1.0e-10);
	      if(rel > maxrel)
		maxrel = rel;
	    }
	  
	  delta_old = delta_new;
	  delta_new = radtransfer_vector_multiply(Zvec, Residue);
	  
	  sum = radtransfer_vector_sum(XVec);
	  res = radtransfer_vector_sum(Residue);
	  
	  if(delta_old)
	    beta = delta_new / delta_old;
	  else
	    beta = 0;
	  
	  for(j = 0; j < N_gas; j++)
	    DVec[j] = Zvec[j] + beta * DVec[j];
	  
	  MPI_Allreduce(&maxrel, &glob_maxrel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	  
	  if(ThisTask == 0)
	    {
	      printf("radtransfer: iter=%3d  |res|/|x|=%12.6g  maxrel=%12.6g  |x|=%12.6g | res|=%12.6g\n",
		     iter, res/sum, glob_maxrel, sum, res);
	      fflush(stdout);
	    }
	  
	  iter++;
	}
      while((res > ACCURACY * sum && glob_maxrel > ACCURACY  && iter < MAX_ITER) || iter < 2);
      
      if(ThisTask == 0)
	printf("\n");
      
      /* update the intensity */
      for(j = 0; j < N_gas; j++)
	if(P[j].Type == 0)
#ifdef RT_VAR_DT
	  if(SphP[j].rt_flag == 1)
#endif
	    {
	      if(XVec[j] < 0 && XVec[j] > -EPSILON)
		XVec[j] = 0;
	      
	      SphP[j].n_gamma = XVec[j];
	    }
      
      /* update the chemistry */
      radtransfer_update_chemistry();
    }
  
  myfree(Diag2);
  myfree(Zvec);
  myfree(Diag);
  myfree(Lambda);
  myfree(Kappa);
  myfree(Residue);
  myfree(DVec);
  myfree(QVec);
  myfree(XVec);
  myfree(NextActiveParticleRT);

}

/* internal product of two vectors */
double radtransfer_vector_multiply(double *a, double *b)
{
  int i;
  double sum, sumall;

  for(i = 0, sum = 0; i < N_gas; i++)
    if(P[i].Type == 0)
#ifdef RT_VAR_DT
      if(SphP[i].rt_flag == 1)
#endif
	sum += a[i] * b[i];
  
  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sumall;
}


double radtransfer_vector_sum(double *a)
{
  int i;
  double sum, sumall;

  for(i = 0, sum = 0; i < N_gas; i++)
    if(P[i].Type == 0)
#ifdef RT_VAR_DT
      if(SphP[i].rt_flag == 1)
#endif
	sum += fabs(a[i]);
  
  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  return sumall;
}


/* this function computes the vector b(out) given the vector x(in) such as Ax = b, where A is a matrix */
void radtransfer_matrix_multiply(double *in, double *out, double *sum)
{
  int i, j, k, ngrp, dummy, ndone, ndone_flag;
  int sendTask, recvTask, nexport, nimport, place;
  double a;

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct radtransferdata_in) +
					     sizeof(struct radtransferdata_out) +
					     sizemax(sizeof(struct radtransferdata_in),
						     sizeof(struct radtransferdata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));
  
  if(All.ComovingIntegrationOn)
    a = All.Time;
  else
    a = 1.0;
  
  i = 0;

  do				/* communication loop */
    {

      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      for(nexport = 0; i < N_gas; i++)
	{
	  if(P[i].Type == 0)
#ifdef RT_VAR_DT
	    if(SphP[i].rt_flag == 1)
#endif
	      if(radtransfer_evaluate(i, 0, in, out, sum, &nexport, Send_count) < 0)
		break;
	  
	}
      
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

      RadTransferDataGet =
	(struct radtransferdata_in *) mymalloc(nimport * sizeof(struct radtransferdata_in));
      RadTransferDataIn = (struct radtransferdata_in *) mymalloc(nexport * sizeof(struct radtransferdata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  for(k = 0; k < 3; k++)
	    {
	      RadTransferDataIn[j].Pos[k] = P[place].Pos[k];
	      RadTransferDataIn[j].ET[k] = SphP[place].ET[k];
	      RadTransferDataIn[j].ET[k + 3] = SphP[place].ET[k + 3];
	    }
	  RadTransferDataIn[j].Hsml = PPP[place].Hsml;
	  RadTransferDataIn[j].Kappa = Kappa[place];
	  RadTransferDataIn[j].Lambda = Lambda[place];
	  RadTransferDataIn[j].Mass = P[place].Mass;
	  RadTransferDataIn[j].Density = SphP[place].d.Density;

	  memcpy(RadTransferDataIn[j].NodeList,
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
		  MPI_Sendrecv(&RadTransferDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct radtransferdata_in), MPI_BYTE,
			       recvTask, TAG_RT_A,
			       &RadTransferDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct radtransferdata_in), MPI_BYTE,
			       recvTask, TAG_RT_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      myfree(RadTransferDataIn);
      RadTransferDataResult =
	(struct radtransferdata_out *) mymalloc(nimport * sizeof(struct radtransferdata_out));
      RadTransferDataOut =
	(struct radtransferdata_out *) mymalloc(nexport * sizeof(struct radtransferdata_out));

      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
	radtransfer_evaluate(j, 1, in, out, sum, &dummy, &dummy);
      
      if(i < N_gas)
	ndone_flag = 0;
      else
	ndone_flag = 1;

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
		  MPI_Sendrecv(&RadTransferDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct radtransferdata_out),
			       MPI_BYTE, recvTask, TAG_RT_B,
			       &RadTransferDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct radtransferdata_out),
			       MPI_BYTE, recvTask, TAG_RT_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  out[place] += RadTransferDataOut[j].Out;
	  sum[place] += RadTransferDataOut[j].Sum;
	}

      myfree(RadTransferDataOut);
      myfree(RadTransferDataResult);
      myfree(RadTransferDataGet);

    }
  while(ndone < NTask);

  /* do final operations on results */
  for(i = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      {
#ifdef RT_VAR_DT
	if(SphP[i].rt_flag == 0)
	  {
	    sum[i] = 0.0;
	    out[i] = 0.0;
	  }
	else
#endif
	  {
	    /* divide c_light by a to get comoving speed of light (because kappa is comoving) */
	    if((1 + dt * (c_light / a) * Kappa[i] + sum[i]) < 0)
	      {
		printf("1 + sum + rate= %g   sum=%g rate=%g i =%d\n",
		       1 + dt * (c_light / a) * Kappa[i] + sum[i],
		       sum[i], dt * (c_light / a) * Kappa[i], i);
		endrun(123);
	      }
	    
	    sum[i] += 1.0 + dt * (c_light / a) * Kappa[i];
	    
	    out[i] += in[i] * sum[i];
	  }
      }
  
  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}

/* this function evaluates parts of the matrix A */
int radtransfer_evaluate(int target, int mode, double *in, double *out, double *sum, int *nexport,
			 int *nsend_local)
{
  int startnode, numngb, listindex = 0;
  int j, n, k;
  MyFloat *ET_aux, ET_j[6], ET_i[6], ET_ij[6];
  MyFloat kappa_i, kappa_j, kappa_ij;
#ifdef RADTRANSFER_FLUXLIMITER
  MyFloat lambda_i, lambda_j;
#endif
  MyDouble *pos;
  MyFloat mass, mass_i, rho, rho_i;
  double sum_out = 0, sum_w = 0, a3, fac = 0;

  double dx, dy, dz;
  double h_j, hinv, hinv4, h_i;
  double dwk_i, dwk_j, dwk;
  double r, r2, r3, u;

#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif

  if(All.ComovingIntegrationOn)
    a3 = All.Time * All.Time * All.Time;
  else
    a3 = 1.0;

  if(mode == 0)
    {
      ET_aux = SphP[target].ET;
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      kappa_i = Kappa[target];
#ifdef RADTRANSFER_FLUXLIMITER
      lambda_i = Lambda[target];
#endif
      mass_i = P[target].Mass;
      rho_i = SphP[target].d.Density;
    }
  else
    {
      ET_aux = RadTransferDataGet[target].ET;
      pos = RadTransferDataGet[target].Pos;
      h_i = RadTransferDataGet[target].Hsml;
      kappa_i = RadTransferDataGet[target].Kappa;
#ifdef RADTRANSFER_FLUXLIMITER
      lambda_i = RadTransferDataGet[target].Lambda;
#endif
      mass_i = RadTransferDataGet[target].Mass;
      rho_i = RadTransferDataGet[target].Density;
    }

#ifdef RADTRANSFER_MODIFY_EDDINGTON_TENSOR
  /*modify Eddington tensor */
  ET_i[0] = 2 * ET_aux[0] - 0.5 * ET_aux[1] - 0.5 * ET_aux[2];
  ET_i[1] = 2 * ET_aux[1] - 0.5 * ET_aux[2] - 0.5 * ET_aux[0];
  ET_i[2] = 2 * ET_aux[2] - 0.5 * ET_aux[0] - 0.5 * ET_aux[1];

  for(k = 3; k < 6; k++)
    ET_i[k] = 2.5 * ET_aux[k];
#else
  for(k = 0; k < 6; k++)
    ET_i[k] = ET_aux[k];
#endif

  if(mode == 0)
    {
      startnode = All.MaxPart;
    }
  else
    {
      startnode = RadTransferDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_pairs(pos, h_i, target, &startnode, mode, nexport, nsend_local);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];
#ifdef RT_VAR_DT
	      if(SphP[j].rt_flag == 1)
#endif
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
		  r = sqrt(r2);
		  r3 = r2 * r;
		  h_j = PPP[j].Hsml;
		  
		  if(r > 0 && (r < h_i || r < h_j))
		    {
		      mass = P[j].Mass;
		      rho = SphP[j].d.Density;
		      kappa_j = Kappa[j];
#ifdef RADTRANSFER_FLUXLIMITER
		      lambda_j = Lambda[j];
#endif
		      
#ifdef RADTRANSFER_MODIFY_EDDINGTON_TENSOR
		      ET_aux = SphP[j].ET;
		      
		      /*modify Eddington tensor */
		      ET_j[0] = 2 * ET_aux[0] - 0.5 * ET_aux[1] - 0.5 * ET_aux[2];
		      ET_j[1] = 2 * ET_aux[1] - 0.5 * ET_aux[2] - 0.5 * ET_aux[0];
		      ET_j[2] = 2 * ET_aux[2] - 0.5 * ET_aux[0] - 0.5 * ET_aux[1];
		      
		      for(k = 3; k < 6; k++)
			ET_j[k] = 2.5 * ET_aux[k];
#else
		      for(k = 0; k < 6; k++)
			ET_j[k] = SphP[j].ET[k];
#endif
		      
		      for(k = 0; k < 6; k++)
			ET_ij[k] = 0.5 * (ET_i[k] + ET_j[k]);
		      
		      if(r < h_i)
			{
			  hinv = 1.0 / h_i;
			  hinv4 = hinv * hinv * hinv * hinv;
			  u = r * hinv;
			  
			  if(u < 0.5)
			    dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
			  else
			    dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
			}
		      else
			dwk_i = 0;
		      
		      if(r < h_j)
			{
			  hinv = 1.0 / h_j;
			  hinv4 = hinv * hinv * hinv * hinv;
			  u = r * hinv;
			  
			  if(u < 0.5)
			    dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
			  else
			    dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
			}
		      else
			dwk_j = 0;
		      
		      kappa_ij = 0.5 * (1/kappa_i + 1/kappa_j);
		      dwk = 0.5 * (dwk_i + dwk_j);	 
		      mass = 0.5 * (mass + mass_i);
		      rho = 0.5 * (rho + rho_i);
		      
		      double tensor = (ET_ij[0] * dx * dx + ET_ij[1] * dy * dy + ET_ij[2] * dz * dz
				       + 2.0 * ET_ij[3] * dx * dy + 2.0 * ET_ij[4] * dy * dz + 2.0 * ET_ij[5] * dz * dx);
		      
		      if(tensor > 0)
			{
			  /* divide c_light by a because kappa is comoving */
			  /* all variables are comoving!!! */
			  fac = -2.0 * dt * (c_light / a3) * (mass / rho) * kappa_ij * dwk / r3 * tensor;
			  
#ifdef RADTRANSFER_FLUXLIMITER
			  fac *= 0.5 * (lambda_i + lambda_j);
#endif
			  
			  sum_out -= fac * in[j];
			  
			  sum_w += fac;
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
	      startnode = RadTransferDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;
	    }
	}
    }
  
  if(mode == 0)
    {
      out[target] = sum_out;
      sum[target] = sum_w;
    }
  else
    {
      RadTransferDataResult[target].Out = sum_out;
      RadTransferDataResult[target].Sum = sum_w;
    }
  
  return 0;
}


 /* this function calculates the mean photon fraction */
void radtransfer_mean(void)
{
  int i;
  double n_gamma, n_gamma_all = 0;

  for(i = 0, n_gamma = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      n_gamma += SphP[i].n_gamma;

  MPI_Allreduce(&n_gamma, &n_gamma_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  n_gamma_all /= All.TotN_gas;

  if(ThisTask == 0)
    {
      printf("n_gamma_all is: %g\n", n_gamma_all);
      fflush(stdout);
    }

}

/* this function sets up simple initial conditions for a single source in a uniform field of gas with constant density*/
void radtransfer_set_simple_inits(void)
{
  int i;

  for(i = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      {
	/* in code units */
	SphP[i].nHII = 0.0;
	SphP[i].nHI = 1 - SphP[i].nHII;
	SphP[i].n_elec = SphP[i].nHII;
	SphP[i].n_gamma = 0.0;
      }
}

/* produces a simple output - particle positions x, y and z, overintensity and neutral fraction*/
void radtransfer_simple_output(void)
{
  char buf[100];
  int i;

  sprintf(buf, "%s%s%i%s%f%s", All.OutputDir, "radtransfer_", ThisTask, "_", All.Time, ".txt");
  FdRadtransfer = fopen(buf, "wa");
  for(i = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      fprintf(FdRadtransfer, "%f %f %f %g %f \n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], SphP[i].n_gamma,
	      SphP[i].nHII);

  fclose(FdRadtransfer);
}

#ifdef RT_VAR_DT
/* produces a simple output - particle positions x, y and z, overintensity and neutral fraction*/
void radtransfer_get_active(void)
{
  int i, k;

  for(i = 0, k = 0; i < N_gas; i++)
    if(SphP[i].rt_flag == 1)
      {
	NextActiveParticleRT[k] = i;
	k++;
      }
  NextActiveParticleRT[k] = -1;
}
#endif


/* solves the recombination equation and updates the HI and HII abundances */
void radtransfer_update_chemistry(void)
{
  int i, j, n;
  double inter_dt, nH, temp, alpha, a3inv, gamma, molecular_weight;
#ifdef RT_PHOTOHEATING
  double dtemp, e1, sigma1, de1, rate1;
  double eV_to_erg = 1.60184e-12;
#endif
#ifdef RT_COOLING
  double de2, de3, de4, de5, rate2, rate3, rate4, rate5;
#endif 
  double a, b, c, s = 1, q;
  
  n = 1;
  inter_dt = dt / FLT(n);

  if(All.ComovingIntegrationOn)   
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1.0;
  

  /* begin substepping */
  for(j = 0; j < n; j++)
    for(i = 0; i < N_gas; i++)
      if(P[i].Type == 0)
#ifdef RT_VAR_DT
	if(SphP[i].rt_flag == 1)
#endif
	  {  
	    SphP[i].n_gamma /= P[i].Mass;
	      
	    /* number of photons should be positive */
	    if(SphP[i].n_gamma < 0)
	      {
		printf("NEGATIVE n_gamma: %g %d %d \n", SphP[i].n_gamma, i, ThisTask);
		endrun(111);
	      }
	    
	    nH = (SphP[i].d.Density * a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam); //physical
	    
	    molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC); //needed to obtain correct temperatures
	    /* molecular_weight should equal 1 in this case, but then HYDROGEN_MASSFRAC has to be set to 1 in allvars.h */
	    
	    temp = SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) * molecular_weight * (PROTONMASS / All.UnitMass_in_g) /
	      (BOLTZMANN / All.UnitEnergy_in_cgs);
	    //printf("temp %g %g %g\n", temp, molecular_weight, HYDROGEN_MASSFRAC);
	    
	    /* collisional ionization rate */
#ifdef RT_COLLISIONAL_IONIZATION
	    gamma = 5.85e-11 * pow(temp,0.5) * exp(-157809.1/temp) / (1.0 + pow(temp/1e5,0.5));
	    gamma *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
	    gamma *= All.HubbleParam * All.HubbleParam;
#else
	    gamma=0;
#endif
	    
	    /* alpha_B recombination coefficient */
	    alpha = 2.59e-13 * pow(temp/1e4,-0.7);
	    alpha *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
	    alpha *= All.HubbleParam * All.HubbleParam;
	    
	    
#ifndef SEMI_IMPLICIT
	    /* implicit scheme for ionization*/
	    a = inter_dt * alpha * nH;
	    b = -2.0 * inter_dt * alpha * nH - inter_dt * c_light * sigma * nH * SphP[i].n_gamma - 1.0;
	    c = SphP[i].nHI +  inter_dt * alpha * nH;
	    
	    if(b > 0)
	      s = 1.0;
	    if(b < 0)
	      s = -1.0;
	    if(b == 0)
	      s = 0.0;
	    
	    q = -0.5 * (b + s * sqrt(b * b - 4.0 * a * c));
	    SphP[i].nHI = c / q;
	    
#ifdef RT_NO_RECOMBINATIONS
	    SphP[i].nHII += inter_dt * c_light * sigma * nH * SphP[i].n_gamma * SphP[i].nHI;
	    SphP[i].nHI = 1.0 -  SphP[i].nHII;
#endif

	    if(SphP[i].nHI < 0 && SphP[i].nHI > -EPSILON)
	      SphP[i].nHI = 0;
	    if(SphP[i].nHI > 1 && SphP[i].nHI < (1+EPSILON))
	      SphP[i].nHI = 1;
	    if(SphP[i].nHI < -EPSILON || SphP[i].nHI > (1+EPSILON))
	      {
		printf("WRONG nHII: %g %d %d \n", SphP[i].nHII, i, ThisTask);
		endrun(222);
	      }
	    
	    SphP[i].nHII = 1 - SphP[i].nHI;
	    SphP[i].n_elec = SphP[i].nHII;
	    
#else
	    /* semi-implicit scheme for ionization*/
	    SphP[i].nHII = (SphP[i].nHII + inter_dt * c_light * sigma * nH * SphP[i].n_gamma + inter_dt * gamma * nH * SphP[i].n_elec) /
	      (1.0 + inter_dt * c_light * sigma * nH * SphP[i].n_gamma +
	       inter_dt * alpha * nH * SphP[i].n_elec +  inter_dt * gamma * nH * SphP[i].n_elec);
	    
	    /* fraction should be between 0 and 1 */
	    if(SphP[i].nHII < 0 || SphP[i].nHII > 1)
	      {
		printf("WRONG nHII: %g %d %d \n", SphP[i].nHII, i, ThisTask);
		endrun(222);
	      }
	    
	    SphP[i].nHI = 1 - SphP[i].nHII;
	    SphP[i].n_elec = SphP[i].nHII;
#endif
	    
#ifdef RT_PHOTOHEATING
	    /*photoheating for 10^5K blackbody*/
	    sigma1 = 1.63e-18 / All.UnitLength_in_cm / All.UnitLength_in_cm; //code units
	    sigma1 *= All.HubbleParam * All.HubbleParam;
	    
	    e1 = 29.65 * eV_to_erg; //real units
	    
	    /* all rates in erg cm^3 s^-1 in code units and energy in ergs*/
	    /* photoheating rate */
	    rate1 = c_light * e1 * sigma1;
	    
	    de1 = SphP[i].nHI * nH * SphP[i].n_gamma * nH * rate1;
	    dtemp = de1 * GAMMA_MINUS1 / nH / BOLTZMANN * inter_dt;
	    SphP[i].Entropy += dtemp * (BOLTZMANN / All.UnitEnergy_in_cgs) / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / 
	      molecular_weight / (PROTONMASS / All.UnitMass_in_g);
#endif
	    
#ifdef RT_COOLING
	    /* recombination cooling rate */
	    rate2 = 8.7e-27 * pow(temp,0.5) * pow(temp/1e3,-0.2) / (1.0 + pow(temp/1e6,0.7));
	    rate2 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm; 
	    rate2 *= All.HubbleParam * All.HubbleParam;
	    
	    /* collisional ionization cooling rate */
	    rate3 = 1.27e-21 * pow(temp,0.5) * exp(-157809.1/temp) / (1.0 + pow(temp/1e5,0.5));
	    rate3 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm; 
	    rate3 *= All.HubbleParam * All.HubbleParam;
	    
	    /* collisional excitation cooling rate */
	    rate4 = 7.5e-19 * (1.0 + pow(temp/1e5,0.5)) * exp(-118348/temp);
	    rate4 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
	    rate4 *= All.HubbleParam * All.HubbleParam;
	    
	    /* Bremsstrahlung cooling rate */
	    rate5 = 1.42e-27 * 1.3 * pow(temp,0.5);
	    rate5 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
	    rate5 *= All.HubbleParam * All.HubbleParam;
	    
	    de2 = SphP[i].nHII * nH * SphP[i].n_elec * nH * rate2;
	    de3 = SphP[i].nHI * nH * SphP[i].n_elec * nH * rate3;
	    de4 = SphP[i].nHI * nH * SphP[i].n_elec * nH * rate4;	  
	    de5 = SphP[i].nHII * nH * SphP[i].n_elec * nH * rate5;
	    
	    dtemp = (de2 + de3 + de4 + de5) * GAMMA_MINUS1 / nH / BOLTZMANN * inter_dt;
	    SphP[i].Entropy -= dtemp * (BOLTZMANN / All.UnitEnergy_in_cgs) / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / 
	      molecular_weight / (PROTONMASS / All.UnitMass_in_g);
#endif
	    SphP[i].n_gamma *= P[i].Mass;
	  }
  
}

#endif
#endif
