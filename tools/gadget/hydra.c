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
#ifdef MACH_NUM
#include "machfinder.h"
#endif
#ifdef CS_MODEL
#include "cs_metals.h"
#endif

#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>

/*! \file hydra.c
*  \brief Computation of SPH forces and rate of entropy generation
*
*  This file contains the "second SPH loop", where the SPH forces are
*  computed, and where the rate of change of entropy due to the shock heating
*  (via artificial viscosity) is computed.
*/


struct hydrodata_in
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat Density;
  MyFloat Pressure;
  MyFloat F1;
  MyFloat DhsmlDensityFactor;
  int Timestep;

#ifdef CS_MODEL
  MyFloat DensityNow;
  MyFloat Entropy;
#endif

#ifdef PARTICLE_DEBUG
  MyIDType ID;			/*!< particle identifier */
#endif

#ifdef MAGNETIC
  MyFloat BPred[3];
#ifdef TIME_DEP_MAGN_DISP
  MyFloat Balpha;
#endif
#ifdef DIVBCLEANING_DEDNER
  MyFloat PhiPred;
#endif
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
  MyFloat RotB[3];
#endif
#endif
#ifdef TIME_DEP_ART_VISC
  MyFloat alpha;
#endif

#if defined(NAVIERSTOKES)
  MyFloat Entropy;
#endif


#ifdef NAVIERSTOKES
  MyFloat stressoffdiag[3];
  MyFloat stressdiag[3];
  MyFloat shear_viscosity;
#endif

#ifdef NAVIERSTOKES_BULK
  MyFloat divvel;
#endif

#ifdef EOS_DEGENERATE
  MyFloat dpdr;
#endif

  int NodeList[NODELISTLENGTH];
}
 *HydroDataIn, *HydroDataGet;


struct hydrodata_out
{
  MyLongDouble Acc[3];
  MyLongDouble DtEntropy;
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  MyFloat MinViscousDt;
#else
  MyFloat MaxSignalVel;
#endif
#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
  MyFloat DtB[3];
#ifdef DIVBCLEANING_DEDNER
  MyFloat DtPhi;
#endif
#endif

#if  defined(CR_SHOCK)
  MyFloat CR_EnergyChange[NUMCRPOP];
  MyFloat CR_BaryonFractionChange[NUMCRPOP];
#endif

#ifdef HYDRO_COST_FACTOR
  int Ninteractions;
#endif
}
 *HydroDataResult, *HydroDataOut;






#ifdef MACHNUM
double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;
#else
static double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;
#endif

/*! This function is the driver routine for the calculation of hydrodynamical
*  force and rate of change of entropy due to shock heating for all active
*  particles .
*/
void hydro_force(void)
{
  int i, j, k, ngrp, ndone, ndone_flag, dummy;
  int sendTask, recvTask, nexport, nimport, place;
  double soundspeed_i;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0, timenetwork = 0;
  double timecomp, timecomm, timewait, tstart, tend, t0, t1;

#if defined(WINDS) || defined(TIME_DEP_ART_VISC) || defined(MAGNETIC)
  double dmax1, dmax2;
#endif
#ifdef NAVIERSTOKES
  double fac;
#endif

#if (!defined(COOLING) && !defined(CR_SHOCK) && (defined(CR_DISSIPATION) || defined(CR_THERMALIZATION)))
  double utherm;
  double dt;
  int CRpop;
#endif

#if defined(CR_SHOCK)
  double rShockEnergy;
  double rNonRethermalizedEnergy;

#ifndef COOLING
  double utherm, CRpop;
#endif
#endif

#ifdef WINDS
  double windspeed, hsml_c;
#endif


#ifdef TIME_DEP_ART_VISC
  double f, cs_h;
#endif
#if defined(MAGNETIC) && defined(MAGFORCE)
#ifdef TIME_DEP_MAGN_DISP
  double mu0 = 1;
#endif
#endif
#ifdef DIVBCLEANING_DEDNER
  double phiphi, tmpb;
#endif
#ifdef HEALPIX
  double r_new, t[3];
  long ipix;
  int count = 0;
  int total_count = 0;
#endif

#ifdef NUCLEAR_NETWORK
  double dedt_nuc;
  int nuc_particles = 0;
  int nuc_particles_sum;
#endif

#ifdef WAKEUP
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	SphP[i].wakeup = 0;
    }
#endif

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      hubble_a = hubble_function(All.Time);
      hubble_a2 = All.Time * All.Time * hubble_a;

      fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;

      fac_egy = pow(All.Time, 3 * (GAMMA - 1));

      fac_vsic_fix = hubble_a * pow(All.Time, 3 * GAMMA_MINUS1);

      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = fac_egy = 1.0;

#if defined(MAGFORCE) && defined(TIME_DEP_MAGN_DISP)
#ifndef MU0_UNITY
  mu0 *= (4 * M_PI);
  mu0 /= All.UnitTime_in_s * All.UnitTime_in_s *
    All.UnitLength_in_cm / (All.UnitMass_in_g * All.HubbleParam * All.HubbleParam);
#endif
#endif


  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc(NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct hydrodata_in) +
					     sizeof(struct hydrodata_out) +
					     sizemax(sizeof(struct hydrodata_in),
						     sizeof(struct hydrodata_out))));
  DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));


  CPU_Step[CPU_HYDMISC] += measure_time();
  t0 = second();

  i = FirstActiveParticle;	/* first particle for this task */

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
	if(P[i].Type == 0)
	  {
	    if(hydro_evaluate(i, 0, &nexport, Send_count) < 0)
	      break;
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

      HydroDataGet = (struct hydrodata_in *) mymalloc(nimport * sizeof(struct hydrodata_in));
      HydroDataIn = (struct hydrodata_in *) mymalloc(nexport * sizeof(struct hydrodata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
	      HydroDataIn[j].Pos[k] = P[place].Pos[k];
	      HydroDataIn[j].Vel[k] = SphP[place].VelPred[k];
	    }
	  HydroDataIn[j].Hsml = PPP[place].Hsml;
	  HydroDataIn[j].Mass = P[place].Mass;
	  HydroDataIn[j].DhsmlDensityFactor = SphP[place].h.DhsmlDensityFactor;
	  HydroDataIn[j].Density = SphP[place].d.Density;
	  HydroDataIn[j].Pressure = SphP[place].Pressure;
	  HydroDataIn[j].Timestep = (P[place].TimeBin ? (1 << P[place].TimeBin) : 0);
#ifdef EOS_DEGENERATE
	  HydroDataIn[j].dpdr = SphP[place].dpdr;
#endif

	  /* calculation of F1 */
#ifndef ALTVISCOSITY
#ifndef EOS_DEGENERATE
	  soundspeed_i = sqrt(GAMMA * SphP[place].Pressure / SphP[place].d.Density);
#else
	  soundspeed_i = sqrt(SphP[place].dpdr);
#endif
#ifndef NAVIERSTOKES
	  HydroDataIn[j].F1 = fabs(SphP[place].v.DivVel) /
	    (fabs(SphP[place].v.DivVel) + SphP[place].r.CurlVel +
	     0.0001 * soundspeed_i / PPP[place].Hsml / fac_mu);
#else
	  HydroDataIn[j].F1 = fabs(SphP[place].v.DivVel) /
	    (fabs(SphP[place].v.DivVel) + SphP[place].u.s.CurlVel +
	     0.0001 * soundspeed_i / PPP[place].Hsml / fac_mu);
#endif

#else
	  HydroDataIn[j].F1 = SphP[place].v.DivVel;
#endif

	  memcpy(HydroDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));

#ifdef CS_MODEL
	  HydroDataIn[j].DensityNow = SphP[place].d.Density;
	  HydroDataIn[j].Entropy = SphP[place].Entropy;
#endif


#ifdef MAGNETIC
	  for(k = 0; k < 3; k++)
	    {
	      HydroDataIn[j].BPred[k] = SphP[place].BPred[k];
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
#ifdef SMOOTH_ROTB
	      HydroDataIn[j].RotB[k] = SphP[place].SmoothedRotB[k];
#else
	      HydroDataIn[j].RotB[k] = SphP[place].RotB[k];
#endif
#endif
	    }
#ifdef DIVBCLEANING_DEDNER
#ifdef SMOOTH_PHI
	  HydroDataIn[j].PhiPred = SphP[place].SmoothPhi;
#else
	  HydroDataIn[j].PhiPred = SphP[place].PhiPred;
#endif
#endif
#endif


#if defined(NAVIERSTOKES)
	  HydroDataIn[j].Entropy = SphP[place].Entropy;
#endif

#ifdef TIME_DEP_ART_VISC
	  HydroDataIn[j].alpha = SphP[place].alpha;
#endif


#ifdef PARTICLE_DEBUG
	  HydroDataIn[j].ID = P[place].ID;
#endif

#ifdef NAVIERSTOKES
	  for(k = 0; k < 3; k++)
	    {
	      HydroDataIn[j].stressdiag[k] = SphP[i].u.s.StressDiag[k];
	      HydroDataIn[j].stressoffdiag[k] = SphP[i].u.s.StressOffDiag[k];
	    }
	  HydroDataIn[j].shear_viscosity = get_shear_viscosity(i);

#ifdef NAVIERSTOKES_BULK
	  HydroDataIn[j].divvel = SphP[i].u.s.DivVel;
#endif
#endif

#ifdef TIME_DEP_MAGN_DISP
	  HydroDataIn[j].Balpha = SphP[place].Balpha;
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
		  MPI_Sendrecv(&HydroDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A,
			       &HydroDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm1 += timediff(tstart, tend);


      myfree(HydroDataIn);
      HydroDataResult = (struct hydrodata_out *) mymalloc(nimport * sizeof(struct hydrodata_out));
      HydroDataOut = (struct hydrodata_out *) mymalloc(nexport * sizeof(struct hydrodata_out));



      /* now do the particles that were sent to us */

      tstart = second();
      for(j = 0; j < nimport; j++)
	{
	  hydro_evaluate(j, 1, &dummy, &dummy);
	}
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
		  MPI_Sendrecv(&HydroDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct hydrodata_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B,
			       &HydroDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct hydrodata_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

	  for(k = 0; k < 3; k++)
	    {
	      SphP[place].a.dHydroAccel[k] += HydroDataOut[j].Acc[k];
	    }

#ifdef GASRETURN
	  if(All.Time-SphP[place].GasAge > 0.00)
	  {	    
	    SphP[place].e.dDtEntropy += HydroDataOut[j].DtEntropy;
	  }else{
	    SphP[place].e.dDtEntropy += 0.0;
	  }
#else
	    SphP[place].e.dDtEntropy += HydroDataOut[j].DtEntropy;
#endif

#ifdef HYDRO_COST_FACTOR
	  P[place].GravCost += HYDRO_COST_FACTOR * HydroDataOut[j].Ninteractions;
#endif

#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
	  if(SphP[place].MinViscousDt > HydroDataOut[j].MinViscousDt)
	    SphP[place].MinViscousDt = HydroDataOut[j].MinViscousDt;
#else
	  if(SphP[place].MaxSignalVel < HydroDataOut[j].MaxSignalVel)
	    SphP[place].MaxSignalVel = HydroDataOut[j].MaxSignalVel;
#endif

#ifdef OUTPUTCOOLRATE
	  SphP[place].CondRate += HydroDataOut[j].CondRate;
#endif

#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
	  for(k = 0; k < 3; k++)
	    SphP[place].DtB[k] += HydroDataOut[j].DtB[k];
#ifdef DIVBCLEANING_DEDNER
	  SphP[place].DtPhi += HydroDataOut[j].DtPhi;
#endif
#endif
	}
      tend = second();
      timecomp1 += timediff(tstart, tend);

      myfree(HydroDataOut);
      myfree(HydroDataResult);
      myfree(HydroDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);


  /* do final operations on results */


#ifdef FLTROUNDOFFREDUCTION
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	SphP[i].e.DtEntropy = FLT(SphP[i].e.dDtEntropy);

	for(j = 0; j < 3; j++)
	  SphP[i].a.HydroAccel[j] = FLT(SphP[i].a.dHydroAccel[j]);
      }
#endif



  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
#ifdef CR_SHOCK
	/* state right here:
	 *
	 * _c denotes comoving quantities
	 * _p denotes physical quantities
	 *
	 *
	 * Delta u_p = rho_p^(gamma-1)/(gamma-1) Delta A
	 *
	 * Delta A = dA/dloga * Delta loga
	 *
	 * dA/dloga = DtE * (gamma-1) / ( H(a) a^2 rho_c^(gamma-1)
	 *
	 * => Delta u_p = DtE * dloga / ( H(a) a^2 a^(3(gamma-1)) )
	 */

	if(SphP[i].e.DtEntropy > 0.0)
	  {
	    rShockEnergy = SphP[i].e.DtEntropy *
	      (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / hubble_a2 / fac_egy;
	  }
	else
	  {
	    rShockEnergy = 0.0;
	  }

#endif /* CR_SHOCK */

#if !defined(EOS_DEGENERATE)

#ifndef TRADITIONAL_SPH_FORMULATION
	/* Translate energy change rate into entropy change rate */
	SphP[i].e.DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1));
#endif

#else
	/* DtEntropy stores the energy change rate in internal units */
	SphP[i].e.DtEntropy *= All.UnitEnergy_in_cgs / All.UnitTime_in_s;
#endif

#ifdef MACHNUM

	/* Estimates the Mach number of particle i for non-radiative runs,
	 * or the Mach number, density jump and specific energy jump
	 * in case of cosmic rays!
	 */
#if (CR_SHOCK == 2)
	GetMachNumberCR(SphP + i);
#else

#ifndef CS_MODEL
	GetMachNumber(SphP + i);
#else
	GetMachNumber(SphP + i, P + i);
#endif /* CS_MODEL */
#endif /* COSMIC_RAYS */
#endif /* MACHNUM */
#ifdef MACHSTATISTIC
	GetShock_DtEnergy(SphP + i);
#endif

#ifdef CR_SHOCK
	if(rShockEnergy > 0.0)
	  {
	    /* Feed fraction "All.CR_ShockEfficiency" into CR and see what
	     * amount of energy instantly gets rethermalized
	     *
	     * for this, we need the physical time step, which is
	     * Delta t_p = Delta t_c / hubble_a
	     */

	    /* The  CR_find_alpha_InjectTo induces an error in the density jump since it can set 
	     *  Particle->Shock_DensityJump = 1.0 + 1.0e-6 which is used in ShockInject as the input DensityJump
	     *  if (NUMCRPOP > 1)
	     *  {
	     *  #if ( CR_SHOCK == 1 )
	     *  InjPopulation = CR_Find_Alpha_to_InjectTo(All.CR_ShockAlpha);
	     *  #else          
	     *  InjPopulation = CR_find_alpha_InjectTo(SphP + i);
	     *  #endif
	     *  }
	     *  else 
	     *  InjPopulation = 0;
	     */

	    rNonRethermalizedEnergy =
	      CR_Particle_ShockInject(SphP + i,
				      rShockEnergy,
				      (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval /
				      hubble_a);

	    /* Fraction of total energy that went and remained in CR is
	     * rNonRethermalizedEnergy / rShockEnergy,
	     * hence, we conserve energy if we do:
	     */
#ifndef CR_NO_CHANGE
	    SphP[i].e.DtEntropy *= (1.0 - rNonRethermalizedEnergy / rShockEnergy);
#endif /* CR_NO_CHANGE */

	    assert(rNonRethermalizedEnergy >= 0.0);

	    assert(rNonRethermalizedEnergy <= (rShockEnergy * All.CR_ShockEfficiency));


#if (!defined(COOLING) && (defined(CR_DISSIPATION) || defined(CR_THERMALIZATION)))
	    utherm = 0.0;
	    for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	      utherm +=
		CR_Particle_ThermalizeAndDissipate(SphP + i,
						   (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) *
						   All.Timebase_interval / hubble_a, CRpop);

	    /* we need to add this thermalized energy to the internal energy */

	    SphP[i].e.DtEntropy += GAMMA_MINUS1 * utherm * fac_egy / pow(SphP[i].d.Density, GAMMA_MINUS1) /
	      ((P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval);
#endif

	  }
#endif /* CR_SHOCK */


#if (!defined(COOLING) && !defined(CR_SHOCK) && (defined(CR_DISSIPATION) || defined(CR_THERMALIZATION)))
	double utherm;
	double dt;
	int CRpop;

	dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / hubble_a;

	if(P[i].TimeBin)	/* upon start-up, we need to protect against dt==0 */
	  {
	    if(dt > 0)
	      {
		for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
		  {
		    utherm = CR_Particle_ThermalizeAndDissipate(SphP + i, dt, CRpop);

		    SphP[i].e.DtEntropy +=
		      GAMMA_MINUS1 * utherm * fac_egy / pow(SphP[i].d.Density,
							    GAMMA_MINUS1) / (dt * hubble_a);
		  }
	      }
	  }
#endif

#ifdef NAVIERSTOKES
	/* sigma_ab * sigma_ab */
	for(k = 0, fac = 0; k < 3; k++)
	  {
	    fac += SphP[i].u.s.StressDiag[k] * SphP[i].u.s.StressDiag[k] +
	      2 * SphP[i].u.s.StressOffDiag[k] * SphP[i].u.s.StressOffDiag[k];
	  }

#ifndef NAVIERSTOKES_CONSTANT	/*entropy increase due to the shear viscosity */
#ifdef NS_TIMESTEP
	SphP[i].ViscEntropyChange = 0.5 * GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  get_shear_viscosity(i) / SphP[i].d.Density * fac *
	  pow((SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);

	SphP[i].e.DtEntropy += SphP[i].ViscEntropyChange;
#else
	SphP[i].e.DtEntropy += 0.5 * GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  get_shear_viscosity(i) / SphP[i].d.Density * fac *
	  pow((SphP[i].Entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);
#endif

#else
	SphP[i].e.DtEntropy += 0.5 * GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  get_shear_viscosity(i) / SphP[i].d.Density * fac;

#ifdef NS_TIMESTEP
	SphP[i].ViscEntropyChange = 0.5 * GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  get_shear_viscosity(i) / SphP[i].d.Density * fac;
#endif

#endif

#ifdef NAVIERSTOKES_BULK	/*entropy increase due to the bulk viscosity */
	SphP[i].e.DtEntropy += GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  All.NavierStokes_BulkViscosity / SphP[i].d.Density * pow(SphP[i].u.s.a4.DivVel, 2);

#ifdef NS_TIMESTEP
	SphP[i].ViscEntropyChange = GAMMA_MINUS1 /
	  (hubble_a2 * pow(SphP[i].d.Density, GAMMA_MINUS1)) *
	  All.NavierStokes_BulkViscosity / SphP[i].d.Density * pow(SphP[i].u.s.a4.DivVel, 2);
#endif

#endif

#endif /* these entropy increases directly follow from the general heat transfer equation */


#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
	/* take care of cosmological dilution */
	if(All.ComovingIntegrationOn)
	  for(k = 0; k < 3; k++)
	    SphP[i].DtB[k] -= 2.0 * SphP[i].BPred[k];
#endif

#ifdef GASRETURN
      if(All.Time-SphP[i].GasAge < 0.00)
	{
	  SphP[i].e.DtEntropy = 0.0;
	}
#endif

#ifdef WINDS
	/* if we have winds, we decouple particles briefly if delaytime>0 */

	if(SphP[i].DelayTime > 0)
	  {
	    for(k = 0; k < 3; k++)
	      SphP[i].a.HydroAccel[k] = 0;

	    SphP[i].e.DtEntropy = 0;

#ifdef NOWINDTIMESTEPPING
	    SphP[i].MaxSignalVel = 2 * sqrt(GAMMA * SphP[i].Pressure / SphP[i].d.Density);
#else
	    windspeed = sqrt(2 * All.WindEnergyFraction * All.FactorSN *
			     All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency) * All.Time;
	    windspeed *= fac_mu;
	    hsml_c = pow(All.WindFreeTravelDensFac * All.PhysDensThresh /
			 (SphP[i].d.Density * a3inv), (1. / 3.));
	    SphP[i].MaxSignalVel = hsml_c * DMAX((2 * windspeed), SphP[i].MaxSignalVel);
#endif
	  }
#endif

#ifdef SPH_BND_PARTICLES
	if(P[i].ID == 0)
	  {
	    SphP[i].e.DtEntropy = 0;
#ifdef NS_TIMESTEP
	    SphP[i].ViscEntropyChange = 0;
#endif
	    for(k = 0; k < 3; k++)
	      SphP[i].a.HydroAccel[k] = 0;
	  }
#endif

#ifdef HEALPIX
	r_new = 0;
	for(k = 0; k < 3; k++)
	  {
	    t[k] = P[i].Pos[k] - SysState.CenterOfMassComp[1][k];
	    r_new = +t[k] * t[k];
	  }
	r_new = sqrt(r_new);
	vec2pix_nest((long) All.Nside, t, &ipix);
	if(All.healpixmap[ipix] * 0.975 < r_new)
	  {
	    SphP[i].e.DtEntropy = 0;
	    for(k = 0; k < 3; k++)
	      SphP[i].a.HydroAccel[k] = 0;
	    SphP[i].v.DivVel = 0.0;
	    for(k = 0; k < 3; k++)
	      {
		SphP[i].VelPred[k] = 0.0;
		P[i].Vel[k] = 0.0;
	      }
	    count++;
	  }
#endif

#ifdef TIME_DEP_ART_VISC
#if !defined(EOS_DEGENERATE)
	cs_h = sqrt(GAMMA * SphP[i].Pressure / SphP[i].d.Density) / PPP[i].Hsml;
#else
	cs_h = sqrt(SphP[i].dpdr) / PPP[i].Hsml;
#endif
	f = fabs(SphP[i].v.DivVel) / (fabs(SphP[i].v.DivVel) + SphP[i].r.CurlVel + 0.0001 * cs_h / fac_mu);
	SphP[i].Dtalpha = -(SphP[i].alpha - All.AlphaMin) * All.DecayTime *
	  0.5 * SphP[i].MaxSignalVel / (PPP[i].Hsml * fac_mu)
	  + f * All.ViscSource * DMAX(0.0, -SphP[i].v.DivVel);
	if(All.ComovingIntegrationOn)
	  SphP[i].Dtalpha /= (hubble_a * All.Time * All.Time);
#endif
#ifdef MAGNETIC
#ifdef TIME_DEP_MAGN_DISP
	SphP[i].DtBalpha = -(SphP[i].Balpha - All.ArtMagDispMin) * All.ArtMagDispTime *
	  0.5 * SphP[i].MaxSignalVel / (PPP[i].Hsml * fac_mu)
#ifndef ROT_IN_MAG_DIS
	  + All.ArtMagDispSource * fabs(SphP[i].divB) / sqrt(mu0 * SphP[i].d.Density);
#else
#ifdef SMOOTH_ROTB
	  + All.ArtMagDispSource / sqrt(mu0 * SphP[i].d.Density) *
	  DMAX(fabs(SphP[i].divB), fabs(sqrt(SphP[i].SmoothedRotB[0] * SphP[i].SmoothedRotB[0] +
					     SphP[i].SmoothedRotB[1] * SphP[i].SmoothedRotB[1] +
					     SphP[i].SmoothedRotB[2] * SphP[i].SmoothedRotB[2])));
#else
	  + All.ArtMagDispSource / sqrt(mu0 * SphP[i].d.Density) *
	  DMAX(fabs(SphP[i].divB), fabs(sqrt(SphP[i].RotB[0] * SphP[i].RotB[0] +
					     SphP[i].RotB[1] * SphP[i].RotB[1] +
					     SphP[i].RotB[2] * SphP[i].RotB[2])));
#endif
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
	phiphi =
	  SphP[i].PhiPred * All.DivBcleanParabolicSigma * 0.5 * SphP[i].MaxSignalVel / (PPP[i].Hsml * fac_mu);
	phiphi += All.DivBcleanHyperbolicSigma * 0.25 * SphP[i].MaxSignalVel * SphP[i].MaxSignalVel
#ifdef SMOOTH_PHI
	  * SphP[i].SmoothDivB / (fac_mu * fac_mu);
#else
	  * SphP[i].divB / (fac_mu * fac_mu);
#endif
	SphP[i].DtPhi -= phiphi;
#endif

#endif
      }


#if defined(CS_MODEL) && defined(CS_FEEDBACK)
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0 && (SphP[i].TempPromotion > 0 || SphP[i].DensPromotion > 0))
      {
	SphP[i].TempPromotion = 0;
	SphP[i].DensPromotion = 0;
      }
#endif
#ifdef HEALPIX
  MPI_Allreduce(&count, &total_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(total_count > 0)
    {
      if(ThisTask == 0)
	printf("hey %i particles where freeezed\n", total_count);
      if(total_count > 150)
	{
	  if(ThisTask == 0)
	    printf(" Next calculation of Healpix\n");
	  healpix_halo(All.healpixmap);
	}
      total_count = 0;
      fflush(stdout);
    }
#endif


#ifdef NUCLEAR_NETWORK
  tstart = second();

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
	/* evaluate network here, but do it only for temperatures > 10^7 K */
	if(SphP[i].temp > 1e7)
	  {
	    nuc_particles++;
	    network_integrate(SphP[i].temp, SphP[i].d.Density * All.UnitDensity_in_cgs, SphP[i].xnuc,
			      SphP[i].dxnuc,
			      (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval *
			      All.UnitTime_in_s, &dedt_nuc);
	    SphP[i].e.DtEntropy += dedt_nuc * All.UnitEnergy_in_cgs / All.UnitTime_in_s;
	  }
	else
	  {
	    for(k = 0; k < EOS_NSPECIES; k++)
	      {
		SphP[i].dxnuc[k] = 0;
	      }
	  }
      }

  tend = second();
  timenetwork += timediff(tstart, tend);

  MPI_Allreduce(&nuc_particles, &nuc_particles_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      printf("Nuclear network done for %d particles.\n", nuc_particles_sum);
    }

  timewait1 += timediff(tend, second());
#endif

  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_HYDCOMPUTE] += timecomp;
  CPU_Step[CPU_HYDWAIT] += timewait;
  CPU_Step[CPU_HYDCOMM] += timecomm;
  CPU_Step[CPU_HYDNETWORK] += timenetwork;
  CPU_Step[CPU_HYDMISC] += timeall - (timecomp + timewait + timecomm + timenetwork);
}




/*! This function is the 'core' of the SPH force computation. A target
*  particle is specified which may either be local, or reside in the
*  communication buffer.
*/
int hydro_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int startnode, numngb, listindex = 0;
  int j, k, n, timestep;
  MyDouble *pos;
  MyFloat *vel;
  MyFloat mass, h_i, dhsmlDensityFactor, rho, pressure, f1, f2;
  MyLongDouble acc[3], dtEntropy;

#ifdef HYDRO_COST_FACTOR
  int ninteractions = 0;
#endif


#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  MyFloat minViscousDt;
#else
  MyFloat maxSignalVel;
#endif
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, hinv, hinv4;
  double p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
  double hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
  double h_j, dwk_j;
  double r, r2, u;
  double hfc_visc;
  double dmin1, dmin2;

#ifdef TRADITIONAL_SPH_FORMULATION
  double hfc_egy;
#endif

#if defined(MAGFORCE)
  double dmax1, dmax2;
#endif
  double BulkVisc_ij;
  int imax1, imax2;

#ifdef NAVIERSTOKES
  double faci, facj;
  MyFloat *stressdiag;
  MyFloat *stressoffdiag;
  MyFloat shear_viscosity;

#ifdef VISCOSITY_SATURATION
  double VelLengthScale_i, VelLengthScale_j;
  double IonMeanFreePath_i, IonMeanFreePath_j;
#endif
#ifdef NAVIERSTOKES_BULK
  double facbi, facbj;
  MyFloat divvel;
#endif
#endif

#if defined(NAVIERSTOKES)
  double Entropy;
#endif


#ifdef TIME_DEP_ART_VISC
  MyFloat alpha;
#endif

#ifdef ALTVISCOSITY
  double mu_i, mu_j;
#endif

#ifndef NOVISCOSITYLIMITER
  double dt;
#endif

#ifdef MAGNETIC
  MyFloat *bpred;

#ifndef EULERPOTENTIALS
  double dtB[3];
#endif
  double dBx, dBy, dBz;
  double magfac, magfac_i, magfac_j, magfac_i_base;
  double mu0_1;

#ifdef MAGNETIC_DIFFUSION
  double wk_i, hinv3;
#endif

#ifdef MAGFORCE
  double mm_i[3][3], mm_j[3][3];
  double b2_i, b2_j;
  int l;
#endif

#if defined(MAGNETIC_DISSIPATION) || defined(DIVBCLEANING_DEDNER)
  double magfac_sym;
#endif
#ifdef MAGNETIC_DISSIPATION
  double dTu_diss_b, Balpha_ij;

#ifdef MAGDISSIPATION_PERPEN
  double mft, mvt[3];
#endif
#ifdef TIME_DEP_MAGN_DISP
  double Balpha;
#endif
#endif
#ifdef DIVBCLEANING_DEDNER
  double PhiPred, DtPhi, phifac;
#endif
#ifdef MAGNETIC_SIGNALVEL
  double magneticspeed_i, magneticspeed_j, vcsa2_i, vcsa2_j, Bpro2_i, Bpro2_j;
#endif
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
  MyFloat *rotb;
#endif
#endif

#ifdef PARTICLE_DEBUG
  MyIDType ID;			/*!< particle identifier */
#endif

#ifdef CONVENTIONAL_VISCOSITY
  double c_ij, h_ij;
#endif

#ifdef CS_MODEL
  double density, entropy;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = PPP[target].Hsml;
      mass = P[target].Mass;
      dhsmlDensityFactor = SphP[target].h.DhsmlDensityFactor;
      rho = SphP[target].d.Density;
      pressure = SphP[target].Pressure;
      timestep = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0);

#ifndef EOS_DEGENERATE
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#else
      soundspeed_i = sqrt(SphP[target].dpdr);
#endif

#ifdef CS_MODEL
      density = SphP[target].d.Density;
      entropy = SphP[target].Entropy;
#endif

#ifndef ALTVISCOSITY
#ifndef NAVIERSTOKES
      f1 = fabs(SphP[target].v.DivVel) /
	(fabs(SphP[target].v.DivVel) + SphP[target].r.CurlVel +
	 0.0001 * soundspeed_i / PPP[target].Hsml / fac_mu);
#else
      f1 = fabs(SphP[target].v.DivVel) /
	(fabs(SphP[target].v.DivVel) + SphP[target].u.s.CurlVel +
	 0.0001 * soundspeed_i / PPP[target].Hsml / fac_mu);
#endif
#else
      f1 = SphP[target].v.DivVel;
#endif

#ifdef MAGNETIC
      bpred = SphP[target].BPred;
#ifdef DIVBCLEANING_DEDNER
#ifdef SMOOTH_PHI
      PhiPred = SphP[target].SmoothPhi;
#else
      PhiPred = SphP[target].PhiPred;
#endif
#endif
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
#ifdef SMOOTH_ROTB
      rotb = SphP[target].SmoothedRotB;
#else
      rotb = SphP[target].RotB;
#endif
#endif
#ifdef TIME_DEP_MAGN_DISP
      Balpha = SphP[target].Balpha;
#endif
#endif /*  MAGNETIC  */

#ifdef TIME_DEP_ART_VISC
      alpha = SphP[target].alpha;
#endif

#if defined(NAVIERSTOKES)
      Entropy = SphP[target].Entropy;
#endif


#ifdef PARTICLE_DEBUG
      ID = P[target].ID;
#endif

#ifdef NAVIERSTOKES
      stressdiag = SphP[target].u.s.StressDiag;
      stressoffdiag = SphP[target].u.s.StressOffDiag;
      shear_viscosity = get_shear_viscosity(target);
#ifdef NAVIERSTOKES_BULK
      divvel = SphP[target].u.s.a4.DivVel;
#endif
#endif

    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;
      mass = HydroDataGet[target].Mass;
      dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
      rho = HydroDataGet[target].Density;
      pressure = HydroDataGet[target].Pressure;
      timestep = HydroDataGet[target].Timestep;
#ifndef EOS_DEGENERATE
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#else
      soundspeed_i = sqrt(HydroDataGet[target].dpdr);
#endif
      f1 = HydroDataGet[target].F1;

#ifdef CS_MODEL
      density = HydroDataGet[target].DensityNow;
      entropy = HydroDataGet[target].Entropy;
#endif

#ifdef MAGNETIC
      bpred = HydroDataGet[target].BPred;
#ifdef DIVBCLEANING_DEDNER
      PhiPred = HydroDataGet[target].PhiPred;
#endif
#if defined(MAGNETIC_DIFFUSION) || defined(ROT_IN_MAG_DIS)
      rotb = HydroDataGet[target].RotB;
#endif
#ifdef TIME_DEP_MAGN_DISP
      Balpha = HydroDataGet[target].Balpha;
#endif
#endif /* MAGNETIC */

#ifdef TIME_DEP_ART_VISC
      alpha = HydroDataGet[target].alpha;
#endif

#if defined(NAVIERSTOKES)
      Entropy = HydroDataGet[target].Entropy;
#endif


#ifdef PARTICLE_DEBUG
      ID = HydroDataGet[target].ID;
#endif


#ifdef NAVIERSTOKES
      stressdiag = HydroDataGet[target].stressdiag;
      stressoffdiag = HydroDataGet[target].stressoffdiag;
      shear_viscosity = HydroDataGet[target].shear_viscosity;
#endif
#ifdef NAVIERSTOKES
      stressdiag = HydroDataGet[target].stressdiag;
      stressoffdiag = HydroDataGet[target].stressoffdiag;
      shear_viscosity = HydroDataGet[target].shear_viscosity;
#ifdef NAVIERSTOKES_BULK
      divvel = HydroDataGet[target].divvel;
#endif
#endif
    }


  /* initialize variables before SPH loop is started */

  acc[0] = acc[1] = acc[2] = dtEntropy = 0;



#ifdef MAGNETIC
#ifndef EULERPOTENTIALS
  for(k = 0; k < 3; k++)
    dtB[k] = 0;
#endif
  mu0_1 = 1;
#ifndef MU0_UNITY
  mu0_1 /= (4 * M_PI);
  mu0_1 *= All.UnitTime_in_s * All.UnitTime_in_s *
    All.UnitLength_in_cm / (All.UnitMass_in_g * All.HubbleParam * All.HubbleParam);
#endif
#ifdef DIVBCLEANING_DEDNER
  DtPhi = 0;
#endif
#ifdef MAGFORCE
  magfac_i_base = 1 / (rho * rho);
#ifndef MU0_UNITY
  magfac_i_base /= (4 * M_PI);
#endif
#ifdef CORRECTBFRC
  magfac_i_base *= dhsmlDensityFactor;
#endif
  for(k = 0, b2_i = 0; k < 3; k++)
    {
      b2_i += bpred[k] * bpred[k];
      for(l = 0; l < 3; l++)
	mm_i[k][l] = bpred[k] * bpred[l];
    }
  for(k = 0; k < 3; k++)
    mm_i[k][k] -= 0.5 * b2_i;
#ifdef MAGNETIC_SIGNALVEL
#ifdef ALFVEN_VEL_LIMITER
  vcsa2_i = soundspeed_i * soundspeed_i +
    DMIN(mu0_1 * b2_i / rho, ALFVEN_VEL_LIMITER * soundspeed_i * soundspeed_i);
#else
  vcsa2_i = soundspeed_i * soundspeed_i + mu0_1 * b2_i / rho;
#endif
#endif
#endif /* end of MAGFORCE */
#endif /* end of MAGNETIC */

  p_over_rho2_i = pressure / (rho * rho);
#ifndef TRADITIONAL_SPH_FORMULATION
  p_over_rho2_i *= dhsmlDensityFactor;
#endif
  h_i2 = h_i * h_i;


#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
  minViscousDt = 1.0e32;
#else
  maxSignalVel = soundspeed_i;
#endif


  /* Now start the actual SPH computation for this particle */

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = HydroDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
#ifdef CS_MODEL
	  numngb =
	    cs_ngb_treefind_pairs(pos, h_i, target, &startnode, density, entropy, &vel[0], mode, nexport,
				  nsend_local);
#else
	  numngb = ngb_treefind_pairs(pos, h_i, target, &startnode, mode, nexport, nsend_local);
#endif

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

#ifdef HYDRO_COST_FACTOR
	      ninteractions++;
#endif

#ifdef BLACK_HOLES
	      if(P[j].Mass == 0)
		continue;
#endif

#ifdef NOWINDTIMESTEPPING
#ifdef WINDS
	      if(P[j].Type == 0)
		if(SphP[j].DelayTime > 0)	/* ignore the wind particles */
		  continue;
#endif
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
	      h_j = PPP[j].Hsml;
	      if(r2 < h_i2 || r2 < h_j * h_j)
		{
		  r = sqrt(r2);
		  if(r > 0)
		    {
		      p_over_rho2_j = SphP[j].Pressure / (SphP[j].d.Density * SphP[j].d.Density);
#ifndef EOS_DEGENERATE
		      soundspeed_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].d.Density);
#else
		      soundspeed_j = sqrt(SphP[j].dpdr);
#endif
		      dvx = vel[0] - SphP[j].VelPred[0];
		      dvy = vel[1] - SphP[j].VelPred[1];
		      dvz = vel[2] - SphP[j].VelPred[2];
		      vdotr = dx * dvx + dy * dvy + dz * dvz;
		      rho_ij = 0.5 * (rho + SphP[j].d.Density);

		      if(All.ComovingIntegrationOn)
			vdotr2 = vdotr + hubble_a2 * r2;
		      else
			vdotr2 = vdotr;

		      if(r2 < h_i2)
			{
			  hinv = 1.0 / h_i;
#ifndef  TWODIMS
			  hinv4 = hinv * hinv * hinv * hinv;
#else
			  hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
			  u = r * hinv;
			  if(u < 0.5)
			    dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
			  else
			    dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
#ifdef MAGNETIC_DIFFUSION
#ifndef  TWODIMS
			  hinv3 = hinv * hinv * hinv;
#else
			  hinv3 = hinv * hinv / boxSize_Z;
#endif
			  if(u <= 0.5)
			    wk_i = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			  else
			    wk_i = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
#endif
			}
		      else
			{
			  dwk_i = 0;
#ifdef MAGNETIC_DIFFUSION
			  wk_i = 0;
#endif
			}

		      if(r2 < h_j * h_j)
			{
			  hinv = 1.0 / h_j;
#ifndef  TWODIMS
			  hinv4 = hinv * hinv * hinv * hinv;
#else
			  hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
			  u = r * hinv;
			  if(u < 0.5)
			    dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
			  else
			    dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
			}
		      else
			{
			  dwk_j = 0;
			}

#ifdef MAGNETIC
		      dBx = bpred[0] - SphP[j].BPred[0];
		      dBy = bpred[1] - SphP[j].BPred[1];
		      dBz = bpred[2] - SphP[j].BPred[2];
		      magfac = P[j].Mass / r;	/* we moved 'dwk_i / rho' down ! */
		      if(All.ComovingIntegrationOn)
			magfac *= 1. / (hubble_a * All.Time * All.Time);
		      /* last factor takes care of all cosmological prefactor */
#ifdef CORRECTDB
		      magfac *= dhsmlDensityFactor;
#endif
#if defined(MAGNETIC_DISSIPATION) || defined(DIVBCLEANING_DEDNER)
		      magfac_sym = magfac * (dwk_i + dwk_j) * 0.5;
#endif
#ifdef MAGNETIC_DISSIPATION
#ifdef TIME_DEP_MAGN_DISP
		      Balpha_ij = 0.5 * (Balpha + SphP[j].Balpha);
#else
		      Balpha_ij = All.ArtMagDispConst;
#endif
#endif
		      magfac *= dwk_i / rho;
#ifndef EULERPOTENTIALS
		      dtB[0] +=
			magfac * ((bpred[0] * dvy - bpred[1] * dvx) * dy +
				  (bpred[0] * dvz - bpred[2] * dvx) * dz);
		      dtB[1] +=
			magfac * ((bpred[1] * dvz - bpred[2] * dvy) * dz +
				  (bpred[1] * dvx - bpred[0] * dvy) * dx);
		      dtB[2] +=
			magfac * ((bpred[2] * dvx - bpred[0] * dvz) * dx +
				  (bpred[2] * dvy - bpred[1] * dvz) * dy);
#endif
#ifdef MAGNETIC_DIFFUSION
		      magfac *= All.MagneticEta;
		      dtB[0] += magfac * (rotb[1] * dz - rotb[2] * dy);
		      dtB[1] += magfac * (rotb[2] * dx - rotb[0] * dz);
		      dtB[2] += magfac * (rotb[0] * dy - rotb[1] * dx);
		      magfac *= (r * wk_i / (dwk_i * mu0_1));
		      dtEntropy += magfac * (rotb[0] * rotb[0] + rotb[1] * rotb[1] + rotb[2] * rotb[2]);
#endif
#ifdef DIVBCLEANING_DEDNER
#ifdef SMOOTH_PHI
		      phifac = magfac_sym * (PhiPred - SphP[j].SmoothPhi) / rho;
#else
		      phifac = magfac_sym * (PhiPred - SphP[j].PhiPred) / rho;
#endif
		      dtB[0] -= phifac * dx;
		      dtB[1] -= phifac * dy;
		      dtB[2] -= phifac * dz;
#endif
#ifdef MAGFORCE
		      magfac_j = 1 / (SphP[j].d.Density * SphP[j].d.Density);
#ifndef MU0_UNITY
		      magfac_j /= (4 * M_PI);
#endif
#ifdef CORRECTBFRC
		      magfac_j *= dwk_j * SphP[j].h.DhsmlDensityFactor;
		      magfac_i = dwk_i * magfac_i_base;
#else
		      magfac_i = magfac_i_base;
#endif
		      for(k = 0, b2_j = 0; k < 3; k++)
			{
			  b2_j += SphP[j].BPred[k] * SphP[j].BPred[k];
			  for(l = 0; l < 3; l++)
			    mm_j[k][l] = SphP[j].BPred[k] * SphP[j].BPred[l];
			}
		      for(k = 0; k < 3; k++)
			mm_j[k][k] -= 0.5 * b2_j;
#ifdef MAGNETIC_SIGNALVEL
#ifdef ALFVEN_VEL_LIMITER
		      vcsa2_j = soundspeed_j * soundspeed_j +
			DMIN(mu0_1 * b2_j / SphP[j].d.Density,
			     ALFVEN_VEL_LIMITER * soundspeed_j * soundspeed_j);
#else
		      vcsa2_j = soundspeed_j * soundspeed_j + mu0_1 * b2_j / SphP[j].d.Density;
#endif
		      Bpro2_j = (SphP[j].BPred[0] * dx + SphP[j].BPred[1] * dy + SphP[j].BPred[2] * dz) / r;
		      Bpro2_j *= Bpro2_j;
		      magneticspeed_j = sqrt(vcsa2_j +
					     sqrt(DMAX((vcsa2_j * vcsa2_j -
							4 * soundspeed_j * soundspeed_j * Bpro2_j
							* mu0_1 / SphP[j].d.Density), 0))) / 1.4142136;
		      Bpro2_i = (bpred[0] * dx + bpred[1] * dy + bpred[2] * dz) / r;
		      Bpro2_i *= Bpro2_i;
		      magneticspeed_i = sqrt(vcsa2_i +
					     sqrt(DMAX((vcsa2_i * vcsa2_i -
							4 * soundspeed_i * soundspeed_i * Bpro2_i
							* mu0_1 / rho), 0))) / 1.4142136;
#endif
#ifdef MAGNETIC_DISSIPATION
		      dTu_diss_b = -magfac_sym * Balpha_ij * (dBx * dBx + dBy * dBy + dBz * dBz);
#endif
#ifdef CORRECTBFRC
		      magfac = P[j].Mass / r;
#else
		      magfac = P[j].Mass * 0.5 * (dwk_i + dwk_j) / r;
#endif
		      if(All.ComovingIntegrationOn)
			magfac *= pow(All.Time, 3 * GAMMA);
		      /* last factor takes care of all cosmological prefactor */
#ifndef MU0_UNITY
		      magfac *= All.UnitTime_in_s * All.UnitTime_in_s *
			All.UnitLength_in_cm / (All.UnitMass_in_g * All.HubbleParam * All.HubbleParam);
		      /* take care of B unit conversion into GADGET units ! */
#endif
		      for(k = 0; k < 3; k++)
			acc[k] +=
			  magfac * ((mm_i[k][0] * magfac_i + mm_j[k][0] * magfac_j) * dx +
				    (mm_i[k][1] * magfac_i + mm_j[k][1] * magfac_j) * dy +
				    (mm_i[k][2] * magfac_i + mm_j[k][2] * magfac_j) * dz);
#if defined(DIVBFORCE) && !defined(EULERPOTENTIALS)
		      for(k = 0; k < 3; k++)
			acc[k] -=
			  magfac * (((bpred[k] * bpred[0]) * magfac_i
				     + (bpred[k] * SphP[j].BPred[0]) * magfac_j) * dx
				    + ((bpred[k] * bpred[1]) * magfac_i
				       + (bpred[k] * SphP[j].BPred[1]) * magfac_j) * dy
				    + ((bpred[k] * bpred[2]) * magfac_i +
				       (bpred[k] * SphP[j].BPred[2]) * magfac_j) * dz);
#endif
#endif
#endif /* end of MAGNETIC */


#ifndef MAGNETIC_SIGNALVEL
		      vsig = soundspeed_i + soundspeed_j;
#else
		      vsig = magneticspeed_i + magneticspeed_j;
#endif


#ifndef ALTERNATIVE_VISCOUS_TIMESTEP
		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;
#endif
		      if(vdotr2 < 0)	/* ... artificial viscosity */
			{
#ifndef ALTVISCOSITY
#ifndef CONVENTIONAL_VISCOSITY
			  mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */
#else
			  c_ij = 0.5 * (soundspeed_i + soundspeed_j);
			  h_ij = 0.5 * (h_i + h_j);
			  mu_ij = fac_mu * h_ij * vdotr2 / (r2 + 0.0001 * h_ij * h_ij);
#endif
#ifdef MAGNETIC
			  vsig -= 1.5 * mu_ij;
#else
			  vsig -= 3 * mu_ij;
#endif


#ifndef ALTERNATIVE_VISCOUS_TIMESTEP
			  if(vsig > maxSignalVel)
			    maxSignalVel = vsig;
#endif

#ifndef NAVIERSTOKES
			  f2 =
			    fabs(SphP[j].v.DivVel) / (fabs(SphP[j].v.DivVel) + SphP[j].r.CurlVel +
						      0.0001 * soundspeed_j / fac_mu / PPP[j].Hsml);
#else
			  f2 =
			    fabs(SphP[j].v.DivVel) / (fabs(SphP[j].v.DivVel) + SphP[j].u.s.CurlVel +
						      0.0001 * soundspeed_j / fac_mu / PPP[j].Hsml);
#endif

#ifdef NO_SHEAR_VISCOSITY_LIMITER
			  f1 = f2 = 1;
#endif
#ifdef TIME_DEP_ART_VISC
			  BulkVisc_ij = 0.5 * (alpha + SphP[j].alpha);
#else
			  BulkVisc_ij = All.ArtBulkViscConst;
#endif
#ifndef CONVENTIONAL_VISCOSITY
			  visc = 0.25 * BulkVisc_ij * vsig * (-mu_ij) / rho_ij * (f1 + f2);
#else
			  visc =
			    (-BulkVisc_ij * mu_ij * c_ij + 2 * BulkVisc_ij * mu_ij * mu_ij) /
			    rho_ij * (f1 + f2) * 0.5;
#endif

#else /* start of ALTVISCOSITY block */
			  if(f1 < 0)
			    mu_i = h_i * fabs(f1);	/* f1 hold here the velocity divergence of particle i */
			  else
			    mu_i = 0;
			  if(SphP[j].u.s.a4.DivVel < 0)
			    mu_j = h_j * fabs(SphP[j].u.s.a4.DivVel);
			  else
			    mu_j = 0;
			  visc = All.ArtBulkViscConst * ((soundspeed_i + mu_i) * mu_i / rho +
							 (soundspeed_j + mu_j) * mu_j / SphP[j].d.Density);
#endif /* end of ALTVISCOSITY block */


			  /* .... end artificial viscosity evaluation */
			  /* now make sure that viscous acceleration is not too large */
#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
			  if(visc > 0)
			    {
			      dt = fac_vsic_fix * vdotr2 /
				(0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * visc);

			      dt /= hubble_a;

			      if(dt < minViscousDt)
				minViscousDt = dt;
			    }
#endif

#ifndef NOVISCOSITYLIMITER
			  dt =
			    2 * IMAX(timestep,
				     (P[j].TimeBin ? (1 << P[j].TimeBin) : 0)) * All.Timebase_interval;
			  if(dt > 0 && (dwk_i + dwk_j) < 0)
			    {
#ifdef BLACK_HOLES
			      if((mass + P[j].Mass) > 0)
#endif
				visc = DMIN(visc, 0.5 * fac_vsic_fix * vdotr2 /
					    (0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
			    }
#endif
			}
		      else
			{
			  visc = 0;
			}
#ifdef TRADITIONAL_SPH_FORMULATION
		      hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;


		      hfc = hfc_visc +
			0.5 * P[j].Mass * (dwk_i + dwk_j) / r * (p_over_rho2_i + p_over_rho2_j);

		      /* hfc_egy = 0.5 * P[j].Mass * (dwk_i + dwk_j) / r * (p_over_rho2_i + p_over_rho2_j); */
		      hfc_egy = P[j].Mass * (dwk_i + dwk_j) / r * (p_over_rho2_i);
#else
		      p_over_rho2_j *= SphP[j].h.DhsmlDensityFactor;

		      hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;
		      /* Formulation derived from the Lagrangian */
		      hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
#endif

#ifdef WINDS
		      if(P[j].Type == 0)
			if(SphP[j].DelayTime > 0)	/* No force by wind particles */
			  {
			    hfc = hfc_visc = 0;
			  }
#endif

#ifndef NOACCEL
		      acc[0] += FLT(-hfc * dx);
		      acc[1] += FLT(-hfc * dy);
		      acc[2] += FLT(-hfc * dz);
#endif

#if !defined(EOS_DEGENERATE) && !defined(TRADITIONAL_SPH_FORMULATION)
		      dtEntropy += FLT(0.5 * hfc_visc * vdotr2);
#else

#ifdef TRADITIONAL_SPH_FORMULATION
		      dtEntropy += FLT(0.5 * (hfc_visc + hfc_egy) * vdotr2);
#else
		      dtEntropy += FLT(0.5 * hfc * vdotr2);
#endif
#endif


#ifdef NAVIERSTOKES
		      faci = mass * shear_viscosity / (rho * rho) * dwk_i / r;

#ifndef NAVIERSTOKES_CONSTANT
		      faci *= pow((Entropy * pow(rho * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);	/*multiplied by E^5/2 */
#endif
		      facj = P[j].Mass * get_shear_viscosity(j) /
			(SphP[j].d.Density * SphP[j].d.Density) * dwk_j / r;

#ifndef NAVIERSTOKES_CONSTANT
		      facj *= pow((SphP[j].Entropy * pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.5);	/*multiplied by E^5/2 */
#endif

#ifdef NAVIERSTOKES_BULK
		      facbi = mass * All.NavierStokes_BulkViscosity / (rho * rho) * dwk_i / r;
		      facbj = P[j].Mass * All.NavierStokes_BulkViscosity /
			(SphP[j].d.Density * SphP[j].d.Density) * dwk_j / r;
#endif

#ifdef WINDS
		      if(P[j].Type == 0)
			if(SphP[j].DelayTime > 0)	/* No visc for wind particles */
			  {
			    faci = facj = 0;
#ifdef NAVIERSTOKES_BULK
			    facbi = facbj = 0;
#endif
			  }
#endif

#ifdef VISCOSITY_SATURATION
		      IonMeanFreePath_i = All.IonMeanFreePath * pow((Entropy * pow(rho * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.0) / rho;	/* u^2/rho */

		      IonMeanFreePath_j = All.IonMeanFreePath * pow((SphP[j].Entropy * pow(SphP[j].d.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.0) / SphP[j].d.Density;	/* u^2/rho */

		      for(k = 0, VelLengthScale_i = 0, VelLengthScale_j = 0; k < 3; k++)
			{
			  if(fabs(stressdiag[k]) > 0)
			    {
			      VelLengthScale_i = 2 * soundspeed_i / fabs(stressdiag[k]);

			      if(VelLengthScale_i < IonMeanFreePath_i && VelLengthScale_i > 0)
				{
				  stressdiag[k] = stressdiag[k] * (VelLengthScale_i / IonMeanFreePath_i);

				}
			    }
			  if(fabs(SphP[j].u.s.StressDiag[k]) > 0)
			    {
			      VelLengthScale_j = 2 * soundspeed_j / fabs(SphP[j].u.s.StressDiag[k]);

			      if(VelLengthScale_j < IonMeanFreePath_j && VelLengthScale_j > 0)
				{
				  SphP[j].u.s.StressDiag[k] = SphP[j].u.s.StressDiag[k] *
				    (VelLengthScale_j / IonMeanFreePath_j);

				}
			    }
			  if(fabs(stressoffdiag[k]) > 0)
			    {
			      VelLengthScale_i = 2 * soundspeed_i / fabs(stressoffdiag[k]);

			      if(VelLengthScale_i < IonMeanFreePath_i && VelLengthScale_i > 0)
				{
				  stressoffdiag[k] =
				    stressoffdiag[k] * (VelLengthScale_i / IonMeanFreePath_i);
				}
			    }
			  if(fabs(SphP[j].u.s.StressOffDiag[k]) > 0)
			    {
			      VelLengthScale_j = 2 * soundspeed_j / fabs(SphP[j].u.s.StressOffDiag[k]);

			      if(VelLengthScale_j < IonMeanFreePath_j && VelLengthScale_j > 0)
				{
				  SphP[j].u.s.StressOffDiag[k] = SphP[j].u.s.StressOffDiag[k] *
				    (VelLengthScale_j / IonMeanFreePath_j);
				}
			    }
			}
#endif

		      /* Acceleration due to the shear viscosity */
		      acc[0] += faci * (stressdiag[0] * dx + stressoffdiag[0] * dy + stressoffdiag[1] * dz)
			+ facj * (SphP[j].u.s.StressDiag[0] * dx + SphP[j].u.s.StressOffDiag[0] * dy +
				  SphP[j].u.s.StressOffDiag[1] * dz);

		      acc[1] += faci * (stressoffdiag[0] * dx + stressdiag[1] * dy + stressoffdiag[2] * dz)
			+ facj * (SphP[j].u.s.StressOffDiag[0] * dx + SphP[j].u.s.StressDiag[1] * dy +
				  SphP[j].u.s.StressOffDiag[2] * dz);

		      acc[2] += faci * (stressoffdiag[1] * dx + stressoffdiag[2] * dy + stressdiag[2] * dz)
			+ facj * (SphP[j].u.s.StressOffDiag[1] * dx + SphP[j].u.s.StressOffDiag[2] * dy +
				  SphP[j].u.s.StressDiag[2] * dz);

		      /*Acceleration due to the bulk viscosity */
#ifdef NAVIERSTOKES_BULK
#ifdef VISCOSITY_SATURATION
		      VelLengthScale_i = 0;
		      VelLengthScale_j = 0;

		      if(fabs(divvel) > 0)
			{
			  VelLengthScale_i = 3 * soundspeed_i / fabs(divvel);

			  if(VelLengthScale_i < IonMeanFreePath_i && VelLengthScale_i > 0)
			    {
			      divvel = divvel * (VelLengthScale_i / IonMeanFreePath_i);
			    }
			}

		      if(fabs(SphP[j].u.s.a4.DivVel) > 0)
			{
			  VelLengthScale_j = 3 * soundspeed_j / fabs(SphP[j].u.s.a4.DivVel);

			  if(VelLengthScale_j < IonMeanFreePath_j && VelLengthScale_j > 0)
			    {
			      SphP[j].u.s.a4.DivVel = SphP[j].u.s.a4.DivVel *
				(VelLengthScale_j / IonMeanFreePath_j);

			    }
			}
#endif


		      acc[0] += facbi * divvel * dx + facbj * SphP[j].u.s.a4.DivVel * dx;
		      acc[1] += facbi * divvel * dy + facbj * SphP[j].u.s.a4.DivVel * dy;
		      acc[2] += facbi * divvel * dz + facbj * SphP[j].u.s.a4.DivVel * dz;
#endif
#endif /* end NAVIERSTOKES */


#ifdef MAGNETIC
#ifdef MAGNETIC_DISSIPATION
		      magfac_sym *= vsig * 0.5 * Balpha_ij * r * rho / (rho_ij * rho_ij);
		      dtEntropy += dTu_diss_b * 0.25 * vsig * mu0_1 * r / (rho_ij * rho_ij);
		      dtB[0] += magfac_sym * dBx;
		      dtB[1] += magfac_sym * dBy;
		      dtB[2] += magfac_sym * dBz;
#endif
#endif

#ifdef WAKEUP
		      if(maxSignalVel > 1.1 * SphP[j].MaxSignalVel)
			{
			  SphP[j].wakeup = 1;
			}
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
	      startnode = HydroDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].a.dHydroAccel[k] = acc[k];
      SphP[target].e.dDtEntropy = dtEntropy;
#ifdef HYDRO_COST_FACTOR
      P[target].GravCost += HYDRO_COST_FACTOR * ninteractions;
#endif

#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
      SphP[target].MinViscousDt = minViscousDt;
#else
      SphP[target].MaxSignalVel = maxSignalVel;
#endif
#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
      for(k = 0; k < 3; k++)
	SphP[target].DtB[k] = dtB[k];
#ifdef DIVBCLEANING_DEDNER
      SphP[target].DtPhi = DtPhi;
#endif
#endif
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k];
      HydroDataResult[target].DtEntropy = dtEntropy;
#ifdef HYDRO_COST_FACTOR
      HydroDataResult[target].Ninteractions = ninteractions;
#endif

#ifdef ALTERNATIVE_VISCOUS_TIMESTEP
      HydroDataResult[target].MinViscousDt = minViscousDt;
#else
      HydroDataResult[target].MaxSignalVel = maxSignalVel;
#endif
#if defined(MAGNETIC) && !defined(EULERPOTENTIALS)
      for(k = 0; k < 3; k++)
	HydroDataResult[target].DtB[k] = dtB[k];
#ifdef DIVBCLEANING_DEDNER
      HydroDataResult[target].DtPhi = DtPhi;
#endif
#endif
    }

  return 0;
}
