#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>


#include "allvars.h"
#include "proto.h"
#if defined(SUBFIND_RESHUFFLE_CATALOGUE)
#include "subfind.h"
#endif

/* This function reads initial conditions that are in the default file format
 * of Gadget, i.e. snapshot files can be used as input files.  However, when a
 * snapshot file is used as input, not all the information in the header is
 * used: THE STARTING TIME NEEDS TO BE SET IN THE PARAMETERFILE.
 * Alternatively, the code can be started with restartflag==2, then snapshots
 * from the code can be used as initial conditions-files without having to
 * change the parameterfile.  For gas particles, only the internal energy is
 * read, the density and mean molecular weight will be recomputed by the code.
 * When InitGasTemp>0 is given, the gas temperature will be initialzed to this
 * value assuming a mean colecular weight either corresponding to complete
 * neutrality, or full ionization.
 */

#ifdef AUTO_SWAP_ENDIAN_READIC
int swap_file = 8;
#endif

#if defined(SAVE_HSML_IN_IC_ORDER) || defined(SUBFIND_RESHUFFLE_CATALOGUE)
static unsigned long FileNr;
static long long *NumPartPerFile;
#endif

void read_ic(char *fname)
{
  int i, num_files, rest_files, ngroups, gr, filenr, masterTask, lastTask, groupMaster;
  double u_init, molecular_weight, dmax1, dmax2;
  char buf[500];

  CPU_Step[CPU_MISC] += measure_time();

#ifdef RESCALEVINI
  if(ThisTask == 0 && RestartFlag == 0)
    {
      fprintf(stdout, "\nRescaling v_ini !\n\n");
      fflush(stdout);
    }
#endif

  NumPart = 0;
  N_gas = 0;
  All.TotNumPart = 0;

  num_files = find_files(fname);

#if defined(SAVE_HSML_IN_IC_ORDER) || defined(SUBFIND_RESHUFFLE_CATALOGUE)
  NumPartPerFile = (long long *) mymalloc(num_files * sizeof(long long));

  if(ThisTask == 0)
    get_particle_numbers(fname, num_files);

  MPI_Bcast(NumPartPerFile, num_files * sizeof(long long), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  rest_files = num_files;

  while(rest_files > NTask)
    {
      sprintf(buf, "%s.%d", fname, ThisTask + (rest_files - NTask));
      if(All.ICFormat == 3)
	sprintf(buf, "%s.%d.hdf5", fname, ThisTask + (rest_files - NTask));
#if defined(SAVE_HSML_IN_IC_ORDER) || defined(SUBFIND_RESHUFFLE_CATALOGUE)
      FileNr = ThisTask + (rest_files - NTask);
#endif

      ngroups = NTask / All.NumFilesWrittenInParallel;
      if((NTask % All.NumFilesWrittenInParallel))
	ngroups++;
      groupMaster = (ThisTask / ngroups) * ngroups;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if(ThisTask == (groupMaster + gr))	/* ok, it's this processor's turn */
	    read_file(buf, ThisTask, ThisTask);
	  MPI_Barrier(MPI_COMM_WORLD);
	}

      rest_files -= NTask;
    }


  if(rest_files > 0)
    {
      distribute_file(rest_files, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

      if(num_files > 1)
	{
	  sprintf(buf, "%s.%d", fname, filenr);
	  if(All.ICFormat == 3)
	    sprintf(buf, "%s.%d.hdf5", fname, filenr);
#if defined(SAVE_HSML_IN_IC_ORDER) || defined(SUBFIND_RESHUFFLE_CATALOGUE)
	  FileNr = filenr;
#endif
	}
      else
	{
	  sprintf(buf, "%s", fname);
	  if(All.ICFormat == 3)
	    sprintf(buf, "%s.hdf5", fname);
#if defined(SAVE_HSML_IN_IC_ORDER) || defined(SUBFIND_RESHUFFLE_CATALOGUE)
	  FileNr = 0;
#endif
	}

      ngroups = rest_files / All.NumFilesWrittenInParallel;
      if((rest_files % All.NumFilesWrittenInParallel))
	ngroups++;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
	    read_file(buf, masterTask, lastTask);
	  MPI_Barrier(MPI_COMM_WORLD);
	}
    }

#if defined(SUBFIND_RESHUFFLE_CATALOGUE)
  subfind_reshuffle_free();
#endif

  myfree_msg(CommBuffer, "CommBuffer");


  if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
    {
      /* this makes sure that masses are initialized in the case that the mass-block
         is empty for this particle type */
      for(i = 0; i < NumPart; i++)
	{
	  if(All.MassTable[P[i].Type] != 0)
	    P[i].Mass = All.MassTable[P[i].Type];
	}
    }


#ifdef GENERATE_GAS_IN_ICS
  int count, j;
  double fac, d, a, b, rho;

  if(RestartFlag == 0)
    {
      header.flag_entropy_instead_u = 0;

      for(i = 0, count = 0; i < NumPart; i++)
	if(P[i].Type == 1)
	  count++;

      memmove(P + count, P, sizeof(struct particle_data) * NumPart);

      NumPart += count;
      N_gas += count;

      if(N_gas > All.MaxPartSph)
        {
          printf("Task=%d ends up getting more SPH particles (%d) than allowed (%d)\n",
                 ThisTask, N_gas, All.MaxPartSph);
          endrun(111);
       }

      fac = All.OmegaBaryon / All.Omega0;
      rho = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      for(i = count, j = 0; i < NumPart; i++)
	if(P[i].Type == 1)
	  {
	    P[j] = P[i];

	    d = pow(P[i].Mass / rho, 1.0 / 3);
	    a = 0.5 * All.OmegaBaryon / All.Omega0 * d;
	    b = 0.5 * (All.Omega0 - All.OmegaBaryon) / All.Omega0 * d;

	    P[j].Mass *= fac;
	    P[i].Mass *= (1 - fac);
	    P[j].Type = 0;
	    P[j].ID += 1000000000;

	    P[i].Pos[0] += a;
	    P[i].Pos[1] += a;
	    P[i].Pos[2] += a;
	    P[j].Pos[0] -= b;
	    P[j].Pos[1] -= b;
	    P[j].Pos[2] -= b;

	    j++;
	  }

      All.MassTable[0] = fac * All.MassTable[1];
      All.MassTable[1] *= (1 - fac);
    }
#endif



#if defined(BLACK_HOLES) && defined(SWALLOWGAS)
  if(RestartFlag == 0)
    {
      All.MassTable[5] = 0;
    }
#endif

#ifdef SFR
  if(RestartFlag == 0)
    {
      if(All.MassTable[4] == 0 && All.MassTable[0] > 0)
	{
	  All.MassTable[0] = 0;
	  All.MassTable[4] = 0;
	}
    }
#endif


  u_init = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.InitGasTemp;
  u_init *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;	/* unit conversion */

  if(All.InitGasTemp > 1.0e4)	/* assuming FULL ionization */
    molecular_weight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  else				/* assuming NEUTRAL GAS */
    molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);

  u_init /= molecular_weight;

  All.InitGasU = u_init;



  if(RestartFlag == 0)
    {
      if(All.InitGasTemp > 0)
	{
	  for(i = 0; i < N_gas; i++)
	    {
	      if(ThisTask == 0 && i == 0 && SphP[i].Entropy == 0)
		printf("Initializing u from InitGasTemp !\n");

	      if(SphP[i].Entropy == 0)
		SphP[i].Entropy = All.InitGasU;

	      /* Note: the coversion to entropy will be done in the function init(),
	         after the densities have been computed */
	    }
	}
    }

  for(i = 0; i < N_gas; i++)
    SphP[i].Entropy = DMAX(All.MinEgySpec, SphP[i].Entropy);

#ifdef EOS_DEGENERATE
  for(i = 0; i < N_gas; i++)
    SphP[i].u = 0;
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("reading done.\n");
      fflush(stdout);
    }

  if(ThisTask == 0)
    {
      printf("Total number of particles :  %d%09d\n\n",
	     (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
      fflush(stdout);
    }

  CPU_Step[CPU_SNAPSHOT] += measure_time();
}


/*! This function reads out the buffer that was filled with particle data.
 */
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type)
{
  int n, k;
  MyInputFloat *fp;
  MyIDType *ip;
  float *fp_single;

#if defined(DISTORTIONTENSORPS) && defined(DISTORTION_READALL) && !defined(COSMIC_DISTORTION)
  int alpha, beta;
#endif

#ifdef AUTO_SWAP_ENDIAN_READIC
  int vt, vpb;
  char *cp;
#endif

  fp = (MyInputFloat *) CommBuffer;
  fp_single = (float *) CommBuffer;
  ip = (MyIDType *) CommBuffer;

#ifdef AUTO_SWAP_ENDIAN_READIC
  if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP && blocknr != IO_DMHSML_V && blocknr != IO_DMDENSITY_V)
    {
      cp = (char *) CommBuffer;
      vt = get_datatype_in_block(blocknr);
      vpb = get_values_per_blockelement(blocknr);
      if(vt == 2)
	swap_Nbyte(cp, pc * vpb, 8);
      else
	{
#ifdef INPUT_IN_DOUBLEPRECISION
	  if(vt == 1)
	    swap_Nbyte(cp, pc * vpb, 8);
	  else
#endif
	    swap_Nbyte(cp, pc * vpb, 4);
	}
    }
#endif

#ifdef COSMIC_RAYS
  int CRpop;
#endif

  switch (blocknr)
    {
    case IO_POS:		/* positions */
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  P[offset + n].Pos[k] = *fp++;

      for(n = 0; n < pc; n++)
	P[offset + n].Type = type;	/* initialize type here as well */
      break;

    case IO_VEL:		/* velocities */
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
#ifdef RESCALEVINI
	  /* scaling v to use same IC's for different cosmologies */
	  if(RestartFlag == 0)
	  P[offset + n].Vel[k] = (*fp++) * All.VelIniScale;
	  else
	    P[offset + n].Vel[k] = *fp++;
#else
	  P[offset + n].Vel[k] = *fp++;
#endif
      break;

    case IO_ID:		/* particle ID */
      for(n = 0; n < pc; n++)
	P[offset + n].ID = *ip++;
      break;

    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; n++)
	P[offset + n].Mass = *fp++;
      break;


    case IO_SHEET_ORIENTATION:	/* initial particle sheet orientation */
#if defined(DISTORTIONTENSORPS) && !defined(COSMIC_DISTORTION)
      for(n = 0; n < pc; n++)
	{
	  P[offset + n].V_matrix[0][0] = *fp++;
	  P[offset + n].V_matrix[0][1] = *fp++;
	  P[offset + n].V_matrix[0][2] = *fp++;
	  P[offset + n].V_matrix[1][0] = *fp++;
	  P[offset + n].V_matrix[1][1] = *fp++;
	  P[offset + n].V_matrix[1][2] = *fp++;
	  P[offset + n].V_matrix[2][0] = *fp++;
	  P[offset + n].V_matrix[2][1] = *fp++;
	  P[offset + n].V_matrix[2][2] = *fp++;
	}
#endif
      break;

    case IO_INIT_DENSITY:	/* initial stream density */
#if defined(DISTORTIONTENSORPS) && !defined(COSMIC_DISTORTION)
      for(n = 0; n < pc; n++)
	P[offset + n].init_density = *fp++;
      break;
#endif

    case IO_CAUSTIC_COUNTER:	/* initial caustic counter */
#if defined(DISTORTIONTENSORPS) && !defined(COSMIC_DISTORTION)
      for(n = 0; n < pc; n++)
	P[offset + n].caustic_counter = *fp++;
      break;
#endif

    case IO_DISTORTIONTENSORPS:	/* phase-space distortion tensor */
#if defined(DISTORTIONTENSORPS) && defined(DISTORTION_READALL) && !defined(COSMIC_DISTORTION)
      for(n = 0; n < pc; n++)
	{
         for (alpha = 0; alpha < 6; alpha++)
          for (beta = 0; beta < 6; beta++)
 	   P[offset + n].distortion_tensorps[alpha][beta] = *fp++;
	}

#endif
      break;

    case IO_SECONDORDERMASS:
      for(n = 0; n < pc; n++)
	{
	  P[offset + n].OldAcc = P[offset + n].Mass;	/* use this to temporarily store the masses in the 2plt IC case */
	  P[offset + n].Mass = *fp++;
	}
      break;

    case IO_U:			/* temperature */
      for(n = 0; n < pc; n++)
	SphP[offset + n].Entropy = *fp++;
      break;

    case IO_RHO:		/* density */
      for(n = 0; n < pc; n++)
	SphP[offset + n].d.Density = *fp++;
      break;

    case IO_NE:		/* electron abundance */
#ifdef COOLING
      for(n = 0; n < pc; n++)
	SphP[offset + n].Ne = *fp++;
#endif
      break;

#ifdef CHEMISTRY
    case IO_ELECT:		/* electron abundance */
      for(n = 0; n < pc; n++)
	SphP[offset + n].elec = *fp++;
      break;

    case IO_HI:		/* neutral hydrogen abundance */
      for(n = 0; n < pc; n++)
	SphP[offset + n].HI = *fp++;
      break;

    case IO_HII:		/* ionized hydrogen abundance */
      for(n = 0; n < pc; n++)
	SphP[offset + n].HII = *fp++;
      break;

    case IO_HeI:		/* neutral Helium */
      for(n = 0; n < pc; n++)
	SphP[offset + n].HeI = *fp++;
      break;

    case IO_HeII:		/* ionized Heluum */
      for(n = 0; n < pc; n++)
	SphP[offset + n].HeII = *fp++;

    case IO_HeIII:		/* double ionised Helium */
      for(n = 0; n < pc; n++)
	SphP[offset + n].HeIII = *fp++;
      break;

    case IO_H2I:		/* H2 molecule */
      for(n = 0; n < pc; n++)
	SphP[offset + n].H2I = *fp++;
      break;

    case IO_H2II:		/* ionised H2 molecule */
      for(n = 0; n < pc; n++)
	SphP[offset + n].H2II = *fp++;

    case IO_HM:		/* H minus */
      for(n = 0; n < pc; n++)
	SphP[offset + n].HM = *fp++;
      break;
#else
    case IO_ELECT:		/* electron abundance */
    case IO_HI:		/* neutral hydrogen abundance */
    case IO_HII:		/* ionized hydrogen abundance */
    case IO_HeI:		/* neutral Helium */
    case IO_HeII:		/* ionized Heluum */
    case IO_HeIII:		/* double ionised Helium */
    case IO_H2I:		/* H2 molecule */
    case IO_H2II:		/* ionised H2 molecule */
    case IO_HM:		/* H minus */
      break;
#endif

    case IO_HSML:		/* SPH smoothing length */
      for(n = 0; n < pc; n++)
	PPP[offset + n].Hsml = *fp++;
      break;


    case IO_AGE:		/* Age of stars */
#ifdef STELLARAGE
      for(n = 0; n < pc; n++)
	P[offset + n].StellarAge = *fp++;
#endif
      break;

    case IO_Z:			/* Gas and star metallicity */
#ifdef METALS
#ifndef CS_MODEL
      for(n = 0; n < pc; n++)
	P[offset + n].Metallicity = *fp++;
#else
      for(n = 0; n < pc; n++)
	for(k = 0; k < 12; k++)
	  P[offset + n].Zm[k] = *fp++;
#endif
#endif
      break;

    case IO_EGYPROM:		/* SN Energy Reservoir */
#ifdef CS_FEEDBACK
      for(n = 0; n < pc; n++)
	P[offset + n].EnergySN = *fp++;
#endif
      break;

    case IO_EGYCOLD:		/* Cold  SN Energy Reservoir */
#ifdef CS_FEEDBACK
      for(n = 0; n < pc; n++)
	P[offset + n].EnergySNCold = *fp++;
#endif
      break;

    case IO_BFLD:		/* Magnetic field */
#ifdef MAGNETIC
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  SphP[offset + n].BPred[k] = *fp++;
#ifdef TRACEDIVB
	  SphP[offset + n].divB=0;
#endif
#endif
      break;

    case IO_CR_C0:		/* Adiabatic invariant for cosmic rays */
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; n++)
	for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	  SphP[offset + n].CR_C0[CRpop] = *fp++;
#endif
      break;

    case IO_CR_Q0:		/* Adiabatic invariant for cosmic rays */
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; n++)
	for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	  SphP[offset + n].CR_q0[CRpop] = *fp++;
#endif
      break;

    case IO_CR_P0:
      break;

    case IO_CR_E0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; n++)
	for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	  SphP[offset + n].CR_E0[CRpop] = *fp++;
#endif
      break;

    case IO_CR_n0:
#ifdef COSMIC_RAYS
      for(n = 0; n < pc; n++)
	for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	  SphP[offset + n].CR_n0[CRpop] = *fp++;
#endif
      break;

    case IO_CR_ThermalizationTime:
    case IO_CR_DissipationTime:
      break;

    case IO_BHMASS:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Mass = *fp++;
#endif
      break;

    case IO_BHMDOT:
#ifdef BLACK_HOLES
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Mdot = *fp++;
#endif
      break;

    case IO_BHPROGS:
#ifdef BH_COUNTPROGS
      for(n = 0; n < pc; n++)
	P[offset + n].BH_CountProgs = *ip++;
#endif
      break;

    case IO_BHMBUB:
#ifdef BH_BUBBLES
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Mass_bubbles = *fp++;
#endif
      break;

    case IO_BHMINI:
#ifdef BH_BUBBLES
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Mass_ini = *fp++;
#endif
      break;

    case IO_BHMRAD:
#ifdef UNIFIED_FEEDBACK
      for(n = 0; n < pc; n++)
	P[offset + n].BH_Mass_radio = *fp++;
#endif
      break;

    case IO_EOSXNUC:
#ifdef EOS_DEGENERATE
      for(n = 0; n < pc; n++)
	for(k = 0; k < EOS_NSPECIES; k++)
	  SphP[offset + n].xnuc[k] = *fp++;
#endif
      break;

    case IO_DMHSML:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].DM_Hsml = *fp_single++;
#endif
      break;

    case IO_DMDENSITY:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].u.DM_Density = *fp_single++;
#endif
      break;

    case IO_DMVELDISP:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].v.DM_VelDisp = *fp_single++;
#endif
      break;

    case IO_DMHSML_V:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].DM_Hsml_V = *fp_single++;
#endif
      break;

    case IO_DMDENSITY_V:
#if defined(SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI) && defined(SUBFIND)
      for(n = 0; n < pc; n++)
	P[offset + n].DM_Density_V = *fp_single++;
#endif
      break;

      /* the other input fields (if present) are not needed to define the 
         initial conditions of the code */

    case IO_NH:
    case IO_SFR:
    case IO_POT:
    case IO_ACCEL:
    case IO_DTENTR:
    case IO_STRESSDIAG:
    case IO_STRESSOFFDIAG:
    case IO_STRESSBULK:
    case IO_SHEARCOEFF:
    case IO_TSTP:
    case IO_DBDT:
    case IO_DIVB:
    case IO_ABVC:
    case IO_COOLRATE:
    case IO_CONDRATE:
    case IO_BSMTH:
    case IO_DENN:
    case IO_MACH:
    case IO_DTENERGY:
    case IO_PRESHOCK_DENSITY:
    case IO_PRESHOCK_ENERGY:
    case IO_PRESHOCK_XCR:
    case IO_DENSITY_JUMP:
    case IO_ENERGY_JUMP:
    case IO_CRINJECT:
    case IO_AMDC:
    case IO_PHI:
    case IO_TIDALTENSORPS:
    case IO_ROTB:
    case IO_SROTB:
    case IO_EULERA:
    case IO_EULERB:
    case IO_FLOW_DETERMINANT:    
    case IO_STREAM_DENSITY:
    case IO_PHASE_SPACE_DETERMINANT:
    case IO_ANNIHILATION_RADIATION:
    case IO_EOSTEMP:
    case IO_PRESSURE:
    case IO_PRESHOCK_CSND:
    case IO_nHII:
    case IO_RADGAMMA:
    case IO_SHELL_INFO:
    case IO_LAST_CAUSTIC:
      break;

    case IO_LASTENTRY:
      endrun(220);
      break;
    }
}



/*! This function reads a snapshot file and distributes the data it contains
 *  to tasks 'readTask' to 'lastTask'.
 */
void read_file(char *fname, int readTask, int lastTask)
{
  int blockmaxlen;
  int i, n_in_file, n_for_this_task, ntask, pc, offset = 0, task;
  int blksize1, blksize2;
  MPI_Status status;
  FILE *fd = 0;
  int nall, nread;
  int type, bnr;
  char label[4], expected_label[4], buf[500];
  int nstart, bytes_per_blockelement, npart, nextblock, typelist[6];
  enum iofields blocknr;
  size_t bytes;

#ifdef HAVE_HDF5
  int rank, pcsum;
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_dataspace_in_file;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_memory, hdf5_dataset;
  hsize_t dims[2], count[2], start[2];
#endif

#if defined(COSMIC_RAYS) && (!defined(CR_IC))
  int CRpop;
#endif

#define SKIP  {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
	{
	  if(!(fd = fopen(fname, "r")))
	    {
	      printf("can't open file `%s' for reading initial conditions.\n", fname);
	      endrun(123);
	    }

	  if(All.ICFormat == 2)
	    {
	      SKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_file = blksize1;
#endif
	      my_fread(&label, sizeof(char), 4, fd);
	      my_fread(&nextblock, sizeof(int), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_Nbyte((char *) &nextblock, 1, 4);
#endif
	      printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3],
		     nextblock);
	      SKIP2;
	    }

	  SKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
	  if(All.ICFormat == 1)
	    {
	      if(blksize1 != 256)
		swap_file = 1;
	    }
#endif
	  my_fread(&header, sizeof(header), 1, fd);
	  SKIP2;
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_Nbyte((char *) &blksize1, 1, 4);
	  swap_Nbyte((char *) &blksize2, 1, 4);
#endif

	  if(blksize1 != 256 || blksize2 != 256)
	    {
	      printf("incorrect header format\n");
	      fflush(stdout);
	      endrun(890);
	      /* Probable error is wrong size of fill[] in header file. Needs to be 256 bytes in total. */
	    }
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_header();
#endif
	}


#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
	{
	  read_header_attributes_in_hdf5(fname);

	  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

	  for(type = 0; type < 6; type++)
	    {
	      if(header.npart[type] > 0)
		{
		  sprintf(buf, "/PartType%d", type);
		  hdf5_grp[type] = H5Gopen(hdf5_file, buf);
		}
	    }
	}
#endif

      for(task = readTask + 1; task <= lastTask; task++)
        {
	  MPI_Ssend(&header, sizeof(header), MPI_BYTE, task, TAG_HEADER, MPI_COMM_WORLD);
#ifdef AUTO_SWAP_ENDIAN_READIC
	  MPI_Ssend(&swap_file, sizeof(int), MPI_INT, task, TAG_SWAP, MPI_COMM_WORLD);
#endif
        }

    }
  else
    {
      MPI_Recv(&header, sizeof(header), MPI_BYTE, readTask, TAG_HEADER, MPI_COMM_WORLD, &status);
#ifdef AUTO_SWAP_ENDIAN_READIC
      MPI_Recv(&swap_file, sizeof(int), MPI_INT, readTask, TAG_SWAP, MPI_COMM_WORLD, &status);
#endif
    }

#ifdef INPUT_IN_DOUBLEPRECISION
  if(header.flag_doubleprecision == 0)
    {
      if(ThisTask == 0)
	printf
	  ("\nProblem: Code compiled with INPUT_IN_DOUBLEPRECISION, but input files are in single precision!\n");
      endrun(11);
    }
#else
  if(header.flag_doubleprecision)
    {
      if(ThisTask == 0)
	printf
	  ("\nProblem: Code not compiled with INPUT_IN_DOUBLEPRECISION, but input files are in double precision!\n");
      endrun(10);
    }
#endif


  if(All.TotNumPart == 0)
    {
      if(header.num_files <= 1)
	for(i = 0; i < 6; i++)
	  {
	    header.npartTotal[i] = header.npart[i];
#ifdef SFR
	    header.npartTotalHighWord[i] = 0;
#endif
	  }

      All.TotN_gas = header.npartTotal[0] + (((long long) header.npartTotalHighWord[0]) << 32);

      for(i = 0, All.TotNumPart = 0; i < 6; i++)
	{
	  All.TotNumPart += header.npartTotal[i];
	  All.TotNumPart += (((long long) header.npartTotalHighWord[i]) << 32);
	}

#ifdef GENERATE_GAS_IN_ICS
      if(RestartFlag == 0)
	{
	  All.TotN_gas += header.npartTotal[1];
	  All.TotNumPart += header.npartTotal[1];
	}
#endif

      for(i = 0; i < 6; i++)
	All.MassTable[i] = header.mass[i];

      All.MaxPart = (int) (All.PartAllocFactor * (All.TotNumPart / NTask));	/* sets the maximum number of particles that may */
#ifdef GASRETURN
      All.MaxPartSph = (int) (All.PartAllocFactor * (All.TotNumPart / NTask));	/* sets the maximum number of particles that may */
#else
      All.MaxPartSph = (int) (All.PartAllocFactor * (All.TotN_gas / NTask));	/* sets the maximum number of particles that may */
#endif

#ifdef INHOMOG_GASDISTR_HINT
      All.MaxPartSph = All.MaxPart;
#endif

      allocate_memory();

      if(!(CommBuffer = mymalloc(bytes = All.BufferSize * 1024 * 1024)))
	{
	  printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(2);
	}

      if(RestartFlag >= 2)
	All.Time = All.TimeBegin = header.time;
    }

  if(ThisTask == readTask)
    {
      for(i = 0, n_in_file = 0; i < 6; i++)
	n_in_file += header.npart[i];

      printf("\nreading file `%s' on task=%d (contains %d particles.)\n"
	     "distributing this file to tasks %d-%d\n"
	     "Type 0 (gas):   %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 1 (halo):  %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 2 (disk):  %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 3 (bulge): %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 4 (stars): %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 5 (bndry): %8d  (tot=%6d%09d) masstab=%g\n\n", fname, ThisTask, n_in_file, readTask,
	     lastTask, header.npart[0], (int) (header.npartTotal[0] / 1000000000),
	     (int) (header.npartTotal[0] % 1000000000), All.MassTable[0], header.npart[1],
	     (int) (header.npartTotal[1] / 1000000000), (int) (header.npartTotal[1] % 1000000000),
	     All.MassTable[1], header.npart[2], (int) (header.npartTotal[2] / 1000000000),
	     (int) (header.npartTotal[2] % 1000000000), All.MassTable[2], header.npart[3],
	     (int) (header.npartTotal[3] / 1000000000), (int) (header.npartTotal[3] % 1000000000),
	     All.MassTable[3], header.npart[4], (int) (header.npartTotal[4] / 1000000000),
	     (int) (header.npartTotal[4] % 1000000000), All.MassTable[4], header.npart[5],
	     (int) (header.npartTotal[5] / 1000000000), (int) (header.npartTotal[5] % 1000000000),
	     All.MassTable[5]);
      fflush(stdout);
    }


  ntask = lastTask - readTask + 1;


  /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */

  for(type = 0, nall = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
	n_for_this_task++;


      if(type == 0)
	{
	  if(N_gas + n_for_this_task > All.MaxPartSph)
	    {
	      printf("Not enough space on task=%d for SPH particles (space for %d, need at least %d)\n",
		     ThisTask, All.MaxPartSph, N_gas + n_for_this_task);
	      fflush(stdout);
	      endrun(172);
	    }
	}

      nall += n_for_this_task;
    }

  if(NumPart + nall > All.MaxPart)
    {
      printf("Not enough space on task=%d (space for %d, need at least %d)\n",
	     ThisTask, All.MaxPart, NumPart + nall);
      fflush(stdout);
      endrun(173);
    }

  memmove(&P[N_gas + nall], &P[N_gas], (NumPart - N_gas) * sizeof(struct particle_data));
  nstart = N_gas;



  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum iofields) bnr;

      if(blocknr == IO_LASTENTRY)
	break;

      if(blockpresent(blocknr))
	{
#ifdef CR_IC
	  if(RestartFlag == 0 && ((blocknr > IO_CR_Q0 && blocknr != IO_BFLD)
				  || (blocknr >= IO_RHO && blocknr <= IO_ACCEL)))
#else
#ifdef EOS_DEGENERATE
	  if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_EOSXNUC))
#else
#ifndef CHEMISTRY
#ifndef READ_HSML
            /* normal */
          if(RestartFlag == 0 && blocknr > IO_U && blocknr != IO_BFLD && blocknr != IO_Z && blocknr != IO_AGE)
#else
	  if(RestartFlag == 0 && blocknr > IO_U && blocknr != IO_BFLD && blocknr != IO_HSML)
#endif
#else
	  if(RestartFlag == 0 && blocknr > IO_HM)
#endif
#endif
#endif
#if defined(DISTORTIONTENSORPS) && !defined(COSMIC_DISTORTION)
	  if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_SHEET_ORIENTATION))
          if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_INIT_DENSITY))
    	  if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_CAUSTIC_COUNTER))
#ifdef DISTORTION_READALL
    	  if(RestartFlag == 0 && (blocknr > IO_U && blocknr != IO_DISTORTIONTENSORPS))
#endif
#endif

		  continue;	/* ignore all other blocks in initial conditions */




#ifdef BINISET
	  if(RestartFlag == 0 && blocknr == IO_BFLD)
	    continue;
#endif
	  if(ThisTask == readTask)
	    {
	      get_dataset_name(blocknr, buf);
	      printf("reading block %d (%s)...\n", blocknr, buf);
	      fflush(stdout);
	    }

	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 1);

	  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / bytes_per_blockelement;

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP && blocknr != IO_DMHSML_V && blocknr != IO_DMDENSITY_V)
		if(ThisTask == readTask)
		  {
		    if(All.ICFormat == 2)
		      {
			SKIP;
			my_fread(&label, sizeof(char), 4, fd);
			my_fread(&nextblock, sizeof(int), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
			swap_Nbyte((char *) &nextblock, 1, 4);
#endif
			printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2],
			       label[3], nextblock);
			SKIP2;

			get_Tab_IO_Label(blocknr, expected_label);
			if(strncmp(label, expected_label, 4) != 0)
			  {
			    printf("incorrect block-structure!\n");
			    printf("expected '%c%c%c%c' but found '%c%c%c%c'\n",
				   label[0], label[1], label[2], label[3],
				   expected_label[0], expected_label[1], expected_label[2],
				   expected_label[3]);
			    fflush(stdout);
			    endrun(1890);
			  }
		      }

		    if(All.ICFormat == 1 || All.ICFormat == 2)
		      SKIP;
		  }

	      for(type = 0, offset = 0, nread = 0; type < 6; type++)
		{
		  n_in_file = header.npart[type];
#ifdef HAVE_HDF5
		  pcsum = 0;
#endif
		  if(typelist[type] == 0)
		    {
		      n_for_this_task = n_in_file / ntask;
		      if((ThisTask - readTask) < (n_in_file % ntask))
			n_for_this_task++;

		      offset += n_for_this_task;
		    }
		  else
		    {
		      for(task = readTask; task <= lastTask; task++)
			{
			  n_for_this_task = n_in_file / ntask;
			  if((task - readTask) < (n_in_file % ntask))
			    n_for_this_task++;

			  if(task == ThisTask)
			    if(NumPart + n_for_this_task > All.MaxPart)
			      {
				printf("too many particles. %d %d %d\n", NumPart, n_for_this_task,
				       All.MaxPart);
				endrun(1313);
			      }

			  do
			    {
			      pc = n_for_this_task;

			      if(pc > blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == readTask)
				{
				  if(All.ICFormat == 1 || All.ICFormat == 2)
				    {
				      if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP && blocknr != IO_DMHSML_V && blocknr != IO_DMDENSITY_V)
					{
					  my_fread(CommBuffer, bytes_per_blockelement, pc, fd);
					  nread += pc;
					}
				      else
					{
#ifdef SUBFIND_RESHUFFLE_CATALOGUE
					  read_hsml_files(CommBuffer, pc, blocknr,
							  NumPartPerFile[FileNr] + nread);
#endif
					  nread += pc;
					}
				    }

#ifdef HAVE_HDF5
				  if(All.ICFormat == 3 && pc > 0)
				    {
				      get_dataset_name(blocknr, buf);
				      hdf5_dataset = H5Dopen(hdf5_grp[type], buf);

				      dims[0] = header.npart[type];
				      dims[1] = get_values_per_blockelement(blocknr);
				      if(dims[1] == 1)
					rank = 1;
				      else
					rank = 2;

				      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);

				      dims[0] = pc;
				      hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);

				      start[0] = pcsum;
				      start[1] = 0;

				      count[0] = pc;
				      count[1] = get_values_per_blockelement(blocknr);
				      pcsum += pc;

				      H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
							  start, NULL, count, NULL);

				      switch (get_datatype_in_block(blocknr))
					{
					case 0:
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
					  break;
					case 1:
#ifdef INPUT_IN_DOUBLEPRECISION
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
					  break;
					case 2:
					  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
					  break;
					}

				      H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory,
					      hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

				      H5Tclose(hdf5_datatype);
				      H5Sclose(hdf5_dataspace_in_memory);
				      H5Sclose(hdf5_dataspace_in_file);
				      H5Dclose(hdf5_dataset);
				    }
#endif
				}

			      if(ThisTask == readTask && task != readTask && pc > 0)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA,
					  MPI_COMM_WORLD);

			      if(ThisTask != readTask && task == ThisTask && pc > 0)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, readTask,
					 TAG_PDATA, MPI_COMM_WORLD, &status);

			      if(ThisTask == task)
				{
				  empty_read_buffer(blocknr, nstart + offset, pc, type);

				  offset += pc;
				}

			      n_for_this_task -= pc;
			    }
			  while(n_for_this_task > 0);
			}
		    }
		}
	      if(ThisTask == readTask)
		{
		  if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP && blocknr != IO_DMHSML_V && blocknr != IO_DMDENSITY_V)
		    if(All.ICFormat == 1 || All.ICFormat == 2)
		      {
			SKIP2;
#ifdef AUTO_SWAP_ENDIAN_READIC
			swap_Nbyte((char *) &blksize1, 1, 4);
			swap_Nbyte((char *) &blksize2, 1, 4);
#endif
			if(blksize1 != blksize2)
			  {
			    printf("incorrect block-sizes detected!\n");
			    printf("Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask, blocknr,
				   blksize1, blksize2);
			    if(blocknr == IO_ID)
			      {
				printf
				  ("Possible mismatch of 32bit and 64bit ID's in IC file and GADGET compilation !\n");
			      }
			    fflush(stdout);
			    endrun(1889);
			  }
		      }
		}
	    }
	}
    }


#ifdef SAVE_HSML_IN_IC_ORDER
  MyIDType IdCount = 0;

  for(type = 0, offset = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      for(task = readTask; task <= lastTask; task++)
	{
	  n_for_this_task = n_in_file / ntask;
	  if((task - readTask) < (n_in_file % ntask))
	    n_for_this_task++;

	  if(ThisTask == task)
	    {
	      int i;

	      for(i = 0; i < n_for_this_task; i++)
		P[nstart + offset + i].ID_ic_order = NumPartPerFile[FileNr] + IdCount + i;

	      offset += n_for_this_task;
	    }

	  IdCount += n_for_this_task;
	}
    }
#endif


  for(type = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
	n_for_this_task++;

      NumPart += n_for_this_task;

      if(type == 0)
	N_gas += n_for_this_task;
    }

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
	fclose(fd);
#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
	{
	  for(type = 5; type >= 0; type--)
	    if(header.npart[type] > 0)
	      H5Gclose(hdf5_grp[type]);
	  H5Fclose(hdf5_file);
	}
#endif
    }

#if defined(COSMIC_RAYS) && (!defined(CR_IC))
  for(i = 0; i < n_for_this_task; i++)
    {
      if(P[i].Type != 0)
	{
	  break;
	}

      for(CRpop = 0; CRpop < NUMCRPOP; CRpop++)
	{
	  SphP[i].CR_C0[CRpop] = 0.0;
	  SphP[i].CR_q0[CRpop] = 1.0e10;
	}
    }
#endif

}



/*! This function determines on how many files a given snapshot is distributed.
 */
int find_files(char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;

  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if(All.ICFormat == 3)
    {
      sprintf(buf, "%s.%d.hdf5", fname, 0);
      sprintf(buf1, "%s.hdf5", fname);
    }

#ifndef  HAVE_HDF5
  if(All.ICFormat == 3)
    {
      if(ThisTask == 0)
	printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun(0);
    }
#endif

  header.num_files = 0;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf, "r")))
	{
	  if(All.ICFormat == 1 || All.ICFormat == 2)
	    {
	      if(All.ICFormat == 2)
		{
		  my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		  swap_file = dummy;
#endif
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		}

	      my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      if(All.SnapFormat == 1)
		{
		  if(dummy == 256)
		    swap_file = 8;
		  else
		    swap_file = dummy;
		}
#endif
	      my_fread(&header, sizeof(header), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_header();
#endif
	      my_fread(&dummy, sizeof(dummy), 1, fd);
	    }
	  fclose(fd);

#ifdef HAVE_HDF5
	  if(All.ICFormat == 3)
	    read_header_attributes_in_hdf5(buf);
#endif
	}
    }

#ifdef AUTO_SWAP_ENDIAN_READIC
  MPI_Bcast(&swap_file, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
#endif
  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf1, "r")))
	{
	  if(All.ICFormat == 1 || All.ICFormat == 2)
	    {
	      if(All.ICFormat == 2)
		{
		  my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		  swap_file = dummy;
#endif
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		}

	      my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      if(All.SnapFormat == 1)
		{
		  if(dummy == 256)
		    swap_file = 8;
		  else
		    swap_file = dummy;
		}
#endif
	      my_fread(&header, sizeof(header), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_header();
#endif
	      my_fread(&dummy, sizeof(dummy), 1, fd);
	    }
	  fclose(fd);

#ifdef HAVE_HDF5
	  if(All.ICFormat == 3)
	    read_header_attributes_in_hdf5(buf1);
#endif

	  header.num_files = 1;
	}
    }

#ifdef AUTO_SWAP_ENDIAN_READIC
  MPI_Bcast(&swap_file, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
#endif
  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      printf("\nCan't find initial conditions file.");
      printf("neither as '%s'\nnor as '%s'\n", buf, buf1);
      fflush(stdout);
    }

  endrun(0);
  return 0;
}

#if defined(SAVE_HSML_IN_IC_ORDER) || defined(SUBFIND_RESHUFFLE_CATALOGUE)
void get_particle_numbers(char *fname, int num_files)
{
  char buf[1000];
  int blksize1, blksize2;
  char label[4];
  int nextblock;
  int i, j;

  printf("num_files=%d\n", num_files);

  for(i = 0; i < num_files; i++)
    {
      if(num_files > 1)
	{
	  sprintf(buf, "%s.%d", fname, i);
	  if(All.ICFormat == 3)
	    sprintf(buf, "%s.%d.hdf5", fname, i);
	}
      else
	{
	  sprintf(buf, "%s", fname);
	  if(All.ICFormat == 3)
	    sprintf(buf, "%s.hdf5", fname);
	}

#define SKIP  {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}

      if(All.ICFormat == 1 || All.ICFormat == 2)
	{
	  FILE *fd;

	  if(!(fd = fopen(buf, "r")))
	    {
	      printf("can't open file `%s' for reading initial conditions.\n", buf);
	      endrun(1239);
	    }

	  if(All.ICFormat == 2)
	    {
	      SKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_file = blksize1;
#endif
	      my_fread(&label, sizeof(char), 4, fd);
	      my_fread(&nextblock, sizeof(int), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_Nbyte((char *) &nextblock, 1, 4);
#endif
	      SKIP2;
	    }

	  SKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
	  if(All.ICFormat == 1)
	    {
	      if(blksize1 != 256)
		swap_file = 1;
	    }
#endif
	  my_fread(&header, sizeof(header), 1, fd);
	  SKIP2;
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_Nbyte((char *) &blksize1, 1, 4);
	  swap_Nbyte((char *) &blksize2, 1, 4);
#endif

	  if(blksize1 != 256 || blksize2 != 256)
	    {
	      printf("incorrect header format\n");
	      fflush(stdout);
	      endrun(890);
	    }
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_header();
#endif
	  fclose(fd);
	}

#ifdef HAVE_HDF5
      if(All.ICFormat == 3)
	{
	  read_header_attributes_in_hdf5(buf);
	}
#endif

      NumPartPerFile[i] = 0;

      for(j = 0; j < 6; j++)
	{
#if defined(SUBFIND_RESHUFFLE_CATALOGUE)
	  if(((1 << j) & (FOF_PRIMARY_LINK_TYPES)))
#endif
	    NumPartPerFile[i] += header.npart[j];
	}

      printf("File=%4d:  NumPart= %d\n", i, (int) (NumPartPerFile[i]));
    }


  long long n, sum;

  for(i = 0, sum = 0; i < num_files; i++)
    {
      n = NumPartPerFile[i];

      NumPartPerFile[i] = sum;

      sum += n;
    }
}
#endif




/*! This function assigns a certain number of files to processors, such that
 *  each processor is exactly assigned to one file, and the number of cpus per
 *  file is as homogenous as possible. The number of files may at most be
 *  equal to the number of processors.
 */
void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last)
{
  int ntask, filesleft, filesright, tasksleft, tasksright;

  if(nfiles > 1)
    {
      ntask = lasttask - firsttask + 1;

      filesleft = (int) ((((double) (ntask / 2)) / ntask) * nfiles);
      if(filesleft <= 0)
	filesleft = 1;
      if(filesleft >= nfiles)
	filesleft = nfiles - 1;

      filesright = nfiles - filesleft;

      tasksleft = ntask / 2;
      tasksright = ntask - tasksleft;

      distribute_file(filesleft, firstfile, firsttask, firsttask + tasksleft - 1, filenr, master, last);
      distribute_file(filesright, firstfile + filesleft, firsttask + tasksleft, lasttask, filenr, master,
		      last);
    }
  else
    {
      if(ThisTask >= firsttask && ThisTask <= lasttask)
	{
	  *filenr = firstfile;
	  *master = firsttask;
	  *last = lasttask;
	}
    }
}



#ifdef HAVE_HDF5
void read_header_attributes_in_hdf5(char *fname)
{
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;

  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header");

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npart);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);

#ifdef LT_STELLAREVOLUTION
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "BoxSize");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
  H5Aclose(hdf5_attribute);
#endif

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_IC_Info");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_ic_info);
  H5Aclose(hdf5_attribute);


  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);
}
#endif





#ifdef AUTO_SWAP_ENDIAN_READIC
/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_Nbyte(char *data, int n, int m)
{
  int i, j;
  char old_data[16];

  if(swap_file != 8)
    {
      for(j = 0; j < n; j++)
	{
	  memcpy(&old_data[0], &data[j * m], m);
	  for(i = 0; i < m; i++)
	    {
	      data[j * m + i] = old_data[m - i - 1];
	    }
	}
    }
}

/*------------------------------------------------------------------*/
/*----------- procedure to swap header if needed -------------------*/
/*------------------------------------------------------------------*/

void swap_header()
{
  swap_Nbyte((char *) &header.npart, 6, 4);
  swap_Nbyte((char *) &header.mass, 6, 8);
  swap_Nbyte((char *) &header.time, 1, 8);
  swap_Nbyte((char *) &header.redshift, 1, 8);
  swap_Nbyte((char *) &header.flag_sfr, 1, 4);
  swap_Nbyte((char *) &header.flag_feedback, 1, 4);
  swap_Nbyte((char *) &header.npartTotal, 6, 4);
  swap_Nbyte((char *) &header.flag_cooling, 1, 4);
  swap_Nbyte((char *) &header.num_files, 1, 4);
  swap_Nbyte((char *) &header.BoxSize, 1, 8);
  swap_Nbyte((char *) &header.Omega0, 1, 8);
  swap_Nbyte((char *) &header.OmegaLambda, 1, 8);
  swap_Nbyte((char *) &header.HubbleParam, 1, 8);
  swap_Nbyte((char *) &header.flag_stellarage, 1, 4);
  swap_Nbyte((char *) &header.flag_metals, 1, 4);
  swap_Nbyte((char *) &header.npartTotalHighWord, 6, 4);
  swap_Nbyte((char *) &header.flag_entropy_instead_u, 1, 4);
  swap_Nbyte((char *) &header.flag_doubleprecision, 1, 4);
#ifdef COSMIC_RAYS 
  swap_Nbyte((char *) &header.SpectralIndex_CR_Pop, NUMCRPOP, 8);
#endif
}

#endif
