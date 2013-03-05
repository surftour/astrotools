#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*! \file pm_periodic.c
 *  \brief routines for periodic PM-force computation
 */

#ifdef PMGRID
#ifdef PERIODIC

#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif


#include "allvars.h"
#include "proto.h"

#define  PMGRID2 (2*(PMGRID/2 + 1))

#if (PMGRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif

#ifdef FLTROUNDOFFREDUCTION
#define d_fftw_real MyLongDouble
#else
#define d_fftw_real fftw_real
#endif

static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;

static int slab_to_task[PMGRID];
static int *slabs_per_task;
static int *first_slab_of_task;

static int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

static int fftsize, maxfftsize;

static fftw_real *rhogrid, *forcegrid, *workspace;
static d_fftw_real *d_rhogrid, *d_forcegrid, *d_workspace;

static fftw_complex *fft_of_rhogrid;


static MyFloat to_slab_fac;

void pm_periodic_transposeA(fftw_real * field, fftw_real * scratch);
void pm_periodic_transposeB(fftw_real * field, fftw_real * scratch);
int pm_periodic_compare_sortindex(const void *a, const void *b);

static struct part_slab_data
{
  large_array_offset globalindex;
  int partindex;
  int localindex;
} *part;

static int *part_sortindex;


/*! This routines generates the FFTW-plans to carry out the parallel FFTs
 *  later on. Some auxiliary variables are also initialized.
 */
void pm_init_periodic(void)
{
  int i;
  int slab_to_task_local[PMGRID];

  All.Asmth[0] = ASMTH * All.BoxSize / PMGRID;
  All.Rcut[0] = RCUT * All.Asmth[0];

  /* Set up the FFTW plan files. */

  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
					     FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
					     FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  /* Workspace out the ranges on each processor. */

  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);

  for(i = 0; i < PMGRID; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, slab_to_task, PMGRID, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&nslab_x, &smallest_slab, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  slabs_per_task = (int *) mymalloc(NTask * sizeof(int));
  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  first_slab_of_task = (int *) mymalloc(NTask * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  to_slab_fac = PMGRID / All.BoxSize;

  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
}


/*! This function allocates the memory neeed to compute the long-range PM
 *  force. Three fields are used, one to hold the density (and its FFT, and
 *  then the real-space potential), one to hold the force field obtained by
 *  finite differencing, and finally a workspace field, which is used both as
 *  workspace for the parallel FFT, and as buffer for the communication
 *  algorithm used in the force computation.
 */
void pm_init_periodic_allocate(void)
{
  double bytes_tot = 0;
  size_t bytes;

  /* allocate the memory to hold the FFT fields */

  if(!(rhogrid = (fftw_real *) mymalloc(bytes = maxfftsize * sizeof(d_fftw_real))))
    {
      printf("failed to allocate memory for `FFT-rhogrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!(forcegrid = (fftw_real *) mymalloc(bytes = maxfftsize * sizeof(d_fftw_real))))
    {
      printf("failed to allocate memory for `FFT-forcegrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!(part = (struct part_slab_data *) mymalloc(bytes = 8 * NumPart * sizeof(struct part_slab_data))))
    {
      printf("failed to allocate memory for `part' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!(part_sortindex = (int *) mymalloc(bytes = 8 * NumPart * sizeof(int))))
    {
      printf("failed to allocate memory for `part_sortindex' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(ThisTask == 0)
    printf("Using %g MByte for periodic FFT computation. (presently allocated=%g MB)\n",
	   bytes_tot / (1024.0 * 1024.0), AllocatedBytes / (1024.0 * 1024.0));

  workspace = forcegrid;

  fft_of_rhogrid = (fftw_complex *) & rhogrid[0];

  d_rhogrid = (d_fftw_real *) rhogrid;
  d_forcegrid = (d_fftw_real *) forcegrid;
  d_workspace = (d_fftw_real *) workspace;
}



/*! This routine frees the space allocated for the parallel FFT algorithm.
 */
void pm_init_periodic_free(void)
{
  /* allocate the memory to hold the FFT fields */
  myfree(part_sortindex);
  myfree(part);
  myfree(forcegrid);
  myfree(rhogrid);
}



/*! Calculates the long-range periodic force given the particle positions
 *  using the PM method.  The force is Gaussian filtered with Asmth, given in
 *  mesh-cell units. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The potential is finite differenced
 *  using a 4-point finite differencing formula, and the forces are
 *  interpolated tri-linearly to the particle positions. The CIC kernel is
 *  deconvolved. Note that the particle distribution is not in the slab
 *  decomposition that is used for the FFT. Instead, overlapping patches
 *  between local domains and FFT slabs are communicated as needed.
 */
void pmforce_periodic(void)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, acc_dim;
  int i, j, slab, level, sendTask, recvTask, task;
  int x, y, z, yl, zl, yr, zr, yll, zll, yrr, zrr, ip, dim;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int num_on_grid, num_field_points, pindex, xx, yy, zz;
  MPI_Status status;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;

#ifdef SCALARFIELD
  int phase;
  double kscreening2;

  kscreening2 = pow(All.BoxSize / All.ScalarScreeningLength / (2 * M_PI), 2);
#endif


  force_treefree();

  if(ThisTask == 0)
    {
      printf("Starting periodic PM calculation.  (presently allocated=%g MB)\n",
	     AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */
  fac *= 1 / (2 * All.BoxSize / PMGRID);	/* for finite differencing */


  pm_init_periodic_allocate();

#ifdef SCALARFIELD
  for(phase = 0; phase < 2; phase++)
    {
#endif

      /* determine the cells each particles accesses */
      for(i = 0, num_on_grid = 0; i < NumPart; i++)
	{
#ifdef SCALARFIELD
	  if(phase == 1)
	    if(P[i].Type == 0)	/* don't bin baryonic mass in this phase */
	      continue;
#endif
	  slab_x = (int) (to_slab_fac * P[i].Pos[0]);
	  slab_y = (int) (to_slab_fac * P[i].Pos[1]);
	  slab_z = (int) (to_slab_fac * P[i].Pos[2]);

	  if(slab_x >= PMGRID)
	    slab_x = PMGRID - 1;
	  if(slab_y >= PMGRID)
	    slab_y = PMGRID - 1;
	  if(slab_z >= PMGRID)
	    slab_z = PMGRID - 1;

	  for(xx = 0; xx < 2; xx++)
	    for(yy = 0; yy < 2; yy++)
	      for(zz = 0; zz < 2; zz++)
		{
		  slab_xx = slab_x + xx;
		  slab_yy = slab_y + yy;
		  slab_zz = slab_z + zz;

		  if(slab_xx >= PMGRID)
		    slab_xx -= PMGRID;
		  if(slab_yy >= PMGRID)
		    slab_yy -= PMGRID;
		  if(slab_zz >= PMGRID)
		    slab_zz -= PMGRID;

		  offset = ((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_zz;

		  part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
		  part[num_on_grid].globalindex = offset;
		  part_sortindex[num_on_grid] = num_on_grid;
		  num_on_grid++;
		}
	}

      /* note: num_on_grid will be  8 times larger than the particle number,  
         but num_field_points will generally be much smaller */

      /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
      mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#else
      qsort(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#endif

      /* determine the number of unique field points */
      for(i = 0, num_field_points = 0; i < num_on_grid; i++)
	{
	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	      continue;

	  num_field_points++;
	}

      /* allocate the local field */
      localfield_globalindex = (large_array_offset *) mymalloc(num_field_points * sizeof(large_array_offset));
      localfield_d_data = (d_fftw_real *) mymalloc(num_field_points * sizeof(d_fftw_real));
      localfield_data = (fftw_real *) localfield_d_data;
      localfield_first = (int *) mymalloc(NTask * sizeof(int));
      localfield_count = (int *) mymalloc(NTask * sizeof(int));
      localfield_offset = (int *) mymalloc(NTask * sizeof(int));
      localfield_togo = (int *) mymalloc(NTask * NTask * sizeof(int));

      for(i = 0; i < NTask; i++)
	{
	  localfield_first[i] = 0;
	  localfield_count[i] = 0;
	}

      /* establish the cross link between the part[] array and the local list of 
         mesh points. Also, count on which CPU how many of the needed field points are stored */
      for(i = 0, num_field_points = 0; i < num_on_grid; i++)
	{
	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
	      num_field_points++;

	  part[part_sortindex[i]].localindex = num_field_points;

	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	      continue;

	  localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

	  slab = part[part_sortindex[i]].globalindex / (PMGRID * PMGRID2);
	  task = slab_to_task[slab];
	  if(localfield_count[task] == 0)
	    localfield_first[task] = num_field_points;
	  localfield_count[task]++;
	}
      num_field_points++;

      for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
	localfield_offset[i] = localfield_offset[i - 1] + localfield_count[i - 1];

      /* now bin the local particle data onto the mesh list */

      for(i = 0; i < num_field_points; i++)
	localfield_d_data[i] = 0;

      for(i = 0; i < num_on_grid; i += 8)
	{
	  pindex = (part[i].partindex >> 3);

	  slab_x = (int) (to_slab_fac * P[pindex].Pos[0]);
	  slab_y = (int) (to_slab_fac * P[pindex].Pos[1]);
	  slab_z = (int) (to_slab_fac * P[pindex].Pos[2]);

	  dx = to_slab_fac * P[pindex].Pos[0] - slab_x;
	  dy = to_slab_fac * P[pindex].Pos[1] - slab_y;
	  dz = to_slab_fac * P[pindex].Pos[2] - slab_z;

	  localfield_d_data[part[i + 0].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	  localfield_d_data[part[i + 1].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * dz;
	  localfield_d_data[part[i + 2].localindex] += P[pindex].Mass * (1.0 - dx) * dy * (1.0 - dz);
	  localfield_d_data[part[i + 3].localindex] += P[pindex].Mass * (1.0 - dx) * dy * dz;
	  localfield_d_data[part[i + 4].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
	  localfield_d_data[part[i + 5].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * dz;
	  localfield_d_data[part[i + 6].localindex] += P[pindex].Mass * (dx) * dy * (1.0 - dz);
	  localfield_d_data[part[i + 7].localindex] += P[pindex].Mass * (dx) * dy * dz;
	}

      /* clear local FFT-mesh density field */
      for(i = 0; i < fftsize; i++)
	d_rhogrid[i] = 0;

      /* exchange data and add contributions to the local mesh-path */

      MPI_Allgather(localfield_count, NTask, MPI_INT, localfield_togo, NTask, MPI_INT, MPI_COMM_WORLD);

      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      if(level > 0)
		{
		  import_d_data =
		    (d_fftw_real *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] *
					     sizeof(d_fftw_real));
		  import_globalindex =
		    (large_array_offset *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] *
						    sizeof(large_array_offset));

		  if(localfield_togo[sendTask * NTask + recvTask] > 0
		     || localfield_togo[recvTask * NTask + sendTask] > 0)
		    {
		      MPI_Sendrecv(localfield_d_data + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(d_fftw_real),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_A, import_d_data,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

		      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_B, import_globalindex,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_B, MPI_COMM_WORLD, &status);
		    }
		}
	      else
		{
		  import_d_data = localfield_d_data + localfield_offset[ThisTask];
		  import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
		}

	      for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
		{
		  /* determine offset in local FFT slab */
		  offset =
		    import_globalindex[i] -
		    first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);

		  d_rhogrid[offset] += import_d_data[i];
		}

	      if(level > 0)
		{
		  myfree(import_globalindex);
		  myfree(import_d_data);
		}
	    }
	}

#ifdef FLTROUNDOFFREDUCTION
      for(i = 0; i < fftsize; i++)	/* clear local density field */
	rhogrid[i] = FLT(d_rhogrid[i]);
#endif

      /* Do the FFT of the density field */

      rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

      /* multiply with Green's function for the potential */

      for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
	for(x = 0; x < PMGRID; x++)
	  for(z = 0; z < PMGRID / 2 + 1; z++)
	    {
	      if(x > PMGRID / 2)
		kx = x - PMGRID;
	      else
		kx = x;
	      if(y > PMGRID / 2)
		ky = y - PMGRID;
	      else
		ky = y;
	      if(z > PMGRID / 2)
		kz = z - PMGRID;
	      else
		kz = z;

	      k2 = kx * kx + ky * ky + kz * kz;

	      if(k2 > 0)
		{
#ifdef SCALARFIELD
		  if(phase == 1)
		    smth = -All.ScalarBeta * exp(-k2 * asmth2) / (k2 + kscreening2);
		  else
#endif
		    smth = -exp(-k2 * asmth2) / k2;

		  /* do deconvolution */

		  fx = fy = fz = 1;
		  if(kx != 0)
		    {
		      fx = (M_PI * kx) / PMGRID;
		      fx = sin(fx) / fx;
		    }
		  if(ky != 0)
		    {
		      fy = (M_PI * ky) / PMGRID;
		      fy = sin(fy) / fy;
		    }
		  if(kz != 0)
		    {
		      fz = (M_PI * kz) / PMGRID;
		      fz = sin(fz) / fz;
		    }
		  ff = 1 / (fx * fy * fz);
		  smth *= ff * ff * ff * ff;

		  /* end deconvolution */

		  ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
		  fft_of_rhogrid[ip].re *= smth;
		  fft_of_rhogrid[ip].im *= smth;
		}
	    }

      if(slabstart_y == 0)
	fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;

      /* Do the inverse FFT to get the potential */

      rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

      /* Now rhogrid holds the potential */


#ifdef EVALPOTENTIAL		/* now read out the potential */

      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      if(level > 0)
		{
		  import_data =
		    (fftw_real *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] * sizeof(fftw_real));
		  import_globalindex =
		    (large_array_offset *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] *
						    sizeof(large_array_offset));

		  if(localfield_togo[sendTask * NTask + recvTask] > 0
		     || localfield_togo[recvTask * NTask + sendTask] > 0)
		    {
		      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_C, import_globalindex,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
		    }
		}
	      else
		{
		  import_data = localfield_data + localfield_offset[ThisTask];
		  import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
		}

	      for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
		{
		  offset =
		    import_globalindex[i] -
		    first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
		  import_data[i] = rhogrid[offset];
		}

	      if(level > 0)
		{
		  MPI_Sendrecv(import_data,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
			       recvTask, TAG_NONPERIOD_A,
			       localfield_data + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
			       recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

		  myfree(import_globalindex);
		  myfree(import_data);
		}
	    }
	}

      /* read out the potential values, which all have been assembled in localfield_data */

      double pot;

      for(i = 0, j = 0; i < NumPart; i++)
	{
	  while(j < num_on_grid && (part[j].partindex >> 3) != i)
	    j++;

	  slab_x = (int) (to_slab_fac * P[i].Pos[0]);
	  dx = to_slab_fac * P[i].Pos[0] - slab_x;

	  slab_y = (int) (to_slab_fac * P[i].Pos[1]);
	  dy = to_slab_fac * P[i].Pos[1] - slab_y;

	  slab_z = (int) (to_slab_fac * P[i].Pos[2]);
	  dz = to_slab_fac * P[i].Pos[2] - slab_z;

	  pot =
	    +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
	    + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
	    + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
	    + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
	    + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
	    + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
	    + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
	    + localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

	  P[i].PM_Potential += pot * fac * (2 * All.BoxSize / PMGRID);
	  /* compensate the finite differencing factor */ ;
	}
#endif


      /* get the force components by finite differencing the potential for each dimension, 
         and send back the results to the right CPUs */

      for(dim = 2; dim >= 0; dim--)	/* Calculate each component of the force. */
	{			/* we do the x component last, because for differencing the potential in the x-direction, we need to contruct the transpose */
	  if(dim == 0)
	    pm_periodic_transposeA(rhogrid, forcegrid);	/* compute the transpose of the potential field */

	  for(xx = slabstart_x; xx < (slabstart_x + nslab_x); xx++)
	    for(y = 0; y < PMGRID; y++)
	      for(z = 0; z < PMGRID; z++)
		{
		  x = xx - slabstart_x;

		  yrr = yll = yr = yl = y;
		  zrr = zll = zr = zl = z;

		  switch (dim)
		    {
		    case 0:	/* note: for the x-direction, we difference the transposed direction (y) */
		    case 1:
		      yr = y + 1;
		      yl = y - 1;
		      yrr = y + 2;
		      yll = y - 2;
		      if(yr >= PMGRID)
			yr -= PMGRID;
		      if(yrr >= PMGRID)
			yrr -= PMGRID;
		      if(yl < 0)
			yl += PMGRID;
		      if(yll < 0)
			yll += PMGRID;
		      break;
		    case 2:
		      zr = z + 1;
		      zl = z - 1;
		      zrr = z + 2;
		      zll = z - 2;
		      if(zr >= PMGRID)
			zr -= PMGRID;
		      if(zrr >= PMGRID)
			zrr -= PMGRID;
		      if(zl < 0)
			zl += PMGRID;
		      if(zll < 0)
			zll += PMGRID;
		      break;
		    }

		  if(dim == 0)
		    {
		      forcegrid[PMGRID * (x + y * nslab_x) + z]
			=
			fac * ((4.0 / 3) *
			       (rhogrid[PMGRID * (x + yl * nslab_x) + zl] -
				rhogrid[PMGRID * (x + yr * nslab_x) + zr]) -
			       (1.0 / 6) * (rhogrid[PMGRID * (x + yll * nslab_x) + zll] -
					    rhogrid[PMGRID * (x + yrr * nslab_x) + zrr]));
		    }
		  else
		    forcegrid[PMGRID2 * (PMGRID * x + y) + z]
		      =
		      fac * ((4.0 / 3) *
			     (rhogrid[PMGRID2 * (PMGRID * x + yl) + zl] -
			      rhogrid[PMGRID2 * (PMGRID * x + yr) + zr]) -
			     (1.0 / 6) * (rhogrid[PMGRID2 * (PMGRID * x + yll) + zll] -
					  rhogrid[PMGRID2 * (PMGRID * x + yrr) + zrr]));
		}

	  if(dim == 0)
	    pm_periodic_transposeB(forcegrid, rhogrid);	/* compute the transpose of the potential field */

	  /* send the force components to the right processors */

	  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ level;

	      if(recvTask < NTask)
		{
		  if(level > 0)
		    {
		      import_data =
			(fftw_real *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] *
					       sizeof(fftw_real));
		      import_globalindex =
			(large_array_offset *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] *
							sizeof(large_array_offset));

		      if(localfield_togo[sendTask * NTask + recvTask] > 0
			 || localfield_togo[recvTask * NTask + sendTask] > 0)
			{
			  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				       localfield_togo[sendTask * NTask +
						       recvTask] * sizeof(large_array_offset), MPI_BYTE,
				       recvTask, TAG_NONPERIOD_C, import_globalindex,
				       localfield_togo[recvTask * NTask +
						       sendTask] * sizeof(large_array_offset), MPI_BYTE,
				       recvTask, TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
			}
		    }
		  else
		    {
		      import_data = localfield_data + localfield_offset[ThisTask];
		      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
		    }

		  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
		    {
		      /* determine offset in local FFT slab */
		      offset =
			import_globalindex[i] -
			first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
		      import_data[i] = forcegrid[offset];
		    }

		  if(level > 0)
		    {
		      MPI_Sendrecv(import_data,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_NONPERIOD_A,
				   localfield_data + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

		      myfree(import_globalindex);
		      myfree(import_data);
		    }
		}
	    }

	  /* read out the forces, which all have been assembled in localfield_data */

	  for(i = 0, j = 0; i < NumPart; i++)
	    {
#ifdef SCALARFIELD
	      if(phase == 1)
		if(P[i].Type == 0)	/* baryons don't get an extra scalar force */
		  continue;
#endif
	      while(j < num_on_grid && (part[j].partindex >> 3) != i)
		j++;

	      slab_x = (int) (to_slab_fac * P[i].Pos[0]);
	      dx = to_slab_fac * P[i].Pos[0] - slab_x;

	      slab_y = (int) (to_slab_fac * P[i].Pos[1]);
	      dy = to_slab_fac * P[i].Pos[1] - slab_y;

	      slab_z = (int) (to_slab_fac * P[i].Pos[2]);
	      dz = to_slab_fac * P[i].Pos[2] - slab_z;

	      acc_dim =
		+localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
		+ localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
		+ localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
		+ localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
		+ localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
		+ localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
		+ localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
		+ localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

	      P[i].GravPM[dim] += acc_dim;
	    }
	}

      /* free locallist */
      myfree(localfield_togo);
      myfree(localfield_offset);
      myfree(localfield_count);
      myfree(localfield_first);
      myfree(localfield_d_data);
      myfree(localfield_globalindex);
#ifdef SCALARFIELD
    }
#endif

  pm_init_periodic_free();
  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);

  if(ThisTask == 0)
    {
      printf("done PM.\n");
      fflush(stdout);
    }
}


/*! Calculates the long-range potential using the PM method.  The potential is
 *  Gaussian filtered with Asmth, given in mesh-cell units. We carry out a CIC
 *  charge assignment, and compute the potenial by Fourier transform
 *  methods. The CIC kernel is deconvolved.
 */
void pmpotential_periodic(void)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, pot;
  int i, j, slab, level, sendTask, recvTask, task;
  int x, y, z, ip;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int num_on_grid, num_field_points, pindex, xx, yy, zz;
  MPI_Status status;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;

  force_treefree();

  if(ThisTask == 0)
    {
      printf("Starting periodic PM-potential calculation.  (presently allocated=%g MB)\n",
	     AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */

  pm_init_periodic_allocate();


  /* determine the cells each particles accesses */
  for(i = 0, num_on_grid = 0; i < NumPart; i++)
    {
      slab_x = (int) (to_slab_fac * P[i].Pos[0]);
      slab_y = (int) (to_slab_fac * P[i].Pos[1]);
      slab_z = (int) (to_slab_fac * P[i].Pos[2]);

      if(slab_x >= PMGRID)
	slab_x = PMGRID - 1;
      if(slab_y >= PMGRID)
	slab_y = PMGRID - 1;
      if(slab_z >= PMGRID)
	slab_z = PMGRID - 1;

      for(xx = 0; xx < 2; xx++)
	for(yy = 0; yy < 2; yy++)
	  for(zz = 0; zz < 2; zz++)
	    {
	      slab_xx = slab_x + xx;
	      slab_yy = slab_y + yy;
	      slab_zz = slab_z + zz;

	      if(slab_xx >= PMGRID)
		slab_xx -= PMGRID;
	      if(slab_yy >= PMGRID)
		slab_yy -= PMGRID;
	      if(slab_zz >= PMGRID)
		slab_zz -= PMGRID;

	      offset = ((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_zz;

	      part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
	      part[num_on_grid].globalindex = offset;
	      part_sortindex[num_on_grid] = num_on_grid;
	      num_on_grid++;
	    }
    }

  /* note: num_on_grid will be  8 times larger than the particle number,  
     but num_field_points will generally be much smaller */

  /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
  mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#else
  qsort(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#endif

  /* determine the number of unique field points */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
	if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	  continue;

      num_field_points++;
    }

  /* allocate the local field */
  localfield_globalindex = (large_array_offset *) mymalloc(num_field_points * sizeof(large_array_offset));
  localfield_d_data = (d_fftw_real *) mymalloc(num_field_points * sizeof(d_fftw_real));
  localfield_data = (fftw_real *) localfield_d_data;
  localfield_first = (int *) mymalloc(NTask * sizeof(int));
  localfield_count = (int *) mymalloc(NTask * sizeof(int));
  localfield_offset = (int *) mymalloc(NTask * sizeof(int));
  localfield_togo = (int *) mymalloc(NTask * NTask * sizeof(int));

  for(i = 0; i < NTask; i++)
    {
      localfield_first[i] = 0;
      localfield_count[i] = 0;
    }

  /* establish the cross link between the part[] array and the local list of 
     mesh points. Also, count on which CPU how many of the needed field points are stored */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
	if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
	  num_field_points++;

      part[part_sortindex[i]].localindex = num_field_points;

      if(i > 0)
	if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	  continue;

      localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

      slab = part[part_sortindex[i]].globalindex / (PMGRID * PMGRID2);
      task = slab_to_task[slab];
      if(localfield_count[task] == 0)
	localfield_first[task] = num_field_points;
      localfield_count[task]++;
    }
  num_field_points++;

  for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
    localfield_offset[i] = localfield_offset[i - 1] + localfield_count[i - 1];

  /* now bin the local particle data onto the mesh list */

  for(i = 0; i < num_field_points; i++)
    localfield_d_data[i] = 0;

  for(i = 0; i < num_on_grid; i += 8)
    {
      pindex = (part[i].partindex >> 3);

      slab_x = (int) (to_slab_fac * P[pindex].Pos[0]);
      slab_y = (int) (to_slab_fac * P[pindex].Pos[1]);
      slab_z = (int) (to_slab_fac * P[pindex].Pos[2]);

      dx = to_slab_fac * P[pindex].Pos[0] - slab_x;
      dy = to_slab_fac * P[pindex].Pos[1] - slab_y;
      dz = to_slab_fac * P[pindex].Pos[2] - slab_z;

      localfield_d_data[part[i + 0].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      localfield_d_data[part[i + 1].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * dz;
      localfield_d_data[part[i + 2].localindex] += P[pindex].Mass * (1.0 - dx) * dy * (1.0 - dz);
      localfield_d_data[part[i + 3].localindex] += P[pindex].Mass * (1.0 - dx) * dy * dz;
      localfield_d_data[part[i + 4].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
      localfield_d_data[part[i + 5].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * dz;
      localfield_d_data[part[i + 6].localindex] += P[pindex].Mass * (dx) * dy * (1.0 - dz);
      localfield_d_data[part[i + 7].localindex] += P[pindex].Mass * (dx) * dy * dz;
    }

  /* clear local FFT-mesh density field */
  for(i = 0; i < fftsize; i++)
    d_rhogrid[i] = 0;

  /* exchange data and add contributions to the local mesh-path */

  MPI_Allgather(localfield_count, NTask, MPI_INT, localfield_togo, NTask, MPI_INT, MPI_COMM_WORLD);

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{
	  if(level > 0)
	    {
	      import_d_data =
		(d_fftw_real *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] * sizeof(d_fftw_real));
	      import_globalindex =
		(large_array_offset *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] *
						sizeof(large_array_offset));

	      if(localfield_togo[sendTask * NTask + recvTask] > 0
		 || localfield_togo[recvTask * NTask + sendTask] > 0)
		{
		  MPI_Sendrecv(localfield_d_data + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_NONPERIOD_A,
			       import_d_data,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_B, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_B, MPI_COMM_WORLD, &status);
		}
	    }
	  else
	    {
	      import_d_data = localfield_d_data + localfield_offset[ThisTask];
	      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
	    }

	  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
	    {
	      /* determine offset in local FFT slab */
	      offset =
		import_globalindex[i] -
		first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);

	      d_rhogrid[offset] += import_d_data[i];
	    }

	  if(level > 0)
	    {
	      myfree(import_globalindex);
	      myfree(import_d_data);
	    }
	}
    }

#ifdef FLTROUNDOFFREDUCTION
  for(i = 0; i < fftsize; i++)	/* clear local density field */
    rhogrid[i] = FLT(d_rhogrid[i]);
#endif

  /* Do the FFT of the density field */

  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* multiply with Green's function for the potential */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
	{
	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      smth = -exp(-k2 * asmth2) / k2 * fac;

	      /* do deconvolution */

	      fx = fy = fz = 1;
	      if(kx != 0)
		{
		  fx = (M_PI * kx) / PMGRID;
		  fx = sin(fx) / fx;
		}
	      if(ky != 0)
		{
		  fy = (M_PI * ky) / PMGRID;
		  fy = sin(fy) / fy;
		}
	      if(kz != 0)
		{
		  fz = (M_PI * kz) / PMGRID;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth *= ff * ff * ff * ff;

	      /* end deconvolution */

	      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
	      fft_of_rhogrid[ip].re *= smth;
	      fft_of_rhogrid[ip].im *= smth;
	    }
	}

  if(slabstart_y == 0)
    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;

  /* Do the inverse FFT to get the potential */

  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* Now rhogrid holds the potential */


  /* now read out the potential */

  /* send the force components to the right processors */

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{
	  if(level > 0)
	    {
	      import_data =
		(fftw_real *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] * sizeof(fftw_real));
	      import_globalindex =
		(large_array_offset *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] *
						sizeof(large_array_offset));

	      if(localfield_togo[sendTask * NTask + recvTask] > 0
		 || localfield_togo[recvTask * NTask + sendTask] > 0)
		{
		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_C, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
		}
	    }
	  else
	    {
	      import_data = localfield_data + localfield_offset[ThisTask];
	      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
	    }

	  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
	    {
	      /* determine offset in local FFT slab */
	      offset =
		import_globalindex[i] -
		first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
	      import_data[i] = rhogrid[offset];
	    }

	  if(level > 0)
	    {
	      MPI_Sendrecv(import_data,
			   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_NONPERIOD_A,
			   localfield_data + localfield_offset[recvTask],
			   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

	      myfree(import_globalindex);
	      myfree(import_data);
	    }
	}
    }

  /* read out the potential values, which all have been assembled in localfield_data */

  for(i = 0, j = 0; i < NumPart; i++)
    {
      while(j < num_on_grid && (part[j].partindex >> 3) != i)
	j++;

      slab_x = (int) (to_slab_fac * P[i].Pos[0]);
      dx = to_slab_fac * P[i].Pos[0] - slab_x;

      slab_y = (int) (to_slab_fac * P[i].Pos[1]);
      dy = to_slab_fac * P[i].Pos[1] - slab_y;

      slab_z = (int) (to_slab_fac * P[i].Pos[2]);
      dz = to_slab_fac * P[i].Pos[2] - slab_z;

      pot =
	+localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
	+ localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
	+ localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
	+ localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
	+ localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
	+ localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
	+ localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
	+ localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
      P[i].p.Potential += pot;
#endif
    }

  /* free locallist */
  myfree(localfield_togo);
  myfree(localfield_offset);
  myfree(localfield_count);
  myfree(localfield_first);
  myfree(localfield_d_data);
  myfree(localfield_globalindex);

  pm_init_periodic_free();
  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);

  if(ThisTask == 0)
    {
      printf("done PM-Potential.\n");
      fflush(stdout);
    }
}



int pm_periodic_compare_sortindex(const void *a, const void *b)
{
  if(part[*(int *) a].globalindex < part[*(int *) b].globalindex)
    return -1;

  if(part[*(int *) a].globalindex > part[*(int *) b].globalindex)
    return +1;

  return 0;
}

static void msort_pmperiodic_with_tmp(int *b, size_t n, int *t)
{
  int *tmp;
  int *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_pmperiodic_with_tmp(b1, n1, t);
  msort_pmperiodic_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(part[*b1].globalindex <= part[*b2].globalindex)
	{
	  --n1;
	  *tmp++ = *b1++;
	}
      else
	{
	  --n2;
	  *tmp++ = *b2++;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(int));

  memcpy(b, t, (n - n2) * sizeof(int));
}

void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *))
{
  const size_t size = n * s;

  int *tmp = (int *) mymalloc(size);

  msort_pmperiodic_with_tmp((int *) b, n, tmp);

  myfree(tmp);
}

void pm_periodic_transposeA(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	for(z = 0; z < PMGRID; z++)
	  {
	    scratch[PMGRID * (first_slab_of_task[task] * nslab_x +
			      x * slabs_per_task[task] + (y - first_slab_of_task[task])) + z] =
	      field[PMGRID2 * (PMGRID * x + y) + z];
	  }

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc(2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(field + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }

  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);
#else
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
	{
	  MPI_Sendrecv(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       field + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif
}



void pm_periodic_transposeB(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc(2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(field + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }


  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);

#else
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
	{
	  MPI_Sendrecv(field + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	for(z = 0; z < PMGRID; z++)
	  {
	    field[PMGRID2 * (PMGRID * x + y) + z] =
	      scratch[PMGRID * (first_slab_of_task[task] * nslab_x +
				x * slabs_per_task[task] + (y - first_slab_of_task[task])) + z];
	  }

}


#endif
#endif







#ifdef PMGRID
#ifdef PERIODIC



#ifdef DISTORTIONTENSORPS
/*! Calculates the long-range tidal field using the PM method.  The potential is
 *  Gaussian filtered with Asmth, given in mesh-cell units. We carry out a CIC
 *  charge assignment, and compute the potenial by Fourier transform
 *  methods. The CIC kernel is deconvolved.
 *  Note that the k's need a pre-factor of 2 M_PI / All.BoxSize.
 *  The procedure calculates the second derivates of the gravitational potential by "pulling" down k's in fourier space.
 *  Component specifies the entry in the tidal field tensor that should be calculated:
 *  0=xx 1=xy 2=xz 3=yy 4=yz 5=zz
 */
/* NOTES:
   this procedure is not yet  optimized:
   ->force_treefree and allocate are called again and again
   ->the FFT is done in both directions for all components
 */
void pmtidaltensor_periodic(int component)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, tidal;
  int i, j, slab, level, sendTask, recvTask, task;
  int x, y, z, ip;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int num_on_grid, num_field_points, pindex, xx, yy, zz;
  MPI_Status status;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;

  force_treefree();

  if(ThisTask == 0)
    {
      printf("Starting periodic PM-Tidaltensor (component=%d) calculation.  (presently allocated=%g MB)\n",
	     component, AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */

  pm_init_periodic_allocate();


  /* determine the cells each particles accesses */
  for(i = 0, num_on_grid = 0; i < NumPart; i++)
    {
      slab_x = (int) (to_slab_fac * P[i].Pos[0]);
      slab_y = (int) (to_slab_fac * P[i].Pos[1]);
      slab_z = (int) (to_slab_fac * P[i].Pos[2]);

      if(slab_x >= PMGRID)
	slab_x = PMGRID - 1;
      if(slab_y >= PMGRID)
	slab_y = PMGRID - 1;
      if(slab_z >= PMGRID)
	slab_z = PMGRID - 1;

      for(xx = 0; xx < 2; xx++)
	for(yy = 0; yy < 2; yy++)
	  for(zz = 0; zz < 2; zz++)
	    {
	      slab_xx = slab_x + xx;
	      slab_yy = slab_y + yy;
	      slab_zz = slab_z + zz;

	      if(slab_xx >= PMGRID)
		slab_xx -= PMGRID;
	      if(slab_yy >= PMGRID)
		slab_yy -= PMGRID;
	      if(slab_zz >= PMGRID)
		slab_zz -= PMGRID;

	      offset = ((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_zz;

	      part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
	      part[num_on_grid].globalindex = offset;
	      part_sortindex[num_on_grid] = num_on_grid;
	      num_on_grid++;
	    }
    }

  /* note: num_on_grid will be  8 times larger than the particle number,  
     but num_field_points will generally be much smaller */

  /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
  mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#else
  qsort(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#endif

  /* determine the number of unique field points */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
	if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	  continue;

      num_field_points++;
    }

  /* allocate the local field */
  localfield_globalindex = (large_array_offset *) mymalloc(num_field_points * sizeof(large_array_offset));
  localfield_d_data = (d_fftw_real *) mymalloc(num_field_points * sizeof(d_fftw_real));
  localfield_data = (fftw_real *) localfield_d_data;
  localfield_first = (int *) mymalloc(NTask * sizeof(int));
  localfield_count = (int *) mymalloc(NTask * sizeof(int));
  localfield_offset = (int *) mymalloc(NTask * sizeof(int));
  localfield_togo = (int *) mymalloc(NTask * NTask * sizeof(int));

  for(i = 0; i < NTask; i++)
    {
      localfield_first[i] = 0;
      localfield_count[i] = 0;
    }

  /* establish the cross link between the part[] array and the local list of 
     mesh points. Also, count on which CPU how many of the needed field points are stored */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
	if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
	  num_field_points++;

      part[part_sortindex[i]].localindex = num_field_points;

      if(i > 0)
	if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	  continue;

      localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

      slab = part[part_sortindex[i]].globalindex / (PMGRID * PMGRID2);
      task = slab_to_task[slab];
      if(localfield_count[task] == 0)
	localfield_first[task] = num_field_points;
      localfield_count[task]++;
    }
  num_field_points++;

  for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
    localfield_offset[i] = localfield_offset[i - 1] + localfield_count[i - 1];

  /* now bin the local particle data onto the mesh list */

  for(i = 0; i < num_field_points; i++)
    localfield_d_data[i] = 0;

  for(i = 0; i < num_on_grid; i += 8)
    {
      pindex = (part[i].partindex >> 3);

      slab_x = (int) (to_slab_fac * P[pindex].Pos[0]);
      slab_y = (int) (to_slab_fac * P[pindex].Pos[1]);
      slab_z = (int) (to_slab_fac * P[pindex].Pos[2]);

      dx = to_slab_fac * P[pindex].Pos[0] - slab_x;
      dy = to_slab_fac * P[pindex].Pos[1] - slab_y;
      dz = to_slab_fac * P[pindex].Pos[2] - slab_z;

      localfield_d_data[part[i + 0].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      localfield_d_data[part[i + 1].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * dz;
      localfield_d_data[part[i + 2].localindex] += P[pindex].Mass * (1.0 - dx) * dy * (1.0 - dz);
      localfield_d_data[part[i + 3].localindex] += P[pindex].Mass * (1.0 - dx) * dy * dz;
      localfield_d_data[part[i + 4].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
      localfield_d_data[part[i + 5].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * dz;
      localfield_d_data[part[i + 6].localindex] += P[pindex].Mass * (dx) * dy * (1.0 - dz);
      localfield_d_data[part[i + 7].localindex] += P[pindex].Mass * (dx) * dy * dz;
    }

  /* clear local FFT-mesh density field */
  for(i = 0; i < fftsize; i++)
    d_rhogrid[i] = 0;

  /* exchange data and add contributions to the local mesh-path */

  MPI_Allgather(localfield_count, NTask, MPI_INT, localfield_togo, NTask, MPI_INT, MPI_COMM_WORLD);

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{
	  if(level > 0)
	    {
	      import_d_data =
		(d_fftw_real *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] * sizeof(d_fftw_real));
	      import_globalindex =
		(large_array_offset *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] *
						sizeof(large_array_offset));

	      if(localfield_togo[sendTask * NTask + recvTask] > 0
		 || localfield_togo[recvTask * NTask + sendTask] > 0)
		{
		  MPI_Sendrecv(localfield_d_data + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_NONPERIOD_A,
			       import_d_data,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_B, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_B, MPI_COMM_WORLD, &status);
		}
	    }
	  else
	    {
	      import_d_data = localfield_d_data + localfield_offset[ThisTask];
	      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
	    }

	  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
	    {
	      /* determine offset in local FFT slab */
	      offset =
		import_globalindex[i] -
		first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);

	      d_rhogrid[offset] += import_d_data[i];
	    }

	  if(level > 0)
	    {
	      myfree(import_globalindex);
	      myfree(import_d_data);
	    }
	}
    }

#ifdef FLTROUNDOFFREDUCTION
  for(i = 0; i < fftsize; i++)	/* clear local density field */
    rhogrid[i] = FLT(d_rhogrid[i]);
#endif

  /* Do the FFT of the density field */

  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* multiply with Green's function for the potential */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
	{
	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      smth = -exp(-k2 * asmth2) / k2;

	      /* do deconvolution */

	      fx = fy = fz = 1;
	      if(kx != 0)
		{
		  fx = (M_PI * kx) / PMGRID;
		  fx = sin(fx) / fx;
		}
	      if(ky != 0)
		{
		  fy = (M_PI * ky) / PMGRID;
		  fy = sin(fy) / fy;
		}
	      if(kz != 0)
		{
		  fz = (M_PI * kz) / PMGRID;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth *= ff * ff * ff * ff;


	      /* end deconvolution */

	      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;

	      /* modify greens function to get second derivatives of potential ("pulling" down k's) */
	      if(component == 0)
		{
		  fft_of_rhogrid[ip].re *= smth * kx * kx;
		  fft_of_rhogrid[ip].im *= smth * kx * kx;
		}
	      if(component == 1)
		{
		  fft_of_rhogrid[ip].re *= smth * kx * ky;
		  fft_of_rhogrid[ip].im *= smth * kx * ky;
		}
	      if(component == 2)
		{
		  fft_of_rhogrid[ip].re *= smth * kx * kz;
		  fft_of_rhogrid[ip].im *= smth * kx * kz;
		}
	      if(component == 3)
		{
		  fft_of_rhogrid[ip].re *= smth * ky * ky;
		  fft_of_rhogrid[ip].im *= smth * ky * ky;
		}
	      if(component == 4)
		{
		  fft_of_rhogrid[ip].re *= smth * ky * kz;
		  fft_of_rhogrid[ip].im *= smth * ky * kz;
		}
	      if(component == 5)
		{
		  /*
                    FORCE TEST:
                    this calculates F_z by pulling down -i k_z, later on this is compared to the trilinear interpolation
		    i k_z comes from the FFTW backward transformation and -1 because the force is given by the negative gradient 
                    
                    the second derivative that is needed for the tidalfield can be calculated in the same way by pulling down.
                  */

		  /*
		     double rep, imp;
		     rep = fft_of_rhogrid[ip].re;
		     imp = fft_of_rhogrid[ip].im;               

		     fft_of_rhogrid[ip].re = smth*imp*kz * (2*M_PI) / All.BoxSize;
		     fft_of_rhogrid[ip].im = -smth*rep*kz * (2*M_PI) / All.BoxSize;
		   */

		  fft_of_rhogrid[ip].re *= smth * kz * kz;
		  fft_of_rhogrid[ip].im *= smth * kz * kz;

		}

	      /* prefactor = (2*M_PI) / All.BoxSize */
	      /* note: tidal tensor = - d^2 Phi/ dx_i dx_j  IS THE SIGN CORRECT ?!?! */
	      fft_of_rhogrid[ip].re *= (2 * M_PI) * (2 * M_PI) / (All.BoxSize * All.BoxSize);
	      fft_of_rhogrid[ip].im *= (2 * M_PI) * (2 * M_PI) / (All.BoxSize * All.BoxSize);

	    }
	}

  if(slabstart_y == 0)
    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;

  /* Do the inverse FFT to get the tidal tensor component */

  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* Now rhogrid holds the tidal tensor componet */


  /* now read out the tidal tensor component */

  /* send the tidal tensor component to the right processors */

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{
	  if(level > 0)
	    {
	      import_data =
		(fftw_real *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] * sizeof(fftw_real));
	      import_globalindex =
		(large_array_offset *) mymalloc(localfield_togo[recvTask * NTask + ThisTask] *
						sizeof(large_array_offset));

	      if(localfield_togo[sendTask * NTask + recvTask] > 0
		 || localfield_togo[recvTask * NTask + sendTask] > 0)
		{
		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_C, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
		}
	    }
	  else
	    {
	      import_data = localfield_data + localfield_offset[ThisTask];
	      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
	    }

	  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
	    {
	      /* determine offset in local FFT slab */
	      offset =
		import_globalindex[i] -
		first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
	      import_data[i] = rhogrid[offset];
	    }

	  if(level > 0)
	    {
	      MPI_Sendrecv(import_data,
			   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_NONPERIOD_A,
			   localfield_data + localfield_offset[recvTask],
			   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

	      myfree(import_globalindex);
	      myfree(import_data);
	    }
	}
    }

  /* read out the tidal field values, which all have been assembled in localfield_data */

  for(i = 0, j = 0; i < NumPart; i++)
    {
      while(j < num_on_grid && (part[j].partindex >> 3) != i)
	j++;

      slab_x = (int) (to_slab_fac * P[i].Pos[0]);
      dx = to_slab_fac * P[i].Pos[0] - slab_x;

      slab_y = (int) (to_slab_fac * P[i].Pos[1]);
      dy = to_slab_fac * P[i].Pos[1] - slab_y;

      slab_z = (int) (to_slab_fac * P[i].Pos[2]);
      dz = to_slab_fac * P[i].Pos[2] - slab_z;

      tidal =
	+localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
	+ localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
	+ localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
	+ localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
	+ localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
	+ localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
	+ localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
	+ localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

      tidal *= fac;

#ifdef DISTORTIONTENSORPS
      if(component == 0)
	{
	  P[i].tidal_tensorpsPM[0][0] += tidal;
	}
      if(component == 1)
	{
	  P[i].tidal_tensorpsPM[0][1] += tidal;
	  P[i].tidal_tensorpsPM[1][0] += tidal;
	}
      if(component == 2)
	{
	  P[i].tidal_tensorpsPM[0][2] += tidal;
	  P[i].tidal_tensorpsPM[2][0] += tidal;
	}
      if(component == 3)
	{
	  P[i].tidal_tensorpsPM[1][1] += tidal;
	}
      if(component == 4)
	{
	  P[i].tidal_tensorpsPM[1][2] += tidal;
	  P[i].tidal_tensorpsPM[2][1] += tidal;
	}
      if(component == 5)
	{
	  P[i].tidal_tensorpsPM[2][2] += tidal;
	}
#endif

    }

  /* free locallist */
  myfree(localfield_togo);
  myfree(localfield_offset);
  myfree(localfield_count);
  myfree(localfield_first);
  myfree(localfield_d_data);
  myfree(localfield_globalindex);

  pm_init_periodic_free();
  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);

  if(ThisTask == 0)
    {
      printf("done PM-Tidaltensor (component=%d).\n", component);
      fflush(stdout);
    }
}

void check_tidaltensor(int particle_ID)
{
  int i;

  for(i = 0; i < NumPart; i++)
    {

      if(P[i].ID == particle_ID)
	{

	  FdTidaltensor = fopen("Tidaltensor.txt", "a");
#ifdef DISTORTIONTENSORPS
	  fprintf(FdTidaltensor, "Mesh-Force: %f %f %f\n", P[i].GravPM[0], P[i].GravPM[1], P[i].GravPM[2]);
	  fprintf(FdTidaltensor, "Tree-Force: %f %f %f\n", P[i].g.GravAccel[0], P[i].g.GravAccel[1],
		  P[i].g.GravAccel[2]);
	  fprintf(FdTidaltensor, "Sum-Force : %f %f %f\n", P[i].g.GravAccel[0] + P[i].GravPM[0],
		  P[i].g.GravAccel[1] + P[i].GravPM[1], P[i].g.GravAccel[2] + P[i].GravPM[2]);

	  fprintf(FdTidaltensor, "Mesh-Tidal: %f %f %f %f %f %f\n", P[i].tidal_tensorpsPM[0][0],
		  P[i].tidal_tensorpsPM[0][1], P[i].tidal_tensorpsPM[0][2], P[i].tidal_tensorpsPM[1][1],
		  P[i].tidal_tensorpsPM[1][2], P[i].tidal_tensorpsPM[2][2]);
	  fprintf(FdTidaltensor, "Tree-Tidal: %f %f %f %f %f %f\n", P[i].tidal_tensorps[0][0],
		  P[i].tidal_tensorps[0][1], P[i].tidal_tensorps[0][2], P[i].tidal_tensorps[1][1],
		  P[i].tidal_tensorps[1][2], P[i].tidal_tensorps[2][2]);
	  fprintf(FdTidaltensor, "Sum-Tidal: %f %f %f %f %f %f\n",
		  P[i].tidal_tensorpsPM[0][0] + P[i].tidal_tensorps[0][0],
		  P[i].tidal_tensorpsPM[0][1] + P[i].tidal_tensorps[0][1],
		  P[i].tidal_tensorpsPM[0][2] + P[i].tidal_tensorps[0][2],
		  P[i].tidal_tensorpsPM[1][1] + P[i].tidal_tensorps[1][1],
		  P[i].tidal_tensorpsPM[1][2] + P[i].tidal_tensorps[1][2],
		  P[i].tidal_tensorpsPM[2][2] + P[i].tidal_tensorps[2][2]);
	  fprintf(FdTidaltensor, "----------\n");

/*     FORCE TEST:
       fprintf(FdTidaltensor,"Test:\n");
       fprintf(FdTidaltensor,"FORCE: TRILINEAR = %f\n", P[i].GravPM[2]);
       fprintf(FdTidaltensor,"FORCE: FOURIER   = %f\n", P[i].tidal_tensorpsPM[2][2]);       
       fprintf(FdTidaltensor,"------------------------\n");
       
*/

#endif
	  fclose(FdTidaltensor);
	}
    }
}

#endif

#endif
#endif
