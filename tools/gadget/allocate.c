#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"







/* This routine allocates memory for 
 * particle storage, both the collisionless and the SPH particles.
 * The memory for the ordered binary tree of the timeline
 * is also allocated.
 */
void allocate_memory(void)
{
  size_t bytes;
  double bytes_tot = 0;
  int NTaskTimesThreads;

  NTaskTimesThreads = NTask;
#ifdef NUM_THREADS
  NTaskTimesThreads = NUM_THREADS * NTask;
#endif

  Exportflag = (int *) mymalloc(NTaskTimesThreads * sizeof(int));
  Exportindex = (int *) mymalloc(NTaskTimesThreads * sizeof(int));
  Exportnodecount = (int *) mymalloc(NTaskTimesThreads * sizeof(int));

  Send_count = (int *) mymalloc(sizeof(int) * NTask);
  Send_offset = (int *) mymalloc(sizeof(int) * NTask);
  Recv_count = (int *) mymalloc(sizeof(int) * NTask);
  Recv_offset = (int *) mymalloc(sizeof(int) * NTask);
  Sendcount_matrix = (int *) mymalloc(sizeof(int) * NTask * NTask);

  ProcessedFlag = (unsigned char *) mymalloc(bytes = All.MaxPart * sizeof(unsigned char));
  bytes_tot += bytes;

  NextActiveParticle = (int *) mymalloc(bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

  NextInTimeBin = (int *) mymalloc(bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

  PrevInTimeBin = (int *) mymalloc(bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;


  if(All.MaxPart > 0)
    {
      if(!(P = (struct particle_data *) mymalloc(bytes = All.MaxPart * sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory for `P' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("\nAllocated %g MByte for particle storage.\n\n", bytes_tot / (1024.0 * 1024.0));
    }

  if(All.MaxPart > 0)
    {
      bytes_tot = 0;

      if(!
#ifdef GASRETURN
	 (SphP =
	  (struct sph_particle_data *) mymalloc(bytes = All.MaxPart * sizeof(struct sph_particle_data))))
#else
	 (SphP =
	  (struct sph_particle_data *) mymalloc(bytes = All.MaxPartSph * sizeof(struct sph_particle_data))))
#endif
	{
	  printf("failed to allocate memory for `SphP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of SPH data.\n\n", bytes_tot / (1024.0 * 1024.0));

    }
}
