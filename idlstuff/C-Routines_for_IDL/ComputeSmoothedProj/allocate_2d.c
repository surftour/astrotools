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
void allocate_memory_2d(void)
{
  int bytes,bytes_tot=0;

  printf("MaxPart %d\n",All.MaxPart);
  if(All.MaxPart>0)
    {
      if(!(Pn_data=malloc(bytes=All.MaxPart*sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory for `Pn_data' (%d bytes).\n",bytes);
	  endrun(1);
	}
      bytes_tot+=bytes;


      if(!(PTimeTree=malloc(bytes=All.MaxPart*sizeof(struct timetree_data))))
	{
	  printf("failed to allocate memory for `PTimeTree' (%d bytes).\n",bytes);
	  endrun(1);
	}
      bytes_tot+=bytes;


      Pn = Pn_data-1;   /* start with offset 1 */
      PTimeTree--;

      printf("\nAllocated %g MByte for particle storage.\n\n",bytes_tot/(1024.0*1024.0));
    }

  if(All.MaxPartSph>0)
    {
      bytes_tot=0;

      if(!(SphPn_data=malloc(bytes=All.MaxPartSph*sizeof(struct  sph_particle_data))))
	{
	  printf("failed to allocate memory for `SphPn_data' (%d bytes).\n",bytes);
	  endrun(1);
	}
      bytes_tot+=bytes;

      SphPn= SphPn_data-1;   /* start with offset 1 */

      printf("Allocated %g MByte for storage of SPH data.\n\n",bytes_tot/(1024.0*1024.0));
    }
}



/* This routine frees the memory for the particle storage,
 * but we don't actually call it in the code. 
 * When the program terminats, the memory will be automatically
 * freed by the operating system.
 */
void free_memory_2d(void)
{
  if(All.MaxPart>0)
    free(Pn_data);
  
  if(All.MaxPartSph>0)
    free(SphPn_data);
}
