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
  int bytes,bytes_tot=0;

  //printf("MaxPart %d\n",All.MaxPart);
  if(All.MaxPart>0)
    {
      if(!(P_data=malloc(bytes=All.MaxPart*sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory for `P_data' (%d bytes).\n",bytes);
	  endrun(1);
	}
      bytes_tot+=bytes;


      if(!(PTimeTree=malloc(bytes=All.MaxPart*sizeof(struct timetree_data))))
	{
	  printf("failed to allocate memory for `PTimeTree' (%d bytes).\n",bytes);
	  endrun(1);
	}
      bytes_tot+=bytes;


      P= P_data-1;   /* start with offset 1 */
      PTimeTree--;

      //printf("\nAllocated %g MByte for particle storage.\n\n",bytes_tot/(1024.0*1024.0));
    }

  if(All.MaxPartSph>0)
    {
      bytes_tot=0;

      if(!(SphP_data=malloc(bytes=All.MaxPartSph*sizeof(struct  sph_particle_data))))
	{
	  printf("failed to allocate memory for `SphP_data' (%d bytes).\n",bytes);
	  endrun(1);
	}
      bytes_tot+=bytes;

      SphP= SphP_data-1;   /* start with offset 1 */

      //printf("Allocated %g MByte for storage of SPH data.\n\n",bytes_tot/(1024.0*1024.0));
    }
}



/* This routine frees the memory for the particle storage,
 * but we don't actually call it in the code. 
 * When the program terminats, the memory will be automatically
 * freed by the operating system.
 */
void free_memory(void)
{
  if(All.MaxPart>0)
    free(P_data);
  
  if(All.MaxPartSph>0)
    free(SphP_data);
}
