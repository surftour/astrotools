#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "overhead.h" 
#include "proto.h"


/* And now Paul's more straightforward allocation routines for the
 *   direct particle structures (Paul Martini, Sept. 2004)
 */
int allocate_gas(void)
{
  extern struct particle_gas *PG; 
  if(!(PG=malloc(All.N_gas*sizeof(struct particle_gas))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
}
int allocate_halo(void)
{
  extern struct particle_halo *PH; 
  if(!(PH=malloc(All.N_halo*sizeof(struct particle_halo))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
}
int allocate_disk(void)
{
  extern struct particle_disk *PD; 
  if(!(PD=malloc(All.N_disk*sizeof(struct particle_disk))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
}
int allocate_bulge(void)
{
  extern struct particle_bulge *PB; 
  if(!(PB=malloc(All.N_bulge*sizeof(struct particle_bulge))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
}
int allocate_star(void)
{
  extern struct particle_star *PS; 
  if(!(PS=malloc(All.N_star*sizeof(struct particle_star))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
}
int allocate_bh(void)
{
  extern struct particle_bh *PBH; 
  if(!(PBH=malloc(All.N_bh*sizeof(struct particle_bh))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
}




/* ========================================================================== */
/* ========================================================================== */

/* Some array utilities taken/modified from Jim Stone & Tom Gardiner's 
     ATHENA source code */

/* These functions construct and destruct 2D arrays. */
/* Elements in these arrays can be accessed in the standard fashion of
   array[in][im] (where in = 0 -> nr-1 and im = 0 -> nc-1) as if it
   were an array of fixed size at compile time with the statement:

   "Type" array[nr][nc];

   It is the responsibility of the user to remember and respect the
   array bounds.  -- Tom Gardiner -- 8/2/2001 */

/* ========================================================================== */
/* ========================================================================== */


/* Free the 2-D array: array[nr][nc]
   Usage: free_2d_array((void **)array);
*/
void free_2d_array(void **array){

  free(array[0]);
  free(array);

  return;
}


/* Construct a 2-D array: array[nr][nc] 
   Usage: array = (double **)calloc_2d_array(nr,nc,sizeof(double));
   NOTE --- the 2D and 3D arrays created thus run from index 0 to n-1
*/
void** calloc_2d_array(size_t nr, size_t nc, size_t size){

  void **array;
  size_t i;

  if((array = (void **)calloc(nr,sizeof(void*))) == NULL){
    return NULL;
  }

  if((array[0] = (void *)calloc(nr*nc,size)) == NULL){
    free((void *)array);
    return NULL;
  }

  for(i=0; i<nr; i++){
    array[i] = (void *)((unsigned char *)array[0] + i*nc*size);
  }

  return array;
}

/* Free the 2-D array: array[nr][nc]
   Usage: free_2d_array((void **)array);
*/
void free_3d_array(void ***array){

  free(array[0]);
  free(array);

  return;
}


/* Construct a 3-D array: array[nr][nc][nz] 
   Usage: array = (double ***)calloc_3d_array(nr,nc,nz,sizeof(double));
*/
void*** calloc_3d_array(size_t nr, size_t nc, size_t nz, size_t size){
  ////printf("\n CHECK 3.1 ... \n");

  void ***array;
  size_t i;

  if((array = (void ***)calloc(nr,sizeof(void*))) == NULL){
    return NULL;
  }

  if((array[0] = (void **)calloc(nr*nz,sizeof(void*))) == NULL){
    return NULL;
  }

  if((array[0][0] = (void *)calloc(nr*nc*nz,size)) == NULL){
    free((void *)array);
    return NULL;
  }

  ////printf("\n CHECK 3.2 ... \n");
  for(i=0; i<nr; i++){
    array[i] = (void **)calloc_2d_array(nc,nz,size);
  }
  ////printf("\n CHECK 3.3 ... \n");

  return array;
}



