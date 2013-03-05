
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



#include "allvars.h"
#include "proto.h"





allocate_memory()
{
  int i;

  fprintf(stderr,"allocating memory...\n");


  if(NumPart>0)
    {

      if(!(P=malloc(NumPart*sizeof(struct particle_data *))))
	{
	  fprintf(stderr,"failed to allocate memory. (A)\n");
	  exit(0);
	}

      P--;   /* start with offset 1 */
      
      if(!(P[1]=malloc(NumPart*sizeof(struct particle_data))))
	{
	  fprintf(stderr,"failed to allocate memory. (B)\n");
	  exit(0);
	}

      for(i=2;i<=NumPart;i++)   /* initiliaze pointer table */
	P[i]=P[i-1]+1;


      if(!(PartZone=malloc((NumPart+1)*sizeof(int))))
	{
	  fprintf(stderr,"failed to allocate memory. (E)\n");
	  exit(0);
	}

      if(!(IndOrdMassZone=malloc((NumPart+1)*sizeof(int))))
	{
	  fprintf(stderr,"failed to allocate memory. (F)\n");
	  exit(0);
	}
    }


  if(NumPartBary>0)
    {
      if(!(SphP=malloc(NumPartBary*sizeof(struct  sph_particle_data *))))
	{
	  fprintf(stderr,"failed to allocate memory. (C)\n");
	  exit(0);
	}
      
      SphP--;   /* start with offset 1 */
      
      if(!(SphP[1]=malloc(NumPartBary*sizeof(struct  sph_particle_data))))
	{
	  fprintf(stderr,"failed to allocate memory. (D)\n");
	  exit(0);
	}
      
      for(i=2;i<=NumPartBary;i++)   /* initiliaze pointer table */
	SphP[i]=SphP[i-1]+1;


    }

  fprintf(stderr,"allocating memory...done\n");

}






free_memory()
{
  int i;

  if(NumPart>0)
    {
      free(P[1]);
      P++;
      free(P);
    }
      
  if(NumPartBary>0)
    {
      free(SphP[1]);
      SphP++;
      free(SphP);
    }
}
