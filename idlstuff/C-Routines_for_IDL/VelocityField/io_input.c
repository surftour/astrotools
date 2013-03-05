
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "allvars.h"
#include "proto.h"






int loadpositions(char *fname,int allocflag)
{
  FILE *fd;

  float dummy[3];
  int i,k;
  double Redshift;

  if(fd=fopen(fname,"r"))
    {
      fread(&N_gas,sizeof(int),1,fd);
      fread(&N_halo,sizeof(int),1,fd);
      fread(&N_disk,sizeof(int),1,fd);
      fread(&N_bulge,sizeof(int),1,fd);
      fread(&N_stars,sizeof(int),1,fd);
      /*    fread(&i,sizeof(int),1,fd); */
      fread(&Time,sizeof(double),1,fd);

      /*    fread(&Redshift,sizeof(double),1,fd); */


      if(allocflag)
	{
	  NumPartBary=N_gas;
	  NumPart =   N_gas+N_halo+N_disk+N_bulge+N_stars;
	  allocate_memory();
	}

      for(i=1;i<=NumPart;i++)
	{
	  fread(dummy,sizeof(float),3,fd);
	  for(k=0;k<3;k++)
	    P[i]->Pos[k]=dummy[k];
	}

      fclose(fd);
    }
  else
    {
      fprintf(stderr,"Error. Can't read in file '%s'\n",fname);
      exit(0);
    }
}

