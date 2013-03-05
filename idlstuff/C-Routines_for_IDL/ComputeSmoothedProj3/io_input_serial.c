
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "allvars.h"
#include "proto.h"




struct io_header_1
{
  int npart[6];
  double   mass[6];
  double   time;
  double   redshift;

  char     fill[256 - 6*4 - 6*8 - 2*8];

} header1;
 


#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

int loadpositions_serial(char *fname,int files,int allocflag)
{
  FILE *fd;

  char buf[200];
  int i,k,dummy;
  double Redshift;
  int t,j,n,off,pc;
  float xyz[3];

  NumPart=0;


  sprintf(buf,"%s",fname);
  if(!(fd=fopen(buf,"r")))
    {
      printf("can't open file `%s`\n",buf);
      exit(0);
    }
  
  fread(&dummy, sizeof(dummy), 1, fd);
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);
  
  for(k=0;k<5;k++)
    NumPart+=header1.npart[k];
  
  fclose(fd);

  
  printf("number of particles: %d\n",NumPart);

  allocate_memory();


  printf("time %g \n",header1.time);
  printf("redshift %g\n",header1.redshift);



  for(i=0,pc=1;i<files;i++)
    {
      sprintf(buf,"%s",fname);
      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file `%s`\n",buf);
	  exit(0);
	}

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      n=0;

      for(k=0;k<5;k++)
	n+=header1.npart[k];

      N_gas=header1.npart[0];
      N_halo=header1.npart[1];
      N_disk=header1.npart[2];
      N_bulge=header1.npart[3];
      N_stars=header1.npart[4];


      for(t=0,off=0;t<5;t++)
	for(k=1;k<=header1.npart[t];k++)
	  {
	    P[pc+off]->Mass=header1.mass[t];
	    /*
	    P[pc+off]->Type=t;
	    */
	    off++;
	  }


      fread(&dummy, sizeof(dummy), 1, fd);
      for(k=1,off=0;k<=n;k++)
	{
	  fread(&xyz[0],sizeof(float),3,fd);
	  for(j=0;j<3;j++)
	    P[pc+off]->Pos[j]=xyz[j];
	  off++;
	}
      fread(&dummy, sizeof(dummy), 1, fd);



      SKIP;
      fseek(fd, n*sizeof(float)*3, SEEK_CUR);  /* velocities */
      SKIP;

      
      SKIP;
      for(k=1,off=0;k<=n;k++)
	{
	  /*
	  fread(&P[pc+off]->ID,sizeof(int),1,fd);
	  */
	  off++;
	}
      SKIP;
      

      pc+=off;
	  
      fclose(fd);
    }

}

