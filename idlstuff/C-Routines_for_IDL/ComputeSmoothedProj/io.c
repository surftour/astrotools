#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"



/* This wrapper function select the desired output 
 * routine for snapshot files.
 */
void savepositions(int num)
{
  double t0, t1;

  t0=second();

  printf("\nwriting snapshot file... \n");
 
  savepositions_ioformat1(num);

  printf("done with snapshot.\n");
  
  t1=second();

  All.CPU_Snapshot+= timediff(t0,t1);
}




/* This function writes a snapshot of the particle ditribution to
 * one file using Gadget's default file format.
 * Each snapshot file contains a header first, then particle positions, 
 * velocities and ID's.
 * Then particle masses are written for those particle types with zero entry in
 * MassTable.
 * After that, first the internal energies u, and then the density is written
 * for the SPH particles.
 * Finally, if cooling is enabled, the mean molecular weight is written for the gas
 * particles. 
 */
void savepositions_ioformat1(int num)
{
  FILE *fd;
  char buf[100];
  float dummy[3];
  int i,k;
  int   blklen,masscount;
  double a3inv;
#ifdef COOLING
  double ne, nh0;
#endif

#define BLKLEN my_fwrite(&blklen, sizeof(blklen), 1, fd);


  if(All.ComovingIntegrationOn)
    a3inv=  1/(All.Time*All.Time*All.Time);
  else
    a3inv=  1.0;

  sprintf(buf,"%s%s_%03d",All.OutputDir,All.SnapshotFileBase,num);

  if((fd=fopen(buf,"w")))
    {
      header1.npart[0]= header1.npartTotal[0]= All.TotN_gas;
      header1.npart[1]= header1.npartTotal[1]= All.TotN_halo;
      header1.npart[2]= header1.npartTotal[2]= All.TotN_disk;
      header1.npart[3]= header1.npartTotal[3]= All.TotN_bulge;
      header1.npart[4]= header1.npartTotal[4]= All.TotN_stars;
      header1.npart[5]= header1.npartTotal[5]= 0;

      for(i=0;i<6;i++)
	header1.mass[i]=0;

      for(i=0, masscount=0; i<5; i++)
	{
	  header1.mass[i]= All.MassTable[i];
	  if(All.MassTable[i]==0 && header1.npart[i]>0)
	    masscount+= header1.npart[i];
	}

      header1.time= All.Time;

      if(All.ComovingIntegrationOn)
	header1.redshift=1.0/All.Time - 1.0;
      else
	header1.redshift=0;  

      
      header1.flag_sfr=0;
      header1.flag_feedback=0;
      header1.flag_cooling= 0;
#ifdef COOLING
      header1.flag_cooling= 1;
#endif
      header1.num_files= 1;
      header1.BoxSize= All.BoxSize;
      header1.Omega0=  All.Omega0;
      header1.OmegaLambda= All.OmegaLambda;
      header1.HubbleParam= All.HubbleParam;
      
      blklen=sizeof(header1);
      BLKLEN;
      my_fwrite(&header1, sizeof(header1), 1, fd);
      BLKLEN;


      blklen=NumPart*3*sizeof(float);

      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummy[k]=P[i].PosPred[k];
	  my_fwrite(dummy,sizeof(float),3,fd);
	}
      BLKLEN;


      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummy[k]=P[i].VelPred[k];
	  my_fwrite(dummy,sizeof(float),3,fd);
	}
      BLKLEN;
 

      blklen=NumPart*sizeof(int);
      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  my_fwrite(&P[i].ID,sizeof(int),1,fd);
	}
      BLKLEN;

      blklen=masscount*sizeof(float);
      if(masscount)
	BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  dummy[0]= P[i].Mass;
	  if(All.MassTable[P[i].Type]==0)
	    my_fwrite(dummy,sizeof(float),1,fd);
	}
      if(masscount)
	BLKLEN;

      if(N_gas)
	{
	  blklen=N_gas*sizeof(float);
	  BLKLEN;
	  for(i=1;i<=N_gas;i++)
	    {
	      dummy[0]=SphP[i].EgySpecPred;
	      my_fwrite(dummy,sizeof(float),1,fd);
	    }
	  BLKLEN;


	  blklen=N_gas*sizeof(float);  /* added density  */
	  BLKLEN;
	  for(i=1;i<=N_gas;i++)
	    {
	      dummy[0]=SphP[i].DensityPred;
	      my_fwrite(dummy,sizeof(float),1,fd);
	    }
	  BLKLEN;

#ifdef COOLING
	  blklen=N_gas*sizeof(float);  /* electron abundance */
	  BLKLEN;
	  for(i=1;i<=N_gas;i++)
	    {
	      dummy[0]= SphP[i].Ne;
	      my_fwrite(dummy,sizeof(float),1,fd);
	    }
	  BLKLEN;


	  blklen=N_gas*sizeof(float);  /* neutral hydrogen */
	  BLKLEN;
	  for(i=1;i<=N_gas;i++)
	    {
	      ne= SphP[i].Ne;

	      AbundanceRatios(SphP[i].EgySpecPred, SphP[i].DensityPred*a3inv,
			      &ne, &nh0);
	      dummy[0]= nh0;
	      my_fwrite(dummy,sizeof(float),1,fd);
	    }
	  BLKLEN;
#endif
	  blklen=N_gas*sizeof(float);  /* hsml  */
	  BLKLEN;
	  for(i=1;i<=N_gas;i++)
	    {
	      dummy[0]=SphP[i].Hsml;
	      my_fwrite(dummy,sizeof(float),1,fd);
	    }
	  BLKLEN;
	}
      fclose(fd);
    }
  else
    {
      fprintf(stdout,"Error. Can't write in file '%s'\n", buf);
      endrun(10);
    }
}



/* This catches I/O errors occuring for my_my_fwrite(). In this case we better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nwritten;

  if((nwritten=fwrite(ptr, size, nmemb, stream))!=nmemb)
    {
      printf("I/O error (fwrite) on has occured.\n");
      fflush(stdout);
      endrun(777);
    }
  return nwritten;
}


/* This catches I/O errors occuring for fread(). In this case we better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nread;

  if((nread=fread(ptr, size, nmemb, stream))!=nmemb)
    {
      printf("I/O error (fread) has occured.\n");
      fflush(stdout);
      endrun(778);
    }
  return nread;
}








