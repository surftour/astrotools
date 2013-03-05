#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"



static FILE *fd;

static int modus;   /* modus>0 read, modus==0 write */

static void in(int *x);
static void byten(void *x,int n);



/* This function reads or writes the restart files.
 *
 * If modus>0  the restart()-routine reads, 
 * if modus==0 it writes a restart file. 
 */
void restart(int mod)
{
  char buf[200],buf_bak[200],buf_mv[500];
  double save_PartAllocFactor;

  sprintf(buf,"%s%s",All.OutputDir,All.RestartFile);
  sprintf(buf_bak,"%s%s.bak",All.OutputDir,All.RestartFile); 
  sprintf(buf_mv,"mv %s %s", buf, buf_bak);


  modus=mod;

  if(modus) 
    {
      if(!(fd=fopen(buf,"r")))
	{
	  fprintf(stdout,"Restart file '%s' not found.\n",buf);
	  endrun(7870);
	}
    }
  else
    {
      /* system(buf_mv); */   /* move old restart files to .bak files */

      if(!(fd=fopen(buf,"w")))
	{
	  fprintf(stdout,"Restart file '%s' cannot be opened.\n",buf);
	  endrun(7878);
	}
    }
  

  save_PartAllocFactor= All.PartAllocFactor;

  byten(&All,sizeof(struct global_data_all_processes));   /* common data  */
    

  if(modus) /* read */
    {
      All.PartAllocFactor= save_PartAllocFactor;
      All.MaxPart =    All.PartAllocFactor * (All.TotNumPart); 
      All.MaxPartSph=  All.PartAllocFactor * (All.TotN_gas);

      allocate_memory();
    }


  in(&NumPart); 
  if(NumPart)
    byten(&P[1], NumPart*sizeof(struct particle_data));  /* Particle data  */

  in(&N_gas);
  if(N_gas>0)
    byten(&SphP[1], N_gas*sizeof(struct  sph_particle_data));   /* Quantities for SPH particles only */
 
  fclose(fd);
}


/* reads/writes n bytes 
 */
void byten(void *x,int n)
{
  if(modus)
    {
      if(fread(x, 1, n*sizeof(char), fd) != n*sizeof(char))
	{
	  printf("read error (restart file appears to be truncated)\n");
	  fflush(stdout);
	  endrun(7); 
	}
    }
  else
    {
      if(fwrite(x, 1, n*sizeof(char), fd) != n*sizeof(char))
	{
	  printf("write error upon writing restart file\n");
	  fflush(stdout);
	  endrun(7); 
	}
    }
} 


/* reads/writes one int 
 */
void in(int *x)
{
  size_t ret;

  if(modus)
    {
      if(fread(x, sizeof(int), 1, fd) != 1)
	{
	  printf("read error (restart file appears to be truncated)\n");
	  fflush(stdout);
	  endrun(8);
	}
    }
  else
    {
      if(fwrite(x, sizeof(int), 1, fd) != 1)
	{
	  printf("write error upon writing restart file\n");
  	  fflush(stdout);
	  endrun(8);
	}
    }
}





