#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


/* This is a special-purpose read-in routine for initial conditions
 * produced with Bepi Tormen's initial conditions generator ZIC.
 * Using this routine requires that you know pretty well what you're
 * doing...
 * Be aware that there are unit conversion factor `massfac', `posfac',
 * and `velfac', that have to be set appropriately.
 * Also note that there is a boundary for the intermediate resolution
 * zone set by hand below (to a value of 24000.0 in this example).
 */
void read_ic_cluster(char *fname)  
{
  FILE *fd;
  int   i,j,massflag;
  double  sqr_a;
  float dummy[3],pmhr;
  int   type;
  int4byte  blklen;
  float  a0;
  double massfac, posfac, velfac; 
  int    nhr_blocks, nlr_blocks;
  int    blocks;
  int    npart, pc_here;
  int    nhr, nlr;
  int    counttype2=0, counttype3=0;
  double r2;

#define SKIP my_fread(&blklen,sizeof(int4byte),1,fd);
  
  massfac= 0.3 * 3*0.1*0.1/ (8*PI*All.G) * pow(141300.0/760, 3);
  posfac=  141300.0;
  velfac=  14130.0;

  /* Below, Bepi's new format is assumed !!!!!! */
  /* for the old one, the HR particle mass 'pmhr' has to be set
     by hand!
  */



  for(i=0, massflag=0;i<5;i++)
    All.MassTable[i]= 0;


  if((fd=fopen(fname,"r")))
    {
      fprintf(stdout,"Reading file '%s'...\n",fname); fflush(stdout);

      SKIP;
      my_fread(&nhr,sizeof(int4byte),1,fd);
      my_fread(&nlr,sizeof(int4byte),1,fd);
      my_fread(&a0,sizeof(float),1,fd);
      if(blklen==16)
	my_fread(&pmhr,sizeof(float),1,fd);
      else
	{
	  pmhr=  1.0;    /* here set by hand, if necessary */
	}
      SKIP;

      All.MassTable[1]=  pmhr * massfac;  /* high-res particles */

      printf("All.MassTable[1]=%g\n", All.MassTable[1]);

      printf("reading `%s': contains %d HR and %d LR particles.\n",
	     fname,nhr,nlr);
      
      nhr_blocks=nhr/1000000 + 1;
      nlr_blocks=nlr/1000000 + 1;

      All.TotN_gas  =N_gas  = 0;
      All.TotN_halo = nhr;
      All.TotN_disk = nlr;
      All.TotN_bulge= 0;
      All.TotN_stars= 0;

      printf("\nN_sph: %d\nN_halo: %d\nN_disk: %d\nN_bulge: %d\nN_stars: %d\n",
	     All.TotN_gas, All.TotN_halo, All.TotN_disk, All.TotN_bulge, All.TotN_stars);
           
      NumPart = All.TotNumPart =    All.TotN_gas  + All.TotN_halo 
                  + All.TotN_disk + All.TotN_bulge + All.TotN_stars;


      All.MaxPart =  All.PartAllocFactor *  All.TotNumPart;    
      All.MaxPartSph=  All.PartAllocFactor * All.TotN_gas;

      allocate_memory();


      for(blocks=0, pc_here=1; blocks<(nhr_blocks+nlr_blocks); blocks++)
	{
	  SKIP;
	  my_fread(&npart,sizeof(int),1,fd);
	  SKIP;

	  if(blocks < nhr_blocks)
	    type= 1;
	  else
	    type= 2;

	  SKIP;
	  for(i=0; i<npart; i++)
	    {
	      my_fread(&dummy[0],sizeof(float), 3, fd);

	      for(j=0;j<3;j++)
		P[pc_here+i].Pos[j]= dummy[j]; 

	    }
	  SKIP;
	  

	  SKIP;
	  for(i=0; i<npart; i++)
	    {
	      my_fread(&dummy[0], sizeof(float), 3, fd);

	      for(j=0;j<3;j++)
		P[pc_here+i].Vel[j]= dummy[j]; 
	    }
	  SKIP;
	
	  if(type==2)
	    {
	      SKIP;
	      for(i=0; i<npart; i++)
		{
		  my_fread(&dummy[0], sizeof(float), 1, fd);
		  P[pc_here+i].Mass= dummy[0]; 
		}
	      SKIP;
	    }

	  for(i=0; i<npart; i++)
	    {
	      P[pc_here+i].Type= type;
	      P[pc_here+i].ID  = pc_here+i; 
	    }

	  printf("block=%d n_in_file=%d type=%d\n",blocks,npart,type);

	  pc_here+= npart;
	}


      fclose(fd);
      fprintf(stdout,"done with reading.\n"); fflush(stdout);


      sqr_a=sqrt(All.Time);
      printf("sqr_a= %g\n", sqr_a);

      counttype3= counttype2= 0;

      for(i=1; i<=NumPart; i++)
	{
	  for(j=0;j<3;j++)
	    {
	      P[i].Pos[j] = P[i].Pos[j] * posfac; /* here in units of kpc/h */
	      
	      P[i].Vel[j] = P[i].Vel[j] * velfac; /* comoving velocity xdot on km/sec */

	      P[i].Vel[j] *= sqr_a;  /* transform to velocity variable u */

	    }

	  r2= P[i].Pos[0]*P[i].Pos[0] 
	    + P[i].Pos[1]*P[i].Pos[1] 
	    + P[i].Pos[2]*P[i].Pos[2];
          
	  if(P[i].Type==2)
	    {
	      if(sqrt(r2) > 24000.0)   /* boundary of inner LR zone */ 
		{                     
		  P[i].Type=3;
		  counttype3++;
		}
	      else
		counttype2++;
	    }
	    

	  if(P[i].Type==1)
	    {
	      P[i].Mass= All.MassTable[1];
	    }
	  else
	    {
	      P[i].Mass *= massfac;
	    }
	}
    }
  else
    {
      fprintf(stdout,"File %s not found.\n",fname); 
      endrun(7);
    }

  printf("%d particles, %d of type 2, %d of type 3\n", NumPart, counttype2, counttype3);

  All.TotN_disk = counttype2;
  All.TotN_bulge= counttype3;

  fprintf(stdout,"Baryonic particles        :  %d\n", N_gas);
  fprintf(stdout,"Collisionless particles   :  %d\n", NumPart-N_gas);
  fprintf(stdout,"                          ----------\n");
  fprintf(stdout,"Total number of particles :  %d\n\n", NumPart);
}
