#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef GRAPE
/* This function computes the gravitational potential for ALL the particles.
 * It uses GRAPE for that.
 */

void compute_potential_grape() 
{
  double  t0,t1;
  double  fac;
  int     i, j, k, pc, n, partleft, type;
  int     count,count_dummy;
  double  length_scale,mass_scale;
  double  eps, eps2;
  int     ntype[6];
  int     nchp;
  double  xyz[3], masstot[6];
  double  chipAcc[MAX_GRCHIPS][3], chipPos[MAX_GRCHIPS][3], chipPot[MAX_GRCHIPS];
  double  r2;

  t0=second();

  if(All.ComovingIntegrationOn)
    set_softenings();

  printf("Start potential computation..."); fflush(stdout);

  for(i=0; i<6; i++)
    {
      ntype[i]=0;
      masstot[i]= 0;
    }

  for(i=1; i<=NumPart; i++)
    {
      ntype[P[i].Type]++;
      masstot[P[i].Type]+= P[i].Mass;
      P[i].Potential = P[i].Mass/All.SofteningTable[P[i].Type];  /* removes self energy */
    }


   if(All.GrapeAquired==0)
     {
      initialize_grape3();  All.GrapeAquired=1;    /* ... acquire GRAPE board*/
     }
  
   for(type=0, pc=1; type<=4; type++)
     {
       if(ntype[type]>0)
	 {
	   length_scale= All.SofteningTable[type]/4;
	   mass_scale=   masstot[type]/ntype[type]/10;
	   eps=          All.SofteningTable[type];
	
	   set_scales(length_scale,mass_scale);
	   set_mass_correction_factor(1.0);
	   
	   partleft = ntype[type];
      
	   while(partleft>0)
	    {
	      for(count=0; count<All.GrMaxPart && partleft>0; count++)  /* load particles on grape */
		{
		  for(k=0; k<3; k++)
		    xyz[k]= P[pc].PosPred[k];
		  
		  set_xj(count, &xyz[0]); 
		  set_mj(count, P[pc].Mass);  
		  
		  pc++;
		  partleft--;
		}

	      set_n(count);

	      for(i=1, count=0; count<NumPart; ) 
		{
		  count_dummy=count;
		  for(k=i,nchp=0; (count_dummy<NumPart) && (nchp<All.GrNumChips); k++,count_dummy++) 
		    {
		      eps2= dmax(All.SofteningTable[P[k].Type], eps);
		      eps2*=eps2;
		      set_eps2_to_chip(eps2,nchp);
			    
		      chipPos[nchp][0] = P[k].PosPred[0];
		      chipPos[nchp][1] = P[k].PosPred[1];
		      chipPos[nchp][2] = P[k].PosPred[2];
		      
		      set_h2(0,nchp);
		      
		      nchp++;
		    }
			
		  for(n=nchp;n<All.GrNumChips;n++)
		    set_h2(0,n);
		 
		  calculate_force_on_x(chipPos, chipAcc, chipPot, nchp);   /* CALL GRAPE */
		 
		  count_dummy=count;
		  for(k=i, nchp=0; (count_dummy<NumPart) && (nchp<All.GrNumChips); k++,count_dummy++) 
		    {
		      P[k].Potential+= chipPot[nchp];
			    
		      nchp++;
		    }
			
		  i=k;
		  count=count_dummy;
		}
	    }
	  
	}
    }
  
  
  if(All.FreeGrape)
    {
      free_grape3(); All.GrapeAquired=0;   /* ... release GRAPE board */
    }


  if(All.ComovingIntegrationOn)
    {
      fac=0.5*All.Omega0*All.Hubble*All.Hubble;

      for(i=1;i<=NumPart;i++)
	{
	  for(j=0, r2=0; j<3; j++)
	    r2 += P[i].PosPred[j]*P[i].PosPred[j];

	  P[i].Potential = All.G*P[i].Potential - fac*r2;
	}
    }
  else
    {
      fac= -0.5*All.OmegaLambda*All.Hubble*All.Hubble;

      for(i=1;i<=NumPart;i++)
	{
	  P[i].Potential *= All.G;

	  if(fac!=0)
	     {
	       for(j=0,r2=0;j<3;j++)
		 r2 += P[i].PosPred[j]*P[i].PosPred[j];
	       
	       P[i].Potential += fac*r2;
	     }
	}
    }

  printf("done.\n"); fflush(stdout);

  t1=second();

  All.CPU_Potential+= timediff(t0,t1);
}

#endif

















