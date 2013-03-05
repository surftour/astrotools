#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "allvars.h"
#include "proto.h"


#ifdef GRAPE

/*  This function computes the gravitational forces for all active particles.
 *  The GRAPE is used for this purpose.
 *  Be aware of restrictions of the dynamical range of GRAPE. The code tries
 *  to use sensible value for length- and massscales.
 */
void gravity_grape(void)
{
  int     i, j, k, pc, n, partleft, type;
  int     count,count_dummy;
  double  length_scale, max_length, mass_scale;
  double  eps, eps2;
  int     ntype[6];
  int     nchp;
  double  xyz[3], masstot[6];
  double  chipAcc[MAX_GRCHIPS][3], chipPos[MAX_GRCHIPS][3], chipPot[MAX_GRCHIPS];
  double  r2;
  double  s_a,fac1,fac2,fac3;
  double  tf0, tf1, ti, tstart, tend, timewalk=0;


  if(All.ComovingIntegrationOn)
    set_softenings();

  predict_collisionless_only(All.Time);  /* gas particles have already been predicted */

  if(N_gas)  
    {
      /* if one has also SPH, the tree is nevertheless reconstructed here because 
       * it will be needed in the neighbour finding.
       * Note that this version of the code contains only a simplistic grape interface,
       * which is not optimized for GRAPE/SPH simulations. 
       */
      All.NumForcesSinceLastTreeConstruction+= NumForceUpdate;

      tstart=second();
      if(All.NumForcesSinceLastTreeConstruction >= All.TreeUpdateFrequency*NumPart)
	{
	  force_treebuild();
	  All.NumForcesSinceLastTreeConstruction=0;
	}
      else
	ngb_update_nodes();  /* need to make sure that we have outer boundaries for all nodes */
      
      tend=second();
      All.CPU_TreeConstruction += timediff(tstart,tend);
    }

  tf0=second();

  for(i=0; i<6; i++)
    {
      ntype[i]=0;
      masstot[i]= 0;
    }

  for(i=1; i<=NumPart; i++)
    {
      ntype[P[i].Type]++;
      masstot[P[i].Type]+= P[i].Mass;
    }

   for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag,count++) 
     {
       P[i].Potential = P[i].Mass/All.SofteningTable[P[i].Type];  /* removes self energy */
       for(j=0;j<3;j++)
	 P[i].Accel[j] = 0;
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
	   max_length= length_scale*262144;


	   partleft = ntype[type];
      
	   while(partleft>0)
	    {
	      for(count=0; count<All.GrMaxPart && partleft>0; count++)  /* load particles on grape */
		{
		  for(k=0; k<3; k++)
		    {
		      xyz[k]= P[pc].PosPred[k];
		      if(fabs(xyz[k]) > max_length)
			{
			  printf("outside dynamic range!\n");
			  printf("Particle i=%d  ID=%d  Pos=(%g|%g|%g)\n",
				 pc, P[pc].ID, P[pc].PosPred[0], P[pc].PosPred[1], P[pc].PosPred[2]);
			  printf("length_scale=%g  max_length=%g\n",
				 length_scale, max_length);
			  free_grape3();
			  exit(0);
			}
		    }

		  set_xj(count, &xyz[0]); 
		  set_mj(count, P[pc].Mass);  
		  
		  pc++;
		  partleft--;
		}

	      set_n(count);

	      for(i=IndFirstUpdate, count=0; count<NumForceUpdate; ) 
		{
		  count_dummy=count;
		  for(k=i,nchp=0; (count_dummy<NumForceUpdate) && (nchp<All.GrNumChips); k=P[k].ForceFlag,count_dummy++) 
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
			

		  tstart= second();
		  calculate_force_on_x(chipPos, chipAcc, chipPot, nchp);   /* CALL GRAPE */
		  tend= second();
		  timewalk+= timediff(tstart,tend);

		  count_dummy=count;
		  for(k=i, nchp=0; (count_dummy<NumForceUpdate) && (nchp<All.GrNumChips); k=P[k].ForceFlag,count_dummy++) 
		    {
		      P[k].Potential+= chipPot[nchp];
		      
		      P[k].Accel[0] += chipAcc[nchp][0];
		      P[k].Accel[1] += chipAcc[nchp][1];
		      P[k].Accel[2] += chipAcc[nchp][2];	
			    
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
      s_a=sqrt(All.Omega0 + All.Time*(1-All.Omega0-All.OmegaLambda)+pow(All.Time,3)*All.OmegaLambda);

      fac1= All.G/ ( All.Hubble * All.Time * All.Time* s_a );

      fac2= -1.5/All.Time;

      fac3= 0.5*All.Hubble * All.Omega0 / ( All.Time * All.Time * s_a );
	
      for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag,count++)   
	{
	  for(j=0,r2=0;j<3;j++)
	    r2 += P[i].PosPred[j]*P[i].PosPred[j];
	  P[i].Potential = All.G*P[i].Potential - 0.5 * All.Omega0 *All.Hubble*All.Hubble*r2 ;

	  for(j=0;j<3;j++)
	    P[i].Accel[j] = fac1*P[i].Accel[j] 
	                      + fac2*P[i].VelPred[j] 
		              + fac3*P[i].PosPred[j] ; 
	}
    }
  else 
    {
      fac1 = All.OmegaLambda*All.Hubble*All.Hubble;

      if(All.G!=1)
	for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag,count++)   
	  {
	    P[i].Potential *= All.G;
            
            for(j=0,r2=0;j<3;j++)
                r2 += P[i].PosPred[j]*P[i].PosPred[j];
	    
            P[i].Potential += -0.5*fac1*r2;
	    
	    for(j=0;j<3;j++)
	      {
		P[i].Accel[j] *= All.G; 
		P[i].Accel[j] += fac1*P[i].PosPred[j];

              }
	  }
    }

  tf1=second();
  All.CPU_Gravity+= (ti=timediff(tf0,tf1));
  All.CPU_TreeWalk+= timewalk;

  All.TotNumOfForces+= NumForceUpdate;


#ifdef DIAG
  fprintf(FdTimings,"Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
  fprintf(FdTimings,"Nf= %d  total-Nf= %d\n", NumForceUpdate, All.TotNumOfForces);
  fprintf(FdTimings,"part/sec= %g | %g\n", NumForceUpdate/(timewalk+1.0e-20),   NumForceUpdate/(ti+1.0e-20));
  fprintf(FdTimings,"\n");
  fflush(FdTimings);
#endif
}








/*  This function queries basic parameters of the installed
 *  GRAPE system
 */
void inquire_about_grape_system(void)
{
  /* ... read system specifications of the grape system    */

  initialize_grape3(); All.GrapeAquired=1;

  All.GrMaxPart = max_particle_number();
  All.GrNumChips  = number_of_available_chips();
  All.GrNumChipsPerBoard = number_of_chips_per_board();
  All.GrMaxNbr = max_neighbor_number();

  All.GrNumBoards = All.GrNumChips/All.GrNumChipsPerBoard;

  free_grape3(); All.GrapeAquired=0;   /* ... release GRAPE board */
  
  fprintf(stdout,"\n\nSystem specifications\n\n");
  fprintf(stdout,"\tNumber of boards ................... : %6d\n",All.GrNumBoards);
  fprintf(stdout,"\tNumber of chips .................... : %6d\n",All.GrNumChips);
  fprintf(stdout,"\tNumber of chips per board .......... : %6d\n",All.GrNumChipsPerBoard);
  fprintf(stdout,"\tMaximum number of particles ........ : %6d\n",All.GrMaxPart);
  fprintf(stdout,"\tMaximum length of neighbour list ... : %6d\n",All.GrMaxNbr);
  fprintf(stdout,"\n\n");

  printf("Chips used: %d \n",All.GrNumChips);
    
  if(All.GrNumChips>MAX_GRCHIPS)
    {
      fprintf(stdout,"\n\nWARNING: more chips available than dimensioned in the code\n");
      fprintf(stdout,"reset MAX_GRCHIPS or GRAPE performance will be reduced\n");

      All.GrNumChips = MAX_GRCHIPS;
    }

  if(All.GrMaxNbr>MAX_GRNBR)
    {
      fprintf(stdout,"\n\nWARNING: more neighbours possible than dimensioned\n");
      fprintf(stdout,"\nreset MAX_GRNBR or GRAPE performance will be reduced\n");

      All.GrMaxNbr = MAX_GRNBR;
    }
}


#endif















