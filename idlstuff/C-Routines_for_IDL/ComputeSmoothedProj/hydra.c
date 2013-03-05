#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


/* Calculate hydro force and rate of change of internal energy.
 * Also, we do the cooling here in an isochoric approximation,
 * with an implicit integration.
 */
void hydro_force(void)
{
  double h_i, h_i2, hinv, hinv4;
  double p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
  double r, r2, u, mass_j;
  double dx,dy,dz;
  double dwk_i, vdotr, vdotr2, visc, mu_ij, c_ij, rho_ij, h_ij, f1, f2;
  double hfc;
  double dt;	  
  int i,j,ii,n, numngb;
  double a3inv=1.0;
  double rmax;
  double h_j,dwk_j;
  int    jj,count;
  double hubble_a=0, s_a_inverse, sqrt_of_a=0, prefac=0, hfc_egy, fac_vsic_fix=0; 
  float  *r2list;
  int     *ngblist;
#ifdef COOLING
  double  unew, deltat, ne;
#endif


  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */

      sqrt_of_a= sqrt(All.Time);
      a3inv=  1/(All.Time*All.Time*All.Time);
      hubble_a= All.Hubble*sqrt(All.Omega0/(All.Time*All.Time*All.Time) 
			     + (1-All.Omega0-All.OmegaLambda)/(All.Time*All.Time) + All.OmegaLambda);
      s_a_inverse=1/(All.Hubble*sqrt(All.Omega0 + All.Time*(1-All.Omega0-All.OmegaLambda)+
				     (All.Time*All.Time*All.Time)*All.OmegaLambda));
      prefac= s_a_inverse/All.Time;
      fac_vsic_fix= hubble_a * All.Time*All.Time*All.Time;
    }



  /* Loop over all active particles */
  /* Main loop for hydro accel and adiabatic du/dt  */
  for(i=IndFirstUpdate, count=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
    {
      if(P[i].Type==0)  /* only gas particles */
	{
	  SphP[i].DtEgySpec = 0;   /* adiabatic rate will be redetermined */

	  if(SphP[i].DensityPred>0)  
	    {
	      p_over_rho2_i = SphP[i].Pressure/(SphP[i].DensityPred * SphP[i].DensityPred);
	      soundspeed_i  = sqrt(GAMMA*p_over_rho2_i*SphP[i].DensityPred);
	    }
	  else
	    p_over_rho2_i=soundspeed_i=0;

	  h_i = SphP[i].Hsml;
	  h_i2= h_i*h_i;

 	  numngb= ngb_treefind_pairs(P[i].PosPred, SphP[i].Hsml, &ngblist, &r2list);

	  for(n=0, rmax=0; n < numngb; n++)
	    {
	      j = ngblist[n]+1; 
	      r2= r2list[n];
	      h_j = SphP[j].Hsml;

	      if(r2<h_i2 || r2<h_j*h_j)
		{
		  if(i!=j)
		    {		  
		      r = sqrt(r2);		      

		      if(r>0)
			{
			  if(SphP[j].DensityPred>0)
			    {
			      p_over_rho2_j = SphP[j].Pressure/(SphP[j].DensityPred*SphP[j].DensityPred);
			      soundspeed_j  = sqrt(GAMMA*p_over_rho2_j*SphP[j].DensityPred);
			    }
			  else
			    p_over_rho2_j=soundspeed_j=0;
		  
			  dx = P[i].PosPred[0] - P[j].PosPred[0];
			  dy = P[i].PosPred[1] - P[j].PosPred[1];
			  dz = P[i].PosPred[2] - P[j].PosPred[2];
#ifdef PERIODIC
			  dx= periodic(dx);
			  dy= periodic(dy);
			  dz= periodic(dz);
#endif
			  vdotr = ( dx * (P[i].VelPred[0] - P[j].VelPred[0])
			          + dy * (P[i].VelPred[1] - P[j].VelPred[1])
			          + dz * (P[i].VelPred[2] - P[j].VelPred[2]) );

			  if(All.ComovingIntegrationOn) 
			    {
			      vdotr2 = vdotr/sqrt_of_a + hubble_a*r2;
			    }
			  else
			    vdotr2 = vdotr;
	      	  
			  if(r2<h_i2)
			    {
			      hinv = 1.0/h_i;
			      hinv4 = hinv*hinv*hinv*hinv;
			      
			      u = r*hinv;
			      ii = (int)(u*KERNEL_TABLE);
			      dwk_i= hinv4*( KernelDer[ii] + (KernelDer[ii+1]-KernelDer[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
			    }
			  else
			    dwk_i=0;
			  
			  if(r2<h_j*h_j)
			    {
			      hinv = 1.0/h_j;
			      hinv4 = hinv*hinv*hinv*hinv;
			      
			      u = r*hinv;
			      jj = (int)(u*KERNEL_TABLE);
			      dwk_j= hinv4*( KernelDer[jj] + (KernelDer[jj+1]-KernelDer[jj])*(u-KernelRad[jj])*KERNEL_TABLE);
			    }
			  else
			    dwk_j=0;
		      
			  
			  if(vdotr2<0)  /* ... artificial viscosity */		      
			    {
			      c_ij  = 0.5*(soundspeed_i + soundspeed_j);
			      h_ij  = 0.5*(h_i + h_j);
			      
			      if(All.ComovingIntegrationOn) /* comoving variables */
				mu_ij  = All.Time * h_ij*vdotr2/(r2+0.01*h_ij*h_ij);
			      else
				mu_ij  = h_ij*vdotr2/(r2+0.01*h_ij*h_ij);
			  
			      rho_ij = 0.5*(SphP[i].DensityPred+SphP[j].DensityPred);
			  
			      f1 = fabs(SphP[i].DivVel)/(fabs(SphP[i].DivVel) + SphP[i].CurlVel + 0.0001*soundspeed_i/SphP[i].Hsml);
			      f2 = fabs(SphP[j].DivVel)/(fabs(SphP[j].DivVel) + SphP[j].CurlVel + 0.0001*soundspeed_j/SphP[j].Hsml);
			  
			      if(rho_ij>0)
				visc = (-All.ArtBulkViscConst*mu_ij*c_ij + 2*All.ArtBulkViscConst*mu_ij*mu_ij)/rho_ij*(f1+f2)*0.5;
			      else
				visc=0;
			  
			      /* .... end artificial viscosity evaluation */
			      /* now make sure that viscous acceleration is not too large */
			  
			      dt = 2*(All.Time - P[i].CurrentTime); 

			      if(dt>0 && (dwk_i + dwk_j)<0)
				{
				  if(All.ComovingIntegrationOn) /* comoving variables */
				    visc = dmin(visc, fac_vsic_fix* vdotr2/
					        (0.5*(P[i].Mass+P[j].Mass)*(dwk_i + dwk_j)*r*dt)); 
				  else
				    visc = dmin(visc,  vdotr2/
						(0.5*(P[i].Mass+P[j].Mass)*(dwk_i + dwk_j)*r*dt)); 
				}
			    
			    }
			  else
			    visc=0;
		    
		      /* calculate final acceleration and rate of change of th. en. */
		      
			  mass_j= P[j].Mass;
				  
			  if(All.ComovingIntegrationOn) /* comoving variables */
			    {
			      /* arithemtic mean for symmetrization */
				/*
				  hfc  = prefac* 0.5*mass_j*(p_over_rho2_i + p_over_rho2_j + visc)*(dwk_i+dwk_j)/r;
				*/
				/* or geometric mean. Both can be chosen */
			      hfc  = prefac* 0.5*mass_j*(2*sqrt(p_over_rho2_i*p_over_rho2_j) + visc)*(dwk_i+dwk_j)/r;
				
			      hfc_egy  = hfc*All.Time*sqrt_of_a;
			    }
			  else
			    {
			      /* arithemtic mean for symmetrization */
			      /*
				hfc      = 0.5*mass_j*(p_over_rho2_i + p_over_rho2_j + visc)*(dwk_i+dwk_j)/r;
			      */
			      /* or geometric mean. Both can be chosen */
			      hfc      = 0.5*mass_j*(2*sqrt(p_over_rho2_i*p_over_rho2_j) + visc)*(dwk_i+dwk_j)/r;
			      
			      hfc_egy  = hfc; 
			    }
		  
			  P[i].Accel[0] -= hfc*dx;
			  P[i].Accel[1] -= hfc*dy;
			  P[i].Accel[2] -= hfc*dz;
		      
			  SphP[i].DtEgySpec += 0.5 * hfc_egy * vdotr2;
			}
		    }
		}
	    }
	}
    }
  
  
  
#ifdef COOLING
  /* do implicit isochoric cooling */

  for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag,count++)   
    {
      if(P[i].Type==0)  /* only gas particles */
	{
	  dt = 2*(All.Time - P[i].CurrentTime);  /*  the actual time-step */

	  ne= SphP[i].Ne;  /* electron abundance (gives ionization state and mean molecular weight) */
	  
	  if(All.ComovingIntegrationOn) /* comoving variables */
	    {
	      deltat= dt/(All.Time*hubble_a);
	      unew= DoCooling(dmax(All.MinEgySpec, SphP[i].EgySpec+SphP[i].DtEgySpec*dt), 
			      SphP[i].Density * a3inv, deltat, &ne);
	    }
	  else
	    unew= DoCooling(dmax(All.MinEgySpec, SphP[i].EgySpec+SphP[i].DtEgySpec*dt), 
			    SphP[i].Density, dt, &ne);
	  
	  
	  if(P[i].MaxPredTime > P[i].CurrentTime)  /* upon start-up, we need to protect against dt==0 */
	    SphP[i].DtEgySpec= (unew - SphP[i].EgySpec)/dt;
	  
	  SphP[i].Ne= ne;
	}
    }
#endif
}




