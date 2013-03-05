#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"



/* This function computes the local density for each active SPH particle,
 * It sets the smoothing radius to a specified number of neighbours,
 * and the divergence and rotation of the velocity field are computed.
 * The pressure is updated as well.
 */
void density(void)
{
  double rotv[3];
  double h, hinv, hinv3, hinv4;
  double rho, divv, wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  double dt;
  int    i,j,ii,n,count;
  double hubble_a, prefac=0;
  float  *r2list;
  int     *ngblist;

 
  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      hubble_a=All.Hubble*sqrt(All.Omega0/(All.Time*All.Time*All.Time) 
			       + (1-All.Omega0-All.OmegaLambda)/(All.Time*All.Time) + All.OmegaLambda);
      
      prefac= 1.0/(hubble_a*pow(All.Time, 1.5));
    }

 

  /* only active gas particles */

  for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag,count++)   
    {
      if(P[i].Type==0) 
	{
	  rho = divv = rotv[0] = rotv[1] = rotv[2] = 0;
	  
	  if(SphP[i].Hsml < All.MinGasHsml)
	    SphP[i].Hsml= All.MinGasHsml;
	  
	  if(SphP[i].Hsml < 1.01*All.MinGasHsml)
	    {
	      SphP[i].NumNgb= ngb_treefind_variable(P[i].PosPred, SphP[i].Hsml, 0, &ngblist, &r2list);   

	      if(SphP[i].NumNgb<All.DesNumNgb)
		{
		  SphP[i].Hsml=ngb_treefind(P[i].PosPred, All.DesNumNgb , 1.1*SphP[i].Hsml, 0, &ngblist, &r2list);   
		  SphP[i].NumNgb=All.DesNumNgb; 
		  SphP[i].Hsml = sqrt(SphP[i].Hsml);
		}
	    }
	  else
	    {
	      SphP[i].Hsml=ngb_treefind(P[i].PosPred, All.DesNumNgb , 1.1*SphP[i].Hsml, 0, &ngblist, &r2list);   
	      SphP[i].NumNgb=All.DesNumNgb; 
	      SphP[i].Hsml = sqrt(SphP[i].Hsml); 

	      if(SphP[i].Hsml < All.MinGasHsml)
		{
		  SphP[i].Hsml= All.MinGasHsml;
		  SphP[i].NumNgb= ngb_treefind_variable(P[i].PosPred, SphP[i].Hsml, 0, &ngblist, &r2list);   
		}
	    }

	  h = SphP[i].Hsml;
	  hinv  = 1.0/h;
	  hinv3 = hinv*hinv*hinv;
	  hinv4 = hinv3*hinv;
  
	  for(n=0; n<SphP[i].NumNgb; n++)
	    {
	      j  = ngblist[n]+1; 

	      dx = P[i].PosPred[0] - P[j].PosPred[0];
	      dy = P[i].PosPred[1] - P[j].PosPred[1];
	      dz = P[i].PosPred[2] - P[j].PosPred[2];
#ifdef PERIODIC
	      dx= periodic(dx);
	      dy= periodic(dy);
	      dz= periodic(dz);
#endif
	      r2 = dx*dx + dy*dy + dz*dz;
	      
	      r = sqrt(r2);
	      
	      if(r<h)
		{
		  u = r*hinv;
		  
		  ii = (int)(u*KERNEL_TABLE);
		  
		  wk =hinv3*( Kernel[ii]  + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
		  dwk=hinv4*( KernelDer[ii] + (KernelDer[ii+1]-KernelDer[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
		  
		  mass_j = P[j].Mass;
		  
		  rho += mass_j * wk;
		  
		  if(i!=j)
		    {
		      divv -= mass_j * dwk/r *
			( dx * (P[i].VelPred[0] - P[j].VelPred[0])
		        + dy * (P[i].VelPred[1] - P[j].VelPred[1])
		        + dz * (P[i].VelPred[2] - P[j].VelPred[2]) );
		  
		      rotv[0] += mass_j * dwk/r *
			  (  dz * (P[i].VelPred[1] - P[j].VelPred[1])
			   - dy * (P[i].VelPred[2] - P[j].VelPred[2]) );

		      rotv[1] += mass_j * dwk/r *
		           (  dx * (P[i].VelPred[2] - P[j].VelPred[2])
      	                    - dz * (P[i].VelPred[0] - P[j].VelPred[0]) );

		      rotv[2] += mass_j * dwk/r *
		           ( dy * (P[i].VelPred[0] - P[j].VelPred[0])
			   - dx * (P[i].VelPred[1] - P[j].VelPred[1]) );
		    }
		}
	      
	    }      

	  SphP[i].CurlVel= sqrt(rotv[0]*rotv[0] + rotv[1]*rotv[1] + rotv[2]*rotv[2])/rho;
	  SphP[i].Density= SphP[i].DensityPred = rho;
	  SphP[i].DivVel= divv/rho;
	  if(All.ComovingIntegrationOn) /* comoving variables */
	    {
	      SphP[i].DtDensity= - prefac*divv;
	      SphP[i].DtHsml= - SphP[i].Hsml*SphP[i].DtDensity/(3*SphP[i].Density);
	    }
	  else
	    {  /* physical */
	      SphP[i].DtDensity= -divv;
	      SphP[i].DtHsml= SphP[i].Hsml*SphP[i].DivVel/3;
	    }
    
	  /* make sure that predicted values of density and smoothing length
             cannot take on negative values */

	  dt = 2*(All.Time - P[i].CurrentTime);  /*  the actual time-step */
	  if(dt>0)
	    { 
	      SphP[i].DtHsml+= SphP[i].Hsml/(2*dt)*(pow(((double)All.DesNumNgb)/SphP[i].NumNgb, 1.0/3)-1);

	      SphP[i].DtDensity= dmax(-0.9*SphP[i].Density/dt , SphP[i].DtDensity);
	      SphP[i].DtHsml   = dmax(-0.9*SphP[i].Hsml/dt ,    SphP[i].DtHsml);
	    }


	  SphP[i].Pressure = GAMMA_MINUS1*(SphP[i].EgySpecPred)*SphP[i].Density;
      
	}
    }

}


/*  this function wraps the distance x to the closest image
 *  for the given box size
 */
double INLINE_FUNC periodic(double x)
{
  while(x > All.BoxHalf)
    x -=All.BoxSize;

  while(x < -All.BoxHalf)
      x+=All.BoxSize;

  return x;
}

































