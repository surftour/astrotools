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
void grid_density(struct particle_data *GridP, struct sph_particle_data *GridSphP, int ngridp)
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

  for(i=0; i<ngridp;i++)
    {
      if(GridP[i].Type==0) 
	{
	  rho = divv = rotv[0] = rotv[1] = rotv[2] = 0;
	  
	  if(GridSphP[i].Hsml < All.MinGasHsml)
	    GridSphP[i].Hsml= All.MinGasHsml;
	  //printf("All.hsml %e GridSphP Hsml %e\n",All.MinGasHsml,GridSphP[i].Hsml);
	  
	  if(GridSphP[i].Hsml < 1.01*All.MinGasHsml)
	    {
	      GridSphP[i].NumNgb= ngb_treefind_variable_2d(GridP[i].PosPred, GridSphP[i].Hsml, 0, &ngblist, &r2list);   
	  	//printf("here now MGH %e Hsml %e %d %d\n",All.MinGasHsml,GridSphP[i].Hsml,GridSphP[i].NumNgb,All.DesNumNgb);

	      if(GridSphP[i].NumNgb<All.DesNumNgb)
		{
		  GridSphP[i].Hsml=ngb_treefind_2d(GridP[i].PosPred, All.DesNumNgb , 1.1*GridSphP[i].Hsml, 0, &ngblist, &r2list);   
		  GridSphP[i].NumNgb=All.DesNumNgb; 
		  GridSphP[i].Hsml = sqrt(GridSphP[i].Hsml);
		}
	    }
	  else
	    {
	      GridSphP[i].Hsml=ngb_treefind_2d(GridP[i].PosPred, All.DesNumNgb , 1.1*GridSphP[i].Hsml, 0, &ngblist, &r2list);   
	      GridSphP[i].NumNgb=All.DesNumNgb; 
	      GridSphP[i].Hsml = sqrt(GridSphP[i].Hsml); 

	      if(GridSphP[i].Hsml < All.MinGasHsml)
		{
		  GridSphP[i].Hsml= All.MinGasHsml;
		  GridSphP[i].NumNgb= ngb_treefind_variable_2d(GridP[i].PosPred, GridSphP[i].Hsml, 0, &ngblist, &r2list);   
		}
	    }

	  h = GridSphP[i].Hsml;
	  hinv  = 1.0/h;
	  hinv3 = hinv*hinv*hinv;
	  hinv4 = hinv3*hinv;
  
	  for(n=0; n<GridSphP[i].NumNgb; n++)
	    {
	      j  = ngblist[n]+1; 

	      dx = GridP[i].PosPred[0] - Pn[j].PosPred[0];
	      dy = GridP[i].PosPred[1] - Pn[j].PosPred[1];
	      dz = GridP[i].PosPred[2] - Pn[j].PosPred[2];
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
		  wk*=(40.0/56.0);
		  //dwk=hinv4*( KernelDer[ii] + (KernelDer[ii+1]-KernelDer[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
		  
		  mass_j = Pn[j].Mass;
		  
		  rho += mass_j * wk;
		  
		}
	      
	    }      

	  GridSphP[i].Density= GridSphP[i].DensityPred = rho;
	 //printf("i %d here now MGH %e Hsml %e %d %d %e\n",i,All.MinGasHsml,GridSphP[i].Hsml,GridSphP[i].NumNgb,All.DesNumNgb,GridSphP[i].Density);
      
	}
    }

}




















