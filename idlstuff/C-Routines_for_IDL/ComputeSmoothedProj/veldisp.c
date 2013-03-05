#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"



#ifdef VELDISP

/* This function computes the local velocity dispersion
 * of each particle (among particles of the same type)
 */
void veldisp(void)
{
  double vsum[3], v2sum[3];
  double h, hinv, hinv3;
  double rho, wk;
  double dx, dy, dz, r, r2, u;
  int    i,j,k,ii,n,count;
  float  *r2list;
  int    *ngblist;
 

  /* only active particles */


  for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
    {
      if(P[i].Type>0)
	{
	  rho= vsum[0] = vsum[1] = vsum[2] = v2sum[0] = v2sum[1] = v2sum[2] = 0;
  
	  P[i].HsmlVelDisp= sqrt(ngb_treefind(P[i].PosPred, All.DesNumNgb, 
					      1.1*P[i].HsmlVelDisp, 
					      P[i].Type, &ngblist, &r2list));   
	  
	  h = P[i].HsmlVelDisp;
	  hinv  = 1.0/h;
	  hinv3 = hinv*hinv*hinv;
	  
	  for(n=0; n<All.DesNumNgb; n++)
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
		  
		  rho += P[j].Mass * wk;
		}
	      
	      for(k=0; k<3; k++)
		{
		  vsum[k] +=  P[j].VelPred[k];
		  v2sum[k] += P[j].VelPred[k]*P[j].VelPred[k];
		}
	    }
	  
	  P[i].DensVelDisp = rho;
	  P[i].VelDisp = 0;
	  
	  for(k=0; k<3; k++)
	    {
	      vsum[k] /= All.DesNumNgb;
	      v2sum[k]/= All.DesNumNgb;
	      
	      P[i].VelDisp += v2sum[k] - vsum[k]*vsum[k];
	    }
	  
	  if(P[i].VelDisp > 0)
	    P[i].VelDisp= sqrt(P[i].VelDisp);
	  else
	    P[i].VelDisp= 0;
	}
    }
}
#endif
































