#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



#include "proto.h"
#include "ngbtree2d.h"


#define KERNEL_TABLE 10000
#define  PI               3.1415927

float   Kernel[KERNEL_TABLE+1],
        KernelDer[KERNEL_TABLE+1],
        KernelRad[KERNEL_TABLE+1];











int   N;

float *Pos,*Velocity_x,*Velocity_y,*Mass;

float Xmin,Ymin,Xmax,Ymax;

int   Xpixels,Ypixels;

int   DesNgb;

int   Axis1, Axis2; 

float Hmax; 



float *Value;
float *VelFd_x;
float *VelFd_y;






struct particle_2d 
{
  float Pos[2];
} **P2d;




float sb(float x,float y);




int project_and_velfield(int argc,void *argv[])
{
  int i,j,k,n,ii;
  float dummy[2],xy[2];
  float *r2list;
  int   *ngblist;
  float r,r2,h,h2,hinv2,wk,wk_tot,u,*hmaxfp;
  int   ngbfound;
  float *vv,*vx,*vy;
  float *dens,*mass,masstot,m;


  if(argc!=18)
    {
      fprintf(stderr,"\n\nwrong number of arguments ! (found %d) \n\n",argc);
      exit(0);
    }


  /*******************************************/


  N=*(int *)argv[0];
  
  Pos =(float *)argv[1];
  Velocity_x =(float *)argv[2];
  Velocity_y =(float *)argv[3];
  Mass=(float *)argv[4];

  Xmin=*(float *)argv[5];
  Xmax=*(float *)argv[6];
  Ymin=*(float *)argv[7];
  Ymax=*(float *)argv[8];

  Xpixels=*(int *)argv[9];
  Ypixels=*(int *)argv[10];

  DesNgb= *(int *)argv[11];

  Axis1= *(int *)argv[12];
  Axis2= *(int *)argv[13];
  
  Hmax=  *(float *)argv[14];
  hmaxfp= (float *)argv[14];

  Value= (float *)argv[15];
  VelFd_x= (float *)argv[16];
  VelFd_y= (float *)argv[17];

  
  printf("N=%d\n",N);

  printf("Xmin=%f\n",Xmin);
  printf("Xmax=%f\n",Xmax);

  printf("Ymin=%f\n",Ymin);
  printf("Ymax=%f\n",Ymax);

  printf("Axis1=%d\n",Axis1);


  set_sph_kernel();

  allocate_2d();

  ngb2d_treeallocate(N,10*N);

  project();

  ngb2d_treebuild((float **)&P2d[1], N, 0,dummy,dummy);



  for(i=0;i<Xpixels;i++)
    {
      /* printf("%d\n",i); */
      
      for(j=0;j<Ypixels;j++)
	{
	  xy[0]=(Xmax-Xmin)*(i+0.5)/Xpixels + Xmin;
	  xy[1]=(Ymax-Ymin)*(j+0.5)/Ypixels + Ymin;
	  
	  if(j==0)
	    h2=ngb2d_treefind( xy, DesNgb ,0,&ngblist,&r2list, Hmax,&ngbfound); 
	  else
	    h2=ngb2d_treefind( xy, DesNgb ,1.04*h,&ngblist,&r2list, Hmax,&ngbfound); 
	    
	  h=sqrt(h2);


	  if(h>Hmax)
	    h=Hmax;

	  hinv2 = 1.0/(h2);
	  

	  vv= Value+j*Xpixels + i;
	  vx= VelFd_x+j*Xpixels + i;
	  vy= VelFd_y+j*Xpixels + i;

	  *vv=0;
	  *vx=0;
	  *vy=0;
	  wk_tot=0;

	  for(k=0;k<ngbfound;k++)
	    {
	      r=sqrt(r2list[k]);
	      
	      if(r<h)
		{
		  u = r/h;
		  ii = (int)(u*KERNEL_TABLE);
		  wk =hinv2*( Kernel[ii]    + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);

		/* Average the Mass as we usually do */
		  *vv = *vv + Mass[ngblist[k]] * wk;

		/* now average x,y - velocities */
		  *vx = *vx + Velocity_x[ngblist[k]]*wk;
		  *vy = *vy + Velocity_y[ngblist[k]]*wk;
		  wk_tot = wk_tot + wk;
		}      
	    }

	  if(wk_tot>0)
	    {
		*vx = *vx/wk_tot;    /* convert from vel. density to velocity (i think???) */
		*vy = *vy/wk_tot;
	    }
	}
    } 


  /*
  dens=vector(1,N);
  mass=vector(1,N);

  for(i=1,masstot=0;i<=N;i++)
    {
      dens[i]= sb( P2d[i]->Pos[0], P2d[i]->Pos[1] ) ;
      mass[i]= Mass[i-1];
      masstot+= mass[i];
    }
  
  printf("---\n"); fflush(stdout);

  sort2_flt(N, dens, mass);

  for(i=1,m=0;i<=N;i++)
    {
      m+= mass[i];
      if(m>=masstot/2)
	{
	  *hmaxfp= dens[i]; 
	  break;
	}	  
    }

  free_vector(mass,1,N);
  free_vector(dens,1,N);
  */




  ngb2d_treefree();
 
  free_memory_2d();
  

  printf("x\n");
  return 0;
}





float sb(float x,float y)  /* returns surface brightness */
{
  int   i,j;
  float u,v;


  if(x<Xmin || x>=Xmax)
    return 0;

  if(y<Ymin || y>=Ymax)
    return 0;

  x= (x-Xmin)/(Xmax-Xmin) * Xpixels;
  y= (y-Ymin)/(Ymax-Ymin) * Ypixels;


  i=(int)x;
  j=(int)y;

  u=x-i;
  v=y-j;

  return ((1-u)*(1-v)*Value[j*Xpixels + i]+
	  (  u)*(1-v)*Value[j*Xpixels + i+1]+
	  (1-u)*(  v)*Value[(j+1)*Xpixels + i]+
	  (  u)*(  v)*Value[(j+1)*Xpixels + i+1] );
}





project()
{
  int i;
  float *pos;

  /* Positions */
  for(i=1,pos=Pos;i<=N;i++)
    {
      
      P2d[i]->Pos[0] = pos[Axis1];
      P2d[i]->Pos[1] = pos[Axis2];

      pos+=3;
    }

}










set_sph_kernel()  /* with 2D normalization */
{
  int i;
  FILE *fd;

  for(i=0;i<=KERNEL_TABLE;i++)
    KernelRad[i] = ((double)i)/KERNEL_TABLE;
      
  for(i=0;i<=KERNEL_TABLE;i++)
    {
      if(KernelRad[i]<=0.5)
	{
	  Kernel[i] = 40/(7.0*PI) *(1-6*KernelRad[i]*KernelRad[i]*(1-KernelRad[i]));
	  KernelDer[i] = 40/(7.0*PI) *( -12*KernelRad[i] + 18*KernelRad[i]*KernelRad[i]);
	}
      else
	{
	  Kernel[i] = 40/(7.0*PI) * 2*(1-KernelRad[i])*(1-KernelRad[i])*(1-KernelRad[i]);
	  KernelDer[i] = 40/(7.0*PI) *( -6*(1-KernelRad[i])*(1-KernelRad[i]));
	}
    }
  

}











allocate_2d()
{
  int i;

  printf("allocating memory...\n");

  if(N>0)
    {
      if(!(P2d=malloc(N*sizeof(struct particle_2d *))))
	{
	  printf("failed to allocate memory. (A)\n");
	  exit(0);
	}

      P2d--;   /* start with offset 1 */
      
      if(!(P2d[1]=malloc(N*sizeof(struct particle_2d))))
	{
	  printf("failed to allocate memory. (B)\n");
	  exit(0);
	}

      for(i=2;i<=N;i++)   /* initiliaze pointer table */
	P2d[i]=P2d[i-1]+1;

    }

  printf("allocating memory...done\n");
}



free_memory_2d()
{
  if(N>0)
    {
      free(P2d[1]);
      P2d++;
      free(P2d);
    }
}
