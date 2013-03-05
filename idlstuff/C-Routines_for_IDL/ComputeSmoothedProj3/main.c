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

float *Pos,*Mass,*Hsml;

float Xmin,Ymin,Xmax,Ymax;

int   Xpixels,Ypixels;

int   DesNgb;

int   Axis1, Axis2; 

float Hmax;

float *Value;






struct particle_2d 
{
  float Pos[2];
} **P2d;




float sb(float x,float y);




int project_and_smooth(int argc,void *argv[])
{
  int i,j,k,n,ii;
  float dummy[2],xy[2];
  float *r2list;
  int   *ngblist;
  float r,r2,h,h2,hinv2,wk,u,*hmaxfp;
  int   ngbfound;
  float *vv;
  float *dens,*mass,masstot,m;
  float *pos;
  float half_grid_size;


  void project(void);
  void allocate_2d(void);
  void set_sph_kernel(void);
  void free_memory_2d(void);


  if(argc!=13)
    {
      fprintf(stderr,"\n\nwrong number of arguments ! (found %d) \n\n",argc);
      exit(0);
    }


  /*******************************************/


  N=*(int *)argv[0];
  
  Pos =(float *)argv[1];
  Mass=(float *)argv[2];
  Hsml=(float *)argv[3];

  Xmin=*(float *)argv[4];
  Xmax=*(float *)argv[5];
  Ymin=*(float *)argv[6];
  Ymax=*(float *)argv[7];

  Xpixels=*(int *)argv[8];
  Ypixels=*(int *)argv[9];

  /* DesNgb= *(int *)argv[9]; */

  Axis1= *(int *)argv[10];
  Axis2= *(int *)argv[11];
  
  /*
  Hmax=  *(float *)argv[12];
  hmaxfp= (float *)argv[12];
  */

  Value= (float *)argv[12];

  
  printf("N=%d\n",N);

  printf("Xmin=%f\n",Xmin);
  printf("Xmax=%f\n",Xmax);

  printf("Ymin=%f\n",Ymin);
  printf("Ymax=%f\n",Ymax);

  printf("Axis1=%d\n",Axis1);
  printf("Axis2=%d\n",Axis2);


  set_sph_kernel();

  allocate_2d();

  ngb2d_treeallocate((Xpixels*Ypixels), 10*(Xpixels*Ypixels));

  project();

  ngb2d_treebuild((float **)&P2d[1], (Xpixels*Ypixels), 0,dummy,dummy);


  DesNgb= 25000.0;  /* about 11 % of the total number of pixels (assuming 480x480) */
  half_grid_size= 0.5*(Xmax-Xmin)/Xpixels;
  printf("half_grid_size= %g\n",half_grid_size); fflush(stdout);

  for(i=0,pos=Pos;i<N;i++)
    {
	  xy[0]= pos[Axis1];
	  xy[1]= pos[Axis2];
	  pos+=3;

	  Hmax= Hsml[i];
	  if(Hmax<half_grid_size)
		Hmax= half_grid_size*1.05;
	  
	  h2=ngb2d_treefind( xy, DesNgb ,0.84*Hmax,&ngblist,&r2list, Hmax,&ngbfound); 
	    
	  h=sqrt(h2);

	  if(!(i%10000))
	    {
		printf("%d..  xy= %g|%g   mass= %g  Hmax= %g   h= %g   ngbfound=%d\n",i,xy[0],xy[1],Mass[i],Hmax,h,ngbfound); fflush(stdout);
	    }

	  if(h>Hmax)
	    h=Hmax;

	  hinv2 = 1.0/(h2);


	  /* Next, we actually distribute it */
	  for(k=0;k<ngbfound;k++)
	    {
	      r=sqrt(r2list[k]);
	      vv= Value + ngblist[k];
	      
	      if(r<h)
		{
		  u = r/h;
		  ii = (int)(u*KERNEL_TABLE);
		  wk =hinv2*( Kernel[ii]    + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
		  *vv = *vv + (Mass[i] * wk);
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





void project(void)
{
  int i,j,gridi;
  float *vv;

  /* Initialize Grid */
  for(i=0;i<Xpixels;i++)
    {
      for(j=0;j<Ypixels;j++)
	{

	  vv= Value+j*Xpixels + i;
	  *vv=0;

	  gridi= j*Xpixels + i;

	  P2d[gridi+1]->Pos[0]= (Xmax-Xmin)*(i+0.5)/Xpixels + Xmin;
	  P2d[gridi+1]->Pos[1]= (Ymax-Ymin)*(j+0.5)/Ypixels + Ymin;
	}
    }
	  
}










void set_sph_kernel(void)  /* with 2D normalization */
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











void allocate_2d(void)
{
  int i;

  printf("allocating memory...\n");

  if((Xpixels*Ypixels)>0)
    {
      if(!(P2d=malloc((Xpixels*Ypixels)*sizeof(struct particle_2d *))))
	{
	  printf("failed to allocate memory. (A)\n");
	  exit(0);
	}

      P2d--;   /* start with offset 1 */
      
      if(!(P2d[1]=malloc((Xpixels*Ypixels)*sizeof(struct particle_2d))))
	{
	  printf("failed to allocate memory. (B)\n");
	  exit(0);
	}

      for(i=2;i<=(Xpixels*Ypixels);i++)   /* initiliaze pointer table */
	P2d[i]=P2d[i-1]+1;

    }

  printf("allocating memory...done\n");
}



void free_memory_2d(void)
{
  if((Xpixels*Ypixels)>0)
    {
      free(P2d[1]);
      P2d++;
      free(P2d);
    }
}
