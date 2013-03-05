#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



#include "proto.h"
#include "ngbtree3d.h"


#define KERNEL_TABLE 10000
#define  PI               3.1415927

float   Kernel[KERNEL_TABLE+1],
        KernelDer[KERNEL_TABLE+1],
        KernelRad[KERNEL_TABLE+1];











int   N;

float *Pos,*Mass;

int   *Id;

int   Xpixels,Ypixels;

int   DesNgb;

float Hmax; 



float *Value;






struct particle_3d 
{
  float Pos[3];
} **P3d;






int get_neighborhood(int argc,void *argv[])
{
  int i,j,k,n,ii;
  float dummy[3],xyz[3];
  float *r2list;
  int   *ngblist;
  float r,r2,h,h2,hinv2,wk,u;
  int   ngbfound;
  float *vv;
  float *dens,*mass,masstot,m;


  void allocate_3d(void);
  void set_particle_pointer(void);
  void set_sph_kernel(void);
  void free_memory_3d(void);


  if(argc!=9)
    {
      fprintf(stderr,"\n\nwrong number of arguments ! (found %d) \n\n",argc);
      exit(0);
    }


  /*******************************************/


  N=*(int *)argv[0];
  
  Pos =(float *)argv[1];
  Mass=(float *)argv[2];
  Id=(int *)argv[3];

  Xpixels=*(int *)argv[4];
  Ypixels=*(int *)argv[5];

  DesNgb= *(int *)argv[6];

  Hmax=  *(float *)argv[7];

  Value= (float *)argv[8];

  
  printf("N=%d\n",N);


  set_sph_kernel();

  allocate_3d();

  ngb3d_treeallocate(N, 10*N);

  set_particle_pointer();

  ngb3d_treebuild((float **)&P3d[1], N, 0,dummy,dummy);



/*
  for(i=0;i<Xpixels;i++)
    {
      
      for(j=0;j<Ypixels;j++)
	{
	  xy[0]=(Xmax-Xmin)*(i+0.5)/Xpixels + Xmin;
	  xy[1]=(Ymax-Ymin)*(j+0.5)/Ypixels + Ymin;
	  
	  if(j==0)
	    h2=ngb3d_treefind( xy, DesNgb ,0,&ngblist,&r2list, Hmax,&ngbfound); 
	  else
	    h2=ngb3d_treefind( xy, DesNgb ,1.04*h,&ngblist,&r2list, Hmax,&ngbfound); 
	    
	  h=sqrt(h2);


	  if(h>Hmax)
	    h=Hmax;

	  hinv2 = 1.0/(h2);
	  

	  vv= Value+j*Xpixels + i;

	  *vv=0;

	  for(k=0;k<ngbfound;k++)
	    {
	      r=sqrt(r2list[k]);
	      
	      if(r<h)
		{
		  u = r/h;
		  ii = (int)(u*KERNEL_TABLE);
		  wk =hinv2*( Kernel[ii]    + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
		  *vv = *vv + Mass[ngblist[k]] * wk;
		}      
	    }
	}
    } 
*/


printf("Test 1: particles closest to 0,0,0\n"); fflush(stdout);
xyz[0]= 0.0;
xyz[1]= 0.0;
xyz[2]= 0.0;
h2=ngb3d_treefind( xyz, DesNgb ,1.04*h,&ngblist,&r2list, Hmax,&ngbfound);

printf("h= %g\n",sqrt(h2));

hinv2 = 1.0/(h2);

for(k=0;k<ngbfound;k++)
  { 
    r=sqrt(r2list[k]);
    printf("k= %d  id= %d  mass= %g  r=%g\n",k,Id[ngblist[k]],Mass[ngblist[k]],r);
  }





printf("Test 2: particles closest to 10,10,10\n"); fflush(stdout);
xyz[0]= 10.0;
xyz[1]= 10.0;
xyz[2]= 10.0;
h2=ngb3d_treefind( xyz, DesNgb ,1.04*h,&ngblist,&r2list, Hmax,&ngbfound);

printf("h= %g\n",sqrt(h2));

hinv2 = 1.0/(h2);

for(k=0;k<ngbfound;k++)
  { 
    r=sqrt(r2list[k]);
    printf("k= %d  id= %d  mass= %g  r=%g\n",k,Id[ngblist[k]],Mass[ngblist[k]],r);
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




  ngb3d_treefree();
 
  free_memory_3d();
  

  printf("x\n");
  return 0;
}






void set_particle_pointer(void)
{
  int i;
  float *pos;

  for(i=1,pos=Pos;i<=N;i++)
    {
      
      P3d[i]->Pos[0] = pos[0];
      P3d[i]->Pos[1] = pos[1];
      P3d[i]->Pos[2] = pos[2];

      pos+=3;
    }
}







void set_sph_kernel(void)  /* with appropriate normalization */
{
  int i;

  for(i=0;i<=KERNEL_TABLE;i++)
    KernelRad[i] = ((double)i)/KERNEL_TABLE;

  Kernel[KERNEL_TABLE+1] = KernelDer[KERNEL_TABLE+1]= 0;
      
  for(i=0;i<=KERNEL_TABLE;i++)
    {
      if(KernelRad[i]<=0.5)
	{
	/* 2D normalization 
	  Kernel[i] = 40/(7.0*PI) *(1-6*KernelRad[i]*KernelRad[i]*(1-KernelRad[i]));
	  KernelDer[i] = 40/(7.0*PI) *( -12*KernelRad[i] + 18*KernelRad[i]*KernelRad[i]);
	*/
	/* 3D normalization */
	  Kernel[i] = 8/PI *(1-6*KernelRad[i]*KernelRad[i]*(1-KernelRad[i]));
	  KernelDer[i] = 8/PI *( -12*KernelRad[i] + 18*KernelRad[i]*KernelRad[i]);
	}
      else
	{
	/* 2D normalization
	  Kernel[i] = 40/(7.0*PI) * 2*(1-KernelRad[i])*(1-KernelRad[i])*(1-KernelRad[i]);
	  KernelDer[i] = 40/(7.0*PI) *( -6*(1-KernelRad[i])*(1-KernelRad[i]));
	*/
	/* 2D normalization */
	  Kernel[i] = 8/PI * 2*(1-KernelRad[i])*(1-KernelRad[i])*(1-KernelRad[i]);
	  KernelDer[i] = 8/PI *( -6*(1-KernelRad[i])*(1-KernelRad[i]));
	}
    }
  

}











void allocate_3d(void)
{
  int i;

  printf("allocating memory...\n");

  if(N>0)
    {
      if(!(P3d=malloc(N*sizeof(struct particle_3d *))))
	{
	  printf("failed to allocate memory. (A)\n");
	  exit(0);
	}

      P3d--;   /* start with offset 1 */
      
      if(!(P3d[1]=malloc(N*sizeof(struct particle_3d))))
	{
	  printf("failed to allocate memory. (B)\n");
	  exit(0);
	}

      for(i=2;i<=N;i++)   /* initiliaze pointer table */
	P3d[i]=P3d[i-1]+1;

    }

  printf("allocating memory...done\n");
}



void free_memory_3d(void)
{
  if(N>0)
    {
      free(P3d[1]);
      P3d++;
      free(P3d);
    }
}
