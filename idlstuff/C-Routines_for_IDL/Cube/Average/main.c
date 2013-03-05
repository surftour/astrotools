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











int   N,Ngrid;

float *Pos,*Mass,*QuantityToAverage;

/*  Note: all these fields
    are currently absorbed
    into the Pos structure.
float *Velx,*Vely,*Velz;

int   *Id;
*/

int   Xpixels,Ypixels;

int   DesNgb;

float Hmax; 



float *Value;






struct particle_3d 
{
  float Pos[3];
  float Id;
  float Grid[3];
} **P3d;





int average(int argc,void *argv[])
{
  int i,j,k,n,ii;
  float dummy[3],xyz[3];
  float *r2list;
  int   *ngblist;
  float r,r2,h,h2,hinv,hinv3,wk,u;
  int   ngbfound;
  float *vv;
  float *dens,*mass,masstot,m;
  float Density, DensityWeightedAvg;


  void allocate_3d(void);
  void set_particle_pointer(void);
  void set_sph_kernel(void);
  void free_memory_3d(void);


  if(argc!=8)
    {
      fprintf(stderr,"\n\nwrong number of arguments ! (found %d) \n\n",argc);
      exit(0);
    }


  /*******************************************/


  N=*(int *)argv[0];
  Ngrid=*(int *)argv[1];
  
  Pos =(float *)argv[2];
  Mass=(float *)argv[3];
  DesNgb= *(int *)argv[4];
  Hmax=  *(float *)argv[5];
  Value= (float *)argv[6];
  QuantityToAverage=(float *)argv[7];

  
  printf("N=%d\n",N);
  printf("Ngrid=%d\n",Ngrid);
  printf("Hmax=%g\n",Hmax);
  printf("DesNgb=%d\n",DesNgb);

  set_sph_kernel();

  allocate_3d();

  ngb3d_treeallocate(N, 10*N);

  set_particle_pointer();

  ngb3d_treebuild((float **)&P3d[1], N, 0,dummy,dummy);


/*
  for(i=0;i<N;i++)
  for(i=0;i<2;i++)
*/
  for(i=0;i<Ngrid;i++)
    {

          if(!(i%10000))
	    {
              printf("%d..",i); fflush(stdout);
	    }
      
	  xyz[0]=P3d[i+1]->Grid[0];
	  xyz[1]=P3d[i+1]->Grid[1];
	  xyz[2]=P3d[i+1]->Grid[2];
/*
if((i>25785) && (i<25850))
{
printf("i= %d ...  ",i);
printf("xyz= %g|%g|%g\n",xyz[0],xyz[1],xyz[2]); fflush(stdout);
}
*/

	  /*
	  xyz[0]=i*1.0;
	  xyz[1]=i*1.0;
	  xyz[2]=i*1.0;
	  */

	  h= 1.0;
	  h2=ngb3d_treefind( xyz, DesNgb ,1.04*h,&ngblist,&r2list, Hmax,&ngbfound); 
	    
	  h=sqrt(h2);

	  if(h>Hmax)
	    h=Hmax;

	  hinv = 1.0/h;
	  hinv3 = hinv*hinv*hinv;
	  
	  Density= 0.0;
	  DensityWeightedAvg= 0.0;

/*
printf("xyz= %g|%g|%g  h is currently= %g  ngbfound= %d\n",xyz[0],xyz[1],xyz[2],h,ngbfound); fflush(stdout);
*/

	  /* -----------
	      Averages
	     ----------- */
	  for(k=0;k<ngbfound;k++)
	    {
	      r=sqrt(r2list[k]);
	      
	      if(ngblist[k]==0)
		ngblist[k]= i+1;
/*
printf("r= %g  ngblist[%d]= %d\n",r,k,ngblist[k]); fflush(stdout);
*/

	      if(r<h)
		{
		  u = r/h;
		  ii = (int)(u*KERNEL_TABLE);
		  wk =hinv3*( Kernel[ii]    + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);

/*
printf("particle %d (id=%g) xyz=%g|%g|%g    vel=%g|%g|%g\n",k,P3d[ngblist[k]]->Id,P3d[ngblist[k]]->Pos[0],P3d[ngblist[k]]->Pos[1],P3d[ngblist[k]]->Pos[2],P3d[ngblist[k]]->Vel[0],P3d[ngblist[k]]->Vel[1],P3d[ngblist[k]]->Vel[2]); fflush(stdout);
*/
		  Density += Mass[ngblist[k]] * wk;
		  DensityWeightedAvg += QuantityToAverage[ngblist[k]] * Mass[ngblist[k]] * wk;
		}      
	    }

	  if(Density > 0.0)
		DensityWeightedAvg /= Density;
	  else
		DensityWeightedAvg = 0.0;

	  /* --- Cartesian Coordinates --- */
	  vv = Value + i;  *vv = DensityWeightedAvg;
    } 




  ngb3d_treefree();
 
  free_memory_3d();
  

  printf("done\n");
  return 0;
}






void set_particle_pointer(void)
{
  int i;
  float *pos;
  int Nsize;

  Nsize= N;
  if(Ngrid>N)
	Nsize= Ngrid;

  for(i=1,pos=Pos;i<=Nsize;i++)
    {
      
      P3d[i]->Pos[0] = pos[0];
      P3d[i]->Pos[1] = pos[1];
      P3d[i]->Pos[2] = pos[2];

      P3d[i]->Id = pos[3];

      P3d[i]->Grid[0] = pos[4];
      P3d[i]->Grid[1] = pos[5];
      P3d[i]->Grid[2] = pos[6];

      pos+=7;
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
	/* 3D normalization */
	  Kernel[i] = 8/PI * 2*(1-KernelRad[i])*(1-KernelRad[i])*(1-KernelRad[i]);
	  KernelDer[i] = 8/PI *( -6*(1-KernelRad[i])*(1-KernelRad[i]));
	}
    }
  

}











void allocate_3d(void)
{
  int i;
  int Nsize;

  printf("allocating memory...\n");

  Nsize= N;
  if(Ngrid>N)
	Nsize= Ngrid;

  if(Nsize>0)
    {
      if(!(P3d=malloc(Nsize*sizeof(struct particle_3d *))))
	{
	  printf("failed to allocate memory. (A)\n");
	  exit(0);
	}

      P3d--;   /* start with offset 1 */
      
      if(!(P3d[1]=malloc(Nsize*sizeof(struct particle_3d))))
	{
	  printf("failed to allocate memory. (B)\n");
	  exit(0);
	}

      for(i=2;i<=Nsize;i++)   /* initiliaze pointer table */
	P3d[i]=P3d[i-1]+1;

    }


  printf("allocating memory...done\n");
}



void free_memory_3d(void)
{
  int Nsize;

  Nsize= N;
  if(N>Ngrid)
	Nsize= Ngrid;

  if(Nsize>0)
    {
      free(P3d[1]);
      P3d++;
      free(P3d);
    }
}
