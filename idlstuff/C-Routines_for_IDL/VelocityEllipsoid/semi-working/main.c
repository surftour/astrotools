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

float *Velx,*Vely,*Velz;

int   *Id;

int   Xpixels,Ypixels;

int   DesNgb;

float Hmax; 



float *Value;






struct particle_3d 
{
  float Pos[3];
} **P3d;






int velocityellipsoid(int argc,void *argv[])
{
  int i,j,k,n,ii;
  float dummy[3],xyz[3];
  float *r2list;
  int   *ngblist;
  float r,r2,h,h2,hinv2,wk,u;
  int   ngbfound;
  float *vv;
  float *dens,*mass,masstot,m;
  float Ns,Vel_x,Vel_y,Vel_z,Disp_x,Disp_y,Disp_z;
  int NumVelocity=6;


  void allocate_3d(void);
  void set_particle_pointer(void);
  void set_sph_kernel(void);
  void free_memory_3d(void);


  if(argc!=10)
    {
      fprintf(stderr,"\n\nwrong number of arguments ! (found %d) \n\n",argc);
      exit(0);
    }


  /*******************************************/


  N=*(int *)argv[0];
  
  Pos =(float *)argv[1];
  Mass=(float *)argv[2];
printf("Mass[1]= %g\n",Mass[1]); fflush(stdout);
  Velx=(float *)argv[3];
printf("Velx[1]= %g\n",Velx[1]); fflush(stdout);
  Vely=(float *)argv[4];
printf("Vely[1]= %g\n",Vely[1]); fflush(stdout);
  Velz=(float *)argv[5];
printf("Velz[1]= %g\n",Velz[1]); fflush(stdout);
  Id=(int *)argv[6];
printf("Id[1]= %d\n",Id[1]); fflush(stdout);

  DesNgb= *(int *)argv[7];

  Hmax=  *(float *)argv[8];

  Value= (float *)argv[9];

  
  printf("N=%d\n",N);


  set_sph_kernel();

  allocate_3d();

printf(" ---- step 1 ----\n"); fflush(stdout);
printf("Mass[1]= %g\n",Mass[1]); fflush(stdout);
printf("Id[1]= %d\n",Id[1]); fflush(stdout);
printf("Velx[1]= %g\n",Velx[1]); fflush(stdout);
printf("Vely[1]= %g\n",Vely[1]); fflush(stdout);
printf("Velz[1]= %g\n",Velz[1]); fflush(stdout);

  ngb3d_treeallocate(N, 10*N);

printf(" ---- step 2 ----\n"); fflush(stdout);
printf("Mass[1]= %g\n",Mass[1]); fflush(stdout);
printf("Id[1]= %d\n",Id[1]); fflush(stdout);
printf("Velx[1]= %g\n",Velx[1]); fflush(stdout);
printf("Vely[1]= %g\n",Vely[1]); fflush(stdout);
printf("Velz[1]= %g\n",Velz[1]); fflush(stdout);

  set_particle_pointer();

printf(" ---- step 3 ----\n"); fflush(stdout);
printf("Mass[1]= %g\n",Mass[1]); fflush(stdout);
printf("Id[1]= %d\n",Id[1]); fflush(stdout);
printf("Velx[1]= %g\n",Velx[1]); fflush(stdout);
printf("Vely[1]= %g\n",Vely[1]); fflush(stdout);
printf("Velz[1]= %g\n",Velz[1]); fflush(stdout);

  ngb3d_treebuild((float **)&P3d[1], N, 0,dummy,dummy);

printf(" ---- step 4 ----\n"); fflush(stdout);
printf("Mass[1]= %g\n",Mass[1]); fflush(stdout);
printf("Id[1]= %d\n",Id[1]); fflush(stdout);
printf("Velx[1]= %g\n",Velx[1]); fflush(stdout);
printf("Vely[1]= %g\n",Vely[1]); fflush(stdout);
printf("Velz[1]= %g\n",Velz[1]); fflush(stdout);


  for(i=0;i<N;i++)
    {

          if(!(i%10000))
	    {
              printf("%d..",i); fflush(stdout);
	    }
      
/*
	  xyz[0]=P3d[i+1]->Pos[0]+0.01;
	  xyz[1]=P3d[i+1]->Pos[1]+0.01;
	  xyz[2]=P3d[i+1]->Pos[2]+0.01;
*/
	  xyz[0]=0.0;
	  xyz[1]=0.0;
	  xyz[2]=0.0;

	  h= 1.0;
	  h2=ngb3d_treefind( xyz, DesNgb ,1.04*h,&ngblist,&r2list, Hmax,&ngbfound); 
	    
	  h=sqrt(h2);

	  if(h>Hmax)
	    h=Hmax;

	  hinv2 = 1.0/(h2);
	  
	  Ns= 0.0;
	  Vel_x= Vel_y= Vel_z= 0.0;
	  Disp_x= Disp_y= Disp_z= 0.0;

printf("h is currently= %g  ngbfound= %d\n",h,ngbfound); fflush(stdout);

	  /* -----------
	      Averages
	     ----------- */
	  for(k=0;k<ngbfound;k++)
	    {
	      r=sqrt(r2list[k]);
	      
printf("r= %g  ngblist[%d]= %d\n",r,k,ngblist[k]); fflush(stdout);
	      if(r<h)
		{
		  u = r/h;
printf("u= %g\n",u); fflush(stdout);
		  ii = (int)(u*KERNEL_TABLE);
printf("ii= %g\n",ii); fflush(stdout);
		  wk =hinv2*( Kernel[ii]    + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
printf("wk= %g\n",wk); fflush(stdout);

printf("Mass= %g\n",Mass[ngblist[k]]); fflush(stdout);
printf("Id= %d\n",Id[ngblist[k]]); fflush(stdout);
printf("Velx= %g\n",Velx[ngblist[k]]); fflush(stdout);
printf("Velz= %g\n",Velz[ngblist[k]]); fflush(stdout);
printf("Vely[1]= %g\n",Vely[1]); fflush(stdout);
printf("Vely[2]= %g\n",Vely[2]); fflush(stdout);
printf("Vely[3]= %g\n",Vely[3]); fflush(stdout);
printf("Vely= %g\n",Vely[ngblist[k]]); fflush(stdout);
		  Vel_x += Velx[ngblist[k]] * Mass[ngblist[k]] * wk;
		  Vel_y += Vely[ngblist[k]] * Mass[ngblist[k]] * wk;
		  Vel_z += Velz[ngblist[k]] * Mass[ngblist[k]] * wk;
printf("Cur Vel= %g %g %g   Mass= %g\n",Velx[ngblist[k]], Vely[ngblist[k]], Velz[ngblist[k]], Mass[ngblist[k]]); fflush(stdout);

		  Ns += Mass[ngblist[k]] * wk;
		}      
	    }
	  Vel_x = Vel_x / Ns;
	  Vel_y = Vel_y / Ns;
	  Vel_z = Vel_z / Ns;

printf("Avg Vel= %g %g %g \n",Vel_x, Vel_y, Vel_z); fflush(stdout);


          /* -------------
              Dispersions
             ------------- */
          for(k=0;k<ngbfound;k++)
            {
              r=sqrt(r2list[k]);
              
              if(r<h)
                {
                  u = r/h;
                  ii = (int)(u*KERNEL_TABLE);
                  wk =hinv2*( Kernel[ii]    + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);

                  Disp_x += (Velx[ngblist[k]] - Vel_x)*(Velx[ngblist[k]] - Vel_x) * Mass[ngblist[k]] * wk;
                  Disp_y += (Vely[ngblist[k]] - Vel_y)*(Vely[ngblist[k]] - Vel_y) * Mass[ngblist[k]] * wk;
                  Disp_z += (Velz[ngblist[k]] - Vel_z)*(Velz[ngblist[k]] - Vel_z) * Mass[ngblist[k]] * wk;
                }      
            }
	  Disp_x = Disp_x / Ns;
	  Disp_y = Disp_y / Ns;
	  Disp_z = Disp_z / Ns;



	  /*  NumVelocity = number of velocity fiels passed back in array (6 currently)
          */

	  /* --- Cartesian Coordinates --- */
	  vv = Value + i*NumVelocity + 0;  *vv = Vel_x;
	  vv = Value + i*NumVelocity + 1;  *vv = Disp_x;

	  vv = Value + i*NumVelocity + 2;  *vv = Vel_y;
	  vv = Value + i*NumVelocity + 3;  *vv = Disp_y;

	  vv = Value + i*NumVelocity + 4;  *vv = Vel_z;
	  vv = Value + i*NumVelocity + 5;  *vv = Disp_z;

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
