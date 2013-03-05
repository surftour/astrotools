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
  float Vel[3];
  float Id;
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
  float dumb_x, dumb_y, dumb_z;
  int NumVelocity=6;


  void allocate_3d(void);
  void set_particle_pointer(void);
  void set_sph_kernel(void);
  void free_memory_3d(void);


  if(argc!=6)
    {
      fprintf(stderr,"\n\nwrong number of arguments ! (found %d) \n\n",argc);
      exit(0);
    }


  /*******************************************/


  N=*(int *)argv[0];
  
  Pos =(float *)argv[1];
  Mass=(float *)argv[2];
  DesNgb= *(int *)argv[3];
  Hmax=  *(float *)argv[4];
  Value= (float *)argv[5];

  
  printf("N=%d\n",N);
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
  for(i=0;i<N;i++)
    {

          if(!(i%10000))
	    {
              printf("%d..",i); fflush(stdout);
	    }
      
	  xyz[0]=P3d[i+1]->Pos[0]+0.01;
	  xyz[1]=P3d[i+1]->Pos[1]+0.01;
	  xyz[2]=P3d[i+1]->Pos[2]+0.01;

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

	  hinv2 = 1.0/(h2);
	  
	  Ns= 0.0;
	  Vel_x= Vel_y= Vel_z= 0.0;
	  Disp_x= Disp_y= Disp_z= 0.0;
	  dumb_x= dumb_y= dumb_z= 0.0;     /* this was for checking the weighted versus straight-up calc. */

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
		  wk =hinv2*( Kernel[ii]    + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);

/*
printf("particle %d (id=%g) xyz=%g|%g|%g    vel=%g|%g|%g\n",k,P3d[ngblist[k]]->Id,P3d[ngblist[k]]->Pos[0],P3d[ngblist[k]]->Pos[1],P3d[ngblist[k]]->Pos[2],P3d[ngblist[k]]->Vel[0],P3d[ngblist[k]]->Vel[1],P3d[ngblist[k]]->Vel[2]); fflush(stdout);
*/
		  x= P3d[ngblist[k]]->Pos[0];
		  y= P3d[ngblist[k]]->Pos[1];
		  z= P3d[ngblist[k]]->Pos[2];
		  r_xyz= sqrt(x*x + y*y + z*z);
		  r_xy= sqrt(x*x + y*y);
		  r_t= sqrt(x*x*z*z + y*y*z*z + pow(r_xyz,4));

		  vx= P3d[ngblist[k]]->Vel[0];
		  vy= P3d[ngblist[k]]->Vel[1];
		  vz= P3d[ngblist[k]]->Vel[2];

		  v_tot= sqrt(vx*vx + vy*vy + vz*vz);
		  v_rad= (vx*x + vy*y + vz*z)/r_xyz;
		  v_theta= (-vx*x*z - vy*y*z + z*r_xy*r_xy)/r_t;
		  v_phi= (vx*y - vy*x)/r_xyz;
		  v_tan= sqrt(v_tot - v_rad*v_rad);

		  Vel_tot= += v_tot * Mass[ngblist[k]] * wk;
		  Vel_rad += v_rad * Mass[ngblist[k]] * wk;
		  Vel_theta += v_theta * Mass[ngblist[k]] * wk;
		  Vel_phi += v_phi * Mass[ngblist[k]] * wk;
		  Vel_tan += v_tan * Mass[ngblist[k]] * wk;

		  Ns += Mass[ngblist[k]] * wk;

		  /*
		  dumb_x += P3d[ngblist[k]]->Vel[0];
		  dumb_y += P3d[ngblist[k]]->Vel[1];
		  dumb_z += P3d[ngblist[k]]->Vel[2];
		  */

		}      
	    }
	  Vel_tot = Vel_tot / Ns;
	  Vel_rad = Vel_rad / Ns;
	  Vel_theta = Vel_theta / Ns;
	  Vel_phi = Vel_phi / Ns;
	  Vel_tan = Vel_tan / Ns;

/*
printf("Avg Vel= %g %g %g \n",Vel_x, Vel_y, Vel_z); fflush(stdout);
*/

	  /*
	  dumb_x = dumb_x / ngbfound;
	  dumb_y = dumb_y / ngbfound;
	  dumb_z = dumb_z / ngbfound;
	  printf("Dumb Avg Vel= %g %g %g \n",dumb_x, dumb_y, dumb_z); fflush(stdout);
	  */


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

                  Disp_x += (P3d[ngblist[k]]->Vel[0] - Vel_x)*(P3d[ngblist[k]]->Vel[0] - Vel_x) * Mass[ngblist[k]] * wk;
                  Disp_y += (P3d[ngblist[k]]->Vel[1] - Vel_y)*(P3d[ngblist[k]]->Vel[1] - Vel_y) * Mass[ngblist[k]] * wk;
                  Disp_z += (P3d[ngblist[k]]->Vel[2] - Vel_z)*(P3d[ngblist[k]]->Vel[2] - Vel_z) * Mass[ngblist[k]] * wk;
                }      
            }
	  Disp_tot = Disp_tot / Ns;
	  Disp_rad = Disp_rad / Ns;
	  Disp_theta = Disp_theta / Ns;
	  Disp_phi = Disp_phi / Ns;
	  Disp_tan = Disp_tan / Ns;

/*
printf("Dispersion= %g %g %g \n",Disp_x, Disp_y, Disp_z); fflush(stdout);
*/


	  /*  NumVelocity = number of velocity fields passed back in array (6 currently)
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

      P3d[i]->Vel[0] = pos[3];
      P3d[i]->Vel[1] = pos[4];
      P3d[i]->Vel[2] = pos[5];

      P3d[i]->Id = pos[6];

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
