#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

/*
 *  This function initializes the cpu-time counters to 0.
 *  Then begrun() is called, which sets up the simulation
 *  either from IC's or from restart files.
 *  Finally, run() is started, the main simulation loop,
 *  which iterates over the timesteps.
 */
//void project(void);
//void allocate_3d(void);
//void free_memory_3d(void);

int   N;

float *Pos,*Mass,*Hsml,*Density;

float Xmin,Ymin,Xmax,Ymax;

int   Xpixels,Ypixels;

int   DesNgb;

int   Axis1, Axis2; 




float *Value;

struct particle_3d
{
  float Pos[3];
} **P3d;


int smooth(int argc, void *argv[])
{
	FILE *output;
	double gxyz[3];
	double gdxyz[3];
	int nxyz[2];// grid pixel dimensions
	int ixyz[3];
	int ngridp;

	float *img;
	float *x,*y;
	double r,u,wk;
	int ii;
	int ihux;
	int ihuy;
	int ihlx;
	int ihly;
	int ipx;
	int ipy;

	double dx,dy;
	double hinv,hinv2,hinv3;
	double pdxy[2];
	double minh;
	double *grid;

	int i,j,k,m;
	int RestartFlag =0;
	char fname[200];
	float *vv;
	struct global_data_all_processes_all;


	N=*(int *)argv[0];
	NumPart = N;
	Pos =(float*)argv[1];
	Mass=(float *)argv[2];
	Density=(float *)argv[3];
	Hsml=(float *)argv[4];

	Xmin=*(float *)argv[5];
	Xmax=*(float *)argv[6];
	Ymin=*(float *)argv[7];
	Ymax=*(float *)argv[8];

	Xpixels=*(int *)argv[9];
	Ypixels=*(int *)argv[10];

	//DesNgb= *(int *)argv[11];

	Axis1= *(int *)argv[11];
	Axis2= *(int *)argv[12];
  

	Value= (float *)argv[13];


	gxyz[0] = (Xmax-Xmin)/2.0 + Xmin;
	gxyz[1] = (Ymax-Ymin)/2.0 + Ymin; //chunk y center
	gdxyz[0] = Xmax-Xmin; //chunk x width 
	gdxyz[1] = Ymax-Ymin; //chunk y width
	nxyz[0]  = Ypixels;
	nxyz[1]  = Xpixels;
	pdxy[0] = 0.0;
	pdxy[1] = 0.0;
	ngridp = nxyz[0]*nxyz[1];


	//printf("N %d Numpart %d Xmin %e xmax %e ymin %e ymax %e nxyz[0] %d nxyz[1] %d Axis1 %d Axis2 %d\n",N,NumPart,Xmin,Xmax,Ymin,Ymax,Xpixels,Ypixels,Axis1,Axis2);

	All.PartAllocFactor = 1.05;
	All.TreeAllocFactor = 50.0;
	All.MinGasHsml = 0.1;
	All.DesNumNgb = DesNgb;
	N_gas = NumPart;
	All.TotNumPart = NumPart;
	All.TotN_gas= N_gas;
	All.MaxPart =  All.PartAllocFactor *  All.TotNumPart;
	All.MaxPartSph=  All.PartAllocFactor * All.TotN_gas;
	printf("producing image.... ");
	fflush(stdout);
	//printf("Numpart=%d\n", NumPart);
	//printf("TNP=%d\n",All.TotNumPart);
	//printf("MP=%d\n",All.MaxPart);
	//printf("MSP=%d\n",All.MaxPartSph);


	allocate_3d();
	allocate_memory();
	project();
	//printf("....initializing particles....");
	//fflush(stdout);
	for(i=1;i<=NumPart;i++)
	{
		for(j=0;j<3;j++)
		{
			P[i].Pos[j]    = P3d[i]->Pos[j];
			P[i].PosPred[j]=P[i].Pos[j];
			P[i].VelPred[j]=P[i].Vel[j];
			P[i].Accel[j]=0;
		}

		P[i].Mass = Mass[i];
		P[i].Type = 0;
		P[i].ID = i;

		SphP[i].EgySpecPred  = SphP[i].EgySpec = 0;
		SphP[i].DtEgySpec=0;
		SphP[i].DtDensity = 0;
		SphP[i].Density=Density[i];
		SphP[i].DtHsml = 0;
		SphP[i].Hsml=Hsml[i];
	/*	if(!(i%5000))
		{
			printf("...%d...",i);
			fflush(stdout);
		}*/

	}

	printf("particles initialized ....\n");
	fflush(stdout);
	//force_treeallocate(All.TreeAllocFactor*All.MaxPart, All.MaxPart);
	//ngb_treeallocate(MAX_NGB);


	//setup_smoothinglengths(All.DesNumNgb);


	//Set up image array
	img = (float *) malloc(nxyz[0]*nxyz[1]*sizeof(float));
	x   = (float *) malloc(nxyz[0]*sizeof(float));
	y   = (float *) malloc(nxyz[1]*sizeof(float));
	pdxy[0] = gdxyz[0]/((double) nxyz[0]);
	pdxy[1] = gdxyz[1]/((double) nxyz[1]);

	//printf("nxyz[0] %d nxyz[1] %d\n",nxyz[0],nxyz[1]);
	for(i=0;i<nxyz[0];i++)
	{
		x[i] = (gxyz[0] - gdxyz[0]/2.0) + ((double) i)*pdxy[0] + 0.5*pdxy[0];
		for(j=0;j<nxyz[1];j++)
		{
			img[i*nxyz[1]+j] = 0.0;
			if(i==0)
				y[j] = (gxyz[1] - gdxyz[1]/2.0) + ((double) j)*pdxy[1] + 0.5*pdxy[1];
		}
		/*if(!(i%100))
		{
			printf("...%d...",i);
			fflush(stdout);
		}*/
	}
	printf("array initialized....\n");
	fflush(stdout);

	//set the kernel
  	set_sph_kernel();
	//printf("Calculating grid densities....\n");
	//printf("Part\n");
	for(k=1;k<=N_gas;k++)
	{
		dx = P[k].Pos[Axis1]-gxyz[0];
		dy = P[k].Pos[Axis2]-gxyz[1];
		if(dx==0.0)
			P[k].Pos[Axis1]+=1.0e-6;
		if(dy==0.0)
			P[k].Pos[Axis2]+=1.0e-6;
		if( fabs(dx)<(gdxyz[1]/2.0))
		if( fabs(dy)<(gdxyz[0]/2.0))
		{
			hinv = 1.0/SphP[k].Hsml;
			hinv3 = hinv*hinv*hinv;
			hinv2 = hinv*hinv;
			ihux = (int)(SphP[k].Hsml/pdxy[0]) + 1;
			ihuy = (int)(SphP[k].Hsml/pdxy[1]) + 1;
			ipx = (int)((dx/pdxy[0]) + nxyz[0]/2);
			ipy = (int)((dy/pdxy[1]) + nxyz[1]/2);

			if((ipx-ihux-1)<0)
			{
				ihlx = ipx-1;
			}else{
				ihlx = ihux;
			}
			if((ipy-ihuy-1)<0)
			{
				ihly = ipy-1;
			}else{
				ihly = ihuy;
			}
			if((ipx+ihux+1)>nxyz[0])
			{
				ihux = nxyz[0]-1-ipx;
			}
			if((ipy+ihuy+1)>nxyz[1])
			{
				ihuy = nxyz[1]-1-ipy;
			}

			ii = 0;
	
			for(i=(ipx-ihlx-1);i<(ipx+ihux+1);i++)
				for(j=(ipy-ihly-1);j<(ipy+ihuy+1);j++)
				{
					r = sqrt((P[k].Pos[Axis1]-x[i])*(P[k].Pos[Axis1]-x[i]) + (P[k].Pos[Axis2]-y[j])*(P[k].Pos[Axis2]-y[j]));
					if(r < SphP[k].Hsml/2.0)
					{
						u = r*hinv;
						ii = (int)(u*KERNEL_TABLE);
						//wk =hinv3*( Kernel[ii]  + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
						//wk =hinv2*( Kernel[ii]  + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
						//img[i*nxyz[1]+j]+= (40.0/56.0)*P[k].Mass*wk;
						wk = 2.0*( A_kernel(r,SphP[k].Hsml) + B_kernel(r,SphP[k].Hsml) );
					}else{
						if(r < SphP[k].Hsml)
						{
							u = r*hinv;
							ii = (int)(u*KERNEL_TABLE);
							wk = 2.0*F_kernel(r,SphP[k].Hsml);
						}else{
							wk = 0.0;
						}
					}
					img[i*nxyz[1]+j]+= P[k].Mass*wk;
					if(wk*pow(SphP[k].Hsml,2.0) > 2.0)
						printf("wk > 2.0, wk = %e\n",wk*pow(SphP[k].Hsml,2.0));
				}
		}
		if(!(k%10000))
		{
			//printf("....%d..%e..%e..%e.....",k,r,SphP[k].Hsml,wk);
			printf("d... ",k);
			fflush(stdout);
		}
		
		//	printf("%d hsml %e ipx %d ipy %d ihux %d ihlx %d ihuy %d ihly %d\n",k,SphP[k].Hsml,ipx,ipy,ihux,ihlx,ihuy,ihly);
	}

	printf("\nsaving image ....\t");
	fflush(stdout);
	for(i=0;i<nxyz[0];i++)
		for(j=0;j<nxyz[1];j++)
		{
			vv = Value+j*nxyz[0] + i;
			*vv = img[i*nxyz[1]+j];
		}

	printf("done!\n");
	fflush(stdout);

	free(img);
	free(x);
	free(y);
	free_memory();
	free_memory_3d();

	return 0;
}
//void allocate_3d(void)
allocate_3d()
{
  int i;

  //printf("allocating memory...\n");

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

  //printf("allocating memory...done\n");
}

project()
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

free_memory_3d()
{
  if(N>0)
    {
      free(P3d[1]);
      P3d++;
      free(P3d);
    }
}




double kernel2d(double r, double h)
{

}
double A_kernel(double r, double h)
{
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double d = 0.0;
	double y = 0.0;

	double Az = 0.875352;

	if(r/h < 1.0e-6)
		return Az;

	y =  sqrt(h*h/4.0 - r*r);

	a =  (1.0 / (2.0*pow(h,6.0)*PI));
	b =  11.0*h*h*h*y;
	c = -46.0*h*r*r*y;
	d =  36.0*r*r*r*r*( log(h*(h+2*y)) - log(2*h*r) );

	return a*(b+c+d);
}
double B_kernel(double r, double h)
{
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double d = 0.0;
	double e = 0.0;
	double f = 0.0;
	double g = 0.0;
	double y = 0.0;
	double x = 0.0;

	double Bz = 0.0795775;

	if(r/h < 1.0e-6)
		return Bz;

	y =  sqrt(h*h/4.0 - r*r);
	x =  sqrt(h*h     - r*r);

	a =  (1.0 / (2.0*pow(h,6.0)*PI));
	b = -15.0*h*h*h*y;
	c = -58.0*h*r*r*y;
	d =   8.0*h*h*h*x;
	e =  52.0*h*r*r*x;
	f =  48.0*h*h*r*r*( log(h*(h+2*y)) - log(2*h*(h+x)) );
	g =  12.0*r*r*r*r*( log(h*(h+2*y)) - log(2*h*(h+x)) );


	return a*(b+c+d+e+f+g);
}
double F_kernel(double r, double h)
{
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double d = 0.0;
	double e = 0.0;
	double y = 0.0;
	double x = 0.0;

	if(r>=h) return 0.0;

	y =  sqrt(h*h/4.0 - r*r);
	x =  sqrt(h*h     - r*r);

	a =  (2.0 / (pow(h,6.0)*PI));
	b =   2.0*h*h*h*x;
	c =  13.0*h*r*r*x;
	d =  12.0*h*h*r*r*( log(2*h*r) - log(2*h*(h+x)) );
	e =   3.0*r*r*r*r*( log(2*h*r) - log(2*h*(h+x)) );

	return a*(b+c+d+e);
}
