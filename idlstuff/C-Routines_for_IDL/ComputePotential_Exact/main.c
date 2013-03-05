#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>




#define  PI               3.1415927
#define KERN_LEN   10000  /* length of arrays for tabulated functions  */



int   N;

float *Pos,*Mass;



float *Potential;


static double  knlrad  [KERN_LEN+1],   /* tabulated functions for force/potential softening */
               knlpot  [KERN_LEN+1];






struct particle_data
{
  float Pos[3];
} **P;



void setkernel(void);
double dmax(double x,double y);






int calc_potential(int argc,void *argv[])
{
  int i,j;
  float r,r_x,r_y,r_z;
  float phi_ij= 0.0;

  double epsilon, h, h_inv;
  int ii;
  double u, wp, ff;

  if(argc!=4)
    {
      fprintf(stderr,"\n\nwrong number of arguments ! (found %d) \n\n",argc);
      exit(0);
    }


  /*******************************************/


  N=*(int *)argv[0];
  
  Pos =(float *)argv[1];
  Mass=(float *)argv[2];


  Potential= (float *)argv[3];

  
  printf("N=%d\n",N);



  setkernel();

  allocate();

  extract();


  for(i=1;i<=N;i++)
    {
      if((i%200)==0)
        printf("...%d\n",i);
      
      for(j=1;j<=N;j++)
	{
/*	  
if((i<10) && (j<10)) { printf("i= %d (%g,%g,%g - m=%g) j=%d (%g,%g,%g - m=%g)\n",i,P[i]->Pos[0],P[i]->Pos[1],P[i]->Pos[2],Mass[i],j,P[j]->Pos[0],P[j]->Pos[1],P[j]->Pos[2],Mass[j]); }
*/
	  r_x= P[i]->Pos[0] - P[j]->Pos[0];
	  r_y= P[i]->Pos[1] - P[j]->Pos[1];
	  r_z= P[i]->Pos[2] - P[j]->Pos[2];

	  r= sqrt(r_x*r_x + r_y*r_y + r_z*r_z);

/*	  epsilon= dmax(Smoothing[i],Smoothing[j]);   */
	  epsilon= 0.1;

	  /* these need to be edited based upon the specific snapshot
	   * looking at.  I should change this to pass the smoothing
	   * lengths, but late */
	  if((i>500) && (i<=1500)) { epsilon= 0.4; }
	  if((j>500) && (j<=1500)) { epsilon= 0.4; }
  	  h = 2.8*epsilon;
	  h_inv=1/h;

	  if(r>0.0)
	    {
		  u=r*h_inv;
/*
if((i>400) && (i<410) && (j>400) && (j<450)) { printf("i= %d  j=%d  r=%g  ep=%g  m_i=%g  m_j%g\n",i,j,r,epsilon,Mass[i],Mass[j]); }
*/

	          if(u>=1)
	            {
	              phi_ij = -1.0*Mass[j]/r;
	            }
	          else
	            {
	              ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
	              wp=knlpot[ii]+(knlpot[ii+1]-knlpot[ii])*ff;
/*
printf("it's in here: i,j=%d,%d  wp=%g, h_inv=%g,  u=%g,  r=%g\n",i,j,wp,h_inv,u,r);
*/

	              phi_ij = +1.0 * Mass[j]*h_inv*wp;
	            }
	    }
	  else
		phi_ij= 0.0;


	/* this doesn't have the offset of 1, hence we subtract 1 */
	  Potential[i-1] += phi_ij;

	}
    } 


 
  free_memory();
  

  printf("Done with potential calculation.\n");
  return 0;
}








extract()
{
  int i;
  float *pos;

  printf("extracting...\t");
  for(i=1,pos=Pos;i<=N;i++)
    {
      
      P[i]->Pos[0] = pos[0];
      P[i]->Pos[1] = pos[1];
      P[i]->Pos[2] = pos[2];

      pos+=3;
    }

  printf("          ...done\n"); fflush(stdout);
}












allocate()
{
  int i;

  printf("allocating memory...\t");

  if(N>0)
    {
      if(!(P=malloc(N*sizeof(struct particle_data *))))
	{
	  printf("failed to allocate memory. (A)\n");
	  exit(0);
	}

      P--;   /* start with offset 1 */
      
      if(!(P[1]=malloc(N*sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory. (B)\n");
	  exit(0);
	}

      for(i=2;i<=N;i++)   /* initiliaze pointer table */
	P[i]=P[i-1]+1;

    }

  printf("   ...done\n");
}



free_memory()
{
  if(N>0)
    {
      free(P[1]);
      P++;
      free(P);
    }
}








/* this routine tabulates the kernel used for the softening
 * of the force and potential of monopole and quadrupole
 * moments
 */
void setkernel(void)
{
  int i;
  double u;

  for(i=0;i<=KERN_LEN;i++)
    {
      u=((double)i)/KERN_LEN;

      knlrad[i] = u;

      if(u<=0.5)
        {
          knlpot[i]=16.0/3*pow(u,2)-48.0/5*pow(u,4)+32.0/5*pow(u,5)-14.0/5;
        }
      else
        {
          knlpot[i]=1.0/15/u +32.0/3*pow(u,2)-16.0*pow(u,3)+48.0/5*pow(u,4)-32.0/15*pow(u,5)-16.0/5;
        }
    }
}







/* returns the maximum of two double
 */
double dmax(double x,double y)
{
  if(x>y)
    return x;
  else
    return y;
}

