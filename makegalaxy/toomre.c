#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"


void plot_toomre_stability(FILE *fd);

void plot_circular_speeds(FILE *fd);



double epicyclic_kappa2(double R)
{
  double dR, dphi, dphi_;

  if(R > 0)
    {
      dR = R * 0.01;

      dphi = comp_Dphi_R(R, 0);

      /*comp_Dphi_R_halo(R, 0) + comp_Dphi_R_disk(R, 0);
       */

      dphi_ = comp_Dphi_R(R + dR, 0);

      /*
         comp_Dphi_R_halo(R + dR, 0) + comp_Dphi_R_disk(R + dR, 0);
       */

      return 3 * dphi / R + (dphi_ - dphi) / dR;

    }
  else
    return 0;
}


void write_toomre(void)
{
  FILE *fd;
  char buf[100];

  strcpy(buf, "");
  printf("OutputFile= %s\n",OutputFile);
  if(strcmp(&OutputFile[strlen(OutputFile)-5],".hdf5") == 0)
    strncpy(buf, OutputFile, strlen(OutputFile)-5);
  else
    strcpy(buf, OutputFile);
  strcat(buf, ".vc");
  printf("writing header, vc, toomre Q to file: %s\n",buf);
  if((fd = fopen(buf, "w")))
    {
      printf("writing circular velocity curve + Toomre's Q\n");
      plot_circular_speeds(fd);
      plot_toomre_stability(fd);
      fclose(fd);
    }
  else
    {
      fprintf(stderr, "Can't open file '%s'.\n", buf);
      exit(0);
    }
  printf("done.\n");
}



void plot_toomre_stability(FILE * fd)
{

// GB modified to include gas
  double *Q;
  double *Qg;
  int i;
  double Sigma0;
  double Sigma0g;
  int count;

  for(i = 0, count = 0; i <= RSIZE; i++)
    if(list_R[i] <= 6 * H)
      count++;


  Sigma0 = (M_DISK) / (2 * PI * H * H);
  Sigma0g = (M_GAS) / (2 * PI * H * H);
  Q = dvector(0, RSIZE);
  Qg = dvector(0, RSIZE);


  for(i = 0; i < count; i++)
    Q[i] = RadialDispersionFactor * sqrt(VelDispRz_disk[i][0]) * sqrt(epi_kappa2[i]) / (3.36 * G * Sigma0 * exp(-list_R[i] / H));


/* GB:  assume all gas at U4  (q=0) 
  vsound =  sqrt(gamma * uinternal)  uinternal has units of pressure/density
  not sure how to include the extended gas disk in this .. */  

 for(i=0; i<count; i++)
    Qg[i] = sqrt(epi_kappa2[i])*sqrt(GAMMA*U4*pow(RhoGas[i][0]/PhysDensThresh,-0.5))/(PI* G * Sigma0g * exp(-list_R[i]/H));

  fprintf(fd, "\n%d\n", count);

  for(i = 0; i < count; i++)
    fprintf(fd, "%g\n", list_R[i]);

  for(i = 0; i < count; i++)
    fprintf(fd, "%g\n", Q[i]);

  for(i=0; i < count; i++)
    fprintf(fd, "%g\n", Qg[i]);



// end GB
}







void plot_circular_speeds(FILE * fd)
{
  int i;
  double R;
  double RMAX;
  double vc2;

#define POINTS 1000

  RMAX = R200;

  fprintf(fd, "%d\n", POINTS);

  for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      fprintf(fd, "%f\n", R);
    }
  for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      fprintf(fd, "%f\n", vc2 = R * comp_Dphi_R(R, 0));
    }

  for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      fprintf(fd, "%f\n", vc2 = R * comp_Dphi_R_disk_tree(R, 0));
    }
  for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      fprintf(fd, "%f\n", vc2 = R * comp_Dphi_R_halo(R, 0));
    }

  for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      fprintf(fd, "%f\n", vc2 = R * comp_Dphi_R_bulge(R, 0));
    }

#undef POINTS
}
