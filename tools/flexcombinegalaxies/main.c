/*** code to combine two galaxies such that they collide
     on a parabolic encounter ***/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "globvars.h"




int main(int argc, char *argv[])
{
  char parameter_filename[100];
  char p_gal1[120], p_gal2[120], cmd[120];  /* , temp[120]; */
  FILE *fd;

  if(argc != 2)
    {
      printf("\nwrong number of argument(s).  Specify a parameterfile.\n\n");
      exit(0);
    }

#ifdef ECC_NONZERO
#ifdef ORBITCORRECTION
  printf(" can't have ECC_NONZERO and ORBITCORRECTION turned on at the same time.\n");
  printf(" please fix.\n");
  exit(0);
#endif
#endif

#ifdef ORBITCORRECTION
  if(strstr(combinetype,"manual"))
    {
      printf("WARNING: can't do ORBITCORRECTION and manually set the velocities ....\n");
      exit(0);
    }
#endif

  strcpy(parameter_filename, argv[1]);

  read_parameterfile(parameter_filename);



  /* Do we have the parameter files */
  /* ------------------------------ */
  strcpy(p_gal1,"");
  strcpy(p_gal1, gal_fname1);    /* old method: *.dat.parameters */
  strcat(p_gal1, ".parameters");
  if((fd=fopen(p_gal1,"r")))
        fclose(fd);
  else
  {
        fprintf(stderr,"\nProblem: can't find 1st galaxy IC file '%s'.\n",p_gal1);
        exit(0);
  }
  
  strcpy(p_gal2,"");
  strcpy(p_gal2, gal_fname2);
  strcat(p_gal2, ".parameters");
  if((fd=fopen(p_gal2,"r")))
        fclose(fd);
  else
  {
        fprintf(stderr,"\nProblem: can't find 2nd galaxy IC file '%s'.\n",p_gal2);
        exit(0);
  }



  /* Now do the work */
  /* --------------- */
  load_particles(gal_fname1, &gal[0], &header[0]);
  load_particles(gal_fname2, &gal[1], &header[1]);

  turn_galaxy(&gal[0], theta1, phi1);
  turn_galaxy(&gal[1], theta2, phi2);


  if(strstr(combinetype,"kepler"))
    kepler_move_galaxies(&gal[0], &gal[1]);

  if(strstr(combinetype,"manual"))
    manual_move_galaxies(&gal[0], &gal[1]);


  /* output file */
  if(!(strstr(gal_output,"hdf5")) || (strlen(strstr(gal_output,"hdf5")) > 4))
    strcat(gal_output, ".hdf5");


  save_combined(gal_output, &gal[0], &gal[1]);



  /* Write a parameter file */
  /* ---------------------- */
  strcat(gal_output,".parameters");
  if(!(fd=fopen(gal_output,"w")))
  {
        fprintf(stderr,"\nCan't find file '%s'.\n",gal_output);
        exit(0);
  }

  fprintf(fd,"gal_out \t%s\n",gal_output);
  fprintf(fd,"gal1    \t%s\n",gal_fname1);
  fprintf(fd,"theta1  \t%g\n",theta1);
  fprintf(fd,"phi1    \t%g\n",phi1);
  fprintf(fd,"gal2    \t%s\n",gal_fname2);
  fprintf(fd,"theta2  \t%g\n",theta2);
  fprintf(fd,"phi2    \t%g\n",phi2);
  fprintf(fd,"cmbtype \t%s\n",combinetype);
  fprintf(fd,"rperi   \t%g\n",rperi);
  fprintf(fd,"rstart  \t%g\n",rstart);
  fprintf(fd,"ecc     \t%g\n",ecc);
  fprintf(fd,"g1xstrt \t%g\n",g1_x_start);
  fprintf(fd,"g1ystrt \t%g\n",g1_y_start);
  fprintf(fd,"g1zstrt \t%g\n",g1_z_start);
  fprintf(fd,"g1vxstrt\t%g\n",g1_vx_start);
  fprintf(fd,"g1vystrt\t%g\n",g1_vy_start);
  fprintf(fd,"g1vzstrt\t%g\n",g1_vz_start);
  fprintf(fd,"g2xstrt \t%g\n",g2_x_start);
  fprintf(fd,"g2ystrt \t%g\n",g2_y_start);
  fprintf(fd,"g2zstrt \t%g\n",g2_z_start);
  fprintf(fd,"g2vxstrt\t%g\n",g2_vx_start);
  fprintf(fd,"g2vystrt\t%g\n",g2_vy_start);
  fprintf(fd,"g2vzstrt\t%g\n",g2_vz_start);
  fclose(fd);

  strcpy(cmd, "");
  strcpy(cmd,"echo GAL1 >> ");
  strcat(cmd,gal_output);
  printf("cmd= %s\n",cmd);
  system(cmd);
  
  strcpy(cmd, "");
  strcpy(cmd,"head -n 35 ");
  strcat(cmd,p_gal1);
  strcat(cmd," >> ");
  strcat(cmd,gal_output);
  printf("cmd= %s\n",cmd);
  system(cmd);

  strcpy(cmd, "");
  strcpy(cmd,"echo GAL2 >> ");
  strcat(cmd,gal_output);
  printf("cmd= %s\n",cmd);
  system(cmd);

  strcpy(cmd, "");
  strcpy(cmd,"head -n 35 ");
  strcat(cmd,p_gal2);
  strcat(cmd," >> ");
  strcat(cmd,gal_output);
  printf("cmd= %s\n",cmd);
  system(cmd);



  return 0;
}

