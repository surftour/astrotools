#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "globvars.h"



void move_galaxy(struct galaxy_data *g, double x, double y, double z, double vx, double vy, double vz);


void manual_move_galaxies(struct galaxy_data *g1, struct galaxy_data *g2)
{
  double x, y, z, vx, vy, vz;
  double mu1, mu2;
  double xcom, ycom, zcom;
  double vxcom, vycom, vzcom;

  mu1 = g1->Mtot / (g1->Mtot + g2->Mtot);
  mu2 = g2->Mtot / (g1->Mtot + g2->Mtot);
  printf("g1 Mass= %g   (mu= %g)\n", g1->Mtot, mu1);
  printf("g2 Mass= %g   (mu= %g)\n", g2->Mtot, mu2);
  printf("Total Mass= %g\n", (g1->Mtot + g2->Mtot));

  xcom= mu1 * g1_x_start + mu2 * g2_x_start;
  ycom= mu1 * g1_y_start + mu2 * g2_y_start;
  zcom= mu1 * g1_z_start + mu2 * g2_z_start;
  printf("computed c.o.m.=  %g,  %g,  %g\n", xcom, ycom, zcom);
 
  vxcom= mu1 * g1_vx_start + mu2 * g2_vx_start;
  vycom= mu1 * g1_vy_start + mu2 * g2_vy_start;
  vzcom= mu1 * g1_vz_start + mu2 * g2_vz_start;
  printf("computed c.o.v.=  %g,  %g,  %g\n", vxcom, vycom, vzcom);

  printf("\n");
 
  /*  ----  Galaxy 1  ---- */
  x= g1_x_start - xcom/(mu1 + mu2);
  y= g1_y_start - ycom/(mu1 + mu2);
  z= g1_z_start - zcom/(mu1 + mu2);

  vx= g1_vx_start - vxcom/(mu1 + mu2);
  vy= g1_vy_start - vycom/(mu1 + mu2);
  vz= g1_vz_start - vzcom/(mu1 + mu2);

  printf("moving galaxy 1 to c.o.m and c.o.v. frame...\n");
  move_galaxy(g1, x, y, z, vx, vy, vz);

  /*  ----  Galaxy 2  ---- */
  x= g2_x_start - xcom/(mu1 + mu2);
  y= g2_y_start - ycom/(mu1 + mu2);
  z= g2_z_start - zcom/(mu1 + mu2);

  vx= g2_vx_start - vxcom/(mu1 + mu2);
  vy= g2_vy_start - vycom/(mu1 + mu2);
  vz= g2_vz_start - vzcom/(mu1 + mu2);

  printf("moving galaxy 2 to c.o.m and c.o.v. frame...\n");
  move_galaxy(g2, x, y, z, vx, vy, vz);

}


void move_galaxy(struct galaxy_data *g, double x, double y, double z, double vx, double vy, double vz)
{
  double xyz[3];

  xyz[0]= x;
  xyz[1]= y;
  xyz[2]= z;
  printf("galaxy moved x,y,z= %g, %g, %g\n", xyz[0],xyz[1],xyz[2]);
  translate(g, &xyz[0]);

  xyz[0] = vx;
  xyz[1] = vy;
  xyz[2] = vz;
  printf("galaxy velocity now = %g, %g, %g\n", xyz[0],xyz[1],xyz[2]);
  vel_translate(g, &xyz[0]);

  printf("\n");
}



