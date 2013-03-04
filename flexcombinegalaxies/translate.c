#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "globvars.h"



void translate(struct galaxy_data *g, double vec[3])
{
  int i, j;


  for(i = 1; i <= g->Ntot; i++)
    for(j = 0; j <= 2; j++)
      g->pos[i][j + 1] += vec[j];
}


void vel_translate(struct galaxy_data *g, double vec[3])
{
  int i, j;


  for(i = 1; i <= g->Ntot; i++)
    for(j = 0; j <= 2; j++)
      g->vel[i][j + 1] += vec[j];
}
