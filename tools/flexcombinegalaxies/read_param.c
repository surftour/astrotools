#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "globvars.h"


void read_parameterfile(char *fname)
{
#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  char buf[200], buf1[200], buf2[200], buf3[200];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0;

  /* read parameter file on all processes for simplicty */

  nt = 0;

  strcpy(tag[nt], "gal_fname1");
  addr[nt] = &gal_fname1;
  id[nt++] = STRING;

  strcpy(tag[nt], "gal_fname2");
  addr[nt] = &gal_fname2;
  id[nt++] = STRING;

  strcpy(tag[nt], "gal_output");
  addr[nt] = &gal_output;
  id[nt++] = STRING;

  strcpy(tag[nt], "combinetype");
  addr[nt] = &combinetype;
  id[nt++] = STRING;

  strcpy(tag[nt], "theta1");
  addr[nt] = &theta1;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "phi1");
  addr[nt] = &phi1;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "theta2");
  addr[nt] = &theta2;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "phi2");
  addr[nt] = &phi2;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "rperi");
  addr[nt] = &rperi;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "rstart");
  addr[nt] = &rstart;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ecc");
  addr[nt] = &ecc;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g1_x_start");
  addr[nt] = &g1_x_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g1_y_start");
  addr[nt] = &g1_y_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g1_z_start");
  addr[nt] = &g1_z_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g1_vx_start");
  addr[nt] = &g1_vx_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g1_vy_start");
  addr[nt] = &g1_vy_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g1_vz_start");
  addr[nt] = &g1_vz_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g2_x_start");
  addr[nt] = &g2_x_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g2_y_start");
  addr[nt] = &g2_y_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g2_z_start");
  addr[nt] = &g2_z_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g2_vx_start");
  addr[nt] = &g2_vx_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g2_vy_start");
  addr[nt] = &g2_vy_start;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "g2_vz_start");
  addr[nt] = &g2_vz_start;
  id[nt++] = FLOAT;


  if((fd = fopen(fname, "r")))
    {
      while(!feof(fd))
	{
	  buf[0] = 0;
	  fgets(buf, 200, fd);

	  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
	    continue;

	  if(buf1[0] == '%')
	    continue;

	  for(i = 0, j = -1; i < nt; i++)
	    if(strcmp(buf1, tag[i]) == 0)
	      {
		j = i;
		tag[i][0] = 0;
		break;
	      }

	  if(j >= 0)
	    {
	      switch (id[j])
		{
		case FLOAT:
		  *((double *) addr[j]) = atof(buf2);
		  break;
		case STRING:
		  strcpy(addr[j], buf2);
		  break;
		case INT:
		  *((int *) addr[j]) = atoi(buf2);
		  break;
		}
	    }
	  else
	    {
	      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);
	      errorFlag = 1;
	    }
	}
      fclose(fd);

    }
  else
    {
      fprintf(stdout, "Parameter file %s not found.\n", fname);
      errorFlag = 1;
    }


  for(i = 0; i < nt; i++)
    {
      if(*tag[i])
	{
	  fprintf(stdout, "Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	  errorFlag = 1;
	}
    }

  if(errorFlag)
    {
      exit(1);
    }



#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS
}
