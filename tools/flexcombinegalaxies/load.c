#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>

#include "globvars.h"
#include "nrsrc/nrutil.h"


void load_particles_standard(char *fname, struct galaxy_data *g, struct io_header *header);
void load_particles_hdf5(char *fname, struct galaxy_data *g, struct io_header *header);



void load_particles(char *fname, struct galaxy_data *g, struct io_header *header)
{
  printf("opening: %s\n", fname);
  if(strstr(fname,"hdf5"))
        load_particles_hdf5(fname, g, header);
  else
        load_particles_standard(fname, g, header);
}





/* ---------------------
    std binary format
   --------------------- */
void load_particles_standard(char *fname, struct galaxy_data *g, struct io_header *header)
{
  FILE *fd;
  int i, type, pc, blksize;
  int ntot_withmasses;

#define SKIP  fread(&blksize,sizeof(int),1,fd);

  if(!(fd = fopen(fname, "r")))
    {
      fprintf(stderr, "File '%s' not found.\n", fname);
      exit(0);
    }

  SKIP;
  fread(header, sizeof(struct io_header), 1, fd);
  SKIP;



  g->Ngas = header->npart[0];
  g->Ndm = header->npart[1];

  for(i = 0, g->Ntot = 0; i < 6; i++)
    g->Ntot += header->npart[i];


  g->pos = matrix(1, g->Ntot, 1, 3);
  g->vel = matrix(1, g->Ntot, 1, 3);
  g->id = ivector(1, g->Ntot);
  g->m = vector(1, g->Ntot);
  g->z = vector(1, g->Ntot);

  if(g->Ngas)
    g->u = vector(1, g->Ngas);


  SKIP;
  fread(&g->pos[1][1], sizeof(float), 3 * g->Ntot, fd);
  SKIP;

  SKIP;
  fread(&g->vel[1][1], sizeof(float), 3 * g->Ntot, fd);
  SKIP;

  SKIP;
  fread(&g->id[1], sizeof(int), g->Ntot, fd);
  SKIP;


  for(i = 0, ntot_withmasses = 0; i < 6; i++)
    {
      if(header->mass[i] == 0)
	ntot_withmasses += header->npart[i];
    }

  if(ntot_withmasses)
    SKIP;
  for(type = 0, pc = 1; type < 6; type++)
    {
      for(i = 0; i < header->npart[type]; i++)
	{
	  if(header->mass[type] == 0)
	    fread(&g->m[pc], sizeof(float), 1, fd);
	  else
	    g->m[pc] = header->mass[type];

	  pc++;
	}
    }
  if(ntot_withmasses)
    SKIP;

  if(g->Ngas)
    {
      SKIP;
      fread(&g->u[1], sizeof(float), g->Ngas, fd);
      SKIP;
    }

  fclose(fd);


  for(i = 1; i <= g->Ntot; i++)
    g->z[i] = 0.0;

  for(i = 1, g->Mtot = 0; i <= g->Ntot; i++)
    g->Mtot += g->m[i];

  for(i = 1 + g->Ngas, g->Mdm = 0; i <= (g->Ngas + g->Ndm); i++)
    g->Mdm += g->m[i];
}










/* ---------------------
       HDF5 format
   --------------------- */
void load_particles_hdf5(char *fname, struct galaxy_data *g, struct io_header *header)
{
  char   buf[200];
  int    i,k,ntot_withmasses;
  int    add_idx;
  
  float *xyz, *x;
  int *id;
  
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;
  hid_t hdf5_grp[6];
  hid_t hdf5_dataset;
  

  sprintf(buf,"%s",fname);

  hdf5_file = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header");
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, header->npart);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header->npartTotal);
  H5Aclose(hdf5_attribute);
  
/*
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header->npartTotalHighWord);
  H5Aclose(hdf5_attribute);
*/

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header->mass);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header->time);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Redshift");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header->redshift);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "BoxSize");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header->BoxSize);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header->num_files);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Omega0");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header->Omega0);
  H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "OmegaLambda");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header->OmegaLambda);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "HubbleParam");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header->HubbleParam);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Sfr");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header->flag_sfr);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Cooling");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header->flag_cooling);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_StellarAge");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header->flag_stellarage);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Metals");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header->flag_metals);
  H5Aclose(hdf5_attribute);

/*
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_DoublePrecision");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header->flag_doubleprecision);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_IC_Info");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header->flag_ic_info);
  H5Aclose(hdf5_attribute);
*/

  H5Gclose(hdf5_headergrp);


  g->Ngas = header->npart[0];
  g->Ndm = header->npart[1];
  for(k=0, g->Ntot=0; k<6; k++)
        g->Ntot += header->npart[k];

/*
printf("Ngas= %d\n",g->Ngas);
printf("Ntot= %d\n",g->Ntot);
*/

  for(k=0, ntot_withmasses=0; k<6; k++)
    {
        if(header->mass[k]==0)
          ntot_withmasses+= header->npart[k];
    }


  /* allocate memory for the temporary read-in arrays */
  g->pos = matrix(1, g->Ntot, 1, 3);
  g->vel = matrix(1, g->Ntot, 1, 3);
  g->id = ivector(1, g->Ntot);
  g->m = vector(1, g->Ntot);
  g->z = vector(1, g->Ntot);

  if(g->Ngas)
    g->u = vector(1, g->Ngas);


  /* temporary array's used to write file */
  if(!(id = (int *) malloc((g->Ntot)*sizeof(int))))
    {
        printf("failed to allocate memory (id).\n");
        exit(0);
    }
  if(!(x = (float *) malloc((g->Ntot)*sizeof(float))))
    {
        printf("failed to allocate memory (x).\n");
        exit(0);
    }
  if(!(xyz = (float *) malloc((3*g->Ntot)*sizeof(float))))
    {
        printf("failed to allocate memory (xyz).\n");
        exit(0);
    }





  /* Now, let's open the actual particle information */
  for(k = 0; k < 6; k++)
  {
    if(header->npart[k] > 0)
      {
        sprintf(buf, "/PartType%d", k);
        hdf5_grp[k] = H5Gopen(hdf5_file, buf);
      }
  }


  if(header->npart[0] > 0)
    {
        printf("reading gas\n"); fflush(stdout);

        hdf5_dataset = H5Dopen(hdf5_grp[0], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[0]; i++)
          {
                g->pos[i][1]= xyz[(i-1)*3 + 0];
                g->pos[i][2]= xyz[(i-1)*3 + 1];
                g->pos[i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[0], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[0]; i++)
          {
                g->vel[i][1]= xyz[(i-1)*3 + 0];
                g->vel[i][2]= xyz[(i-1)*3 + 1];
                g->vel[i][3]= xyz[(i-1)*3 + 2];
          }


        hdf5_dataset = H5Dopen(hdf5_grp[0], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[0]; i++)
                g->id[i]= id[i-1];

        if(header->mass[0] > 0)
          for(i=1; i<=header->npart[0]; i++)
                g->m[i]= header->mass[0];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[0], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header->npart[0];i++)
                g->m[i]= x[i-1];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[0], "InternalEnergy");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[0]; i++)
                g->u[i]= x[i-1];

        hdf5_dataset = H5Dopen(hdf5_grp[0], "Metallicity");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[0]; i++)
                g->z[i]= x[i-1];
    }


  if(header->npart[1] > 0)
    {
        printf("reading halo\n"); fflush(stdout);

        add_idx= header->npart[0];

        hdf5_dataset = H5Dopen(hdf5_grp[1], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[1];i++)
          {
                g->pos[add_idx + i][1]= xyz[(i-1)*3 + 0];
                g->pos[add_idx + i][2]= xyz[(i-1)*3 + 1];
                g->pos[add_idx + i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[1], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[1];i++)
          {
                g->vel[add_idx + i][1]= xyz[(i-1)*3 + 0];
                g->vel[add_idx + i][2]= xyz[(i-1)*3 + 1];
                g->vel[add_idx + i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[1], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[1];i++)
                g->id[add_idx + i]= id[i-1];

        if(header->mass[1] > 0)
          for(i=1; i<=header->npart[1]; i++)
                g->m[add_idx + i]= header->mass[1];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[1], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header->npart[1];i++)
                g->m[add_idx + i]= x[i-1];
          }
    }



  if(header->npart[2] > 0)
    {
        printf("reading disk\n"); fflush(stdout);

        add_idx= header->npart[0] + header->npart[1];

        hdf5_dataset = H5Dopen(hdf5_grp[2], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[2];i++)
          {
                g->pos[add_idx + i][1]= xyz[(i-1)*3 + 0];
                g->pos[add_idx + i][2]= xyz[(i-1)*3 + 1];
                g->pos[add_idx + i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[2], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[2];i++)
          {
                g->vel[add_idx + i][1]= xyz[(i-1)*3 + 0];
                g->vel[add_idx + i][2]= xyz[(i-1)*3 + 1];
                g->vel[add_idx + i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[2], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[2];i++)
                g->id[add_idx + i]= id[i-1];

        if(header->mass[2] > 0)
          for(i=1; i<=header->npart[2]; i++)
                g->m[add_idx + i]= header->mass[2];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[2], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header->npart[2];i++)
                g->m[add_idx + i]= x[i-1];
          }
    }



  if(header->npart[3] > 0)
    {
        printf("reading bulge\n"); fflush(stdout);

        add_idx= header->npart[0] + header->npart[1] + header->npart[2];

        hdf5_dataset = H5Dopen(hdf5_grp[3], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[3];i++)
          {
                g->pos[add_idx + i][1]= xyz[(i-1)*3 + 0];
                g->pos[add_idx + i][2]= xyz[(i-1)*3 + 1];
                g->pos[add_idx + i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[3], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[3];i++)
          {
                g->vel[add_idx + i][1]= xyz[(i-1)*3 + 0];
                g->vel[add_idx + i][2]= xyz[(i-1)*3 + 1];
                g->vel[add_idx + i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[3], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[3];i++)
                g->id[add_idx + i]= id[i-1];

        if(header->mass[3] > 0)
          for(i=1; i<=header->npart[3]; i++)
                g->m[add_idx + i]= header->mass[3];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[3], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header->npart[3];i++)
                g->m[add_idx + i]= x[i-1];
          }
    }




  if(header->npart[4] > 0)
    {
        printf("reading new stars\n"); fflush(stdout);

        add_idx= header->npart[0] + header->npart[1] + header->npart[2] + header->npart[3];

        hdf5_dataset = H5Dopen(hdf5_grp[4], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[4];i++)
          {
                g->pos[add_idx + i][1]= xyz[(i-1)*3 + 0];
                g->pos[add_idx + i][2]= xyz[(i-1)*3 + 1];
                g->pos[add_idx + i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[4], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[4];i++)
          {
                g->vel[add_idx + i][1]= xyz[(i-1)*3 + 0];
                g->vel[add_idx + i][2]= xyz[(i-1)*3 + 1];
                g->vel[add_idx + i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[4], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[4];i++)
                g->id[add_idx + i]= id[i-1];

        if(header->mass[4] > 0)
          for(i=1; i<=header->npart[4]; i++)
                g->m[add_idx + i]= header->mass[4];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[4], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header->npart[4];i++)
                g->m[add_idx + i]= x[i-1];
          }
    }


  if(header->npart[5] > 0)
    {
        printf("reading black holes\n"); fflush(stdout);

        add_idx= header->npart[0] + header->npart[1] + header->npart[2] + header->npart[3] + header->npart[4];

        hdf5_dataset = H5Dopen(hdf5_grp[5], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[5];i++)
          {
                g->pos[add_idx + i][1]= xyz[(i-1)*3 + 0];
                g->pos[add_idx + i][2]= xyz[(i-1)*3 + 1];
                g->pos[add_idx + i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[5], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[5];i++)
          {
                g->vel[add_idx + i][1]= xyz[(i-1)*3 + 0];
                g->vel[add_idx + i][2]= xyz[(i-1)*3 + 1];
                g->vel[add_idx + i][3]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[5], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header->npart[5];i++)
                g->id[add_idx + i]= id[i-1];

        if(header->mass[5] > 0)
          for(i=1; i<=header->npart[5]; i++)
                g->m[add_idx + i]= header->mass[5];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[5], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header->npart[5];i++)
                g->m[add_idx + i]= x[i-1];
          }
    }




  /* close this thing */
  for(k = 5; k >= 0; k--)
    if(header->npart[k] > 0)
      H5Gclose(hdf5_grp[k]);
  H5Fclose(hdf5_file);





  for(i = 1, g->Mtot = 0; i <= g->Ntot; i++)
    g->Mtot += g->m[i];

  for(i = 1 + g->Ngas, g->Mdm = 0; i <= (g->Ngas + g->Ndm); i++)
    g->Mdm += g->m[i];



}






