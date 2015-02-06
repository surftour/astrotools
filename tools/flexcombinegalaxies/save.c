#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>

#include "globvars.h"


void save_combined_regular(char *fname, struct galaxy_data *g1, struct galaxy_data *g2);
void save_combined_hdf5(char *fname, struct galaxy_data *g1, struct galaxy_data *g2);



void save_combined(char *fname, struct galaxy_data *g1, struct galaxy_data *g2)
{
  printf("saving to: %s\n", fname);

  if(strstr(fname,"hdf5"))
        save_combined_hdf5(fname, g1, g2);
  else
        save_combined_regular(fname, g1, g2);
}








/*
 *
 * this routine saves to the old binary Gadget format
 *
 */
void save_combined_regular(char *fname, struct galaxy_data *g1, struct galaxy_data *g2)
{
  FILE *fd;
  int i, type, pc1, pc2, blklen, ntot_withmasses;
  struct io_header new_header;

#define BLKLEN fwrite(&blklen, sizeof(blklen), 1, fd);

  printf("\n\n  -----  WARNING: the code is somehow saving the output to the OLD binary format!  ----- \n\n");

  if(!(fd = fopen(fname, "w")))
    {
      fprintf(stderr, "Error in writing to '%s'\n", fname);
      exit(0);
    }

  new_header = header[0];

  for(i = 0; i < 6; i++)
    {
      new_header.npart[i] += header[1].npart[i];
      new_header.npartTotal[i] += header[1].npartTotal[i];

      if(header[0].mass[i] != header[1].mass[i])
	new_header.mass[i] = 0;
    }

  for(i = 1; i <= g2->Ntot; i++)
    g2->id[i] += g1->Ntot;




  blklen = sizeof(struct io_header);
  BLKLEN;
  fwrite(&new_header, sizeof(struct io_header), 1, fd);
  BLKLEN;


  blklen = 3 * (g1->Ntot + g2->Ntot) * sizeof(float);
  BLKLEN;
  for(type = 0, pc1 = pc2 = 1; type < 6; type++)
    {
      for(i = 1; i <= header[0].npart[type]; i++)
	fwrite(&g1->pos[pc1++][1], sizeof(float), 3, fd);

      for(i = 1; i <= header[1].npart[type]; i++)
	fwrite(&g2->pos[pc2++][1], sizeof(float), 3, fd);
    }
  BLKLEN;


  blklen = 3 * (g1->Ntot + g2->Ntot) * sizeof(float);
  BLKLEN;
  for(type = 0, pc1 = pc2 = 1; type < 6; type++)
    {
      for(i = 1; i <= header[0].npart[type]; i++)
	fwrite(&g1->vel[pc1++][1], sizeof(float), 3, fd);

      for(i = 1; i <= header[1].npart[type]; i++)
	fwrite(&g2->vel[pc2++][1], sizeof(float), 3, fd);
    }
  BLKLEN;


  blklen = (g1->Ntot + g2->Ntot) * sizeof(int);
  BLKLEN;
  for(type = 0, pc1 = pc2 = 1; type < 6; type++)
    {
      for(i = 1; i <= header[0].npart[type]; i++)
	fwrite(&g1->id[pc1++], sizeof(int), 1, fd);

      for(i = 1; i <= header[1].npart[type]; i++)
	fwrite(&g2->id[pc2++], sizeof(int), 1, fd);
    }
  BLKLEN;



  for(i = 0, ntot_withmasses = 0; i < 6; i++)
    {
      if(new_header.mass[i] == 0)
	ntot_withmasses += new_header.npart[i];
    }

  blklen = ntot_withmasses * sizeof(float);
  if(ntot_withmasses)
    BLKLEN;
  for(type = 0, pc1 = pc2 = 1; type < 6; type++)
    {
      for(i = 1; i <= header[0].npart[type]; i++)
	if(new_header.mass[type] == 0)
	  fwrite(&g1->m[pc1++], sizeof(float), 1, fd);
	else
	  pc1++;

      for(i = 1; i <= header[1].npart[type]; i++)
	if(new_header.mass[type] == 0)
	  fwrite(&g2->m[pc2++], sizeof(float), 1, fd);
	else
	  pc2++;
    }
  if(ntot_withmasses)
    BLKLEN;


  if(new_header.npart[0] > 0)
    {
      blklen = new_header.npartTotal[0] * sizeof(float);
      BLKLEN;
      for(i = 1; i <= header[0].npart[0]; i++)
	fwrite(&g1->u[i], sizeof(float), 1, fd);
      for(i = 1; i <= header[1].npart[0]; i++)
	fwrite(&g2->u[i], sizeof(float), 1, fd);
      BLKLEN;
    }

  fclose(fd);
}









/*
 *
 * this routine saves to HDF5 format
 *
 */
void save_combined_hdf5(char *fname, struct galaxy_data *g1, struct galaxy_data *g2)
{
  int i, j, type, add_idx;
  struct io_header new_header;
  float *xyz, *x;
  int *id;

  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0;
  hid_t hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  hsize_t dims[2];
  char buf[500];
  int dummyzero=0, dummyone=1, dummyzeroarray[6];

  new_header = header[0];

  for(i = 0; i < 6; i++)
    {
      new_header.npart[i] += header[1].npart[i];
      new_header.npartTotal[i] += header[1].npartTotal[i];

      if(header[0].mass[i] != header[1].mass[i])
        new_header.mass[i] = 0;
    }

  for(i = 1; i <= g2->Ntot; i++)
    g2->id[i] += g1->Ntot;



  /* temporary array's used to write file */
  if(!(id = (int *) malloc((g1->Ntot + g2->Ntot)*sizeof(int))))
    {
        printf("failed to allocate memory (id).\n");
        exit(0);
    }
  if(!(x = (float *) malloc((g1->Ntot + g2->Ntot)*sizeof(float))))
    {
        printf("failed to allocate memory (x).\n");
        exit(0);
    }
  if(!(xyz = (float *) malloc((3*(g1->Ntot + g2->Ntot))*sizeof(float))))
    {
        printf("failed to allocate memory (xyz).\n");
        exit(0);
    }





 /* ------------------------------
    Now, start writing to the file
    ------------------------------ */

  for(i=0; i<6; i++)
    dummyzeroarray[i]=0;


  sprintf(buf, "%s", fname);

  hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

  for(type = 0; type < 6; type++)
    {
      if(new_header.npart[type] > 0)
        {
          sprintf(buf, "/PartType%d", type);
          hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
        }
    }





  /* ---------------- */
  /* Write the Header */
  /* ---------------- */

  printf("writing header ... \n"); fflush(stdout);

  hsize_t adim[1] = { 6 };
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, new_header.npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, new_header.npartTotal);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &dummyzeroarray);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, new_header.mass);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &new_header.time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &new_header.redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &new_header.BoxSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &dummyone);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &new_header.Omega0);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &new_header.OmegaLambda);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &new_header.HubbleParam);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &new_header.flag_sfr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &new_header.flag_cooling);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &dummyzero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &dummyone);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &dummyzero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &dummyzero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Flag_IC_Info", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &dummyzero);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);




  /* -------------------------------- */
  /* Done with header, move onto data */


  if(new_header.npart[0] > 0)
    {
        printf("writing gas ... \n"); fflush(stdout);

        /* Positions */
        for(i = 1, j=1; i <= header[0].npart[0]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g1->pos[i][1];
                xyz[(j-1)*3 + 1] = g1->pos[i][2];
                xyz[(j-1)*3 + 2] = g1->pos[i][3];
          }
        for(i = 1; i <= header[1].npart[0]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g2->pos[i][1];
                xyz[(j-1)*3 + 1] = g2->pos[i][2];
                xyz[(j-1)*3 + 2] = g2->pos[i][3];
          }

        dims[0] = new_header.npart[0];
        dims[1] = 3;                   /* x, y, z */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);    /* 2 -> 1 when dims[1]= 1 */
        hdf5_dataset =  H5Dcreate(hdf5_grp[0], "Coordinates", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Velocities */
        for(i = 1, j=1; i <= header[0].npart[0]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g1->vel[i][1];
                xyz[(j-1)*3 + 1] = g1->vel[i][2];
                xyz[(j-1)*3 + 2] = g1->vel[i][3];
          }
        for(i = 1; i <= header[1].npart[0]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g2->vel[i][1];
                xyz[(j-1)*3 + 1] = g2->vel[i][2];
                xyz[(j-1)*3 + 2] = g2->vel[i][3];
          }

        dims[0] = new_header.npart[0];
        dims[1] = 3;                  /* vx, vy, vz */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[0], "Velocities", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* IDs */
        for(i = 1, j=1; i <= header[0].npart[0]; i++, j++)
          id[j-1] = g1->id[i];
        for(i = 1; i <= header[1].npart[0]; i++, j++)
          id[j-1] = g2->id[i];

        dims[0] = new_header.npart[0];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[0], "ParticleIDs", H5T_NATIVE_UINT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Masses - are in header*/
	/* if(new_header.mass[0] == 0)   - not sure why this wasn't working */
	if((1.0e6 * new_header.mass[0]) == 0.0)
	  {
        	for(i = 1, j=1; i <= header[0].npart[0]; i++, j++)
                  x[j-1] = g1->m[i];
        	for(i = 1; i <= header[1].npart[0]; i++, j++)
                  x[j-1] = g2->m[i];

        	dims[0] = new_header.npart[0];
        	dims[1] = 1;
        	hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        	hdf5_dataset =  H5Dcreate(hdf5_grp[0], "Masses", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        	hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        	H5Dclose(hdf5_dataset);
		H5Sclose(hdf5_dataspace_in_file);
	  }

        /* Internal Energy */
        for(i = 1, j=1; i <= header[0].npart[0]; i++, j++)
          x[j-1] = g1->u[i];
        for(i = 1; i <= header[1].npart[0]; i++, j++)
          x[j-1] = g2->u[i];

        dims[0] = new_header.npart[0];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[0], "InternalEnergy", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Metallicity */
        for(i = 1, j=1; i <= header[0].npart[0]; i++, j++)
          x[j-1] = g1->z[i];
        for(i = 1; i <= header[1].npart[0]; i++, j++)
          x[j-1] = g2->z[i];

        dims[0] = new_header.npart[0];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[0], "Metallicity", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

    }


  if(new_header.npart[1] > 0)
    {
        printf("writing halo ... \n"); fflush(stdout);

        /* Positions */
	add_idx= header[0].npart[0];
        for(i=1, j=1; i <= header[0].npart[1]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g1->pos[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g1->pos[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g1->pos[add_idx + i][3];
          }
	add_idx= header[1].npart[0];
        for(i = 1; i <= header[1].npart[1]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g2->pos[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g2->pos[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g2->pos[add_idx + i][3];
          }

        dims[0] = new_header.npart[1];
        dims[1] = 3;                   /* x, y, z */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);    /* 2 -> 1 when dims[1]= 1 */
        hdf5_dataset =  H5Dcreate(hdf5_grp[1], "Coordinates", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Velocities */
	add_idx= header[0].npart[0];
        for(i = 1, j=1; i <= header[0].npart[1]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g1->vel[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g1->vel[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g1->vel[add_idx + i][3];
          }
	add_idx= header[1].npart[0];
        for(i = 1; i <= header[1].npart[1]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g2->vel[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g2->vel[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g2->vel[add_idx + i][3];
          }

        dims[0] = new_header.npart[1];
        dims[1] = 3;                  /* vx, vy, vz */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[1], "Velocities", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* IDs */
	add_idx= header[0].npart[0];
        for(i = 1, j=1; i <= header[0].npart[1]; i++, j++)
          id[j-1] = g1->id[add_idx + i];
	add_idx= header[1].npart[0];
        for(i = 1; i <= header[1].npart[1]; i++, j++)
          id[j-1] = g2->id[add_idx + i];

        dims[0] = new_header.npart[1];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[1], "ParticleIDs", H5T_NATIVE_UINT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Masses - if not in header*/
	if(new_header.mass[0] == 0)
	  {
		add_idx= header[0].npart[0];
        	for(i = 1, j=1; i <= header[0].npart[1]; i++, j++)
                  x[j-1] = g1->m[add_idx + i];
		add_idx= header[1].npart[0];
        	for(i = 1; i <= header[1].npart[1]; i++, j++)
                  x[j-1] = g2->m[add_idx + i];

        	dims[0] = new_header.npart[1];
        	dims[1] = 1;
        	hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        	hdf5_dataset =  H5Dcreate(hdf5_grp[1], "Masses", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        	hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        	H5Dclose(hdf5_dataset);
		H5Sclose(hdf5_dataspace_in_file);
	  }
    }




  if(new_header.npart[2] > 0)
    {
        printf("writing disk ... \n"); fflush(stdout);

        /* Positions */
	add_idx= header[0].npart[0] + header[0].npart[1];
        for(i = 1, j=1; i <= header[0].npart[2]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g1->pos[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g1->pos[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g1->pos[add_idx + i][3];
          }
	add_idx= header[1].npart[0] + header[1].npart[1];
        for(i = 1; i <= header[1].npart[2]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g2->pos[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g2->pos[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g2->pos[add_idx + i][3];
          }

        dims[0] = new_header.npart[2];
        dims[1] = 3;                   /* x, y, z */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);    /* 2 -> 1 when dims[1]= 1 */
        hdf5_dataset =  H5Dcreate(hdf5_grp[2], "Coordinates", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Velocities */
	add_idx= header[0].npart[0] + header[0].npart[1];
        for(i = 1, j=1; i <= header[0].npart[2]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g1->vel[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g1->vel[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g1->vel[add_idx + i][3];
          }
	add_idx= header[1].npart[0] + header[1].npart[1];
        for(i = 1; i <= header[1].npart[2]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g2->vel[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g2->vel[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g2->vel[add_idx + i][3];
          }

        dims[0] = new_header.npart[2];
        dims[1] = 3;                  /* vx, vy, vz */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[2], "Velocities", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* IDs */
	add_idx= header[0].npart[0] + header[0].npart[1];
        for(i = 1, j=1; i <= header[0].npart[2]; i++, j++)
          id[j-1] = g1->id[add_idx + i];
	add_idx= header[1].npart[0] + header[1].npart[1];
        for(i = 1; i <= header[1].npart[2]; i++, j++)
          id[j-1] = g2->id[add_idx + i];

        dims[0] = new_header.npart[2];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[2], "ParticleIDs", H5T_NATIVE_UINT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Masses - are in header*/
	if(new_header.mass[0] == 0)
	  {
		add_idx= header[0].npart[0] + header[0].npart[1];
        	for(i = 1, j=1; i <= header[0].npart[2]; i++, j++)
                  x[j-1] = g1->m[add_idx + i];
		add_idx= header[1].npart[0] + header[1].npart[1];
        	for(i = 1; i <= header[1].npart[2]; i++, j++)
                  x[j-1] = g2->m[add_idx + i];

        	dims[0] = new_header.npart[2];
        	dims[1] = 1;
        	hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        	hdf5_dataset =  H5Dcreate(hdf5_grp[2], "Masses", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        	hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        	H5Dclose(hdf5_dataset);
		H5Sclose(hdf5_dataspace_in_file);
	  }
    }





  if(new_header.npart[3] > 0)
    {
        printf("writing bulge ... \n"); fflush(stdout);

        /* Positions */
	add_idx= header[0].npart[0] + header[0].npart[1] + header[0].npart[2];
        for(i = 1, j=1; i <= header[0].npart[3]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g1->pos[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g1->pos[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g1->pos[add_idx + i][3];
          }
	add_idx= header[1].npart[0] + header[1].npart[1] + header[1].npart[2];
        for(i = 1; i <= header[1].npart[3]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g2->pos[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g2->pos[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g2->pos[add_idx + i][3];
          }

        dims[0] = new_header.npart[3];
        dims[1] = 3;                   /* x, y, z */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);    /* 2 -> 1 when dims[1]= 1 */
        hdf5_dataset =  H5Dcreate(hdf5_grp[3], "Coordinates", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Velocities */
	add_idx= header[0].npart[0] + header[0].npart[1] + header[0].npart[2];
        for(i = 1, j=1; i <= header[0].npart[3]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g1->vel[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g1->vel[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g1->vel[add_idx + i][3];
          }
	add_idx= header[1].npart[0] + header[1].npart[1] + header[1].npart[2];
        for(i = 1; i <= header[1].npart[3]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g2->vel[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g2->vel[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g2->vel[add_idx + i][3];
          }

        dims[0] = new_header.npart[3];
        dims[1] = 3;                  /* vx, vy, vz */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[3], "Velocities", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* IDs */
	add_idx= header[0].npart[0] + header[0].npart[1] + header[0].npart[2];
        for(i = 1, j=1 ; i <= header[0].npart[3]; i++, j++)
          id[j-1] = g1->id[add_idx + i];
	add_idx= header[1].npart[0] + header[1].npart[1] + header[1].npart[2];
        for(i = 1; i <= header[1].npart[3]; i++, j++)
          id[j-1] = g2->id[add_idx + i];

        dims[0] = new_header.npart[3];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[3], "ParticleIDs", H5T_NATIVE_UINT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Masses - are in header*/
	if(new_header.mass[0] == 0)
	  {
		add_idx= header[0].npart[0] + header[0].npart[1] + header[0].npart[2];
        	for(i = 1, j=1; i <= header[0].npart[3]; i++, j++)
                  x[j-1] = g1->m[add_idx + i];
		add_idx= header[1].npart[0] + header[1].npart[1] + header[1].npart[2];
        	for(i = 1; i <= header[1].npart[3]; i++, j++)
                  x[j-1] = g2->m[add_idx + i];

        	dims[0] = new_header.npart[3];
        	dims[1] = 1;
        	hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        	hdf5_dataset =  H5Dcreate(hdf5_grp[3], "Masses", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        	hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        	H5Dclose(hdf5_dataset);
		H5Sclose(hdf5_dataspace_in_file);
	  }
    }




  if(new_header.npart[5] > 0)
    {
        printf("writing black holes ... \n"); fflush(stdout);

        /* Coordinates */
	add_idx= header[0].npart[0] + header[0].npart[1] + header[0].npart[2] + header[0].npart[3] + header[0].npart[4];
        for(i = 1, j=1; i <= header[0].npart[5]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g1->pos[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g1->pos[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g1->pos[add_idx + i][3];
          }
	add_idx= header[1].npart[0] + header[1].npart[1] + header[1].npart[2] + header[1].npart[3] + header[1].npart[4];
        for(i = 1; i <= header[1].npart[5]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g2->pos[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g2->pos[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g2->pos[add_idx + i][3];
          }

        dims[0] = new_header.npart[5];
        dims[1] = 3;                   /* x, y, z */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);    /* 2 -> 1 when dims[1]= 1 */
        hdf5_dataset =  H5Dcreate(hdf5_grp[5], "Coordinates", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Velocities */
	add_idx= header[0].npart[0] + header[0].npart[1] + header[0].npart[2] + header[0].npart[3] + header[0].npart[4];
        for(i = 1, j=1; i <= header[0].npart[5]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g1->vel[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g1->vel[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g1->vel[add_idx + i][3];
          }
	add_idx= header[1].npart[0] + header[1].npart[1] + header[1].npart[2] + header[1].npart[3] + header[1].npart[4];
        for(i = 1; i <= header[1].npart[5]; i++, j++)
          {
                xyz[(j-1)*3 + 0] = g2->vel[add_idx + i][1];
                xyz[(j-1)*3 + 1] = g2->vel[add_idx + i][2];
                xyz[(j-1)*3 + 2] = g2->vel[add_idx + i][3];
          }
        dims[0] = new_header.npart[5];
        dims[1] = 3;                  /* vx, vy, vz */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[5], "Velocities", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* IDs */
	add_idx= header[0].npart[0] + header[0].npart[1] + header[0].npart[2] + header[0].npart[3] + header[0].npart[4];
        for(i = 1, j=1 ; i <= header[0].npart[5]; i++, j++)
          id[j-1] = g1->id[add_idx + i];
	add_idx= header[1].npart[0] + header[1].npart[1] + header[1].npart[2] + header[1].npart[3] + header[1].npart[4];
        for(i = 1; i <= header[1].npart[5]; i++, j++)
          id[j-1] = g2->id[add_idx + i];

        dims[0] = new_header.npart[5];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[5], "ParticleIDs", H5T_NATIVE_UINT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Masses - are in header*/
	if(new_header.mass[5] == 0)
	  {
		add_idx= header[0].npart[0] + header[0].npart[1] + header[0].npart[2] + header[0].npart[3] + header[0].npart[4];
        	for(i = 1, j=1; i <= header[0].npart[5]; i++, j++)
                  x[j-1] = g1->m[add_idx + i];
		add_idx= header[1].npart[0] + header[1].npart[1] + header[1].npart[2] + header[1].npart[3] + header[1].npart[4];
        	for(i = 1; i <= header[1].npart[5]; i++, j++)
                  x[j-1] = g2->m[add_idx + i];

        	dims[0] = new_header.npart[5];
        	dims[1] = 1;
        	hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        	hdf5_dataset =  H5Dcreate(hdf5_grp[5], "Masses", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        	hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        	H5Dclose(hdf5_dataset);
		H5Sclose(hdf5_dataspace_in_file);
	  }

        /* Added BH_Mdot and BH_mass */
    }



  free(x);
  free(xyz);
  free(id);


  /* ------------------------------------------ */
  /* OK, we're done writing, now close and exit */

  printf("done!\n"); fflush(stdout);

  for(type = 5; type >= 0; type--)
    if(new_header.npart[type] > 0)
      H5Gclose(hdf5_grp[type]);
  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);


}

