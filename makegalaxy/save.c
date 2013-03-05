#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"

#include <hdf5.h>


void save_particles_hdf(char *fname);
void save_particles_old(char *fname);


#ifdef T3E
typedef short int int4byte;	/* Note: int has 8 Bytes on the T3E ! */
#else
typedef int int4byte;
#endif





struct io_header_1
{
  int4byte npart[6];
  double mass[6];
  double time;
  double redshift;
  int4byte flag_sfr;
  int4byte flag_feedback;
  int4byte npartTotal[6];
  int4byte flag_cooling;
  int4byte num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
}
header1;




void save_particles(void)
{
  char buf[100];

  mkdir(OutputDir, 02755);

  /* given as *.dat, but change it to *.hdf5 */
  if(strstr(OutputFile,"dat"))
    {
        strcpy(buf, "");
        strncpy(buf, OutputFile, strlen(OutputFile)-3);
        strcat(buf, "hdf5");
        strcpy(OutputFile, buf);
    }

  /* add .hdf5 if it doesn't have one (no matter what the current one is) */
  if(!(strstr(OutputFile,"hdf5")) || (strlen(strstr(OutputFile,"hdf5")) > 4))
    strcat(OutputFile, ".hdf5");

  sprintf(buf, "%s/%s", OutputDir, OutputFile);

  save_particles_hdf(buf);
}



void save_particles_old(char *fname)
{
  FILE *fd;
  int i;
  float xyz[3];
  int4byte blklen;

#define BLKLEN fwrite(&blklen, sizeof(blklen), 1, fd);


  if(!(fd = fopen(fname, "w")))
    {
      printf("error opening file %s\n", fname);
      exit(0);
    }

  printf("\nsaving initial conditions to file `%s'\n", fname);

  header1.npart[0] = header1.npartTotal[0] = N_GAS;
  header1.npart[1] = header1.npartTotal[1] = N_HALO;
  header1.npart[2] = header1.npartTotal[2] = N_DISK;
  header1.npart[3] = header1.npartTotal[3] = N_BULGE;
  header1.npart[4] = 0;
  header1.npart[5] = header1.npartTotal[5] = N_BLACKHOLE;

  header1.num_files = 0;

  for(i = 0; i < 6; i++)
    header1.mass[i] = 0;

  if(N_GAS)
    header1.mass[0] = mp_gas[1];

  if(N_HALO)
    header1.mass[1] = mp_halo[1];

  if(N_DISK)
    header1.mass[2] = mp_disk[1];

  if(N_BULGE)
    header1.mass[3] = mp_bulge[1];

  if(N_BLACKHOLE)
    header1.mass[5] = M_BLACKHOLE;

  header1.flag_sfr = 0;
  header1.flag_feedback = 0;
  header1.flag_cooling = 0;


  blklen = sizeof(header1);
  BLKLEN;
  fwrite(&header1, sizeof(header1), 1, fd);
  BLKLEN;


  blklen = 3 * (N_GAS + N_HALO + N_DISK + N_BULGE + N_BLACKHOLE) * sizeof(float);
  BLKLEN;
  for(i = 1; i <= N_GAS; i++)
    {
      xyz[0] = xp_gas[i];
      xyz[1] = yp_gas[i];
      xyz[2] = zp_gas[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_HALO; i++)
    {
      xyz[0] = xp_halo[i];
      xyz[1] = yp_halo[i];
      xyz[2] = zp_halo[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_DISK; i++)
    {
      xyz[0] = xp_disk[i];
      xyz[1] = yp_disk[i];
      xyz[2] = zp_disk[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_BULGE; i++)
    {
      xyz[0] = xp_bulge[i];
      xyz[1] = yp_bulge[i];
      xyz[2] = zp_bulge[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }

  for(i = 1; i <= N_BLACKHOLE; i++)
    {
      xyz[0] = xyz[1] = xyz[2] = 0;
      fwrite(xyz, sizeof(float), 3, fd);
    }

  BLKLEN;



  blklen = 3 * (N_GAS + N_HALO + N_DISK + N_BULGE + N_BLACKHOLE) * sizeof(float);
  BLKLEN;
  for(i = 1; i <= N_GAS; i++)
    {
      xyz[0] = vxp_gas[i];
      xyz[1] = vyp_gas[i];
      xyz[2] = vzp_gas[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_HALO; i++)
    {
      xyz[0] = vxp_halo[i];
      xyz[1] = vyp_halo[i];
      xyz[2] = vzp_halo[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_DISK; i++)
    {
      xyz[0] = vxp_disk[i];
      xyz[1] = vyp_disk[i];
      xyz[2] = vzp_disk[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_BULGE; i++)
    {
      xyz[0] = vxp_bulge[i];
      xyz[1] = vyp_bulge[i];
      xyz[2] = vzp_bulge[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_BLACKHOLE; i++)
    {
      xyz[0] = xyz[1] = xyz[2] = 0;
      fwrite(xyz, sizeof(float), 3, fd);
    }
  BLKLEN;


  blklen = (N_GAS + N_HALO + N_DISK + N_BULGE + N_BLACKHOLE) * sizeof(int);
  BLKLEN;
  for(i = 1; i <= (N_GAS + N_HALO + N_DISK + N_BULGE + N_BLACKHOLE); i++)
    {
      fwrite(&i, sizeof(int), 1, fd);	/* ID */
    }
  BLKLEN;

  if(N_GAS)
    {
      blklen = (N_GAS) * sizeof(float);
      BLKLEN;
      for(i = 1; i <= N_GAS; i++)
	{
	  xyz[0] = u_gas[i];
	  fwrite(xyz, sizeof(float), 1, fd);
	}
      BLKLEN;
    }

  fclose(fd);
}





/*  The case were we save everything in HDF5 format */ 
void save_particles_hdf(char *fname)
{
  int i, type;
  float *xyz, *x;
  int *id;
  double R;
  int idnum= 1, dummyzero= 0, dummyone= 1;
  int dummyzeroarray[6];

  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0;
  hid_t hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  hsize_t dims[2];
  char buf[500];


  for(i = 0; i < 6; i++)
    dummyzeroarray[i]= 0;


  /* temporary array's used to write file */
  if(!(id = (int *) malloc((N_TOTAL)*sizeof(int))))
    {
	printf("failed to allocate memory (id).\n");
	exit(0);
    }
  if(!(x = (float *) malloc((N_TOTAL)*sizeof(float))))
    {
	printf("failed to allocate memory (x).\n");
	exit(0);
    }
  if(!(xyz = (float *) malloc((3*N_TOTAL)*sizeof(float))))
    {
	printf("failed to allocate memory (xyz).\n");
	exit(0);
    }





  /* setup the header */

  header1.npart[0] = header1.npartTotal[0] = N_GAS;
  header1.npart[1] = header1.npartTotal[1] = N_HALO;
  header1.npart[2] = header1.npartTotal[2] = N_DISK;
  header1.npart[3] = header1.npartTotal[3] = N_BULGE;
  header1.npart[4] = 0;
  header1.npart[5] = header1.npartTotal[5] = N_BLACKHOLE;

  header1.num_files = 0;

  for(i = 0; i < 6; i++)
    header1.mass[i] = 0;

  if(N_GAS && !WriteDensity)
    header1.mass[0] = mp_gas[1];

  if(N_HALO)
    header1.mass[1] = mp_halo[1];

  if(N_DISK)
    header1.mass[2] = mp_disk[1];

  if(N_BULGE)
    header1.mass[3] = mp_bulge[1];

  if(N_BLACKHOLE)
    header1.mass[5] = M_BLACKHOLE;

  header1.flag_sfr = 0;
  header1.flag_feedback = 0;
  header1.flag_cooling = 0;
  header1.HubbleParam= HUBBLE;



  /* Now, start writing to the file */

  /* sprintf(buf, "%s.hdf5", fname);   - this was already taken care of in main.c  */
  sprintf(buf, "%s", fname);
  printf("\nsaving initial conditions to file `%s'\n", buf);
 
  hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

  for(type = 0; type < 6; type++)
    {
      if(header1.npart[type] > 0)
        {
          sprintf(buf, "/PartType%d", type);
          hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
        }
    }



  /* ---------------- */
  /* Write the Header */
  printf("writing header ... \n"); fflush(stdout);

  hsize_t adim[1] = { 6 };
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header1.npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header1.npartTotal);
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
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header1.mass);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.BoxSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &dummyone);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.Omega0);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.OmegaLambda);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

 hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.HubbleParam);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header1.flag_sfr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(hdf5_headergrp, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header1.flag_cooling);
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


  if(N_GAS > 0)
    {
	printf("writing gas ... \n"); fflush(stdout);

	/* Positions */
	for(i = 1; i <= N_GAS; i++)
	  {
		xyz[(i-1)*3 + 0] = xp_gas[i];
		xyz[(i-1)*3 + 1] = yp_gas[i];
		xyz[(i-1)*3 + 2] = zp_gas[i];
	  }
	dims[0] = header1.npart[0];
	dims[1] = 3;                   /* x, y, z */
	hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);    /* 2 -> 1 when dims[1]= 1 */
	hdf5_dataset =  H5Dcreate(hdf5_grp[0], "Coordinates", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
	hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
	H5Dclose(hdf5_dataset);
	H5Sclose(hdf5_dataspace_in_file);

	/* Velocities */
        for(i = 1; i <= N_GAS; i++)
          {
                xyz[(i-1)*3 + 0] = vxp_gas[i];
                xyz[(i-1)*3 + 1] = vyp_gas[i];
                xyz[(i-1)*3 + 2] = vzp_gas[i];
          }

        dims[0] = header1.npart[0];
        dims[1] = 3;                  /* vx, vy, vz */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[0], "Velocities", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

	/* IDs */
        for(i = 1; i <= N_GAS; i++, idnum++)
                id[i-1] = idnum;

        dims[0] = header1.npart[0];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[0], "ParticleIDs", H5T_NATIVE_UINT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

	/* Masses - might be in header*/
	if(WriteDensity) {
	  for(i = 1; i <= N_GAS; i++)
	    x[i-1]= rho_gas[i];
	  dims[0] = header1.npart[0];
	  dims[1] = 1;
	  hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
	  hdf5_dataset =  H5Dcreate(hdf5_grp[0], "Masses", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
	  hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
	  H5Dclose(hdf5_dataset);
	  H5Sclose(hdf5_dataspace_in_file);
	}

	/* Internal Energy */
	for(i = 1; i <= N_GAS; i++)
		x[i-1]= u_gas[i];
        dims[0] = header1.npart[0];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[0], "InternalEnergy", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

	/* Metallicity */
        for(i = 1; i <= N_GAS; i++)
          {
		R= sqrt(xp_gas[i]*xp_gas[i] + yp_gas[i]*yp_gas[i]);
                x[i-1] = Zmet0 * exp(R * log(10) * ZmetGrad);
          }
        dims[0] = header1.npart[0];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[0], "Metallicity", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

    }




  if(N_HALO > 0)
    {
	printf("writing halo ... \n"); fflush(stdout);

        /* Positions */
        for(i = 1; i <= N_HALO; i++)
          {
		xyz[(i-1)*3 + 0] = xp_halo[i];
		xyz[(i-1)*3 + 1] = yp_halo[i];
		xyz[(i-1)*3 + 2] = zp_halo[i];
          }
  
        dims[0] = header1.npart[1];
        dims[1] = 3;                   /* x, y, z */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);    /* 2 -> 1 when dims[1]= 1 */
        hdf5_dataset =  H5Dcreate(hdf5_grp[1], "Coordinates", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);
  
        /* Velocities */
        for(i = 1; i <= N_HALO; i++)
          {
		xyz[(i-1)*3 + 0] = vxp_halo[i];
		xyz[(i-1)*3 + 1] = vyp_halo[i];
		xyz[(i-1)*3 + 2] = vzp_halo[i];
          }
  
        dims[0] = header1.npart[1];
        dims[1] = 3;                  /* vx, vy, vz */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[1], "Velocities", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);
  
        /* IDs */
        for(i = 1; i <= N_HALO; i++, idnum++)
                id[i-1] = idnum;
  
        dims[0] = header1.npart[1];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[1], "ParticleIDs", H5T_NATIVE_UINT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);
  
        /* Masses - are in header*/
    }
  



  if(N_DISK > 0)
    {
	printf("writing disk ... \n"); fflush(stdout);

        /* Positions */ 
        for(i = 1; i <= N_DISK; i++)
          {
                xyz[(i-1)*3 + 0] = xp_disk[i];
                xyz[(i-1)*3 + 1] = yp_disk[i];
                xyz[(i-1)*3 + 2] = zp_disk[i];
          }

        dims[0] = header1.npart[2];
        dims[1] = 3;                   /* x, y, z */ 
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);    /* 2 -> 1 when dims[1]= 1 */
        hdf5_dataset =  H5Dcreate(hdf5_grp[2], "Coordinates", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Velocities */
        for(i = 1; i <= N_DISK; i++)
          {
                xyz[(i-1)*3 + 0] = vxp_disk[i];
                xyz[(i-1)*3 + 1] = vyp_disk[i];
                xyz[(i-1)*3 + 2] = vzp_disk[i];
          }

        dims[0] = header1.npart[2];
        dims[1] = 3;                  /* vx, vy, vz */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[2], "Velocities", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* IDs */
        for(i = 1; i <= N_DISK; i++, idnum++)
                id[i-1] = idnum;

        dims[0] = header1.npart[2];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[2], "ParticleIDs", H5T_NATIVE_UINT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Masses - are in header*/

        /* Metallicity */
        for(i = 1; i <= N_DISK; i++)
          {
                R= sqrt(xp_disk[i]*xp_disk[i] + yp_disk[i]*yp_disk[i]);
                x[i-1] = Zmet0 * exp(R * log(10) * ZmetGrad);
          }
        dims[0] = header1.npart[2];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[2], "Metallicity", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Formation Time */
	if(strcmp(DiskPopMode,"constant") == 0)
          for(i = 1; i <= N_DISK; i++)
              x[i-1] = -DiskPopAge * drand48();
	else if(strcmp(DiskPopMode,"instantaneous") == 0)
          for(i = 1; i <= N_DISK; i++)
              x[i-1] = -DiskPopAge;
	else if(strcmp(DiskPopMode,"exponential") == 0)
          for(i = 1; i <= N_DISK; i++)
              x[i-1] = -DiskPopTau * log(1. + drand48() * (exp(DiskPopAge/DiskPopTau) - 1.));
	else
	  {
	    printf("Hey, use the correct DiskPopMode, not, %s!\n",DiskPopMode);
	    exit(0);
	  }

        dims[0] = header1.npart[2];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[2], "StellarFormationTime", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);
    }





  if(N_BULGE > 0)
    {
	printf("writing bulge ... \n"); fflush(stdout);

        /* Positions */
        for(i = 1; i <= N_BULGE; i++)
          {
                xyz[(i-1)*3 + 0] = xp_bulge[i];
                xyz[(i-1)*3 + 1] = yp_bulge[i];
                xyz[(i-1)*3 + 2] = zp_bulge[i];
          }

        dims[0] = header1.npart[3];
        dims[1] = 3;                   /* x, y, z */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);    /* 2 -> 1 when dims[1]= 1 */
        hdf5_dataset =  H5Dcreate(hdf5_grp[3], "Coordinates", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Velocities */
        for(i = 1; i <= N_BULGE; i++)
          {
                xyz[(i-1)*3 + 0] = vxp_bulge[i];
                xyz[(i-1)*3 + 1] = vyp_bulge[i];
                xyz[(i-1)*3 + 2] = vzp_bulge[i];
          }

        dims[0] = header1.npart[3];
        dims[1] = 3;                  /* vx, vy, vz */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[3], "Velocities", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* IDs */
        for(i = 1; i <= N_BULGE; i++, idnum++)
                id[i-1] = idnum;

        dims[0] = header1.npart[3];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[3], "ParticleIDs", H5T_NATIVE_UINT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Masses - are in header*/

        /* Metallicity */
        for(i = 1; i <= N_BULGE; i++)
          {
                R= sqrt(xp_bulge[i]*xp_bulge[i] + yp_bulge[i]*yp_bulge[i]);
                x[i-1] = Zmet0 * exp(R * log(10) * ZmetGrad);
          }
        dims[0] = header1.npart[3];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[3], "Metallicity", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Formation Time */
        if(strcmp(BulgePopMode,"constant") == 0)
          for(i = 1; i <= N_BULGE; i++)
              x[i-1] = -BulgePopAge * drand48();
        else if(strcmp(BulgePopMode,"instantaneous") == 0)
          for(i = 1; i <= N_BULGE; i++)
              x[i-1] = -BulgePopAge;
        else if(strcmp(BulgePopMode,"exponential") == 0)
          for(i = 1; i <= N_BULGE; i++)
              x[i-1] = -BulgePopTau * log(1. + drand48() * (exp(BulgePopAge/BulgePopTau) - 1.));
        else
          {
            printf("Hey, use the correct BulgePopMode, not, %s!\n",BulgePopMode);
            exit(0);
          }

        dims[0] = header1.npart[3];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[3], "StellarFormationTime", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);
    }



  if(N_BLACKHOLE > 0)
    {
	printf("writing black holes ... \n"); fflush(stdout);

        /* Coordinates */
        x[0] = x[1] = x[2] = 0.0;
        dims[0] = header1.npart[5];
        dims[1] = 3;                   /* x, y, z */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);    /* 2 -> 1 when dims[1]= 1 */
        hdf5_dataset =  H5Dcreate(hdf5_grp[5], "Coordinates", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Velocities */
        x[0] = x[1] = x[2] = 0.0;
        dims[0] = header1.npart[5];
        dims[1] = 3;                  /* vx, vy, vz */
        hdf5_dataspace_in_file = H5Screate_simple(2, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[5], "Velocities", H5T_NATIVE_FLOAT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* IDs */
        id[0] = idnum;
        dims[0] = header1.npart[5];
        dims[1] = 1;
        hdf5_dataspace_in_file = H5Screate_simple(1, dims, NULL);
        hdf5_dataset =  H5Dcreate(hdf5_grp[5], "ParticleIDs", H5T_NATIVE_UINT, hdf5_dataspace_in_file, H5P_DEFAULT);
        hdf5_status = H5Dwrite(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        H5Sclose(hdf5_dataspace_in_file);

        /* Masses - are in header*/

	/* Added BH_Mdot and BH_mass */
    }



  free(x);
  free(xyz);
  free(id);


  /* ------------------------------------------ */
  /* OK, we're done writing, now close and exit */

  printf("done!\n"); fflush(stdout);

  for(type = 5; type >= 0; type--)
    if(header1.npart[type] > 0)
      H5Gclose(hdf5_grp[type]);
  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);


}


