#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef HDF5
#include <hdf5.h>
#endif

#include "globvars.h"


int load_snapshot_regular(char *fname, int files, int quiteon, int headeronly);
int load_snapshot_hdf5(char *fname, int files, int quiteon, int headeronly);


int load_snapshot(char *fname, int files, int quiteon, int headeronly)
{
  if(quiteon==0) printf("opening: %s\n", fname);
  if(strstr(fname,"hdf5"))
	load_snapshot_hdf5(fname, files, quiteon, headeronly);
  else
	load_snapshot_regular(fname, files, quiteon, headeronly);

}





/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot_regular(char *fname, int files, int quiteon, int headeronly)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;
  int    first_n, dummy_total=0;
  int    dog;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file `%s`\n",buf);
	  exit(0);
	}

      /*printf("reading `%s' ...\n",buf); fflush(stdout); */
      if(quiteon==0) printf("\n");
      fflush(stdout);

      SKIP;
      fread(&header1, sizeof(header1), 1, fd);
      SKIP;
      if(SnapInfo) printf("%10d\tHeader",dummy);  dummy_total+= dummy+8;

      if(files==1)
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<6; k++)
	    NumPart+= header1.npart[k];
	  Ngas= header1.npart[0];
	}
      else
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<6; k++)
	    NumPart+= header1.npartTotal[k];
	  Ngas= header1.npartTotal[0];
	}
      if(SnapInfo) printf(" (NumPart=%10d\tNgas=%10d)\n",NumPart,Ngas);

      for(k=0, ntot_withmasses=0; k<6; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      if(headeronly==1)
	{
        fclose(fd);
	return(0);
	}

      if(i==0)
	allocate_memory();

      SKIP;
      if(SnapInfo) printf("%10d\tPositions\n",dummy); dummy_total+=dummy+8;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      if(SnapInfo) printf("%10d\tVelocities\n",dummy); dummy_total+=dummy+8;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;
    

      SKIP;
      if(SnapInfo) printf("%10d\tIDs\n",dummy); dummy_total+=dummy+8;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses>0)
	{
	  SKIP;
	  if(SnapInfo)
	    if (quiteon==0) printf("%10d\tMasses (ntot_withmasses= %d)\n",dummy,ntot_withmasses);
	  else
	    {printf("ntot_withmasses = %d\n",ntot_withmasses); fflush(stdout);}
	  dummy_total+=dummy+8;
	}
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;

	      if(header1.mass[k]==0)
		fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
	SKIP;
      

	  /* --------- Gas Internal Energy --------- */
	  if(header1.npart[0]>0)
	    {
	     SKIP;
	     if(SnapInfo) printf("%10d\tInteral Energy (Ngas)\n",dummy); dummy_total+=dummy+8;
	     for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	       {
	          fread(&P[pc_sph].U, sizeof(float), 1, fd);
	          pc_sph++;
	        }
	      SKIP;
	     }

	  /* --------- Gas Density ----------------- */
	  if(SnapFile && (header1.npart[0]>0))
	    {
		SKIP;
		if(SnapInfo) printf("%10d\tDensity (Ngas)\n",dummy); dummy_total+=dummy+8;
		for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		  {
		    fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
		    pc_sph++;
		  }
		SKIP;
	    }
	  else
	    {
		for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		   {
			P[pc_sph].Rho= 0;
			pc_sph++;
		   }
	    }

	  /* -------- Electron Abundance --------- */
	  if(header1.flag_cooling && (header1.npart[0]>0))
	    {
	      SKIP;
	      if(SnapInfo) printf("%10d\tNe (Ngas)\n",dummy); dummy_total+=dummy+8;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Ne= 1.0;
		pc_sph++;
	      }

	  /* -------- Neutral Hydrogen Density -------- */
	  if(header1.flag_cooling && (header1.npart[0]>0))
	  {
              SKIP;
	      if(SnapInfo) printf("%10d\tNH0 (Ngas)\n",dummy); dummy_total+=dummy+8;
              for(n=0, pc_sph=pc; n<header1.npart[0];n++)
                {
                  fread(&P[pc_sph].Nh, sizeof(float), 1, fd);
                  pc_sph++;
                }
              SKIP;
          }
          else
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
              {
                P[pc_sph].Nh= 0.0;
                pc_sph++;
	      }

#ifdef TJ_VERSION
	  /* --------- Formed Stellar Mass --------- */
	  if((header1.flag_sfr) && (header1.flag_stargens==0) && (header1.npart[0]>0))
	  {
              SKIP;
	      if(SnapInfo) printf("%10d\tMfs (Ngas)\n",dummy); dummy_total+=dummy+8;
              for(n=0, pc_sph=pc; n<header1.npart[0];n++)
                {
                  fread(&P[pc_sph].Mfs, sizeof(float), 1, fd);
                  pc_sph++;
                }
              SKIP;
          }
          else
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
              {
                P[pc_sph].Mfs = 0.0;
                pc_sph++;
              }

	  /* --------- Mass of Cold Clouds --------- */
	  if(header1.flag_multiphase && (header1.npart[0]>0))
	  {
              SKIP;
	      if(SnapInfo) printf("%10d\tMclouds (Ngas)\n",dummy); dummy_total+=dummy+8;
              for(n=0, pc_sph=pc; n<header1.npart[0];n++)
                {
                  fread(&P[pc_sph].Mclouds, sizeof(float), 1, fd);
                  pc_sph++;
                }
              SKIP;
          }
          else
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
              {
                P[pc_sph].Mclouds = 0.0;
                pc_sph++;
              }
#endif

	  /* --------- Gas Smoothing length --------- */
	  if(SnapFile && (header1.npart[0]>0))
	  {
              SKIP;
	      if(SnapInfo) printf("%10d\tHsml (Ngas)\n",dummy); dummy_total+=dummy+8;
              for(n=0, pc_sph=pc; n<header1.npart[0];n++)
                {
                  fread(&P[pc_sph].hsml, sizeof(float), 1, fd);
                  pc_sph++;
                }
              SKIP;
          }
          else
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
              {
                P[pc_sph].hsml = 0.0;
                pc_sph++;
              }

	  /* --------- Star Formation Rate ---------- */
	  if(header1.flag_sfr && (header1.npart[0]>0))
	  {
              SKIP;
	      if(SnapInfo) printf("%10d\tSfr (Ngas)\n",dummy); dummy_total+=dummy+8;
              for(n=0, pc_sph=pc; n<header1.npart[0];n++)
                {
                  fread(&P[pc_sph].sfr, sizeof(float), 1, fd);
                  pc_sph++;
                }
              SKIP;
          }
          else
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
              {
                P[pc_sph].sfr = 0.0;
                pc_sph++;
              }

	  /* --------- Mean formation redshift ------ */
	  if(header1.flag_stellarage)
	  {

#ifdef TJ_VERSION
	      if((header1.flag_stargens) && (header1.npart[4]>0))
#else
              if(header1.npart[4]>0)
#endif
		{
		  SKIP;
		  if(SnapInfo) printf("%10d\tStellarAge (Nstars)\n",dummy); dummy_total+=dummy+8;
		  first_n= header1.npart[0]+header1.npart[1]+header1.npart[2]+header1.npart[3]+1;
		  for(n=0, pc_sph=first_n; n<header1.npart[4];n++)
			{
			  fread(&P[pc_sph].meanz, sizeof(float), 1, fd);
			  pc_sph++;
			}
		  SKIP;
		}

#ifdef TJ_VERSION
	      if((!(header1.flag_stargens)) && (header1.npart[0]>0))
		{
              	SKIP;
	        if(SnapInfo) printf("%10d\tStellarAge (Ngas)\n",dummy); dummy_total+=dummy+8;
		for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		   {
		     fread(&P[pc_sph].meanz, sizeof(float), 1, fd);
		     pc_sph++;
		   }
		SKIP;
		}
#endif
          }
          else
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
              {
                P[pc_sph].meanz = 0.0;
                pc_sph++;
              }


          /* --------- Turbulent Energy Reservior ------ */
#ifdef TJ_VERSION
          if(header1.flag_feedbacktp && (header1.npart[0]>0))
          {
              SKIP;
	      if(SnapInfo) printf("%10d\tUtp (Ngas)\n",dummy); dummy_total+=dummy+8;
              for(n=0, pc_sph=pc; n<(header1.npart[0]);n++)
                {
                  fread(&P[pc_sph].TpU, sizeof(float), 1, fd);
                  pc_sph++;
                }
              SKIP;
          }
          else
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
              {
                P[pc_sph].TpU = 0.0;
                pc_sph++;
	      }
#endif



#ifdef TJ_VERSION
          /* --------- Potential Energy ------ */
          if(header1.flag_snaphaspot)
          {
              SKIP;
	      if(SnapInfo) printf("%10d\tPotential\n",dummy); dummy_total+=dummy+8;
              for(k=0, pc_sph=pc; k<6;k++)
                {
		  for(n=0;n<header1.npart[k];n++)
		    {
                      fread(&P[pc_sph].Potential, sizeof(float), 1, fd);
                      pc_sph++;
		    }
                }
              SKIP;
          }
          else
	    for(k=0, pc_sph=pc; k<6;k++)
	      {
            	for(n=0;n<header1.npart[k];n++)
           	  {
 		    P[pc_sph].Potential = 0.0;
                    pc_sph++;
                  }
	      }
#endif


              
#ifdef TJ_VERSION
          /* --------- Metallicity - Gas ------ */
	  if(header1.flag_metals && (header1.npart[0]>0))
	  {
	      SKIP;
	      if(SnapInfo) printf("%10d\tMetallicity (Ngas)\n",dummy); dummy_total+=dummy+8;
	      for(n=0, pc_sph=pc; n<(header1.npart[0]);n++)
		{
		  fread(&P[pc_sph].Metallicity, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	  }
	  else
	    for(n=0, pc_sph=pc; n<NumPart;n++)
	      {
		P[pc_sph].Metallicity= 0.0;
		pc_sph++;
	      }


	  /* --------- Metallicity - Stars ------ */
	  if((header1.flag_metals) && (header1.npart[4]>0))
	    {
	      SKIP;
	      if(SnapInfo) printf("%10d\tMetallicity (Nstars)\n",dummy); dummy_total+=dummy+8;
	      first_n= header1.npart[0]+header1.npart[1]+header1.npart[2]+header1.npart[3]+1;
	      for(n=0, pc_sph=first_n; n<(header1.npart[4]);n++)
	        {
	          fread(&P[pc_sph].Metallicity, sizeof(float), 1, fd);
	          pc_sph++;
	        }
	      SKIP;
	    }
#else
          /* --------- Metallicity - Gas+Stars ------ */
          if(header1.flag_metals && ((header1.npart[0]+header1.npart[4])>0))
          {
              SKIP;
              if(SnapInfo) printf("%10d\tMetallicity (Ngas+Nstars)\n",dummy); dummy_total+=dummy+8;
              for(n=0, pc_sph=pc; n<(header1.npart[0]);n++)
                {
                  fread(&P[pc_sph].Metallicity, sizeof(float), 1, fd);
                  pc_sph++;
                }
              first_n= header1.npart[0]+header1.npart[1]+header1.npart[2]+header1.npart[3]+1;
              for(n=0, pc_sph=first_n; n<(header1.npart[4]);n++)
                {
                  fread(&P[pc_sph].Metallicity, sizeof(float), 1, fd);
                  pc_sph++;
                }
              SKIP;
          }
          else
            for(n=0, pc_sph=pc; n<NumPart;n++)
              {
                P[pc_sph].Metallicity= 0.0;
                pc_sph++;
              }
#endif


#ifdef TJ_VERSION
          /* --------- Total Radiated Energy ------ */
	  if(header1.flag_energydetail && (header1.npart[0]>0))
	  {
	    SKIP;
	    if(SnapInfo) printf("%10d\tTotal Radiated Energy (Ngas)\n",dummy); dummy_total+=dummy+8;
	    for(n=0, pc_sph=pc; n<(header1.npart[0]);n++)
	      {
		fread(&P[pc_sph].totrad, sizeof(float), 1, fd);
		pc_sph++;
	      }
	    SKIP;
	  }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].totrad = 0.0;
		pc_sph++;
	      }


	  /* --------- Total Shocked Energy ------ */
	  if(header1.flag_energydetail && (header1.npart[0]>0))
	    {
	      SKIP;
	      if(SnapInfo) printf("%10d\tTotal Shocked Energy (Ngas)\n",dummy); dummy_total+=dummy+8;
	      for(n=0, pc_sph=pc; n<(header1.npart[0]);n++)
	        {
	          fread(&P[pc_sph].totshock, sizeof(float), 1, fd);
	          pc_sph++;
	        }
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
	        P[pc_sph].totshock = 0.0;
	        pc_sph++;
	      }



	  /* --------- Total Fb Energy ------ */
	  if(header1.flag_energydetail && (header1.npart[0]>0))
	    {
	      SKIP;
	      if(SnapInfo) printf("%10d\tTotal Fb Energy (Ngas)\n",dummy); dummy_total+=dummy+8;
	      for(n=0, pc_sph=pc; n<(header1.npart[0]);n++)
	        {
	          fread(&P[pc_sph].totfb, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
	        P[pc_sph].totfb = 0.0;
	        pc_sph++;
	      }



          /* --------- Parent ID ------ */
          if((header1.flag_parentid) && (header1.npart[4]>0))
            {
              SKIP;
              if(SnapInfo) printf("%10d\tParent ID (Nstars)\n",dummy); dummy_total+=dummy+8;
	      first_n= header1.npart[0]+header1.npart[1]+header1.npart[2]+header1.npart[3]+1;
              for(n=0, pc_sph=first_n; n<(header1.npart[4]);n++)
                {
                  fread(&dog, sizeof(int), 1, fd);
		  P[pc_sph].ParentID= dog;
                  pc_sph++;
                }
              SKIP;
            }
/*
          else
            for(n=0, pc_sph=pc; n<header1.npart[0];n++)
              {
                P[pc_sph].ParentID= 0.0;
                pc_sph++;
              }
*/



          /* --------- Original Gas Mass ------ */
          if((header1.flag_starorig) && (header1.npart[4]>0))
            {
              SKIP;
              if(SnapInfo) printf("%10d\tOriginal Mass (Nstars)\n",dummy); dummy_total+=dummy+8;
	      first_n= header1.npart[0]+header1.npart[1]+header1.npart[2]+header1.npart[3]+1;
              for(n=0, pc_sph=first_n; n<(header1.npart[4]);n++)
                {
                  fread(&P[pc_sph].OrigMass, sizeof(float), 1, fd);
                  pc_sph++;
                }
              SKIP;
            }
/*
          else
            for(n=0, pc_sph=pc; n<header1.npart[4];n++)
              {
                P[pc_sph].OrigMass= 0.0;
                pc_sph++;
              }
*/



          /* --------- Original Gas Hsml ------ */
          if((header1.flag_starorig) && (header1.npart[4]>0))
            {
              SKIP;
              if(SnapInfo) printf("%10d\tOriginal Hsml (Nstars)\n",dummy); dummy_total+=dummy+8;
              first_n= header1.npart[0]+header1.npart[1]+header1.npart[2]+header1.npart[3]+1;
              for(n=0, pc_sph=first_n; n<(header1.npart[4]);n++)
                {
                  fread(&P[pc_sph].OrigHsml, sizeof(float), 1, fd);
                  pc_sph++;
                }
              SKIP;
            }
#endif

#ifndef TJ_VERSION
          /* --------- Potential Energy ------ */
          /* if(header1.flag_snaphaspot) */
          /* if(0==1) */
          if(header1.time>0)
          {   
              SKIP;
              if(SnapInfo) printf("%10d\tPotential\n",dummy); dummy_total+=dummy+8;
              for(k=0, pc_sph=pc; k<6;k++)
                { 
                  for(n=0;n<header1.npart[k];n++)
                    { 
                      fread(&P[pc_sph].Potential, sizeof(float), 1, fd);
                      pc_sph++;
                    }
                }
              SKIP;
          }
          else
            for(k=0, pc_sph=pc; k<6;k++)
              { 
                for(n=0;n<header1.npart[k];n++)
                  { 
                    P[pc_sph].Potential = 0.0;
                    pc_sph++;
                  }
              }


          /* --------- Black Hole Information ------ */
          /* if(0==1) */
          if(0==1)
          if((header1.npart[5]>0) && (header1.time>0))
	  {
              SKIP;
              if(SnapInfo) printf("%10d\tBlack Hole Mass\n",dummy); dummy_total+=dummy+8;
              for(n=0, pc_sph=pc; n<header1.npart[5];n++)
                { 
                   /* fread(&P[pc_sph].Potential, sizeof(float), 1, fd); */
                   pc_sph++;
                }
              SKIP;
              SKIP;
              if(SnapInfo) printf("%10d\tBlack Hole Mdot\n",dummy); dummy_total+=dummy+8;
              for(n=0, pc_sph=pc; n<header1.npart[5];n++)
                { 
                   /* fread(&P[pc_sph].Potential, sizeof(float), 1, fd); */
                   pc_sph++;
                }
              SKIP;
	  }
          


#endif



      fclose(fd);
    }

  if(SnapInfo)
    {
      printf("----------\n");
      printf("%10d\tTotal Bytes (included padding)\n\n",dummy_total);
    }


  Time= header1.time;
  Redshift= header1.time;
}








/* this routine loads particle data from a Gadget
 * fill in HDF5 format.
 */
int load_snapshot_hdf5(char *fname, int files, int quiteon, int headeronly)
{
  char   buf[200];
  int    i,k,ntot_withmasses;
  int    add_idx;
  int    dog;

  float *xyz, *x;
  int *id;

  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;
  hid_t hdf5_grp[6], hdf5_dataspace_in_file;
  hid_t hdf5_dataset;
  hsize_t dims[2];


  sprintf(buf,"%s",fname);

  hdf5_file = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header");

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, header1.npart);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header1.npartTotal);
  H5Aclose(hdf5_attribute);

/*
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header1.npartTotalHighWord);
  H5Aclose(hdf5_attribute);
*/

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header1.mass);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.time);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Redshift");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.redshift);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "BoxSize");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.BoxSize);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header1.num_files);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Omega0");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.Omega0);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "OmegaLambda");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.OmegaLambda);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "HubbleParam");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header1.HubbleParam);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Sfr");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header1.flag_sfr);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Cooling");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header1.flag_cooling);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_StellarAge");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header1.flag_stellarage);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Metals");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header1.flag_metals);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Feedback");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header1.flag_feedbacktp);
  H5Aclose(hdf5_attribute);

/*
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_DoublePrecision");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header1.flag_doubleprecision);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_IC_Info");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header1.flag_ic_info);
  H5Aclose(hdf5_attribute);
*/

  H5Gclose(hdf5_headergrp);


  for(k=0, NumPart=0, ntot_withmasses=0; k<6; k++)
	NumPart+= header1.npart[k];
  Ngas= header1.npart[0];

  if(SnapInfo) printf(" (NumPart=%10d\tNgas=%10d)\n",NumPart,Ngas);

  for(k=0, ntot_withmasses=0; k<6; k++)
    {
	if(header1.mass[k]==0)
	  ntot_withmasses+= header1.npart[k];
    }

  if(headeronly==1)
    {
	H5Fclose(hdf5_file);
        return(0);
    }


  /* Allocate memory for the temporary read-in arrays */
  allocate_memory();

  /* temporary array's used to write file */
  if(!(id = (int *) malloc((NumPart)*sizeof(int))))
    {
        printf("failed to allocate memory (id).\n");
        exit(0);
    }
  if(!(x = (float *) malloc((NumPart)*sizeof(float))))
    {
        printf("failed to allocate memory (x).\n");
        exit(0);
    }
  if(!(xyz = (float *) malloc((3*NumPart)*sizeof(float))))
    {
        printf("failed to allocate memory (xyz).\n");
        exit(0);
    }





  /* Now, let's open the actual particle information */
  for(k = 0; k < 6; k++)
  {
    if(header1.npart[k] > 0)
      {
        sprintf(buf, "/PartType%d", k);
        hdf5_grp[k] = H5Gopen(hdf5_file, buf);
      }
  }


  if(header1.npart[0] > 0)
    {
	if(quiteon==0) printf("reading gas\n"); fflush(stdout);

	for(i=1; i<=header1.npart[0]; i++)
		P[i].Type= 0;

	hdf5_dataset = H5Dopen(hdf5_grp[0], "Coordinates");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
	H5Dclose(hdf5_dataset);
	for(i=1; i<=header1.npart[0]; i++)
	  {
		P[i].Pos[0]= xyz[(i-1)*3 + 0];
		P[i].Pos[1]= xyz[(i-1)*3 + 1];
		P[i].Pos[2]= xyz[(i-1)*3 + 2];
	  }

	hdf5_dataset = H5Dopen(hdf5_grp[0], "Velocities");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
	H5Dclose(hdf5_dataset);
	for(i=1; i<=header1.npart[0]; i++)
	  {
		P[i].Vel[0]= xyz[(i-1)*3 + 0];
		P[i].Vel[1]= xyz[(i-1)*3 + 1];
		P[i].Vel[2]= xyz[(i-1)*3 + 2];
	  }


	hdf5_dataset = H5Dopen(hdf5_grp[0], "ParticleIDs");
	H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
	H5Dclose(hdf5_dataset);
	for(i=1; i<=header1.npart[0]; i++)
		Id[i]= id[i-1];

        if(header1.mass[0] > 0)
          for(i=1; i<=header1.npart[0]; i++)
                P[i].Mass= header1.mass[0];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[0], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header1.npart[0];i++)
                P[i].Mass= x[i-1];
          }

	hdf5_dataset = H5Dopen(hdf5_grp[0], "InternalEnergy");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
	H5Dclose(hdf5_dataset);
	for(i=1; i<=header1.npart[0]; i++)
		P[i].U= x[i-1];

	hdf5_dataset = H5Dopen(hdf5_grp[0], "Metallicity");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
	H5Dclose(hdf5_dataset);
    }



  if(header1.npart[1] > 0)
    {
        if(quiteon==0) printf("reading halo\n"); fflush(stdout);

	add_idx= header1.npart[0];

	for(i=1; i<=header1.npart[1]; i++)
		P[add_idx + i].Type= 1;

        hdf5_dataset = H5Dopen(hdf5_grp[1], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[1];i++)
          {
                P[add_idx + i].Pos[0]= xyz[(i-1)*3 + 0];
                P[add_idx + i].Pos[1]= xyz[(i-1)*3 + 1];
                P[add_idx + i].Pos[2]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[1], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[1];i++)
          {
                P[add_idx + i].Vel[0]= xyz[(i-1)*3 + 0];
                P[add_idx + i].Vel[1]= xyz[(i-1)*3 + 1];
                P[add_idx + i].Vel[2]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[1], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[1];i++)
                Id[add_idx + i]= id[i-1];

	if(header1.mass[1] > 0)
	  for(i=1; i<=header1.npart[1]; i++)
                P[add_idx + i].Mass= header1.mass[1];
	else
	  {
	    hdf5_dataset = H5Dopen(hdf5_grp[1], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header1.npart[1];i++)
                P[add_idx + i].Mass= x[i-1];
	  }
    }



  if(header1.npart[2] > 0)
    {
        if(quiteon==0) printf("reading disk\n"); fflush(stdout);

        add_idx= header1.npart[0] + header1.npart[1];

        for(i=1; i<=header1.npart[2]; i++)
                P[add_idx + i].Type= 2;

        hdf5_dataset = H5Dopen(hdf5_grp[2], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[2];i++)
          {
                P[add_idx + i].Pos[0]= xyz[(i-1)*3 + 0];
                P[add_idx + i].Pos[1]= xyz[(i-1)*3 + 1];
                P[add_idx + i].Pos[2]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[2], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[2];i++)
          {
                P[add_idx + i].Vel[0]= xyz[(i-1)*3 + 0];
                P[add_idx + i].Vel[1]= xyz[(i-1)*3 + 1];
                P[add_idx + i].Vel[2]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[2], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[2];i++)
                Id[add_idx + i]= id[i-1];

        if(header1.mass[2] > 0)
          for(i=1; i<=header1.npart[2]; i++)
                P[add_idx + i].Mass= header1.mass[2];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[2], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header1.npart[2];i++)
                P[add_idx + i].Mass= x[i-1];
          }
    }




  if(header1.npart[3] > 0)
    {
        if(quiteon==0) printf("reading bulge\n"); fflush(stdout);

        add_idx= header1.npart[0] + header1.npart[1] + header1.npart[2];

        for(i=1; i<=header1.npart[3]; i++)
                P[add_idx + i].Type= 3;

        hdf5_dataset = H5Dopen(hdf5_grp[3], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[3];i++)
          {
                P[add_idx + i].Pos[0]= xyz[(i-1)*3 + 0];
                P[add_idx + i].Pos[1]= xyz[(i-1)*3 + 1];
                P[add_idx + i].Pos[2]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[3], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[3];i++)
          {
                P[add_idx + i].Vel[0]= xyz[(i-1)*3 + 0];
                P[add_idx + i].Vel[1]= xyz[(i-1)*3 + 1];
                P[add_idx + i].Vel[2]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[3], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[3];i++)
                Id[add_idx + i]= id[i-1];

        if(header1.mass[3] > 0)
          for(i=1; i<=header1.npart[3]; i++)
                P[add_idx + i].Mass= header1.mass[3];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[3], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header1.npart[3];i++)
                P[add_idx + i].Mass= x[i-1];
          }
    }




  if(header1.npart[4] > 0)
    {
        if(quiteon==0) printf("reading new stars\n"); fflush(stdout);

        add_idx= header1.npart[0] + header1.npart[1] + header1.npart[2] + header1.npart[3];

        for(i=1; i<=header1.npart[4]; i++)
                P[add_idx + i].Type= 4;

        hdf5_dataset = H5Dopen(hdf5_grp[4], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[4];i++)
          {
                P[add_idx + i].Pos[0]= xyz[(i-1)*3 + 0];
                P[add_idx + i].Pos[1]= xyz[(i-1)*3 + 1];
                P[add_idx + i].Pos[2]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[4], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[4];i++)
          {
                P[add_idx + i].Vel[0]= xyz[(i-1)*3 + 0];
                P[add_idx + i].Vel[1]= xyz[(i-1)*3 + 1];
                P[add_idx + i].Vel[2]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[4], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[4];i++)
                Id[add_idx + i]= id[i-1];

        if(header1.mass[4] > 0)
          for(i=1; i<=header1.npart[4]; i++)
                P[add_idx + i].Mass= header1.mass[4];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[4], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header1.npart[4];i++)
                P[add_idx + i].Mass= x[i-1];
          }
    }




  if(header1.npart[5] > 0)
    {
        if(quiteon==0) printf("reading black holes\n"); fflush(stdout);

        add_idx= header1.npart[0] + header1.npart[1] + header1.npart[2] + header1.npart[3] + header1.npart[4];

        for(i=1; i<=header1.npart[5]; i++)
                P[add_idx + i].Type= 5;

        hdf5_dataset = H5Dopen(hdf5_grp[5], "Coordinates");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[5];i++)
          {
                P[add_idx + i].Pos[0]= xyz[(i-1)*3 + 0];
                P[add_idx + i].Pos[1]= xyz[(i-1)*3 + 1];
                P[add_idx + i].Pos[2]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[5], "Velocities");
        H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[5];i++)
          {
                P[add_idx + i].Vel[0]= xyz[(i-1)*3 + 0];
                P[add_idx + i].Vel[1]= xyz[(i-1)*3 + 1];
                P[add_idx + i].Vel[2]= xyz[(i-1)*3 + 2];
          }

        hdf5_dataset = H5Dopen(hdf5_grp[5], "ParticleIDs");
        H5Dread(hdf5_dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
        H5Dclose(hdf5_dataset);
        for(i=1; i<=header1.npart[5];i++)
                Id[add_idx + i]= id[i-1];

        if(header1.mass[5] > 0)
          for(i=1; i<=header1.npart[5]; i++)
                P[add_idx + i].Mass= header1.mass[5];
        else
          {
            hdf5_dataset = H5Dopen(hdf5_grp[5], "Masses");
            H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
            H5Dclose(hdf5_dataset);
            for(i=1; i<=header1.npart[5];i++)
                P[add_idx + i].Mass= x[i-1];
          }
    }




  /* close this thing */
  for(k = 5; k >= 0; k--)
    if(header1.npart[k] > 0)
      H5Gclose(hdf5_grp[k]);
  H5Fclose(hdf5_file);


}










/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  /*printf("allocating memory...\n");*/

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      fprintf(stderr,"NumPart = %d\n",NumPart);
      fprintf(stderr,"sizeof = %d\n",sizeof(struct particle_data));
      fprintf(stderr,"\n\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      fprintf(stderr,"NumPart = %d\n",NumPart);
      fprintf(stderr,"sizeof = %d\n",sizeof(int));
      fprintf(stderr,"\n\n");
      exit(0);
    }
  
  Id--;   /* start with offset 1 */

  /*printf("allocating memory...done\n");*/
}




/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(void)
{
  int i,j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i=1; i<=NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource= P[i];
	  idsource=Id[i];
	  dest=Id[i];

	  do
	    {
	      psave= P[dest];
	      idsave=Id[dest];

	      P[dest]= psource;
	      Id[dest]= idsource;
	      
	      if(dest == i) 
		break;

	      psource= psave;
	      idsource=idsave;

	      dest=idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");

  Id++;   
  free(Id);

  printf("space for particle ID freed\n");
}






  











