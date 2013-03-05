#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "overhead.h"
#include "proto.h" 

#define DEBUG NO
#define HEAD YES

void readsnap(char *filename)
{
  int    i,dummy,ntot_withmasses;
  int    t,n,off; 
  char   buf[200];
  extern struct snap_header header; 
  FILE *fd;


  sprintf(buf,"%s",filename);

  if(!(fd=fopen(buf,"r"))) {
    printf("Can't open snapshot file %s\n",buf);
    exit(0);
  }


  /* Read in the header */
  SKIP; 
  fread(&header, sizeof(header), 1, fd);
  SKIP; 

  All.N_gas = header.npart[0]; 
  All.N_halo = header.npart[1]; 
  All.N_disk = header.npart[2]; 
  All.N_bulge = header.npart[3]; 
  All.N_star = header.npart[4]; 
  All.N_bh = header.npart[5];
  for(i=0,All.N_total=0; i<6; i++) All.N_total+=header.npart[i]; 

  if (DEBUG) { 
    fprintf(stderr, "readsnap(): Ngas = %d\n", All.N_gas); 
    fprintf(stderr, "readsnap(): Nstar = %d\n", All.N_star); 
    fprintf(stderr, "readsnap(): Ntotal = %d\n", All.N_total); 
    fprintf(stderr, "readsnap(): Time = %f\n", header.time); 
  } 

  if (HEAD) printhead(); 

  /* allocate memory for the different particle types */ 
  allocate_gas(); 
  allocate_halo(); 
  allocate_disk(); 
  allocate_bulge(); 
  allocate_star();
  allocate_bh();
  int Ngas=All.N_gas,Nhalo=All.N_halo,Ndisk=All.N_disk;
  int Nbulge=All.N_bulge,Nstar=All.N_star,Nbh=All.N_bh,Ntotal=All.N_total; 
  if (DEBUG) fprintf(stderr, "Memory successfully allocated\n"); 

  /* read in the position info for all particle types */
  SKIP;
  if (Ngas > 0) for(n=0;n<Ngas;n++) fread(&PG[n].pos[0], sizeof(float), 3, fd);
  if (Nhalo > 0) for(n=0;n<Nhalo;n++) fread(&PH[n].pos[0], sizeof(float), 3, fd);
  if (Ndisk > 0) for(n=0;n<Ndisk;n++) fread(&PD[n].pos[0], sizeof(float), 3, fd);
  if (Nbulge > 0) for(n=0;n<Nbulge;n++) fread(&PB[n].pos[0], sizeof(float), 3, fd);
  if (Nstar > 0) for(n=0;n<Nstar;n++) fread(&PS[n].pos[0], sizeof(float), 3, fd);
  if (Nbh > 0) for(n=0;n<Nbh;n++) fread(&PBH[n].pos[0], sizeof(float), 3, fd);
  SKIP;

  /* read in the velocity info for all particle types */
  SKIP;
  if (Ngas > 0) for(n=0;n<Ngas;n++) fread(&PG[n].vel[0], sizeof(float), 3, fd);
  if (Nhalo > 0) for(n=0;n<Nhalo;n++) fread(&PH[n].vel[0], sizeof(float), 3, fd);
  if (Ndisk > 0) for(n=0;n<Ndisk;n++) fread(&PD[n].vel[0], sizeof(float), 3, fd);
  if (Nbulge > 0) for(n=0;n<Nbulge;n++) fread(&PB[n].vel[0], sizeof(float), 3, fd);
  if (Nstar > 0) for(n=0;n<Nstar;n++) fread(&PS[n].vel[0], sizeof(float), 3, fd);
  if (Nbh > 0) for(n=0;n<Nbh;n++) fread(&PBH[n].vel[0], sizeof(float), 3, fd);
  SKIP;
  

  if (DEBUG) fprintf(stderr, "Finished reading in position and velocity information\n"); 
    
  /* read in ID info */ 
  SKIP;
  if (Ngas > 0) for(n=0;n<Ngas;n++) fread(&PG[n].ID, sizeof(int), 1, fd);
  if (Nhalo > 0) for(n=0;n<Nhalo;n++) fread(&PH[n].ID, sizeof(int), 1, fd);
  if (Ndisk > 0) for(n=0;n<Ndisk;n++) fread(&PD[n].ID, sizeof(int), 1, fd);
  if (Nbulge > 0) for(n=0;n<Nbulge;n++) fread(&PB[n].ID, sizeof(int), 1, fd);
  if (Nstar > 0) for(n=0;n<Nstar;n++) fread(&PS[n].ID, sizeof(int), 1, fd);
  if (Nbh > 0) for(n=0;n<Nbh;n++) fread(&PBH[n].ID, sizeof(int), 1, fd);
  SKIP;


  /* read in mass data, or assign if particle type all has same mass */
  SKIP;
  if (header.mass[0]==0.) for(n=0;n<Ngas;n++) fread(&PG[n].mass, sizeof(float), 1, fd);
  if (header.mass[1]==0. && Nhalo > 0) fprintf(stderr, "Error, need mass added to PH\n"); 
  if (header.mass[2]==0. && Ndisk > 0) fprintf(stderr, "Error, need mass added to PD\n");
  if (header.mass[3]==0. && Nbulge > 0) fprintf(stderr, "Error, need mass added to PB\n");
  if (header.mass[4]==0.) for(n=0;n<Nstar;n++) fread(&PS[n].mass, sizeof(float), 1, fd);
  if (header.mass[5]==0.) for(n=0;n<Nbh;n++) fread(&PBH[n].mass, sizeof(float), 1, fd);
  SKIP;

  if (Nhalo > 0) All.M_halo   = header.mass[1];
  if (Ndisk > 0) All.M_disk   = header.mass[2];
  if (Nbulge > 0) All.M_bulge = header.mass[3];

  if (DEBUG) fprintf(stderr, "Finished reading in ID and mass information\n"); 

  /* read in data for SPH particles */ 
  if(Ngas>0) {

    /* U (internal energy per particle per unit mass (check)*/ 
    SKIP;
    for(n=0; n<Ngas;n++) fread(&PG[n].u, sizeof(float), 1, fd);
    SKIP;

    /* rho (density in problem units) */ 
    SKIP;
    for(n=0; n<Ngas;n++) fread(&PG[n].rho, sizeof(float), 1, fd);
    SKIP;

  /* 
   *  if flag_cooling, read in additional info for SPH particles 
   */ 

    /* Ne (number of free electrons per hydrogen atom/nucleus) */ 
    if(header.flag_cooling) {
      SKIP;
      for(n=0; n<Ngas;n++) fread(&PG[n].ne, sizeof(float), 1, fd);
      SKIP;
    } else { 
      for(n=0; n<Ngas;n++) PG[n].ne= 0.0;
    }
  
    /* nh (neutral fraction of hydrogen atoms) */ 
    if(header.flag_cooling) {
      SKIP;
      for(n=0; n<Ngas;n++) fread(&PG[n].nh, sizeof(float), 1, fd);
      SKIP;
    } else { 
      for(n=0; n<Ngas;n++) PG[n].nh= 0.0;
    }

    /* hsml (smoothing length for each particle) */ 
    SKIP;
    for(n=0; n<Ngas;n++) fread(&PG[n].hsml, sizeof(float), 1, fd);
    SKIP;

/* 
 *  if flag_sfr, read in additional info for SPH and star particles 
 */ 

    /* sfr (current star formation rate in the SPH particle) */ 
    if(header.flag_sfr) {
      SKIP;
      for(n=0; n<Ngas;n++) fread(&PG[n].sfr, sizeof(float), 1, fd);
      SKIP;
    } else { 
      for(n=0; n<Ngas;n++) PG[n].sfr= 0.0;
    }
  } 

  /* age of stars */ 
  if (Nstar > 0) {
    if(header.flag_sfr) {
      SKIP;
      for(n=0; n<Nstar;n++) fread(&PS[n].age, sizeof(float), 1, fd);
      SKIP;
    } else { 
      for(n=0; n<Nstar;n++) PS[n].age= 0.0;
    }
  } 

  /* Z (metallicity) of gas */ 
  if (Ngas > 0) {
    if(header.flag_sfr) {
      SKIP;
      for(n=0; n<Ngas;n++) fread(&PG[n].z, sizeof(float), 1, fd);
      SKIP;
    } else { 
      for(n=0; n<Ngas;n++) PG[n].z= 0.0;
    }
  }

  /* Z (metallicity) of stars */ 
  if (Nstar > 0) {
    if(header.flag_sfr) {
      SKIP;
      for(n=0; n<Nstar;n++) fread(&PS[n].z, sizeof(float), 1, fd);
      SKIP;
    } else { 
      for(n=0; n<Nstar;n++) PS[n].z= 0.0;
    }
  }
  if (DEBUG) fprintf(stderr, "Finished reading in gas density info information\n"); 

  fclose(fd);

  All.Time= header.time;
  All.TimeBegin= header.time;
  All.Redshift= header.time;

}
