#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "overhead.h"
#include "proto.h"

#define DEBUG NO
#define NAMELEN 60 

int NUMBER_OF_RAYS;

int Ngas,Nstar,Nbh;

float *Pos;

struct particle_3d
{
  float Pos[3];
  float u;
  float rho;
  float hsml;
  float numh;
  float nume;
  float z;
} **P3d;



/* int getnh(int argc, char **argv) */
int getnh(int argc, void *argv[]) 
{

  int i, j, n_gas, ray_num, N_rays, pid, tmpval; 
  int RAY_ORIGIN_CELL; 
  int first_bh_ray, first_star_ray, first_disk_ray, first_bulge_ray; 
  int *oi, *ot; 
  float *OUT_NH; 	/* output NH */ 
  float theta, phi; 
  float nh_prefactor;
  float *on; 
  Ray_struct *ray;
  FILE *op; 

  void allocate_3d(void);
  void set_particle_pointer(void);
  void free_memory_3d(void);

  if (argc != 7) {
    fprintf(stderr, "Expected 7 arguments (found %d)\n",argc); 
    exit(0); 
  }

  All.N_gas = *(int *)argv[0];
  All.N_star = *(int *)argv[1];
  All.N_bh = *(int *)argv[2];
  theta = *(float *)argv[3]; 
  phi = *(float *)argv[4]; 
  Pos = (float *)argv[5];
  /* returned quantities */ 
  OUT_NH = (float *)argv[6]; 


  All.N_halo = 0; 
  All.N_disk = 0; 
  All.N_bulge = 0; 
  All.N_total = All.N_gas + All.N_star + All.N_bh; 
  Ngas = All.N_gas; 
  Nstar = All.N_star; 
  Nbh = All.N_bh; 

  fprintf(stderr, "theta, phi = %.2f %.2f\n", theta, phi); 
  fprintf(stderr, "Ngas = %d Nstar = %d Nbh = %d Ntotal = %d\n", Ngas, Nstar, Nbh, All.N_total);

  if (All.N_bh < 1) { 
    fprintf(stdout,"Warning: No black holes! \n"); 
  }

  /* check that we don't go try to integrate exactly along cell boundaries */

  if (sin(theta) == 0.) theta += 1.e-30; 
  if (sin(phi) == 0.) theta += 1.e-30; 

  allocate_3d();
  set_particle_pointer();
  allocate_gas(); 
  allocate_star(); 
  allocate_bh(); 

  fprintf(stderr, "star, gas, and bh allocation done too\n"); 

  /* --------------------------------------------------

     OK, now we're ready to use whatever variables we
     need.

     We can access the various fields via:

     P3d[i+1]->Pos[0], or P3d[30]->rho= 100.0

     -------------------------------------------------- */

  printf("\n\n");
  printf("Here are some of the P3d variables\n");
  printf("P3d[1]->Pos[0,1,2]= %g|%g|%g\n",P3d[1]->Pos[0],P3d[1]->Pos[1],P3d[1]->Pos[2]);
  printf("P3d[3003]->u= %g\n",P3d[3003]->u);
  printf("P3d[145073]->Pos[0]= %g\n",P3d[145073]->Pos[0]);
  printf("\n\n");


  /* Fill in the structures that are normally done in readsnap */ 
  
  for (i=0; i<All.N_gas; i++) {
    PG[i].pos[0] = P3d[i+1]->Pos[0]; 
    PG[i].pos[1] = P3d[i+1]->Pos[1];
    PG[i].pos[2] = P3d[i+1]->Pos[2];
    PG[i].u = P3d[i+1]->u; 
    PG[i].rho = P3d[i+1]->rho; 
    PG[i].hsml = P3d[i+1]->hsml; 
    PG[i].nh = P3d[i+1]->numh; 
    PG[i].ne = P3d[i+1]->nume; 
    PG[i].z = P3d[i+1]->z; 
  } 
  tmpval = All.N_gas+All.N_star; 
  for (i=All.N_gas,pid=0; i<tmpval; i++,pid++) {
    PS[pid].pos[0] = P3d[i+1]->Pos[0]; 
    PS[pid].pos[1] = P3d[i+1]->Pos[1];
    PS[pid].pos[2] = P3d[i+1]->Pos[2];
  }
  for (i=tmpval,pid=0; i<All.N_total; i++,pid++) {
    PBH[pid].pos[0] = P3d[i+1]->Pos[0]; 
    PBH[pid].pos[1] = P3d[i+1]->Pos[1];
    PBH[pid].pos[2] = P3d[i+1]->Pos[2];
  }

  fprintf(stderr, "\n"); 
  fprintf(stderr, "Quick test: \n"); 
  fprintf(stderr, "PG[0].pos = %g %g %g\n", PG[0].pos[0], PG[0].pos[1], PG[0].pos[2]); 
  fprintf(stderr, "PG[1].pos = %g %g %g\n", PG[1].pos[0], PG[1].pos[1], PG[1].pos[2]); 
  fprintf(stderr, "PG[3003].pos = %g %g %g PG[3003].u = %g\n", PG[3003].pos[0], PG[3003].pos[1], PG[3003].pos[2], PG[3003].u); 
  fprintf(stderr, "PS[0].pos = %g %g %g\n", PS[0].pos[0], PS[0].pos[1], PS[0].pos[2]); 
  fprintf(stderr, "PS[1].pos = %g %g %g\n", PS[1].pos[0], PS[1].pos[1], PS[1].pos[2]); 
  fprintf(stderr, "PBH[0].pos = %g %g %g\n", PBH[0].pos[0], PBH[0].pos[1], PBH[0].pos[2]); 
  fprintf(stderr, "PBH[1].pos = %g %g %g\n", PBH[1].pos[0], PBH[1].pos[1], PBH[1].pos[2]); 

  /* free_memory_3d();  */

  /* 
   *  initialize the tree 
   */

  setup_lineofsight_overhead();
  initialize_the_tree();
  fprintf(stdout, "tree built... \n"); 

  /* 
   *  initialize the rays 
   */ 

  NUMBER_OF_RAYS  = (All.N_bh + All.N_star + All.N_disk + All.N_bulge);
  N_rays = NUMBER_OF_RAYS;
  N_rays = (int)(sqrt(N_rays)+1) * (int)(sqrt(N_rays)+1);		/* Convert to a perfect square */
  fprintf(stderr, "NUMBER_OF_RAYS = %d, N_rays = %d, sizeof(Ray_struct) = %d\n", NUMBER_OF_RAYS, N_rays, sizeof(Ray_struct) ); 
  printf("Initializing (%i) Rays... \n",N_rays);
  /* ray = calloc(N_rays,sizeof(Ray_struct)); */
  ray = calloc(N_rays,sizeof(float));
  if (DEBUG) fprintf(stderr, "Done allocating memory for rays\n"); 
  ALL_CELL_COUNTER = 0;
  nh_prefactor = (All.UnitMass_in_g/PROTONMASS)/(SQR(All.UnitLength_in_cm));

  /* 
   *  set up the rays 
   */ 

  first_bh_ray = 0; 
  first_star_ray = first_bh_ray + All.N_bh; 
  first_disk_ray = first_star_ray + All.N_star; 
  first_bulge_ray = first_disk_ray + All.N_disk; 

/*
 * allocate memory for the IDL arrays 
 */


  /* rays from BH particles */ 
  for (ray_num=first_bh_ray,pid=0; ray_num<first_star_ray; ray_num++,pid++) { 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PBH[pid].pos[j]; 
  } 
  /* rays from star particles */ 
  for (ray_num=first_star_ray,pid=0; ray_num<first_disk_ray; ray_num++,pid++) { 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PS[pid].pos[j]; 
  } 
  /* rays from disk particles */ 
  for (ray_num=first_disk_ray,pid=0; ray_num<first_bulge_ray; ray_num++,pid++) { 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PD[pid].pos[j]; 
  } 
  /* rays from bulge particles */ 
  for (ray_num=first_bulge_ray,pid=0; ray_num<NUMBER_OF_RAYS; ray_num++,pid++) { 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PB[pid].pos[j]; 
  } 
  /* initialize the remaining rays */ 
  for (ray_num=NUMBER_OF_RAYS; ray_num<N_rays; ray_num++) { 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = 0.; 
    ray[ray_num].theta = 0.; 
    ray[ray_num].phi = 0.; 
  } 
  if (DEBUG) fprintf(stderr, "Done setting up the ray IDs and positions\n"); 

  /* 
   *  now perform the calculations 
   */ 

  for (ray_num=0; ray_num<NUMBER_OF_RAYS; ray_num++) { 
    ray[ray_num].theta = theta; 
    ray[ray_num].phi = phi; 
    ray[ray_num].n_hat[0] = sin(ray[ray_num].theta) * cos(ray[ray_num].phi);
    ray[ray_num].n_hat[1] = sin(ray[ray_num].theta) * sin(ray[ray_num].phi);
    ray[ray_num].n_hat[2] = cos(ray[ray_num].theta);
    ray[ray_num].nh=1.e-40;
    ray[ray_num].nh_hot=1.e-40;
    ray[ray_num].Z=0.; 
    ray[ray_num].neutral_frac=0.; 
    RAY_ORIGIN_CELL = find_cell_from_scratch(ray[ray_num].pos); 
    integrate_ray_to_escape(&ray[ray_num],RAY_ORIGIN_CELL);
    ray[ray_num].Z /= ray[ray_num].nh;
    ray[ray_num].neutral_frac /= ray[ray_num].nh;
    ray[ray_num].nh   *= nh_prefactor;
    ray[ray_num].nh_hot *= nh_prefactor;

    if (DEBUG) fprintf(stderr, "ray_num %d origin_cell %d origin_type %d origin_id %d theta %.1f phi %.1f nh %.3f Z %f\n",ray_num, RAY_ORIGIN_CELL, ray[ray_num].ORIGIN_TYPE, ray[ray_num].ORIGIN_ID,ray[ray_num].theta,ray[ray_num].phi,log10(ray[ray_num].nh_hot),ray[ray_num].Z);

/*
    fwrite(&ray[ray_num].nh_hot,sizeof(float),1,op);
*/ 
  
    /* write the results to the input arrays */ 

    on = OUT_NH + ray_num; 
    *on = ray[ray_num].nh_hot; 

/*
    OUT_NH[ray_num] = ray[ray_num].nh_hot; 
*/

/*
    if (DEBUG) fprintf(stderr, "(%d) %d %d %f -- %d %d %f -- %d %d %f\n", ray_num, OUT_ID[ray_num], OUT_TYPE[ray_num], OUT_NH[ray_num], *oi, *ot, *on, ray[ray_num].ORIGIN_ID, ray[ray_num].ORIGIN_TYPE, ray[ray_num].nh_hot); 
*/

  }


  fclose(op); 

  printf("\n Finish call to getnh.so ! \n");

  return 0; 

}


/* grab-bag for random initializations, etc (as 'begrun' from main code) */
void setup_lineofsight_overhead(void)
{
	set_sph_kernel();
	InitCool();
	set_All_struct_terms();
}


/* Generate the ray angles (generated such that they cover solid angle 
 *  equally, trivial to switch around to focus on particular direction)
 */
void generate_ray_angles(Ray_struct *R, int N_rays)
{
  int i0, theta_i, phi_i, n_ang_max = (int)(sqrt(N_rays));
  float theta, costheta, phi;
  float d_costheta = 2.0 / ((float)n_ang_max);
  float d_phi = 2.0 * PI / ((float)n_ang_max);
  
  /* Move in increments of cos(theta) and phi */
  i0 = 0;
  for (theta_i = 0; theta_i < n_ang_max; theta_i++) {
    costheta = 1.0 - (0.5 + (float)theta_i) * d_costheta;
      if (costheta < -1.0) costheta = -1.0;
      if (costheta >  1.0) costheta =  1.0;
    theta = acos(costheta);
    for (phi_i = 0; phi_i < n_ang_max; phi_i++) {
      phi = (0.5 + (float)phi_i) * d_phi;
      
      R[i0].theta = theta;
      R[i0].phi = phi;  
	  i0++;
  }}
}






void set_particle_pointer(void)
{
  int i;
  float *pos;

  for(i=1,pos=Pos;i<=Ngas;i++) 
    {
      
      P3d[i]->Pos[0] = pos[0];
      P3d[i]->Pos[1] = pos[1];
      P3d[i]->Pos[2] = pos[2];
  
      P3d[i]->u = pos[3];
      P3d[i]->rho = pos[4];
      P3d[i]->hsml = pos[5];
      P3d[i]->numh = pos[6];
      P3d[i]->nume = pos[7];
      P3d[i]->z = pos[8];

      pos+=9;
    }

  for(i=Ngas+1,pos=Pos;i<=Nstar+Nbh;i++)
    {

      P3d[i]->Pos[0] = pos[0];
      P3d[i]->Pos[1] = pos[1];
      P3d[i]->Pos[2] = pos[2];

      P3d[i]->u = 0.0;
      P3d[i]->rho = 0.0;
      P3d[i]->hsml = 0.0;
      P3d[i]->numh = 0.0;
      P3d[i]->nume = 0.0;
      P3d[i]->z = 0.0;

      pos+=9;
    }
}

void allocate_3d(void)
{
  int i;
  int Nsize;

  Nsize=Ngas+Nstar+Nbh;

  printf("allocating memory...\n");

  if(Nsize>0)
    {
      if(!(P3d=malloc(Nsize*sizeof(struct particle_3d *))))
        {
          printf("failed to allocate memory. (A)\n");
          exit(0);
        }

      P3d--;   /* start with offset 1 */

      if(!(P3d[1]=malloc(Nsize*sizeof(struct particle_3d)))) 
        {
          printf("failed to allocate memory. (B)\n");
          exit(0);
        }

      for(i=2;i<=Nsize;i++)   /* initiliaze pointer table */
        P3d[i]=P3d[i-1]+1;

    }

  printf("allocating memory...done\n");
}

void free_memory_3d(void)
{ 
  int i=0, Nsize; 
  Nsize=Ngas+Nstar+Nbh;
  while (i < Nsize) {
    free(P3d[i]);
    i++;
  }
  free(P3d);
}
