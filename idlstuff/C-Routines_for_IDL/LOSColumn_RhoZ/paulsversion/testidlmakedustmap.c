#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "overhead.h"
#include "proto.h"

#define DEBUG NO
#define NAMELEN 60 

int NUMBER_OF_RAYS;

/* int getnh(int argc, char **argv) */
int getnh(int argc, char *argv[]) 
{

  int i, j, n_gas, ray_num, N_rays, pid; 
  int RAY_ORIGIN_CELL; 
  int first_bh_ray, first_star_ray, first_disk_ray, first_bulge_ray; 
  int *oi, *ot; 
  char *simfile; 
  char *outfile; 
  int *OUT_ID; 		/* output particle ID */ 
  int *OUT_TYPE;	/* output particle type */ 
  float *OUT_NH; 	/* output NH */ 
  float *coords; 	/* expected input from TJ */ 
  float *particle_nh; 	/* output */ 
  float theta, phi; 
  float nh_prefactor;
  float *on; 
  Ray_struct *ray;

/*
  FILE *op; 
  simfile = (char*) malloc(sizeof(char)*NAMELEN); 
  outfile = (char*) malloc(sizeof(char)*NAMELEN*2); 
*/

  if (argc != 7) {
    fprintf(stderr, "Wrong number of arguments, found %d\n",argc); 
    exit(0); 
  }

/*
  sprintf(argv[0], "/raid4/tcox/vc3vc3b/snapshot_010"); 
  simfile = (char*)argv[5];

  sprintf(simfile, "/raid4/tcox/vc3vc3b/snapshot_010"); 
  printf("This is simfile: %s\n", simfile); 

  printf("This is argv[5]: %s\n", argv[5]); 
  sprintf(simfile, "%s",argv[5]);

  theta = atof(argv[1]); 
  phi = atof(argv[2]); 
*/

  All.N_gas = *(int*)argv[0];  
  All.N_star = *(int*)argv[1];  
  All.N_bh = *(int*)argv[2];  
  theta = *(float *)argv[3]; 
  phi = *(float *)argv[4]; 

  /* returned quantities */ 
/*
  OUT_ID = (int*)argv[2]; 
  OUT_TYPE = (int*)argv[3]; 
  OUT_NH = (float*)argv[4]; 
*/

  printf("Ngas, Nstars, Nbh = %4d %4d %4d\n",All.N_gas, All.N_star, All.N_bh); 
  printf("theta, phi = %4.2f %4.2f\n",theta,phi); 

/* allocate memory */ 

  coords = (float*) f2alloc(9, All.N_gas+All.N_star+All.N_bh); 
  particles_nh = (float*) f1alloc(All.N_star+All.N_bh+All.N_gas); 

  /* check that we don't go try to integrate exactly along cell boundaries */

  if (sin(theta) == 0.) theta += 1.e-30; 
  if (sin(phi) == 0.) phi += 1.e-30; 

  /* parse the snapshot data from dustmap.pro */ 





  /* setup the output file */ 
/*
  sprintf(outfile,"%s_%.1f_%.1f_nh.bin",argv[1], theta, phi);
  if (DEBUG) fprintf(stderr, "outfile = %s\n", outfile); 

  if ((op = fopen(outfile, "wb")) == NULL) 
    fprintf(stderr, "Error: Cannot open %s\n", "outfile");	
*/

/*
  fprintf(stdout, "Loading snapshot %s... \n",simfile);
  readsnap(simfile);
  fprintf(stdout, "Snapshot %s loaded... \n",simfile);
*/

  if (All.N_bh < 1) { 
    fprintf(stdout,"Warning: No black holes in %s! \n", simfile);
  }

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
  printf("NUMBER_OF_RAYS = %d, N_rays = %d\n", NUMBER_OF_RAYS, N_rays); 
  printf("Initializing (%i) Rays... \n",N_rays);
/*
  ray = calloc(N_rays,sizeof(Ray_struct));
*/
  if(!(ray = malloc(N_rays*sizeof(struct Ray_s))))
     {
	printf("failed to allocate memory: N_rays= %d  size= %d\n",N_rays,N_rays*sizeof(struct Ray_s));
	exit(0);
     }
  printf("Done allocating memory for rays\n");
  ALL_CELL_COUNTER = 0;
  nh_prefactor = (All.UnitMass_in_g/PROTONMASS)/(SQR(All.UnitLength_in_cm));
  printf("nh_prefactor = %g\n", nh_prefactor); 

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

/*
  OUT_ID = (int*) calloc(NUMBER_OF_RAYS, sizeof(int)); 
  OUT_TYPE = (int*) calloc(NUMBER_OF_RAYS, sizeof(int)); 
  OUT_NH = (float*) calloc(NUMBER_OF_RAYS, sizeof(float)); 
 */

  /* rays from BH particles */ 
  for (ray_num=first_bh_ray,pid=0; ray_num<first_star_ray; ray_num++,pid++) { 
    ray[ray_num].ORIGIN_TYPE = PTYPE_BH; 
    ray[ray_num].ORIGIN_ID = PBH[pid].ID; 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PBH[pid].pos[j]; 
  } 
  /* rays from star particles */ 
  for (ray_num=first_star_ray,pid=0; ray_num<first_disk_ray; ray_num++,pid++) { 
    ray[ray_num].ORIGIN_TYPE = PTYPE_STAR; 
    ray[ray_num].ORIGIN_ID = PS[pid].ID; 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PS[pid].pos[j]; 
  } 
  /* rays from disk particles */ 
  for (ray_num=first_disk_ray,pid=0; ray_num<first_bulge_ray; ray_num++,pid++) { 
    ray[ray_num].ORIGIN_TYPE = PTYPE_DISK; 
    ray[ray_num].ORIGIN_ID = PD[pid].ID; 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PD[pid].pos[j]; 
  } 
  /* rays from bulge particles */ 
  for (ray_num=first_bulge_ray,pid=0; ray_num<NUMBER_OF_RAYS; ray_num++,pid++) { 
    ray[ray_num].ORIGIN_TYPE = PTYPE_BULGE; 
    ray[ray_num].ORIGIN_ID = PB[pid].ID; 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PB[pid].pos[j]; 
  } 
  /* initialize the remaining rays */ 
  for (ray_num=NUMBER_OF_RAYS; ray_num<N_rays; ray_num++) { 
    ray[ray_num].ORIGIN_TYPE = -1; 
    ray[ray_num].ORIGIN_ID = -1; 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = 0.; 
    ray[ray_num].theta = 0.; 
    ray[ray_num].phi = 0.; 
  } 
  printf("Done setting up the ray IDs and positions\n"); 

  /* 
   *  now perform the calculations 
   */ 

  for (ray_num=0; ray_num<NUMBER_OF_RAYS; ray_num++) { 

    if(!(ray_num%10000))
      {
	printf("%d..",ray_num);
	fflush(stdout);
      }

/*
    if(ray_num>163420)
      {
	printf("starting %d\n",ray_num);
	fflush(stdout);
      }
*/

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

/*
    if(ray_num>163420)
      {
	printf("ray_num %d origin_cell %d origin_type %d origin_id %d theta %.1f phi %.1f nh %.3f Z %f\n",ray_num, RAY_ORIGIN_CELL, ray[ray_num].ORIGIN_TYPE, ray[ray_num].ORIGIN_ID,ray[ray_num].theta,ray[ray_num].phi,log10(ray[ray_num].nh_hot),ray[ray_num].Z);
	fflush(stdout);
      }
*/

/*
    fwrite(&ray[ray_num].ORIGIN_TYPE,sizeof(int),1,op);
    fwrite(&ray[ray_num].ORIGIN_ID,sizeof(int),1,op);
    fwrite(&ray[ray_num].nh_hot,sizeof(float),1,op);
*/ 
  
    /* write the results to the input arrays */ 

    oi = OUT_ID + ray_num; 
    *oi = ray[ray_num].ORIGIN_ID; 
    ot = OUT_TYPE + ray_num; 
    *ot = ray[ray_num].ORIGIN_TYPE; 
    on = OUT_NH + ray_num; 
    *on = ray[ray_num].nh_hot; 

/*
    OUT_ID[ray_num] = ray[ray_num].ORIGIN_ID; 
    OUT_TYPE[ray_num] = ray[ray_num].ORIGIN_TYPE; 
    OUT_NH[ray_num] = ray[ray_num].nh_hot; 
*/

/*
    if(ray_num>163420)
      {
	printf("(%d) %d %d %f -- %d %d %f -- %d %d %f\n", ray_num, OUT_ID[ray_num], OUT_TYPE[ray_num], OUT_NH[ray_num], *oi, *ot, *on, ray[ray_num].ORIGIN_ID, ray[ray_num].ORIGIN_TYPE, ray[ray_num].nh_hot); 
	fflush(stdout);
      }
*/

  }


/*
  fclose(op); 
*/
  printf("Free memory\n");
  free(ray);

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
