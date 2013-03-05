#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "overhead.h"
#include "proto.h"
#include "ngbtree3d.h"

#define NAMELEN 60 
float *Pos;
int Ngas,Nstar,Nbh;
struct particle_3d
{
  float Pos[3];
  float Vel[3];
  float Mass;
  float u;
  float rho;
  float hsml;
  float numh;
  float nume;
  float z;
} **P3d;


/* grab-bag for random initializations, etc (as 'begrun' from main code) */
void setup_lineofsight_overhead(void)
{
	set_sph_kernel();
	InitCool();
	set_All_struct_terms();
}

int main()
{
	int N_rays = (int)NUMBER_OF_RAYS;		/* Number of rays to use to trace the structure */
	N_rays = (int)(sqrt(N_rays)) * (int)(sqrt(N_rays));		/* Convert to a perfect square */

	/* Load the simulation snapshot and check for black holes */	
	char sim_file[100];
	printf("Snapshot filename: ");
	scanf("%s",sim_file);
	printf("Loading snapshot... \n");
	readsnap(sim_file);
	if(All.N_bh < 1)
    {
      fprintf(stderr,"No black holes! \n");
      exit(0);
    }
    free(PH);
    setup_lineofsight_overhead();
    
    
	void allocate_3d(void);
    void set_particle_pointer(void);
    void free_memory_3d(void);
    Ngas = All.N_gas; 
    Nstar = All.N_star; 
    Nbh = All.N_bh; 
    printf("Ngas = %d \n", Ngas);
    allocate_3d();
    set_particle_pointer();  
    float dummy[3];

	/* Build the tree-cell system to sort the gas properties */
	printf("Using: FULL_NEIGHBOR_CALC\n");
	if (USE_CELL_GRADIENTS==0) printf("USE_CELL_GRADIENTS is turned OFF\n");	
	if (USE_CELL_GRADIENTS) printf("USE_CELL_GRADIENTS is turned ON\n");	
	ngb3d_treeallocate(All.N_gas, 10*All.N_gas);
	ngb3d_treebuild((float **)&P3d[1], All.N_gas, 0,dummy,dummy);
	allocate_ngblists();
	initialize_the_tree();


	/* Get the info for the output files */
	int p_file_chk,j;
	float dummf0 = 0.0;
	/*
	char out_file_h[100];
		printf("\nFilename for header information: ");
		scanf("%s",out_file_h);
		printhead_file(out_file_h);	
	*/
	char out_file_2[100];
		printf("\nFilename for integrated values (NH, etc): ");
		scanf("%s",out_file_2);	
		FILE *output_file;
    	if ((output_file = fopen(out_file_2, "wb")) == NULL) 
    		fprintf(stderr, "Cannot open %s\n", "output_file");	
	int dumfcount;


	/* Create & allocate memory for the ray list */
	printf("Initializing (%i) Rays... \n",N_rays);
	int byte,ray_num;
	Ray_struct *ray;
    ray = calloc(N_rays,sizeof(Ray_struct));
	generate_ray_angles(ray,N_rays);
	for (ray_num=0; ray_num < N_rays; ray_num++) 
	{
		ray[ray_num].n_hat[0] = sin(ray[ray_num].theta) * cos(ray[ray_num].phi);
		ray[ray_num].n_hat[1] = sin(ray[ray_num].theta) * sin(ray[ray_num].phi);
		ray[ray_num].n_hat[2] = cos(ray[ray_num].theta);
	}
	ALL_CELL_COUNTER = 0;
	int RAY_ORIGIN_ID, RAY_ORIGIN_CELL, ray_origin_pID;
	int rayorigintype = RAY_ORIGIN_TYPE;
	float origin_pos[3];
	float nh_prefactor = (All.UnitMass_in_g/PROTONMASS)/(SQR(All.UnitLength_in_cm));
	int number_of_ray_origins = header.npart[RAY_ORIGIN_TYPE];
	printf("Beginning to process %i ray sources...\n",number_of_ray_origins);
	
	/* Print some of the overhead info to the binary NH output file */
	printf(" Simulation time at this snapshot = %f \n",header.time);
	float timetemp = header.time;
	fwrite(&timetemp,sizeof(float),1,output_file);
    fwrite(&rayorigintype,sizeof(int),1,output_file);
    fwrite(&number_of_ray_origins,sizeof(int),1,output_file);
    fwrite(&N_rays,sizeof(int),1,output_file);
	for (ray_num=0; ray_num < N_rays; ray_num++) {
    	fwrite(&ray[ray_num].theta,sizeof(float),1,output_file);
   		fwrite(&ray[ray_num].phi,sizeof(float),1,output_file);
	}
    
	for (RAY_ORIGIN_ID=0; RAY_ORIGIN_ID < number_of_ray_origins; RAY_ORIGIN_ID++)
	{
		/*
	    if (RAY_ORIGIN_TYPE == 0) ray_origin_pID = PG[RAY_ORIGIN_ID].ID;
	    if (RAY_ORIGIN_TYPE == 1) ray_origin_pID = PH[RAY_ORIGIN_ID].ID;
	    if (RAY_ORIGIN_TYPE == 2) ray_origin_pID = PD[RAY_ORIGIN_ID].ID;
	    if (RAY_ORIGIN_TYPE == 3) ray_origin_pID = PB[RAY_ORIGIN_ID].ID;
	    if (RAY_ORIGIN_TYPE == 4) ray_origin_pID = PS[RAY_ORIGIN_ID].ID;
	    if (RAY_ORIGIN_TYPE == 5) ray_origin_pID = PBH[RAY_ORIGIN_ID].ID;
	    fwrite(&ray_origin_pID,sizeof(int),1,output_file);
		*/
		
	    fwrite(&RAY_ORIGIN_ID,sizeof(int),1,output_file);
	    if (((float)(RAY_ORIGIN_ID))/500.0 == (float)(RAY_ORIGIN_ID/500)) {
	      printf(" Processing ray (origin %i) ... (N_steps) \n",RAY_ORIGIN_ID);
	    }
	    ALL_CELL_COUNTER++;
	    if (RAY_ORIGIN_TYPE == 0) for(j=0;j<3;j++) origin_pos[j]=PG[RAY_ORIGIN_ID].pos[j];
	    if (RAY_ORIGIN_TYPE == 1) for(j=0;j<3;j++) origin_pos[j]=PH[RAY_ORIGIN_ID].pos[j];
	    if (RAY_ORIGIN_TYPE == 2) for(j=0;j<3;j++) origin_pos[j]=PD[RAY_ORIGIN_ID].pos[j];
	    if (RAY_ORIGIN_TYPE == 3) for(j=0;j<3;j++) origin_pos[j]=PB[RAY_ORIGIN_ID].pos[j];
	    if (RAY_ORIGIN_TYPE == 4) for(j=0;j<3;j++) origin_pos[j]=PS[RAY_ORIGIN_ID].pos[j];
	    if (RAY_ORIGIN_TYPE == 5) for(j=0;j<3;j++) origin_pos[j]=PBH[RAY_ORIGIN_ID].pos[j];
	    RAY_ORIGIN_CELL = find_cell_from_scratch(origin_pos);
	    
	for (ray_num=0; ray_num < N_rays; ray_num++) 
	{
		for (j=0;j<3;j++) ray[ray_num].pos[j] = origin_pos[j];
		integrate_ray_to_escape(&ray[ray_num],RAY_ORIGIN_CELL);

		ray[ray_num].Z /= ray[ray_num].nh;
		ray[ray_num].neutral_frac /= ray[ray_num].nh;
		ray[ray_num].nh 	*= nh_prefactor;
		ray[ray_num].nh_hot *= nh_prefactor;

	    fwrite(&ray[ray_num].nh_hot,sizeof(float),1,output_file);
	    fwrite(&ray[ray_num].Z,sizeof(float),1,output_file);
			printf("%e %e \n",ray[ray_num].nh_hot,ray[ray_num].Z);	    
	}
	}
    fclose(output_file);
    free_ngblists();
	printf("\n Done! \n");
	return 0;
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
  for(i=1;i<=Ngas;i++) 
    {
      P3d[i]->Pos[0] = PG[i-1].pos[0];
      P3d[i]->Pos[1] = PG[i-1].pos[1];
      P3d[i]->Pos[2] = PG[i-1].pos[2];
      P3d[i]->Vel[0] = PG[i-1].vel[0];
      P3d[i]->Vel[1] = PG[i-1].vel[1];
      P3d[i]->Vel[2] = PG[i-1].vel[2];
  
      P3d[i]->Mass = PG[i-1].mass;
      P3d[i]->u = PG[i-1].u;
      P3d[i]->rho = PG[i-1].rho;
      P3d[i]->hsml = PG[i-1].hsml;
      P3d[i]->numh = PG[i-1].nh;
      P3d[i]->nume = PG[i-1].ne;
      P3d[i]->z = PG[i-1].z;
    }
  for(i=1;i<=Nstar;i++)
    {
      P3d[Ngas+i]->Pos[0] = PS[i-1].pos[0];
      P3d[Ngas+i]->Pos[1] = PS[i-1].pos[1];
      P3d[Ngas+i]->Pos[2] = PS[i-1].pos[2];
      P3d[Ngas+i]->Pos[0] = PS[i-1].vel[0];
      P3d[Ngas+i]->Pos[1] = PS[i-1].vel[1];
      P3d[Ngas+i]->Pos[2] = PS[i-1].vel[2];
      P3d[Ngas+i]->Mass = PS[i-1].mass;
      P3d[Ngas+i]->u = 0.0;
      P3d[Ngas+i]->rho = 0.0;
      P3d[Ngas+i]->hsml = 0.0;
      P3d[Ngas+i]->numh = 0.0;
      P3d[Ngas+i]->nume = 0.0;
      P3d[Ngas+i]->z = 0.0;
    }
  for(i=1;i<=Nbh;i++)
    {
      P3d[Ngas+Nstar+i]->Pos[0] = PBH[i-1].pos[0];
      P3d[Ngas+Nstar+i]->Pos[1] = PBH[i-1].pos[1];
      P3d[Ngas+Nstar+i]->Pos[2] = PBH[i-1].pos[2];
      P3d[Ngas+Nstar+i]->Pos[0] = PBH[i-1].vel[0];
      P3d[Ngas+Nstar+i]->Pos[1] = PBH[i-1].vel[1];
      P3d[Ngas+Nstar+i]->Pos[2] = PBH[i-1].vel[2];
      P3d[Ngas+Nstar+i]->Mass = PBH[i-1].mass;
      P3d[Ngas+Nstar+i]->u = 0.0;
      P3d[Ngas+Nstar+i]->rho = 0.0;
      P3d[Ngas+Nstar+i]->hsml = 0.0;
      P3d[Ngas+Nstar+i]->numh = 0.0;
      P3d[Ngas+Nstar+i]->nume = 0.0;
      P3d[Ngas+Nstar+i]->z = 0.0;
    }
}
/// routines to allocate and void memory for the temporary P3d structure 
//     that's used to build the Ngb search tree
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
  int Nsize;
  Nsize=Ngas+Nstar+Nbh;
  if(Nsize>0)
    {
      free(P3d[1]);
      P3d++;
      free(P3d);
    }
}


