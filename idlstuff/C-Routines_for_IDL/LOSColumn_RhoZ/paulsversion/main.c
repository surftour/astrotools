#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "overhead.h"
#include "proto.h"

/* define three modes: */ 
/* eventually make these part of the command-line arguments? */ 

#define QSOMODE 0 	/* original QSO N_H calculation */ 
#define RTMODE 1	/* radiative transfer with N_H between stars and gas */ 
#define VMODE 2		/* calculate N_H from one 'viewing' direction for images */ 

/* set the mode here to one of the three */ 
#define MODE RTMODE 

#define DEBUG YES

/* 
 *  Nomenclature: for RTMODE: 
 *    origin = gas particle 
 *    destination = star or BH particle 
 */ 


int NUMBER_OF_RAYS;

int main(int argc, char **argv)
{

  int n_gas, N_particles; 

	/* Load the simulation snapshot and check for black holes */	
	char sim_file[100];
	if(argc==1)
	{
		printf("Snapshot filename: ");
		scanf("%s",sim_file);
	}else if(argc==2)
	{
		sprintf(sim_file,"%s",argv[1]);
	}else if(argc==3)
	{
		sprintf(sim_file,"%s",argv[1]);
	}else{
		printf("wrong number of arguments\n");
		printf("./sim [snapshot_filename] \n");
		exit(-1);
	}
	printf("Loading snapshot %s... \n",sim_file);
	readsnap(sim_file);
	printf("Snapshot %s loaded... \n",sim_file);

  if(MODE == QSOMODE && All.N_bh < 1) {
    fprintf(stderr,"No black holes! \n");
    exit(0);
  } else if (All.N_bh < 1) { 
    fprintf(stderr,"Warning: No black holes! \n");
  } 

  //free(PH);
  //free(PD);
  //free(PB);

  int N_rays; 
  if (MODE == QSOMODE) { 
    NUMBER_OF_RAYS = 1000; 
  } else if (MODE == RTMODE) { 
    NUMBER_OF_RAYS  = All.N_star+All.N_bh;
  } 
  N_rays = NUMBER_OF_RAYS;		/* Number of rays to use to trace the structure */
  N_rays = (int)(sqrt(N_rays)+1) * (int)(sqrt(N_rays)+1);		/* Convert to a perfect square */
  if (DEBUG) fprintf(stderr, "NUMBER_OF_RAYS = %d, N_rays = %d\n", NUMBER_OF_RAYS, N_rays); 

  /* 
   * initialize the tree 
   */

  setup_lineofsight_overhead();
  /* Build the tree-cell system to sort the gas properties */
  initialize_the_tree();
  printf("tree built... \n",sim_file);


	/* Get the info for the output files */
	int p_file_chk,j;
	float dummf0 = 0.0;
	char out_file_2[100];
	FILE *output_file;

	if((argc==2)||(argc==1))
	{
		printf("\nFilename for integrated values (NH, etc): ");
		scanf("%s",out_file_2);	
	}else{
		sprintf(out_file_2,"%s",argv[2]);
	}
    	if ((output_file = fopen(out_file_2, "wb")) == NULL) 
    		fprintf(stderr, "Cannot open %s\n", "output_file_2");	
	int dumfcount;


  /* Create & allocate memory for the ray list */

  printf("Initializing (%i) Rays... \n",N_rays);
  int byte, ray_num;
  Ray_struct *ray;
  ray = calloc(N_rays,sizeof(Ray_struct));

  int rayorigintype; 
  int RAY_ORIGIN_CELL; 

  ALL_CELL_COUNTER = 0;
  float origin_pos[3];
  float nh_prefactor = (All.UnitMass_in_g/PROTONMASS)/(SQR(All.UnitLength_in_cm));

/* 
 *  QSOMODE 
 *  This is supposed to represent the original purpose of the code, but 
 *  hasn't been actually checked. . . 
 */ 

  if (MODE == QSOMODE) {
    int number_of_ray_origins = 1; 
    rayorigintype = PTYPE_BH; 
    for (j=0;j<3;j++) origin_pos[j]=PBH[0].pos[j]; 
            RAY_ORIGIN_CELL = find_cell_from_scratch(origin_pos);
    generate_ray_angles(ray, N_rays); 
    for (ray_num=0; ray_num < N_rays; ray_num++) {
      ray[ray_num].n_hat[0] = sin(ray[ray_num].theta) * cos(ray[ray_num].phi);
      ray[ray_num].n_hat[1] = sin(ray[ray_num].theta) * sin(ray[ray_num].phi);
      ray[ray_num].n_hat[2] = cos(ray[ray_num].theta);
    }  

    for (ray_num=0; ray_num < N_rays; ray_num++) {
      ray[ray_num].n_hat[0] = sin(ray[ray_num].theta) * cos(ray[ray_num].phi);
      ray[ray_num].n_hat[1] = sin(ray[ray_num].theta) * sin(ray[ray_num].phi);
      ray[ray_num].n_hat[2] = cos(ray[ray_num].theta);
    }

    /* Print some of the overhead info to the binary NH output file */
    fwrite(&rayorigintype,sizeof(int),1,output_file);
    fwrite(&number_of_ray_origins,sizeof(int),1,output_file);
    fwrite(&N_rays,sizeof(int),1,output_file);
 
    for (ray_num=0; ray_num < N_rays; ray_num++) {
      fwrite(&ray[ray_num].theta,sizeof(float),1,output_file);
      fwrite(&ray[ray_num].phi,sizeof(float),1,output_file);
    }

    for (ray_num=0; ray_num < N_rays; ray_num++) {
      for (j=0;j<3;j++) ray[ray_num].pos[j] = PBH[0].pos[j]; 
      integrate_ray_to_escape(&ray[ray_num],RAY_ORIGIN_CELL);
       
      ray[ray_num].Z /= ray[ray_num].nh;
      ray[ray_num].neutral_frac /= ray[ray_num].nh;
      ray[ray_num].nh 	*= nh_prefactor;
      ray[ray_num].nh_hot *= nh_prefactor;
    }
  
    fwrite(&PBH[0].ID,sizeof(int),1,output_file);
 
    fwrite(&ray[ray_num].DEST_ID,sizeof(int),1,output_file);
    fwrite(&ray[ray_num].nh_hot,sizeof(float),1,output_file);
    fwrite(&ray[ray_num].Z,sizeof(float),1,output_file);

    fclose(output_file);

  } /* end of QSOMODE */ 

/* 
 *  RTMODE 
 *  Very simple radiative transfer. Calculated the N_H from all star and BH
 *  particles to each gas particle. 
 */ 

  if (MODE == RTMODE) {
    int number_of_ray_origins = All.N_star + All.N_bh; 
    int number_of_ray_destinations = All.N_gas; 

    printf("In RTMODE\n"); 
    printf("Beginning to process %i ray sources...\n",number_of_ray_origins);
    allocate_gas_sources(); 

    /* loop over gas particles */ 

    /* N_particles = All.N_gas; */
    N_particles = All.N_gas; 
    for (n_gas=0; n_gas<N_particles; n_gas++) {

      fprintf(stderr, "Processing gas particle %d of %d\n", n_gas, N_particles); 
      /* calculate the angles from the current gas particle to the star and BH particles */ 
      set_ray_angles(ray,N_rays,PG[n_gas].pos, PG[n_gas].ID);
       
      /* calculate each of the vectors */ 
      for (ray_num=0; ray_num < N_rays; ray_num++) {
        ray[ray_num].n_hat[0] = sin(ray[ray_num].theta) * cos(ray[ray_num].phi);
        ray[ray_num].n_hat[1] = sin(ray[ray_num].theta) * sin(ray[ray_num].phi);
        ray[ray_num].n_hat[2] = cos(ray[ray_num].theta);
      }

      RAY_ORIGIN_CELL = find_cell_from_scratch(PG[n_gas].pos); 

      /* Calculate N_H along the rays */ 

      for (ray_num=0; ray_num < N_rays; ray_num++) {
        for (j=0;j<3;j++) ray[ray_num].pos[j] = PG[n_gas].pos[j]; 
        /* need to fix this */ 
        /* integrate_ray_to_target(&ray[ray_num],RAY_ORIGIN_CELL); */
        
        ray[ray_num].Z /= ray[ray_num].nh;
        ray[ray_num].neutral_frac /= ray[ray_num].nh;
        ray[ray_num].nh 	*= nh_prefactor;
        ray[ray_num].nh_hot *= nh_prefactor;
  
        PG[n_gas].PTYPE[ray_num] = ray[ray_num].DEST_TYPE; 
        PG[n_gas].PID[ray_num] = ray[ray_num].DEST_ID; 
        PG[n_gas].NH[ray_num] = ray[ray_num].nh_hot;
        PG[n_gas].D[ray_num] = ray[ray_num].D; 

        if (DEBUG) printf("ray_num %d origin cell %d dest type %d dest id %d theta %f phi %f xmax[0] %.3f xmax[1] %.3f xmax[2] %.3f nh %.3f Z %f D %.3f\n",ray_num, RAY_ORIGIN_CELL,PG[n_gas].PTYPE[ray_num], PG[n_gas].PID[ray_num],ray[ray_num].theta,ray[ray_num].phi,ray[ray_num].xmax[0]-origin_pos[0],ray[ray_num].xmax[1]-origin_pos[1],ray[ray_num].xmax[2]-origin_pos[2],log10(PG[n_gas].NH[ray_num]),ray[ray_num].Z,PG[n_gas].D[ray_num]);
 
      } /* end of loop for ray calculations */ 
    }  /* end of loop over destination particles */ 
  }  /* end of RTMODE */ 


  /* 
   *  VMODE: 
   *  Calculate the N_H for all star particles along one theta,phi for a 
   *  quick & dirty treatment of dust in the output images 
   *  will compute an N_H for each star and BH particle 
   *
   */ 

  if (MODE == VMODE) {

    int number_of_sources, n_source; 

    number_of_sources = All.N_star + All.N_bh; 

    for(n_source=0; n_source < number_of_sources; n_source++) { 

    } 

  } /* end of VMODE */ 

  printf("\n Done! \n");
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

/* 
 *  RTMODE: 
 *  initialize rays from origin at gas particles to the star and black hole particles 
 */  

void set_ray_angles(Ray_struct *R, int N_rays, float origin_pos[3], int origin_id)
{
	int i, j;
	double x,y,z;
	for(i=0;i<N_rays;i++) {
		if(i<All.N_star)
		{
			R[i].ORIGIN_TYPE = PTYPE_GAS; 
			R[i].ORIGIN_ID = origin_id; 
			R[i].DEST_TYPE = PTYPE_STAR; 
			R[i].DEST_ID = PS[i].ID; 
			x =(PS[i].pos[0]-origin_pos[0]);
			y =(PS[i].pos[1]-origin_pos[1]);
			z =(PS[i].pos[2]-origin_pos[2]);
			R[i].xmax[0] = PS[i].pos[0];
			R[i].xmax[1] = PS[i].pos[1];
			R[i].xmax[2] = PS[i].pos[2];
			if(y < 0.0)
			{
				if(x < 0.0)
				{
					R[i].phi   = atan(y/x)+PI;
				}else{
					R[i].phi   = atan(y/x)+2.0*PI;
				}
			}else{
				if(x < 0.0)
				{
					R[i].phi   = atan(y/x)+PI;
				}else{
					R[i].phi   = atan(y/x);
				}
			}
			if(z< 0.0)
			{
				R[i].theta = atan(sqrt(x*x+y*y)/z)+PI;
			}else{
				R[i].theta = atan(sqrt(x*x+y*y)/z);
			}
			R[i].D = sqrt(x*x+y*y+z*z); 
                } else if ( (i >= All.N_star) && (i<(All.N_star+All.N_bh)) ) {	/* include the BHs */ 
			j=i-All.N_star;
			R[i].ORIGIN_TYPE = PTYPE_GAS; 
			R[i].ORIGIN_ID = origin_id; 
			R[i].DEST_TYPE = PTYPE_BH; 
			R[i].DEST_ID = PBH[i].ID; 
			x =(PBH[j].pos[0]-origin_pos[0]);
			y =(PBH[j].pos[1]-origin_pos[1]);
			z =(PBH[j].pos[2]-origin_pos[2]);
			R[i].xmax[0] = PBH[j].pos[0];
			R[i].xmax[1] = PBH[j].pos[1];
			R[i].xmax[2] = PBH[j].pos[2];
			if(y < 0.0)
			{
				if(x < 0.0)
				{
					R[i].phi   = atan(y/x)+PI;
				}else{
					R[i].phi   = atan(y/x)+2.0*PI;
				}
			}else{
				if(x < 0.0)
				{
					R[i].phi   = atan(y/x)+PI;
				}else{
					R[i].phi   = atan(y/x);
				}
			}
			if(z< 0.0)
			{
				R[i].theta = atan(sqrt(x*x+y*y)/z)+PI;
			}else{
				R[i].theta = atan(sqrt(x*x+y*y)/z);
			}
			R[i].D = sqrt(x*x+y*y+z*z); 
		}else{
			R[i].DEST_ID = -1;
			R[i].DEST_TYPE = -1;
			R[i].ORIGIN_ID = -1;
			R[i].theta = 0.0;
			R[i].phi= 0.0;
			R[i].xmax[0] = MAX_DISTANCE_FROM0;
			R[i].xmax[1] = MAX_DISTANCE_FROM0;
			R[i].xmax[2] = MAX_DISTANCE_FROM0;
		}
	}
}
