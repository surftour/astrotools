#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "overhead.h"
#include "proto.h"
#include "ngbtree3d.h"

#define DEBUG NO
#define NAMELEN 60 

int NUMBER_OF_RAYS;

int Ngas,Nstar,Nbh;

float *Pos;

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



/* int getnh(int argc, char **argv) */
int getnh(int argc, void *argv[]) 
{

  int i, j, n_gas, ray_num, N_rays, tmpval;
  int RAY_ORIGIN_CELL; 
  int first_bh_ray, first_star_ray, first_disk_ray, first_bulge_ray; 
  int *oi, *ot; 
  float *OUT_NH; 	/* output NH */ 
  float *OUT_Z;
  float *OUT_T;
  float *OUT_ThermalSZ;
  float *OUT_KineticSZ;
  float theta, phi; 
  float nh_prefactor, ne_prefactor, kSZ_prefactor;
  float *on; 
  float *oz;
  Ray_struct *ray;
  int N_side, N_junk, Npix;
  float earth_x, earth_y, earth_z;
  int ia;
  float itheta, iphi, dt, dp;
  float dummy[3];
  FILE *fd;
  char buf[200], buf1[200], buf2[200], buf3[400];


  void allocate_3d(void);
  void set_particle_pointer(void);
  void free_memory_3d(void);

  if (argc != 12) {
    fprintf(stderr, "Expected 12 arguments (found %d)\n",argc); 
    exit(0); 
  }

  All.N_gas = *(int *)argv[0];
  N_side= *(int *)argv[1]; 
  N_junk= *(int *)argv[2]; 
  earth_x = *(float *)argv[3];
  earth_y = *(float *)argv[4];
  earth_z = *(float *)argv[5];
  Pos = (float *)argv[6];
  /* returned quantities */ 
  OUT_NH = (float *)argv[7]; 
  OUT_Z = (float *)argv[8]; 
  OUT_T = (float *)argv[9]; 
  OUT_ThermalSZ = (float *)argv[10]; 
  OUT_KineticSZ = (float *)argv[11]; 


  All.N_star= 0; 
  All.N_bh= 0; 
  All.N_halo = 0; 
  All.N_disk = 0; 
  All.N_bulge = 0; 
  All.N_total = All.N_gas + All.N_star + All.N_bh; 
  Ngas = All.N_gas; 
  Nstar = All.N_star; 
  Nbh = All.N_bh; 

  printf("N_side= %d \n", N_side); 
  printf("N_junk= %d \n", N_junk); 
  printf("Ngas = %d \n", Ngas);


  /* check that we don't go try to integrate exactly along cell boundaries */

  if (sin(theta) == 0.) theta += 1.e-30; 
  if (sin(phi) == 0.) theta += 1.e-30; 

  allocate_3d();
  set_particle_pointer();
  allocate_gas(); 
  allocate_star(); 
  allocate_bh(); 

  printf("star, gas, and bh allocation done too\n"); 

  /* --------------------------------------------------

     OK, now we're ready to use whatever variables we
     need.

     We can access the various fields via:

     P3d[i+1]->Pos[0], or P3d[30]->rho= 100.0


  here's an example:

  printf("\n\n");
  printf("Here are some of the P3d variables\n");
  printf("P3d[1]->Pos[0,1,2]= %g|%g|%g\n",P3d[1]->Pos[0],P3d[1]->Pos[1],P3d[1]->Pos[2]);
  printf("P3d[3003]->u= %g\n",P3d[3003]->u);
  printf("P3d[%d]->Pos[0,1,2]= %g|%g|%g\n",Ngas+Nstar,P3d[Ngas+Nstar]->Pos[0],P3d[Ngas+Nstar]->Pos[1],P3d[Ngas+Nstar]->Pos[2]);
  printf("P3d[%d]->Pos[0,1,2]= %g|%g|%g\n",Ngas+Nstar+1,P3d[Ngas+Nstar+1]->Pos[0],P3d[Ngas+Nstar+1]->Pos[1],P3d[Ngas+Nstar+1]->Pos[2]);
  printf("\n\n");
   */
  printf("P3d[%d]->Vel[0,1,2]= %g|%g|%g\n",1,P3d[1]->Vel[0],P3d[1]->Vel[1],P3d[1]->Vel[2]);
  printf("P3d[%d]->Vel[0,1,2]= %g|%g|%g\n",2,P3d[2]->Vel[0],P3d[2]->Vel[1],P3d[2]->Vel[2]);
  printf("P3d[%d]->Vel[0,1,2]= %g|%g|%g\n",3,P3d[3]->Vel[0],P3d[3]->Vel[1],P3d[3]->Vel[2]);


  /* Fill in the structures that are normally done in readsnap */ 
  
  for (i=0; i<All.N_gas; i++) {
    PG[i].pos[0] = P3d[i+1]->Pos[0]; 
    PG[i].pos[1] = P3d[i+1]->Pos[1];
    PG[i].pos[2] = P3d[i+1]->Pos[2];
    PG[i].vel[0] = P3d[i+1]->Vel[0]; 
    PG[i].vel[1] = P3d[i+1]->Vel[1];
    PG[i].vel[2] = P3d[i+1]->Vel[2];
    PG[i].mass = P3d[i+1]->Mass;
    PG[i].u = P3d[i+1]->u; 
    PG[i].temp = P3d[i+1]->u * 80.766812 * (1.+4.*0.0789474)/(1+0.0789474*P3d[i+1]->nume);
	/* 0.0789474 comes from x_h= 0.76, and this is (1-x_h)/(4*x_h) */
    PG[i].rho = P3d[i+1]->rho; 
    PG[i].hsml = P3d[i+1]->hsml; 
    PG[i].nh = P3d[i+1]->numh; 
    PG[i].ne = P3d[i+1]->nume; 
    PG[i].z = P3d[i+1]->z; 
  } 
  printf("\n");
  /*
  printf("PG[1].u=    %g     PG[1].temp=    %g\n",PG[1].u,PG[1].temp);
  printf("PG[100].u=  %g     PG[100].temp=  %g\n",PG[100].u,PG[100].temp);
  printf("PG[4751].u= %g     PG[4751].temp= %g\n",PG[4751].u,PG[4751].temp);
  */


  /* 
   *  initialize the tree 
   */

  setup_lineofsight_overhead();


#ifndef USE_FULL_NEIGHBOR_CALC
	printf("Using: Faster TREE method\n");
	initialize_the_tree();
	fprintf(stdout, "tree built... \n"); 
#else
	printf("Using: FULL_NEIGHBOR_CALC\n");
	ngb3d_treeallocate(All.N_gas, 10*All.N_gas);
	ngb3d_treebuild((float **)&P3d[1], All.N_gas, 0,dummy,dummy);
	printf("Using New 3D tree.\n");
	allocate_ngblists();
#endif


  /* 
   *  initialize the rays 
   */ 

  NUMBER_OF_RAYS  = (12.0 * N_side * N_side);
  N_rays = NUMBER_OF_RAYS;
  N_rays = (int)(sqrt(N_rays)+1) * (int)(sqrt(N_rays)+1);		/* Convert to a perfect square */
  printf("NUMBER_OF_RAYS = %d, N_rays = %d, sizeof(Ray_struct) = %d\n", NUMBER_OF_RAYS, N_rays, sizeof(Ray_struct) ); 
  printf("Initializing (%i) Rays... \n",N_rays);
  if(!(ray = malloc(N_rays*sizeof(struct Ray_s))))
     {
        printf("failed to allocate memory: N_rays= %d  size= %d\n",N_rays,N_rays*sizeof(struct Ray_s));
        exit(0);
     }
  printf("Done allocating memory for rays\n"); 
  ALL_CELL_COUNTER = 0;
  nh_prefactor = (All.UnitMass_in_g/PROTONMASS)/(SQR(All.UnitLength_in_cm));
  printf("nh_prefactor= %g\n",nh_prefactor);
  ne_prefactor = (All.UnitMass_in_g/PROTONMASS)/pow(All.UnitLength_in_cm,3);
  printf("ne_prefactor= %g\n",ne_prefactor);
  kSZ_prefactor = (All.UnitMass_in_g/PROTONMASS)/All.UnitLength_in_cm/All.UnitTime_in_s;
  printf("kSZ_prefactor= %g\n",kSZ_prefactor);


  /* 
   *  set up the rays 
   */ 


  for (ray_num=0; ray_num<NUMBER_OF_RAYS; ray_num++) { 
    ray[ray_num].pos[0] = earth_x;
    ray[ray_num].pos[1] = earth_y;
    ray[ray_num].pos[2] = earth_z;
  } 

  ray_num= 0;
/*
  dt= 180.0/N_theta;
  dp= 360.0/N_phi;
  printf("dt=%g   dp=%g\n",dt,dp);
  for(it=0; it<N_theta; it++) {
    itheta= dt*0.5 + dt*it;
    for(ip=0; ip<N_phi; ip++) {
	iphi= dp*0.5 + dp*ip;
	iphi= iphi - 180.0;

	if(!(ray_num%10000))
	  {
		printf("theta= %g  phi= %g\n",itheta,iphi);
	  }

	ray[ray_num].theta = itheta * (PI / 180.0);
	ray[ray_num].phi= iphi * (PI / 180.0);
	ray_num++;
    }
  }
  printf("iphi= %g\n",iphi);
  printf("itheta= %g\n",itheta);
  printf("ray_num= %d\n",ray_num);
*/
  

/* ----------------------------------------------------------- */

  printf("opening: /home/tcox/CMBStuff/HealPix_RingPixels/healpixangles.txt");
  if(!(fd=fopen("/home/tcox/CMBStuff/HealPix_RingPixels/healpixangles.txt","r"))) {
    printf("Can't open file: /home/tcox/CMBStuff/HealPix_RingPixels/healpixangles.txt");
    exit(0); 
  }

  Npix= 12.0*N_side*N_side;
  printf("N_side= %g\n",N_side);
  printf("Npix= %g\n",Npix);

  for(ia=0; ia<Npix; ia++) {

 	*buf = 0; 
	fgets(buf, 200, fd);
	sscanf(buf, "%s%s%s", buf1, buf2, buf3);

	ray[ray_num].theta = atof(buf2);
	ray[ray_num].phi= atof(buf3);
	ray_num++;
  }
  fclose(fd);

  printf("DONE setting HealPix angles.\n");

/* ----------------------------------------------------------- */


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

    /* if(!(ray_num%10000)) */
    if(!(ray_num%100))
      {
	printf("%d ....",ray_num); fflush(stdout);
      }

    /* ray[ray_num].theta = theta[ray_num]; */
    /* ray[ray_num].phi = phi[ray_num]; */
    ray[ray_num].n_hat[0] = sin(ray[ray_num].theta) * cos(ray[ray_num].phi);
    ray[ray_num].n_hat[1] = sin(ray[ray_num].theta) * sin(ray[ray_num].phi);
    ray[ray_num].n_hat[2] = cos(ray[ray_num].theta);

#ifndef USE_FULL_NEIGHBOR_CALC
	RAY_ORIGIN_CELL = find_cell_from_scratch(ray[ray_num].pos);
#else
	RAY_ORIGIN_CELL = 0;
#endif
    integrate_ray_to_escape(&ray[ray_num],RAY_ORIGIN_CELL);

    //printf("ray_num %d   nh %.3f Z %f T %f \n",ray_num,log10(ray[ray_num].nh_hot),ray[ray_num].Z,ray[ray_num].Temp);
    ray[ray_num].Z /= ray[ray_num].nh;
    ray[ray_num].neutral_frac /= ray[ray_num].nh;
    ray[ray_num].Temp /= ray[ray_num].nh;

    ray[ray_num].nh   *= nh_prefactor;
    ray[ray_num].nh_hot *= nh_prefactor;

    ray[ray_num].thermalSZ *= nh_prefactor;
    ray[ray_num].kineticSZ *= kSZ_prefactor;

    if (DEBUG) fprintf(stderr, "ray_num %d origin_cell %d origin_type %d origin_id %d theta %.1f phi %.1f nh %.3f Z %f\n",ray_num, RAY_ORIGIN_CELL, ray[ray_num].ORIGIN_TYPE, ray[ray_num].ORIGIN_ID,ray[ray_num].theta,ray[ray_num].phi,log10(ray[ray_num].nh_hot),ray[ray_num].Z);
    //printf("ray_num %d  theta %.2f phi %.2f nh %.3f Z %f T %f \n",ray_num, ray[ray_num].theta,ray[ray_num].phi,log10(ray[ray_num].nh_hot),ray[ray_num].Z,ray[ray_num].Temp);



    /* write the results to the input arrays */ 

    on = OUT_NH + ray_num; 
    *on = ray[ray_num].nh_hot; 

    oz = OUT_Z + ray_num;
    *oz = ray[ray_num].Z;

    oz = OUT_T + ray_num;
    *oz = ray[ray_num].Temp;

    oz = OUT_ThermalSZ + ray_num;
    *oz = ray[ray_num].thermalSZ;

    oz = OUT_KineticSZ + ray_num;
    *oz = ray[ray_num].kineticSZ;


/*
    OUT_NH[ray_num] = ray[ray_num].nh_hot; 
*/

/*
    if (DEBUG) fprintf(stderr, "(%d) %d %d %f -- %d %d %f -- %d %d %f\n", ray_num, OUT_ID[ray_num], OUT_TYPE[ray_num], OUT_NH[ray_num], *oi, *ot, *on, ray[ray_num].ORIGIN_ID, ray[ray_num].ORIGIN_TYPE, ray[ray_num].nh_hot); 
*/

/*
    if(ray_num<10)
	printf("ray_num= %d  OUT_NH= %g  OUT_Z= %g\n",ray_num,ray[ray_num].nh_hot,ray[ray_num].Z);
*/

  }

  printf(" done. \n"); fflush(stdout);

#ifdef USE_FULL_NEIGHBOR_CALC
  ngb3d_treefree();
  //free_ngblists();
#endif


  free(ray);
  free_memory_3d();

  printf("\n Finished call to getnh.so ! \n");

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

      P3d[i]->Vel[0] = pos[3];
      P3d[i]->Vel[1] = pos[4];
      P3d[i]->Vel[2] = pos[5];
  
      P3d[i]->Mass = pos[6];
      P3d[i]->u = pos[7];
      P3d[i]->rho = pos[8];
      P3d[i]->hsml = pos[9];
      P3d[i]->numh = pos[10];
      P3d[i]->nume = pos[11];
      P3d[i]->z = pos[12];

      pos+=13;
    }

  for(i=Ngas+1;i<=(Ngas+Nstar+Nbh);i++)
    {

      P3d[i]->Pos[0] = pos[0];
      P3d[i]->Pos[1] = pos[1];
      P3d[i]->Pos[2] = pos[2];

      P3d[i]->Vel[0] = pos[3];
      P3d[i]->Vel[1] = pos[4];
      P3d[i]->Vel[2] = pos[5];

      P3d[i]->Mass = pos[6];
      P3d[i]->u = 0.0;
      P3d[i]->rho = 0.0;
      P3d[i]->hsml = 0.0;
      P3d[i]->numh = 0.0;
      P3d[i]->nume = 0.0;
      P3d[i]->z = 0.0;

      pos+=13;
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
  int Nsize;
  Nsize=Ngas+Nstar+Nbh;
  if(Nsize>0)
    {
      free(P3d[1]);
      P3d++;
      free(P3d);
    }
}
