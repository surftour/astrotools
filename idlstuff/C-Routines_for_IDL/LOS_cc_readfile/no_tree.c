#include "overhead.h"
#include "proto.h"
#include "ngbtree3d.h"


//float *r_along_ray;
//LOCALVAL *local_quantities;



/* line integral of quantities (such as nh) along a ray, 
 *   for column densities, etc. (just a basic simpsons rule integration, 
 *   i.e. linearly interpolating between each point)
 */
void integrate_ray(Ray_struct *R)
{
	float nh,z,neutral_frac,temp,nh_gadget_units,nh_hot;
	float vrad;
	int i;
	float *vec_to_integrate, *dr;
	dr = (float *)malloc(R->N_steps * sizeof(float));
	vec_to_integrate = (float *)malloc(R->N_steps * sizeof(float));
	
	
	for (i=0; i<R->N_steps; i++) dr[i] = FRACTIONAL_STEPSIZE * R->local[i].h ;
	
	/* First, column density */
	for (i=0; i<R->N_steps; i++) vec_to_integrate[i] = R->local[i].rho ;
	nh_gadget_units = simp(dr, vec_to_integrate, R->N_steps);
	nh  = nh_gadget_units;
	R->nh = nh;
	
	/* Column density of the hot-phase ISM along the ray*/
	for (i=0; i<R->N_steps; i++) vec_to_integrate[i] = R->local[i].rho_hot * R->local[i].rho ;
	nh_hot = simp(dr, vec_to_integrate, R->N_steps);
	R->nh_hot = nh_hot;

	/* Mass-weighted metallicity (so R->Z = total metal fraction) */
	for (i=0; i<R->N_steps; i++) vec_to_integrate[i] = R->local[i].Z * R->local[i].rho ;
	z  = simp(dr, vec_to_integrate, R->N_steps);
	R->Z = z;
	
	/* Mass-weighted neutral fraction (so R->neutral_frac = total neutral fraction) */
	for (i=0; i<R->N_steps; i++) vec_to_integrate[i] = R->local[i].neutral_frac * R->local[i].rho ;
	neutral_frac  = simp(dr, vec_to_integrate, R->N_steps);
	R->neutral_frac = neutral_frac;
	
	/* Mass-weighted temperature (so R->temp = temp fraction) */
	for (i=0; i<R->N_steps; i++) vec_to_integrate[i] = R->local[i].T * R->local[i].rho ;
	temp  = simp(dr, vec_to_integrate, R->N_steps);
	R->Temp = temp;

        /* Thermal SZ: n_e * temperature */
	/*
        for (i=0; i<R->N_steps; i++) vec_to_integrate[i] = R->local[i].rho * R->local[i].T;
        temp  = simp(dr, vec_to_integrate, R->N_steps);
	*/
        R->thermalSZ= temp;    /* same as above */

        /* Kinetic SZ: n_e * radial velocity */
	vrad= R->n_hat[0]*R->local[i].vel[0] +
		R->n_hat[1]*R->local[i].vel[1] +
		  R->n_hat[2]*R->local[i].vel[2];
        for (i=0; i<R->N_steps; i++) vec_to_integrate[i] = R->local[i].rho * vrad;
        temp  = simp(dr, vec_to_integrate, R->N_steps);
        R->kineticSZ= temp;
	
	/* need to free vectors created in this subroutine */
	free(dr);
	free(vec_to_integrate);
} 

float simp(float *dx, float *y, int n)
{
	float z;
	int i;
	z = 0.5 * (y[0]*dx[0] + y[n-1]*dx[n-1]);
	for (i=1; i<n-1; i++) {
	z += y[i]*dx[i];
	}
	return z;
}



/* revised routine to carry a ray out of the box in the good ol' fashioned brute force 
	manner */
void line_of_sight(Ray_struct *ray)
{
	float x_max,x_min,dr_min,x,y,z,dr,rhodr;
	int n_steps, n_steps_max,r_i;
	LOCALVAL local_temp;
	//local_temp = (LOCALVAL *)malloc(sizeof(LOCALVAL));
	
	x_max = 1.0*(float)MAX_DISTANCE_FROM0;	
	x_min = -x_max;
	n_steps_max = (int)LINE_NMAX;
	dr_min = abs(2.0*x_max) / (2000.0 * (float)n_steps_max);

	// Initialize the location to start from 
	n_steps = 0;
	x = ray->pos[0];
	y = ray->pos[1];
	z = ray->pos[2];
	//printf(" x = %e ",x);

	while (x>x_min && x<x_max && y>x_min && y<x_max 
	    && z>x_min && z<x_max && n_steps<n_steps_max) {
	
	  local_vars(x,y,z,&local_temp);
	  dr = (float)FRACTIONAL_STEPSIZE * local_temp.h;
	  	if (dr < dr_min) dr = dr_min;
	  	x += dr * ray->n_hat[0];
	  	y += dr * ray->n_hat[1];
	  	z += dr * ray->n_hat[2];

	  rhodr = local_temp.rho * dr;
	  ray->nh 		+= rhodr;
	  ray->nh_hot	+= local_temp.rho_hot * rhodr;
	  ray->Z		+= local_temp.Z * rhodr;
	  ray->neutral_frac += local_temp.neutral_frac * rhodr;
	  n_steps += 1;
	}
	//free(local_temp);
	return;
}



/* This routine calculates the locally-averages 
 *   values of physical quantities at a given location (x,y,z) (using the brute-force 
 *   calculation to determine the list of nearest neighbors)
 */
void local_vars(float x, float y, float z, LOCALVAL *vals)
{
  int Num_neighbors = (int)NUMBER_OF_NEIGHBORS;	
  	/* Number of local particle neighbors to average over */
  float r0[3];					/* Position vector (from input) */
  	r0[0] = x;
  	r0[1] = y;
  	r0[2] = z;
  int particle_id = 0;			/* Identifies which type of particles we're going
  										to look at for the calculation (could easily be
  										changed, but will almost always be GAS = 0) */

/*
  float *r2_list;
  int   *ngb_list;
printf("Allocating ngb_list and r2_list arrays\n"); fflush(stdout);
    ngb_list= (int *)malloc(Num_neighbors*sizeof(int));
    r2_list= (float *)malloc(Num_neighbors*sizeof(float));
*/
  float h,h2,hsml=10.0;		/* local smoothing length (to be returned) */
  float Hmax= 1000.0;
  int ngbfound;

/* This is Phil's old, brute force method.
  find_neighbors(r0,Num_neighbors,ngb_list,r2_list,&hsml);
  hsml = sqrt(hsml);
  ngbfound= Num_neighbors;
*/
  h= 10.0;
//printf("Heading into ngb3d_treefind - I hope this works.\n"); fflush(stdout);
  h2 = ngb3d_treefind( r0, Num_neighbors,1.04*h,&ngb_list,&r2_list, Hmax,&ngbfound);
//printf("Are we getting this back? h= %g\n",sqrt(h2)); fflush(stdout);

  hsml = sqrt(h2);

  float r,hinv,hinv3,hinv4,u,wk,dwk,rho,pressure,temp,CurlVel,DivVel;
  float divv,rotv[3],mass_j,Z_rho,Ne,neutral_frac,energyd,spec_egy;
  float vx,vy,vz;
  hinv = 1.0/hsml; hinv3 = hinv*hinv*hinv; hinv4 = hinv3*hinv;
  rho = pressure = temp = CurlVel = DivVel = Z_rho = Ne = neutral_frac = energyd = 0.0;
  spec_egy = vx= vy= vz= 0.0;
  int n, ii, j;


  /* Use the kernel to weight each point and sum the mass, etc. about x,y,z */
  for (n=0; n<ngbfound; n++) {
    j = ngb_list[n];
    r = sqrt(r2_list[n]);
    /*
    printf("n= %d    j (ngb_list)= %d  r= %g    h=%g\n",n,j,r,hsml); fflush(stdout);
    */
    if (r<hsml) {
      u = r * hinv;
      ii = (int)(u*KERNEL_TABLE);

	  wk =hinv3*( Kernel[ii]  + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
	  dwk=hinv4*( KernelDer[ii] + (KernelDer[ii+1]-KernelDer[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
	  mass_j = PG[j].mass;

	  rho   += mass_j * wk;
	  /* temp	+= PG[j].temp * wk; */
	  
	  /* Calculate mass-weighted averages of Z, ne, u, and the ionized fraction, 
	   *   noting the definitions of the read-in quantities from the snapshots */
	  Z_rho 	   += 	PG[j].z  *  mass_j * wk;
          Ne    	   += 	PG[j].ne *  mass_j * wk;
          neutral_frac     +=	PG[j].nh *  mass_j * wk;
	  spec_egy	   += 	PG[j].u  *  mass_j * wk;
	  temp		   +=   PG[j].temp * mass_j * wk;
	  vx		   +=   PG[j].vel[0] * mass_j * wk;
	  vy		   +=   PG[j].vel[1] * mass_j * wk;
	  vz		   +=   PG[j].vel[2] * mass_j * wk;

	  
	  /* Div and curl of velocity field, don't need these for now
	  divv -= mass_j * dwk/r *
			   ( dx * (P[i].VelPred[0] - P[j].VelPred[0])
		        + dy * (P[i].VelPred[1] - P[j].VelPred[1])
		        + dz * (P[i].VelPred[2] - P[j].VelPred[2]) );
	  rotv[0] += mass_j * dwk/r *
			      (  dz * (P[i].VelPred[1] - P[j].VelPred[1])
			       - dy * (P[i].VelPred[2] - P[j].VelPred[2]) );
	  rotv[1] += mass_j * dwk/r *
		          (  dx * (P[i].VelPred[2] - P[j].VelPred[2])
      	           - dz * (P[i].VelPred[0] - P[j].VelPred[0]) );
	  rotv[2] += mass_j * dwk/r *
		          (  dy * (P[i].VelPred[0] - P[j].VelPred[0])
			       - dx * (P[i].VelPred[1] - P[j].VelPred[1]) );
	  DivVel = divv/rho;
	  CurlVel = sqrt(rotv[0]*rotv[0] + rotv[1]*rotv[1] + rotv[2]*rotv[2])/rho;

	  SphP[i].Pressure = GAMMA_MINUS1*(SphP[i].EgySpecPred)*SphP[i].Density;
	  */
  }}
  /* Account for mass-weighted averages and convert to LOCALVALS quantities */
  Ne           /= rho;
  Z_rho	       /= rho;
  neutral_frac /= rho;
  spec_egy     /= rho;
  temp         /= rho;
  vx	       /= rho;
  vy	       /= rho;
  vz	       /= rho;
  
  
  vals->h = hsml;		/* returns smoothing length */
  vals->rho = rho;		/* returns density */
  vals->P_eff = GAMMA_MINUS1*(spec_egy*rho);	/* effective pressure */
  vals->u = spec_egy;	/* returns specific internal energy */
  vals->T = temp;		/* returns temperature */
  vals->Z = Z_rho;  	/* mass-weighted metallicity */
  vals->ne = Ne;		/* number of electrons per hydrogen atom/nucleus */
  vals->neutral_frac = neutral_frac;	/* returns fraction of H in neutral atoms */
  vals->vel[0]= vx;
  vals->vel[1]= vy;
  vals->vel[2]= vz;
  
  //float u_phys   = vals->u * All.UnitEnergy_in_cgs/All.UnitMass_in_g;
  //float rho_phys = vals->rho * All.UnitDensity_in_cgs;
  //get_two_phase_breakdown(u_phys, rho_phys, vals->ne, vals);
  get_two_phase_breakdown(vals->u, vals->rho, vals->ne, vals);
  
/*
  free(ngb_list);
  free(r2_list);
*/
  return;
}


/* Brute force routine to find the Ntot nearest neighbors to a given
 *  point and return the indices of each neighbor in nbr_list and the
 *  r-squared distances in r2_list. Searches gas only, but easily modified
 *  to whichever type.
 */
void find_neighbors(float r0[3],int Ntot,int *nbr_list,float *r2_list, float *r2_max)
{
	int n,i,j,k,n_r2max;
	float r2,r2max;
	r2max = 0.0;
	for (n=0; n<All.N_gas; n++) {
	  r2 = 0.0;
	  for (j=0;j<3;j++) r2 += (PG[n].pos[j]-r0[j])*(PG[n].pos[j]-r0[j]);
	  
	  if (n<Ntot) {
	    r2_list[n] = r2;
	    nbr_list[n] = n;
	    if (r2 > r2max) {
	    	r2max = r2;
	    	n_r2max = n;
	    }
	  }
	  
	  if (n>=Ntot) {
	    if (r2 < r2max) {
	      r2_list[n_r2max] = r2;
	      nbr_list[n_r2max] = n;
	      r2max = 0.0;
	   	  for (k=0;k<Ntot;k++) {
	    	if (r2_list[k] > r2max) {
	    		r2max = r2_list[k];
	    		n_r2max = k;
	    	}
	      }
	  }}
	}
	*r2_max = r2max;
}






void allocate_ngblists(void)
{
  int Num_neighbors = (int)NUMBER_OF_NEIGHBORS;	
  int n_steps_max = (int)LINE_NMAX;

  ngb_list= (int *)malloc(Num_neighbors*sizeof(int));
  r2_list= (float *)malloc(Num_neighbors*sizeof(float));

  //local_quantities = (LOCALVAL *) malloc (n_steps_max*sizeof(LOCALVAL));
  //r_along_ray = (float *) malloc (n_steps_max*sizeof(float));    

}


void free_ngblists(void)
{
  free(ngb_list);
  free(r2_list);
  //free(local_quantities);
  //free(r_along_ray);
}
