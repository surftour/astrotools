#include "overhead.h"
#include "proto.h"
/* 
 * Routines for building and using a tree-based grid to compute local
 *   gas properties with order h_sml resolution in a quick manner, and 
 *   to quickly trace rays out of the simulation in calculating column
 *   densities
 *   
 */ 


/* Initialize the tree based on SPH particle list, using the particle
 *   smoothing lengths to determine the depth of the grid anywhere
 */
 void initialize_the_tree(void)
 {
	/* readsnap.c should set up the following:
	 *    Ngas = number of SPH gas particles
	 *    PG[i].hsml = smoothing length of gas particle i (in struct PG)
	 *    PG[i].pos[0,1,2] = x,y,z (respectively) positions of particle i
	 */

    int np,i,j,k,in_box_check,N_SUB_ID,nx,ny,nz,n_ID,PID,n_ID_sub,n_cell;
    int n_subshell_i,n_subshell_f,n_x[3];

	extern struct CELL_STRUCT *CELL; 
	if(!(CELL=(struct CELL_STRUCT *) malloc(All.N_gas*((int)(2.0/TREE_FRAC_HSML))*sizeof(struct CELL_STRUCT))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }


	int CURRENT_CELL_ID = 0;
	ALL_CELL_COUNTER = 1;
    
	for (j=0;j<3;j++) CELL[0].min_x[j] = -1.0*MAX_DISTANCE_FROM0;
	CELL[0].width = 2.0*MAX_DISTANCE_FROM0;
	CELL[0].sub_cell_check = 1;
	CELL[0].parent_ID = -1;
	for (nx=0;nx<N_SUBCELL;nx++) {
	for (ny=0;ny<N_SUBCELL;ny++) {
	for (nz=0;nz<N_SUBCELL;nz++) {

		n_ID_sub = nx*N_SUBCELL_2 + ny*N_SUBCELL + nz;
		n_ID = n_ID_sub + ALL_CELL_COUNTER;
		CELL[n_ID].parent_ID = 0;
		PID = CELL[n_ID].parent_ID;

		CELL[PID].sub_cell_IDs[n_ID_sub] = n_ID;
		CELL[n_ID].width = CELL[PID].width / N_SUBCELL;
		CELL[n_ID].N_particles_contained = 0;

		n_x[0] = nx; n_x[1] = ny; n_x[2] = nz;
		for (j=0;j<3;j++) {
			CELL[n_ID].min_x[j] = CELL[PID].min_x[j] + CELL[n_ID].width * n_x[j];		
		}

	}}}
	CELL[0].N_particles_contained = 0;
	
	/* Alright, so, the first run -- we have box ID 0, the big sucker, and
		we're inside it right now. Do a run of all the particles -- for each, 
		decide which box it's in (dump that to the ticker for the sub-boxes
		counting how many particles go in them), check its h_sml, and also 
		make sure it's inside the giant box itself */
	for (np=0;np<All.N_gas;np++) {
		in_box_check = 1;
		for (j=0;j<3;j++) {
			n_x[j] = (int)((PG[np].pos[j] - CELL[0].min_x[j])/(CELL[0].width/N_SUBCELL));
			if (n_x[j] < 0 || n_x[j] > N_SUBCELL-1) in_box_check = 0;
		}
		if (in_box_check != 0) {
			N_SUB_ID = n_x[0]*N_SUBCELL_2 + n_x[1]*N_SUBCELL + n_x[2];
			CELL[0].N_particles_contained += 1;
			CELL[ALL_CELL_COUNTER + N_SUB_ID].N_particles_contained += 1;
		}
	}

	for (n_cell=ALL_CELL_COUNTER; n_cell<ALL_CELL_COUNTER+N_SUBCELL_3; n_cell++) {
		CELL[n_cell].hsml_min = 1.0e10;
		if (CELL[n_cell].N_particles_contained > 0) {
		if(!(CELL[n_cell].particle_IDs = (int *) malloc(CELL[n_cell].N_particles_contained*sizeof(int))))
		{
			printf("Error allocating particle_IDs for cell %d in tree.c.\n",n_cell);
			exit(-1);
		}
		CELL[n_cell].N_particles_contained = 0;
	}}

	for (np=0;np<All.N_gas;np++) {
		in_box_check = 1;
		for (j=0;j<3;j++) {
			n_x[j] = (int)((PG[np].pos[j] - CELL[0].min_x[j])/(CELL[0].width/N_SUBCELL));
			if (n_x[j] < 0 || n_x[j] > N_SUBCELL-1) in_box_check = 0;
		}
		if (in_box_check != 0) {
			N_SUB_ID = ALL_CELL_COUNTER + n_x[0]*N_SUBCELL_2 + n_x[1]*N_SUBCELL + n_x[2];
			CELL[N_SUB_ID].particle_IDs[CELL[N_SUB_ID].N_particles_contained] = np;
			CELL[N_SUB_ID].N_particles_contained += 1;
			if (PG[np].hsml < CELL[N_SUB_ID].hsml_min) 
				CELL[N_SUB_ID].hsml_min = PG[np].hsml;
		}
	}
	n_subshell_i = ALL_CELL_COUNTER;
	n_subshell_f = ALL_CELL_COUNTER + N_SUBCELL_3;
	ALL_CELL_COUNTER += N_SUBCELL_3;
	
	//printf("Tree opened, building... \n");
	for (n_cell=n_subshell_i; n_cell<n_subshell_f; n_cell++) build_cell(n_cell);
	//printf("Tree constructed, %d total cells \n",ALL_CELL_COUNTER);
	TOTAL_NUMBER_OF_CELLS = ALL_CELL_COUNTER;
	
	return;
}


/* The workhorse of the tree-building phase, this loops over itself to 
	build a recursive tree over the grid and sets the gas values in the 
	tree when a small enough scale is reached */
void build_cell(int cell_ID)
{
	int nx,ny,nz,n_x[3],j,np,p_ID,n_ID_sub,n_ID,n_subshell_i,n_subshell_f,in_box_check;
	
	CELL[cell_ID].sub_cell_check = 1;
	
	if (CELL[cell_ID].N_particles_contained == 0) 
			CELL[cell_ID].sub_cell_check = 0;
	if (CELL[cell_ID].width <= TREE_FRAC_HSML * CELL[cell_ID].hsml_min) 
			CELL[cell_ID].sub_cell_check = 0;
	if (CELL[cell_ID].width <= TREE_MIN_SIZE) 
			CELL[cell_ID].sub_cell_check = 0;
			
	if (CELL[cell_ID].sub_cell_check == 0) {
		get_tree_localvals(cell_ID);
		return;
	}
	
	/* ok, if we've made it this far, we need to build sub-cells for this 
	   particular cell. basically, follow the routine for the first set of sub-cells, 
	   with a few mods */

	/* first set up the basic variables of the new cell */
	for (nx=0;nx<N_SUBCELL;nx++) {
	for (ny=0;ny<N_SUBCELL;ny++) {
	for (nz=0;nz<N_SUBCELL;nz++) {

		n_ID_sub = nx*N_SUBCELL_2 + ny*N_SUBCELL + nz;
		n_ID = n_ID_sub + ALL_CELL_COUNTER;

		CELL[cell_ID].sub_cell_IDs[n_ID_sub] = n_ID;
		CELL[n_ID].parent_ID = cell_ID;
		CELL[n_ID].width = CELL[cell_ID].width / N_SUBCELL;
		CELL[n_ID].N_particles_contained = 0;

		n_x[0] = nx; n_x[1] = ny; n_x[2] = nz;
		for (j=0;j<3;j++)
			CELL[n_ID].min_x[j] = CELL[cell_ID].min_x[j] + CELL[n_ID].width * n_x[j];
	}}}
	/* now loop over the gas particles in the cell and determine which sub-cells 
		they belong to */
	for (np=0; np<CELL[cell_ID].N_particles_contained; np++) {
		//printf("pidchk = %d %d %d %d \n",np,CELL[cell_ID].N_particles_contained,n_ID,CELL[n_ID].N_particles_contained);
		p_ID = CELL[cell_ID].particle_IDs[np];

		in_box_check = 1;
		for (j=0;j<3;j++) {
			n_x[j] = (int)((PG[p_ID].pos[j] - CELL[cell_ID].min_x[j])/(CELL[cell_ID].width/N_SUBCELL));
			if (n_x[j] < 0 || n_x[j] > N_SUBCELL-1) in_box_check = 0;
		}
		if (in_box_check != 0) {
			n_ID_sub = n_x[0]*N_SUBCELL_2 + n_x[1]*N_SUBCELL + n_x[2];
			n_ID = n_ID_sub + ALL_CELL_COUNTER;
			CELL[n_ID].N_particles_contained += 1;
		}
	}

	for (n_ID=ALL_CELL_COUNTER; n_ID<ALL_CELL_COUNTER+N_SUBCELL_3; n_ID++)
	{
		CELL[n_ID].hsml_min = 1.0e10;
		if (CELL[n_ID].N_particles_contained > 0)
		{
			if(!(CELL[n_ID].particle_IDs = (int *) malloc(CELL[n_ID].N_particles_contained*sizeof(int))))
			{
			
				printf("Error allocating particle_IDs for cell %d in tree.c.\n",n_ID);
				exit(-1);
			}
			CELL[n_ID].N_particles_contained = 0;
	
		}
	}

	for (np=0; np<CELL[cell_ID].N_particles_contained; np++) {
		p_ID = CELL[cell_ID].particle_IDs[np];
		
		in_box_check = 1;
		for (j=0;j<3;j++) {
			n_x[j] = (int)((PG[p_ID].pos[j] - CELL[cell_ID].min_x[j])/(CELL[cell_ID].width/N_SUBCELL));
			if (n_x[j] < 0 || n_x[j] > N_SUBCELL-1) in_box_check = 0;
		}
		if (in_box_check != 0) {
			n_ID_sub = n_x[0]*N_SUBCELL_2 + n_x[1]*N_SUBCELL + n_x[2];
			n_ID = n_ID_sub + ALL_CELL_COUNTER;

			CELL[n_ID].particle_IDs[CELL[n_ID].N_particles_contained] = p_ID;
			CELL[n_ID].N_particles_contained += 1;
			if (PG[p_ID].hsml < CELL[n_ID].hsml_min) CELL[n_ID].hsml_min = PG[p_ID].hsml;
	}}

	n_subshell_i = ALL_CELL_COUNTER;
	n_subshell_f = ALL_CELL_COUNTER + N_SUBCELL_3;
	ALL_CELL_COUNTER += N_SUBCELL_3;
	
	//printf("Cell count = %d \n",ALL_CELL_COUNTER);
	for (n_ID=n_subshell_i; n_ID<n_subshell_f; n_ID++) build_cell(n_ID);	
	return;
}



/* Calculate the gas properties in a grid cell, given the cell ID */
void get_tree_localvals(int cell_ID)
{
	if (CELL[cell_ID].N_particles_contained == 0) {
		CELL[cell_ID].U.rho = 0.0;
		return;
	}
	
	/* for now, just a simple average -- easy to extend this to something 
		kernel-weighted later if we decide to do so */
	int i,p_ID;
	float rho,u,Z,ne,nh;

	rho=u=Z=ne=nh=0.0;
	for (i=0; i<CELL[cell_ID].N_particles_contained; i++) {
		p_ID = CELL[cell_ID].particle_IDs[i];

		rho += PG[p_ID].rho;
		u   += PG[p_ID].rho * PG[p_ID].u;
		Z   += PG[p_ID].rho * PG[p_ID].z;
		ne  += PG[p_ID].rho * PG[p_ID].ne;
		nh  += PG[p_ID].rho * PG[p_ID].nh;
	}
	u  /= rho;
	Z  /= rho;
	ne /= rho;
	nh /= rho;
	rho /= (float)(CELL[cell_ID].N_particles_contained);

	CELL[cell_ID].U.rho = rho;
	CELL[cell_ID].U.u   = u;
	CELL[cell_ID].U.Z   = Z;
	CELL[cell_ID].U.ne  = ne;
	CELL[cell_ID].U.neutral_frac = nh;
    get_two_phase_breakdown(u,rho,ne,&CELL[cell_ID].U);

	return;
}

/* Routine to locate the cell_id of the smallest cell a given point is inside */
int find_cell_from_scratch(float pos[3])
{
	// return -1 if outside of the largest cell
	if (fabs(pos[0]) >= MAX_DISTANCE_FROM0 ||
		fabs(pos[1]) >= MAX_DISTANCE_FROM0 ||
		fabs(pos[2]) >= MAX_DISTANCE_FROM0) return -1;

	int j, n_x[3], n_sub_ID, cell_ID;

	cell_ID = 0;
	while (CELL[cell_ID].sub_cell_check == 1) {
		for (j=0; j<3; j++) {
			n_x[j] = (int)((pos[j] - CELL[cell_ID].min_x[j])/(CELL[cell_ID].width/N_SUBCELL));
		}
		n_sub_ID = n_x[0]*N_SUBCELL_2 + n_x[1]*N_SUBCELL + n_x[2];
		cell_ID = CELL[cell_ID].sub_cell_IDs[n_sub_ID];
	}
	return cell_ID;
}

/* Routine to locate the cell_id of a new cell, starting from a given cell */
int find_cell_from_nearby_cell_plim(float cell_ID_0, float pos[3], float xmax[3], float min_pos[3], float D)
{
	double d_epsilon = 1.0e-1;
	double x, y, z;

	// return -1 when you reach the particle
	x = pos[0]-min_pos[0];
	y = pos[1]-min_pos[1];
	z = pos[2]-min_pos[2];
	//printf("cell_ID_0 %5.0f D %8.3f pos[0] %8.3f pos[1] %8.3f pos[2] %8.3f px[0] %8.3f px[1] %8.3f px[2] %8.3f x[0] %8.3f x[1] %8.3f x[2] %8.3f\n",cell_ID_0,D,pos[0],pos[1],pos[2],xmax[0],xmax[1],xmax[2],x,y,z);
	if (sqrt( x*x + y*y + z*z) >= D) return -1;

	int j, n_x[3], n_sub_ID, cell_ID, in_cell_check;


	in_cell_check = 0;
	cell_ID = (int) cell_ID_0; //changed to type cell_ID
	while (in_cell_check == 0) {
		cell_ID = CELL[cell_ID].parent_ID;
		if (cell_ID < 0) return -1;
	
		in_cell_check = 1;
		for (j=0; j<3; j++) {
			if (pos[j] > CELL[cell_ID].min_x[j] + CELL[cell_ID].width ||
				pos[j] < CELL[cell_ID].min_x[j]) in_cell_check = 0;
		}
	}
	
	while (CELL[cell_ID].sub_cell_check == 1) {
		for (j=0; j<3; j++) {
			n_x[j] = (int)((pos[j] - CELL[cell_ID].min_x[j])/(CELL[cell_ID].width/N_SUBCELL));
			if (n_x[j] < 0) n_x[j] = 0;
			if (n_x[j] > N_SUBCELL - 1) n_x[j] = N_SUBCELL - 1;
		}
		n_sub_ID = n_x[0]*N_SUBCELL_2 + n_x[1]*N_SUBCELL + n_x[2];
		cell_ID = CELL[cell_ID].sub_cell_IDs[n_sub_ID];
	}
	return cell_ID;
}

/* Routine to locate the cell_id of a new cell, starting from a given cell */
int find_cell_from_nearby_cell(float cell_ID_0, float pos[3])
{
	// return -1 if outside of the largest cell
	if (fabs(pos[0]) >= MAX_DISTANCE_FROM0 ||
		fabs(pos[1]) >= MAX_DISTANCE_FROM0 ||
		fabs(pos[2]) >= MAX_DISTANCE_FROM0) return -1;

	int j, n_x[3], n_sub_ID, cell_ID, in_cell_check;


	in_cell_check = 0;
	cell_ID = (int) cell_ID_0; //changed to type cell_ID
	while (in_cell_check == 0) {
		cell_ID = CELL[cell_ID].parent_ID;
		if (cell_ID < 0) return -1;
	
		in_cell_check = 1;
		for (j=0; j<3; j++) {
			if (pos[j] > CELL[cell_ID].min_x[j] + CELL[cell_ID].width ||
				pos[j] < CELL[cell_ID].min_x[j]) in_cell_check = 0;
		}
	}
	
	while (CELL[cell_ID].sub_cell_check == 1) {
		for (j=0; j<3; j++) {
			n_x[j] = (int)((pos[j] - CELL[cell_ID].min_x[j])/(CELL[cell_ID].width/N_SUBCELL));
			if (n_x[j] < 0) n_x[j] = 0;
			if (n_x[j] > N_SUBCELL - 1) n_x[j] = N_SUBCELL - 1;
		}
		n_sub_ID = n_x[0]*N_SUBCELL_2 + n_x[1]*N_SUBCELL + n_x[2];
		cell_ID = CELL[cell_ID].sub_cell_IDs[n_sub_ID];
	}
	return cell_ID;
}


/* Routine to carry a ray across a given tree cell, integrating the relevant quantities
	as it moves, GIVEN the cell the ray is inside */
void carry_ray_across_cell(Ray_struct *ray, int cell_id)
{
	int j;
	float dr_i,dr,dx,eps,rhodr;
	
	dx = CELL[cell_id].width;
	dr = 1.0e10;
	eps = 1.0e-2;
	for (j=0; j<3; j++) {
		if (ray->n_hat[j] != 0.0) dr_i=(CELL[cell_id].min_x[j] - ray->pos[j])/ray->n_hat[j];
		if ((dr_i > 0.0) && (dr_i < dr)) dr = dr_i;
		if (ray->n_hat[j] != 0.0) dr_i += dx / ray->n_hat[j];
		if ((dr_i > 0.0) && (dr_i < dr)) dr = dr_i;
	}
	if (dr > 3.0*dx) dr = dx;

	for (j=0;j<3;j++) ray->pos[j] += (dr + eps) * ray->n_hat[j];

	/* finally, integrate the addition to the ray along the line through the cell */
	rhodr = CELL[cell_id].U.rho * dr;
	ray->nh				+= rhodr;
	ray->nh_hot 		+= CELL[cell_id].U.rho_hot * rhodr;
	ray->Z				+= CELL[cell_id].U.Z * rhodr;
	ray->neutral_frac 	+= CELL[cell_id].U.neutral_frac * rhodr;
	return;
}

/* Simple routine to take a given ray and integrate it until it escapes
	from the simulation, given the cell of origin of the ray */
void integrate_ray_to_target(Ray_struct *ray, float min_pos[3], int origin_cell_id)
{

	/* First, make sure the ray is zeroed out */
	//ray->nh			  = 0.0;
	//ray->nh_hot		  = 0.0;
	ray->nh			  = 1.0e-40;
	ray->nh_hot		  = 1.0e-40;
	ray->Z			  = 0.0;
	ray->neutral_frac = 0.0;

	int cell_id = origin_cell_id;
	while (cell_id >= 0) {
		carry_ray_across_cell(ray,cell_id);
		cell_id = find_cell_from_nearby_cell_plim(cell_id,ray->pos,ray->xmax,min_pos,ray->D);
	}
	return;
}

/* Simple routine to take a given ray and integrate it until it escapes
	from the simulation, given the cell of origin of the ray */
void integrate_ray_to_escape(Ray_struct *ray, int origin_cell_id)
{
	/* First, make sure the ray is zeroed out */
	ray->nh			  = 1.0e-40;
	ray->nh_hot		  = 1.0e-40;
	ray->Z			  = 0.0;
	ray->neutral_frac = 0.0;

	int cell_id = origin_cell_id;
	while (cell_id >= 0) {
		carry_ray_across_cell(ray,cell_id);
		cell_id = find_cell_from_nearby_cell(cell_id,ray->pos);
	}
	return;
}




