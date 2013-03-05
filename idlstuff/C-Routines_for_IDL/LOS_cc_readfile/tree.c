#include "overhead.h"
#include "proto.h"
#include "ngbtree3d.h"
/* 
 * Routines for building and using a tree-based grid to compute local
 *   gas properties with order h_sml resolution in a quick manner, and 
 *   to quickly trace rays out of the simulation in calculating column
 *   densities
 *   
 */ 

/*YL add*/
void free_the_tree(void) { 
  extern struct CELL_STRUCT *CELL; 
  int ncell, i; 

  //  ncells = All.N_gas*((int)*(2.0/TREE_FRAC_HSML)); 
  //  for (i = 0; i < ncells; i++) { 
  for (i=1; i < TOTAL_NUMBER_OF_CELLS; i++){
    if (CELL[i].N_particles_contained > 0) free(CELL[i].particle_IDs); 
  }
  free(CELL);

}


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
       	
	np = All.N_gas*(int)(4/TREE_FRAC_HSML);
	if(!(CELL=malloc(np*sizeof(struct CELL_STRUCT))))
	  {
	    fprintf(stderr,"failed to allocate memory.\n");
	    exit(0);
	  }
	printf("Number of cells allocated: %d\n", np);

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
	
	//printf("check 3 \n");
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

	//printf("check 4 \n");
	for (n_cell=ALL_CELL_COUNTER; n_cell<ALL_CELL_COUNTER+N_SUBCELL_3; n_cell++) {
		CELL[n_cell].hsml_min = 1.0e10;
		if (CELL[n_cell].N_particles_contained > 0) {
		//printf("%i \n",CELL[n_cell].N_particles_contained);
		//int *testmalloc;
		//testmalloc = (int *)malloc(3*sizeof(int));
		//printf("test survived \n");
		CELL[n_cell].particle_IDs = malloc(CELL[n_cell].N_particles_contained*sizeof(int));
		CELL[n_cell].N_particles_contained = 0;
	}}

	//printf("check 5 \n");
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
	
	printf("Tree opened, building... \n");
	for (n_cell=n_subshell_i; n_cell<n_subshell_f; n_cell++) build_cell(n_cell);
	printf("Tree constructed, %d total cells \n",ALL_CELL_COUNTER);
	TOTAL_NUMBER_OF_CELLS = ALL_CELL_COUNTER;
		
	for (n_cell=0; n_cell<TOTAL_NUMBER_OF_CELLS; n_cell++) get_tree_localval_derivs(n_cell);
	printf(" ... Gradients calculated \n");
	return;
}


/* The workhorse of the tree-building phase, this loops over itself to 
	build a recursive tree over the grid and sets the gas values in the 
	tree when a small enough scale is reached */
void build_cell(int cell_ID)
{
	//printf("CELL NO = %i \n",ALL_CELL_COUNTER);
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

	for (n_ID=ALL_CELL_COUNTER; n_ID<ALL_CELL_COUNTER+N_SUBCELL_3; n_ID++) {
		CELL[n_ID].hsml_min = 1.0e10;
		if (CELL[n_ID].N_particles_contained > 0) {
		CELL[n_ID].particle_IDs = malloc(CELL[n_ID].N_particles_contained*sizeof(int));
		CELL[n_ID].N_particles_contained = 0;
	}}

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
	int i,p_ID,ii;
	double rho,u,Z,ne,nh,r,hinv,hinv3,hinv4,x,mass_j,Tmp,wk,wt_tot,vx,vy,vz;
  	int   Num_neighbors = (int)NUMBER_OF_NEIGHBORS;
	float *r2_list;
	int *ngb_list;

    float r0[3];
    for (i=0; i<3; i++) r0[i] = CELL[cell_ID].min_x[i] + 0.5*CELL[cell_ID].width;
	/* old brute-force routine to locate the nearest neighbors */
    //find_neighbors(r0,Num_neighbors,ngb_list,r2_list,&hsml);

    float h,h2,hsml=10.0;		/* local smoothing length (to be returned) */
    float Hmax= 10000.0;
    int ngbfound;
    //h= 10.0;
    //Num_neighbors=64;
  	/* use the GADGET tree to find the neighbors */
	//printf("Heading into ngb3d_treefind - I hope this works.\n"); fflush(stdout);
	//printf(" %f %f %f %i %f %f \n",r0[0],r0[1],r0[2],Num_neighbors,h,Hmax);
    h2 = ngb3d_treefind( r0, Num_neighbors,1.04*h,&ngb_list,&r2_list, Hmax,&ngbfound);
/*
    int particle_id = 0;
    ngb_list= (int *)malloc(Num_neighbors*sizeof(int));
    r2_list = (float *)malloc(Num_neighbors*sizeof(float));
    find_neighbors(r0,Num_neighbors,ngb_list,r2_list,&hsml);
    ngbfound= Num_neighbors;
    h2 = hsml;
	printf("Are we getting this back? h= %g\n",sqrt(h2)); fflush(stdout);
    if (Num_neighbors == 0) printf("  neighbor check = %g %i \n",h2,Num_neighbors);
*/
  	hsml = sqrt(h2);
	hinv = 1.0/hsml; hinv3 = hinv*hinv*hinv; hinv4 = hinv3*hinv;
    
	rho=u=Z=ne=nh=Tmp=wt_tot=vx=vy=vz=0.0;
	//for (i=0; i<CELL[cell_ID].N_particles_contained; i++) {
	//	p_ID = CELL[cell_ID].particle_IDs[i];
	
	//printf(" Found %i Neighbors \n",ngbfound);
	
  	for (i=0; i<ngbfound; i++) {
    	p_ID = ngb_list[i];
    	r = sqrt(r2_list[i]);
		if (r>hsml) r=hsml;
    	x = r * hinv;
    	ii = (int)(x*KERNEL_TABLE);
	        wk =hinv3 *( Kernel[ii]  + (Kernel[ii+1]-Kernel[ii])*(x-KernelRad[ii])*KERNEL_TABLE);
	  	//dwk=hinv4*( KernelDer[ii] + (KernelDer[ii+1]-KernelDer[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
	  	mass_j = PG[p_ID].mass;

		//printf(" %e \n",mass_j);
	  	rho += (double)mass_j * (double)wk;
	  	//Tmp	+= PG[p_ID].temp * wk;
	  
	  	/* Calculate mass-weighted averages of Z, ne, u, and the ionized fraction, 
	  	 *   noting the definitions of the read-in quantities from the snapshots */
	  	u	+= 	PG[p_ID].u   *  mass_j * wk;
	  	Z 	+= 	PG[p_ID].z   *  mass_j * wk;
      	ne  += 	PG[p_ID].ne  *  mass_j * wk;
      	nh 	+=	PG[p_ID].nh  *  mass_j * wk;
	  	Tmp +=  PG[p_ID].temp   *  mass_j * wk;
	  	vx	+=  PG[p_ID].vel[0] *  mass_j * wk;
	  	vy	+=  PG[p_ID].vel[1] *  mass_j * wk;
	  	vz	+=  PG[p_ID].vel[2] *  mass_j * wk;
	}
	//free(r2_list);
	//free(ngb_list);
	//printf(" rho = %e \n",rho);
    double epsilon = 1.0e-15;
    rho += epsilon;
	//printf("   rho = %e %e %e %e \n",rho,psi,hinv,hinv3);

	u  /= rho;
	Z  /= rho;
	ne /= rho;
	nh /= rho;
	Tmp/= rho;
	vx /= rho;
	vy /= rho;
	vz /= rho;
	//rho /= (float)(CELL[cell_ID].N_particles_contained);
	//printf("  rho = %e %e %e \n",rho,u,ne);

	CELL[cell_ID].U.h     = hsml;
	CELL[cell_ID].U.rho   = rho;
	CELL[cell_ID].U.P_eff = GAMMA_MINUS1*(u*rho);
	CELL[cell_ID].U.u     = u;
	CELL[cell_ID].U.T     = Tmp;
	CELL[cell_ID].U.Z     = Z;
	CELL[cell_ID].U.ne    = ne;
	CELL[cell_ID].U.neutral_frac = nh;
	CELL[cell_ID].U.vel[0] = vx;
	CELL[cell_ID].U.vel[1] = vy;
	CELL[cell_ID].U.vel[2] = vz;
     
	get_two_phase_breakdown(u,rho,ne,&CELL[cell_ID].U);
	//printf("   --> %e %e %e \n",CELL[cell_ID].U.rho_hot,CELL[cell_ID].U.rho_cold,CELL[cell_ID].U.u_hot);
	return;
}


/* Calculate the local gradient in the gas properties in a grid cell, given the cell ID 
	Note that this method is not strictly correct, as we take the nearest cells in 
	each direction, which may not be perfectly aligned as the cell sizes are adaptive! 
	Still, it should serve as a first-order (actually, second order) correction */
void get_tree_localval_derivs(int cell_ID)
{
	/* first do a subcell check */
	if (CELL[cell_ID].sub_cell_check == 1) return;

	int nbr_cells[2];
	int i,j,p_ID;
	float pos[3];
	float n_nbr_in_grid,dx,wt,dr;
	double eps = 1.0e-10;
	double max_gradient = 10.0;	// maximum allowed gradient in d(log(Y))/dX
	double temp_grad = 0.0;
    CELL[cell_ID].dU = (LOCALVAL *)calloc(3,sizeof(LOCALVAL));
	for (j=0; j<3; j++) {	
		for (i=0; i<3; i++) pos[i] = CELL[cell_ID].min_x[i] + 0.5*CELL[cell_ID].width;
		dr = 0.01;
		pos[j] = CELL[cell_ID].min_x[j] - dr * CELL[cell_ID].width;
			nbr_cells[0] = find_cell_from_nearby_cell(cell_ID,pos);
			while (nbr_cells[0] == cell_ID) {
				dr += dr;
				pos[j] = CELL[cell_ID].min_x[j] - dr * CELL[cell_ID].width;
				nbr_cells[0] = find_cell_from_nearby_cell(cell_ID,pos);		
			}
		dr = 0.01;
		pos[j] = CELL[cell_ID].min_x[j] + (1.0+dr) * CELL[cell_ID].width;
			nbr_cells[1] = find_cell_from_nearby_cell(cell_ID,pos);
			while (nbr_cells[1] == cell_ID) {
				dr += dr;
				pos[j] = CELL[cell_ID].min_x[j] + (1.0+dr) * CELL[cell_ID].width;
				nbr_cells[1] = find_cell_from_nearby_cell(cell_ID,pos);		
			}
		n_nbr_in_grid = 0.0;
		for (i=0; i<2; i++) if (nbr_cells[i] >= 0) n_nbr_in_grid += 1.0;
		pos[j] = CELL[cell_ID].min_x[j] + 0.5*CELL[cell_ID].width;
		//printf(" %i %i \n",nbr_cells[0],nbr_cells[1]);

		if (USE_CELL_GRADIENTS) {	/* USE_CELL_GRADIENTS is turned ON */
		for (i=0; i<2; i++) {
		p_ID = nbr_cells[i];
		if (p_ID >= 0) {
			dx = pos[j] - (CELL[p_ID].min_x[j] + 0.5*CELL[p_ID].width);
			wt = (1./dx) * (1./n_nbr_in_grid);
			//printf(" wt = %e \n",wt);
			/* now set each quantity - do the gradients in log to prevent ugly divergences */   
			temp_grad = (log((CELL[cell_ID].U.rho+eps)/(CELL[p_ID].U.rho+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].rho += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.P_eff+eps) /(CELL[p_ID].U.P_eff+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].P_eff += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.u+eps) /(CELL[p_ID].U.u+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].u += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.T+eps) /(CELL[p_ID].U.T+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].T += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.Z+eps+0.002) /(CELL[p_ID].U.Z+eps+0.002)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].Z += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.ne+eps) /(CELL[p_ID].U.ne+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].ne += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.neutral_frac+eps) /(CELL[p_ID].U.neutral_frac+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].neutral_frac += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.u_hot+eps) /(CELL[p_ID].U.u_hot+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].u_hot += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.cloud_massfrac+eps)/(CELL[p_ID].U.cloud_massfrac+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].cloud_massfrac += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.cloud_fillingfactor+eps)/(CELL[p_ID].U.cloud_fillingfactor+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].cloud_fillingfactor += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.rho_hot+1.0e-8+eps)/(CELL[p_ID].U.rho_hot+1.0e-8+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].rho_hot += temp_grad;
			temp_grad = (log((CELL[cell_ID].U.rho_cold+eps)/(CELL[p_ID].U.rho_cold+eps)));
				temp_grad *= wt;
				if (temp_grad > max_gradient) temp_grad = max_gradient;
				if (temp_grad < -max_gradient) temp_grad = -max_gradient;
				CELL[cell_ID].dU[j].rho_cold += temp_grad;
			//printf("%e %e %f %e \n",dx,wt,n_nbr_in_grid,CELL[cell_ID].dU[j].rho);
		}}
		}
	}
	return;
}




/* Routine to locate the cell_id of the smallest cell a given point is inside */
int find_cell_from_scratch(float pos[3])
{
	// return -1 if outside of the largest cell
	if (abs(pos[0]) >= MAX_DISTANCE_FROM0 ||
		abs(pos[1]) >= MAX_DISTANCE_FROM0 ||
		abs(pos[2]) >= MAX_DISTANCE_FROM0) return -1;

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
int find_cell_from_nearby_cell(float cell_ID_0, float pos[3])
{
	// return -1 if outside of the largest cell
	if (abs(pos[0]) >= MAX_DISTANCE_FROM0 ||
		abs(pos[1]) >= MAX_DISTANCE_FROM0 ||
		abs(pos[2]) >= MAX_DISTANCE_FROM0) return -1;

	int j, n_x[3], n_sub_ID, cell_ID, in_cell_check;


	in_cell_check = 0;
	cell_ID = cell_ID_0;
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
	if (dr > 0.5*dx) dr = dx;

	for (j=0;j<3;j++) ray->pos[j] += (dr + eps) * ray->n_hat[j];


	if (USE_CELL_GRADIENTS==0) {	/* USE_CELL_GRADIENTS is turned OFF */
	/* finally, integrate the addition to the ray along the line through the cell */
	rhodr = CELL[cell_id].U.rho * (dr+eps);
	ray->nh				+= rhodr;
	ray->nh_hot 		+= CELL[cell_id].U.rho_hot * rhodr;
	ray->Z				+= CELL[cell_id].U.Z * rhodr;
	ray->neutral_frac 	+= CELL[cell_id].U.neutral_frac * rhodr;
	}


	if (USE_CELL_GRADIENTS) {	/* USE_CELL_GRADIENTS is turned ON */
	/* replace that with an interpolation using a *ROUGH* estimate of the 
	    local gradients in each quantity within the cell to integrate over the cell */
	   /* simple struct to hold the relevant quantities for each quantity */
	float beta_rho = 0.0;
	float beta = 0.0;
	for (j=0; j<3; j++) beta_rho += CELL[cell_id].dU[j].rho * ray->n_hat[j];
	if (dr == 0. && beta_rho == 0.) return;
	if (beta_rho == 0.) beta_rho = 1.0e-10 / dr;
	//printf(" %e %e %e %e %e \n",beta_rho,((1./beta_rho)*(exp(beta_rho*dr)-1.)),CELL[cell_id].dU[0].rho,CELL[cell_id].dU[1].rho,CELL[cell_id].dU[2].rho);
	ray->nh += CELL[cell_id].U.rho * (1./beta_rho) * (exp(beta_rho*dr) - 1.);

	/* simple calc. for each new quantity, just need to plug it in a couple spots
	 		--- this is to calculate it weighted by density */

	beta = beta_rho;	/* don't change */
		for (j=0; j<3; j++) beta += CELL[cell_id].dU[j].rho_hot * ray->n_hat[j]; /* plug the derivative needed in here */
		if (dr > 0.) {
		if (abs(beta) < 1.0e-10) beta = 1.0e-10 / dr;
		beta = (1./beta) * (exp(beta*dr)-1.) * CELL[cell_id].U.rho; /* don't change */
		ray->nh_hot			+= CELL[cell_id].U.rho_hot * beta; /* plug the value in here */
		}
	beta = beta_rho;
		for (j=0; j<3; j++) beta += CELL[cell_id].dU[j].Z * ray->n_hat[j];
		if (dr > 0.) {
		if (abs(beta) < 1.0e-10) beta = 1.0e-10 / dr;
		beta = (1./beta) * (exp(beta*dr)-1.) * CELL[cell_id].U.rho;
		ray->Z				+= CELL[cell_id].U.Z * beta;
		}
	beta = beta_rho;
		for (j=0; j<3; j++) beta += CELL[cell_id].dU[j].neutral_frac * ray->n_hat[j];
		if (dr > 0.) {
		if (abs(beta) < 1.0e-10) beta = 1.0e-10 / dr;
		beta = (1./beta) * (exp(beta*dr)-1.) * CELL[cell_id].U.rho;
		ray->neutral_frac	+= CELL[cell_id].U.neutral_frac * beta;
		}
	}

	//printf("%d .. ",CELL[cell_id].N_particles_contained);
	return;
}


/* Simple routine to take a given ray and integrate it until it escapes
	from the simulation, given the cell of origin of the ray */
void integrate_ray_to_escape(Ray_struct *ray, int origin_cell_id)
{
	/* First, make sure the ray is zeroed out */
	ray->nh			  = 1.0e-28;
	ray->nh_hot		  = 1.0e-28;
	ray->Z			  = 1.0e-28;
	ray->neutral_frac = 0.0;

	if (USE_FULL_NEIGHBOR_CALC == 0) {
	/* in this case, use the tree-mode to calculate everything */
	
	//int counter = 0;
	int cell_id = origin_cell_id;
	while (cell_id >= 0) {
		carry_ray_across_cell(ray,cell_id);
		cell_id = find_cell_from_nearby_cell(cell_id,ray->pos);
			//printf("%e %e \n",ray->nh_hot,ray->Z);	    
	//counter++;
	}
	//printf("%i \n",counter);

	} else {
	/* this means no tree code initialized - want to use the full (direct) 
	 *    nearest neighbor calculation along the ray */
	    //printf(" Calculating a ray ... \n");
	line_of_sight(ray); 
	//integrate_ray(ray);	
	}
	return;
}




