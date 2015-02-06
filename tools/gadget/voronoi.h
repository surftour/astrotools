#ifndef HAVE_H_VORONOI
#define HAVE_H_VORONOI



#define STACKSIZE_TETRA 1000


#ifndef TWODIMS

#define EDGE_0 1     /* points 0-1 */
#define EDGE_1 2     /* points 0-2 */
#define EDGE_2 4     /* points 0-3 */
#define EDGE_3 8     /* points 1-2 */
#define EDGE_4 16    /* points 1-3 */
#define EDGE_5 32    /* points 2-3 */
#define EDGE_ALL 63

#else

#define EDGE_0 1     /* points 1-2 */
#define EDGE_1 2     /* points 0-2 */
#define EDGE_2 4     /* points 0-1 */
#define EDGE_ALL 7

#endif


#define HSML_INCREASE_FACTOR 1.3


#ifdef TWODIMS    /* will only be compiled in 2D case */
#define DIMS 2
#else
#define DIMS 3
#endif


#define REFL_X_FLAGS 115043766
#define REFL_Y_FLAGS 132379128
#define REFL_Z_FLAGS 115043766


extern struct primexch
{
  MyFloat Pressure;
  MyFloat VelPred[3];
  MyFloat Density;
  MyFloat Mass;
  MyFloat Entropy;
#ifdef VORONOI_MESHRELAX
  MyFloat HydroAccel[3];
  MyFloat Center[3];
#endif
  int task, index;
} 
*PrimExch;


typedef struct 
{
  double x, y, z;
  double xx, yy, zz;
  int task, index;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  unsigned int image_flags;
#endif
  int ix, iy, iz;
}
point;


typedef struct tetra_data
{
  point                 *p[DIMS+1];  /* oriented tetrahedron points */
  struct tetra_data     *t[DIMS+1];  /* adjacent tetrahedrons, always opposite to corresponding point */
  unsigned char          s[DIMS+1];  /* gives the index of the point in the adjacent tetrahedron that
                                        lies opposite to the common face */
  point                 c;  /* describes circumcircle */
  char deleted;
  unsigned char egde_visited;
}
tetra;


typedef struct face_data
{
  point *p1, *p2;
  double area;
  double cx, cy, cz;  /* center-of-mass of face */

#ifdef VORONOI_MESHRELAX
  double mflux;
  double qflux_x, qflux_y, qflux_z;
  double eflux;
#endif
}
face;

#ifdef VORONOI_MESHRELAX
extern struct grad_data
{
  double drho[3];
  double dvel[3][3];
  double dpress[3];
}
*Grad, *GradExch;
#endif

extern struct list_export_data
{
  unsigned int image_bits;
  int origin, index;
  struct list_export_data *nextexport;
}
*ListExports;

extern int Ninlist, MaxNinlist;

extern struct list_P_data
{
  struct list_export_data *firstexport;
  struct list_export_data *currentexport;
} *List_P;

extern int CountInSphereTests, CountInSphereTestsExact;
extern int CountConvexEdgeTest, CountConvexEdgeTestExact;
extern int CountFlips, Count_1_to_3_Flips2d, Count_2_to_4_Flips2d;
extern int Count_1_to_4_Flips, Count_2_to_3_Flips, Count_3_to_2_Flips, Count_4_to_4_Flips;
extern int Count_EdgeSplits, Count_FaceSplits;
extern int Count_InTetra, Count_InTetraExact;

extern int Ninlist, MaxNinlist;
extern int Flag_MaxNinlistReached;

extern int Ndp;			/* number of delaunay points */
extern int MaxNdp;			/* maximum number of delaunay points */
extern point *DP;			/* delaunay points */

extern int Ndt;
extern int MaxNdt;			/* number of delaunary tetrahedra */
extern tetra *DT;			/* Delaunay tetrahedra */


extern int Nvf;			/* number of Voronoi faces */
extern int MaxNvf;			/* maximum number of Voronoi faces */
extern face *VF;			/* Voronoi faces */


extern point *DPinfinity;

extern double CentralOffsetX, CentralOffsetY, CentralOffsetZ, ConversionFac;


void sample_solution_isothermal3d(double S, 
                                  double Rho_L, double Vel_L, double Vel_Ly, double Vel_Lz, 
                                  double Rho_R, double Vel_R, double Vel_Ry, double Vel_Rz,
                                  double Rho, double Vel, double *Rho_flux,
                                  double *Vel_flux, double *Vel_y_flux, double *Vel_z_flux);

void riemann_isotherm(double Rho_L, double Vel_L, double Rho_R, double Vel_R, 
                        double *Rho, double *Vel);
void isothermal_function(double rhostar, double rho, double *F, double *FD);

void voronoi_exchange_ghost_variables(void);
void voronoi_density(void);
void voronoi_hydro_force(void);

void voronoi_meshrelax_update_particles(void);

void limit_gradient(double *d, double phi, double min_phi, double max_phi, double *dphi);
void voronoi_meshrelax_calc_displacements(void);
void voronoi_meshrelax_calculate_gradients(void);
void voronoi_meshrelax_fluxes(void);
void voronoi_meshrelax_drift(void);
void voronoi_meshrelax(void);

void voronoi_mesh(void);
void write_voronoi_mesh(char *fname, int writeTask, int lastTask);
void drift_mesh_generators(void);
void do_hydro_step(void);
void write_voronoi_mesh(char *fname, int writeTask, int lastTask);
void finalize_output_of_voronoi_geometry(void);
void prepare_output_of_voronoi_geometry(void);
void initialize_and_create_first_tetra(void);
void compute_voronoi_faces_and_volumes(void);
void voronoi_exchange_gradients(void);
void process_edge_faces_and_volumes(tetra * t, int nr);
tetra *insert_point(point * p, tetra * tstart);

void make_an_edge_split(tetra *t0, int edge_nr, int count, point *p, tetra **tlist);
void make_a_face_split(tetra *t0, int face_nr, point *p, 
                       tetra *t1, tetra *t2, tetra *q1, tetra *q2);
double calculate_tetra_volume(point * p0, point * p1, point * p2, point * p3);
void make_a_4_to_4_flip(tetra * t, int tip_index, int edge_nr);

void make_a_1_to_4_flip(point * p, tetra * t0, tetra * t1, tetra * t2, tetra * t3);
void make_a_3_to_2_flip(tetra * t0, tetra * t1, tetra * t2, int tip, int edge, int bottom);
void make_a_2_to_3_flip(tetra * t0, int tip, tetra * t1, int bottom, point * qq, tetra * t2);
tetra *get_tetra(point * p, int *moves, tetra * tstart, int *flag, int *edgeface_nr);
int InTetra(tetra *t, point *p, int *edgeface_nr, tetra **nexttetra);
double InSphere(point * p0, point * p1, point * p2, point * p3, point * p);
void update_circumcircle(tetra * t);
int test_tetra_orientation(point * p0, point * p1, point * p2, point * p3);
int test_intersect_triangle(point * p0, point * p1, point * p2, point * q, point * s);
double deter4(point * p0, point * p1, point * p2, point * p3);
double deter3(point * p0, point * p1, point * p2);
double deter4_orient(point * p0, point * p1, point * p2, point * p3);
double determinante3(double *a, double *b, double *c);
double deter_special(double *a, double *b, double *c, double *d);
int voronoi_get_additional_points(void);
int voronoi_exchange_evaluate(int target, int mode, int *nexport, int *nsend_local);
int ngb_treefind_voronoi(MyDouble searchcenter[3], MyFloat hsml, int target, int origin,
			 int *startnode, int mode, int *nexport, int *nsend_local);
void compute_circumcircles(void);
int compute_max_delaunay_radius(void);
void check_for_min_distance(void);
void check_links(void);
void check_orientations(void);
void check_tetras(int npoints);
int voronoi_get_local_particles(void);
void check_for_vertex_crossings(void);

int convex_edge_test(tetra *t, int tip, int *edgenr);
void get_circumcircle_exact(tetra * t, double *x, double *y, double *z);

void do_hydro_calculations(void);
void update_primitive_variables(void);
void update_cells_with_fluxes(void);
void calculate_green_gauss_gradients(void);
void compute_interface_fluxes(void);
void half_step_evolution(void);
void limit_gradient(double *d, double phi, double min_phi, double max_phi, double *dphi);
int analytic_riemann_solution(int argc, char *argv[]);
void godunov_flux(double Rho_L, double Vel_L, double Press_L,
		  double Rho_R, double Vel_R, double Press_R,
		  double *Rho_flux, double *Vel_flux, double *Press_flux);

void godunov_flux_2d(double Rho_L, double Vel_L, double Vel_Ly, double Press_L,
		     double Rho_R, double Vel_R, double Vel_Ry, double Press_R,
		     double *Rho_flux, double *Vel_flux, double *Vel_y_flux, double *Press_flux);

void sample_solution_2d(double S, double Rho_L, double Vel_L, double Vel_Ly, double Press_L, double Csnd_L,
			double Rho_R, double Vel_R, double Vel_Ry, double Press_R, double Csnd_R,
			double Press, double Vel, double *Rho_flux, double *Vel_flux, double *Vel_y_flux,
			double *Press_flux);

void sample_solution(double S, double Rho_L, double Vel_L, double Press_L, double Csnd_L,
		     double Rho_R, double Vel_R, double Press_R, double Csnd_R,
		     double Press, double Vel, double *Rho_flux, double *Vel_flux, double *Press_flux);

void riemann(double Rho_L, double Vel_L, double Press_L, double Csnd_L,
	     double Rho_R, double Vel_R, double Press_R, double Csnd_R, double *Press, double *Vel);

void pressure_function(double P, double Rho, double Press, double Csnd, double *F, double *FD);

double guess_for_pressure(double Rho_L, double Vel_L, double Press_L, double Csnd_L,
			  double Rho_R, double Vel_R, double Press_R, double Csnd_R);

void godunov_flux_3d(double Rho_L, double Vel_L, double Vel_Ly, double Vel_Lz, double Press_L,
		     double Rho_R, double Vel_R, double Vel_Ry, double Vel_Rz, double Press_R,
		     double *Rho_flux, double *Vel_flux, double *Vel_y_flux, double *Vel_z_flux, double *Press_flux);

void sample_solution_3d(double S, double Rho_L, double Vel_L, double Vel_Ly, double Vel_Lz, double Press_L, double Csnd_L,
			double Rho_R, double Vel_R, double Vel_Ry, double Vel_Rz, double Press_R, double Csnd_R,
			double Press, double Vel, double *Rho_flux, double *Vel_flux, double *Vel_y_flux, double *Vel_z_flux,
			double *Press_flux);


void voronoi_exchange_primitive_variables(void);
void voronoi_update_primitive_variables_and_exchange_gradients(void);

int compare_primexch(const void *a, const void *b);



/* 2D voronoi routines */
void check_edge_and_flip_if_needed(point *p, tetra *t);
tetra *get_triangle(point * p, int *moves, int *degenerate_flag, tetra * tstart);
int InTriangle(point * p0, point * p1, point * p2, point * p);
double InCircle(point * p0, point * p1, point * p2, point * p);	
double v2d_deter3(point * p0, point * p1, point * p2);
void make_a_1_to_3_flip(point * p, tetra * t0, tetra * t1, tetra * t2);
double test_triangle_orientation(point * p0, point * p1, point * p2);
void get_circle_center(point * a, point * b, point * c, double *x, double *y);
void do_special_dump(int num);

void make_a_2_to_4_flip(point * p, tetra * t0, tetra * t1, 
                        tetra * t2, tetra * t3, int i0, int j0);

void set_vertex_velocities(void);
void dump_points(void);
int voronoi_bitcount(unsigned int x); 

void set_integers_for_point(point *p);

int solve_linear_equations(double *m, double *res);


void check_triangles(int npoints);


int InCircle_Quick(point * p0, point * p1, point * p2, point * p);
int InCircle_Errorbound(point * p0, point * p1, point * p2, point * p);
int InCircle_Exact(point * p0, point * p1, point * p2, point * p);

int Orient2d_Exact(point * p0, point * p1, point * p2);
int Orient2d_Quick(point * p0, point * p1, point * p2);


int FindTriangle(tetra *t, point *p, int *degnerate_flag, tetra **nexttetra);

int InSphere_Exact(point * p0, point * p1, point * p2, point * p3, point * p);
int InSphere_Quick(point * p0, point * p1, point * p2, point * p3, point * p);
int InSphere_Gauss(point * p0, point * p1, point * p2, point * p3, point * p);
int InSphere_Errorbound(point * p0, point * p1, point * p2, point * p3, point * p);

int Orient3d_Exact(point * p0, point * p1, point * p2, point * p3);
int Orient3d_Quick(point * p0, point * p1, point * p2, point * p3);



#endif
