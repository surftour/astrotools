
#ifndef ALLVARS_H
#include "allvars.h"
#endif
#include "forcetree.h"
#include "cooling.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif
#ifdef LT_STELLAREVOLUTION
#include "lt.h"
#endif

#ifdef HAVE_HDF5
#include <hdf5.h>
void write_header_attributes_in_hdf5(hid_t handle);
void read_header_attributes_in_hdf5(char *fname);
void write_parameters_attributes_in_hdf5(hid_t handle);
void write_units_attributes_in_hdf5(hid_t handle);
void write_constants_attributes_in_hdf5(hid_t handle);
#endif
void report_VmRSS(void);

void read_fof(int num);
int fof_compare_ID_list_ID(const void *a, const void *b);

void myfree_msg(void *p, char *msg);

int snIaheating_evaluate(int target, int mode, int *nexport, int *nSend_local);
void snIa_heating(void);

int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                       int dest, int sendtag, void *recvbufreal, int recvcount,
                       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status);

int MPI_Sizelimited_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
			     int dest, int sendtag, void *recvbuf, int recvcount,
			     MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status);

void compare_partitions(void);
void assign_unique_ids(void);
int permut_data_compare(const void *a, const void *b);
void  generate_permutation_in_active_list(void);
void get_particle_numbers(char *fname, int num_files);

void conduction(void);
void conduction_matrix_multiply(double *in, double *out);
double conduction_vector_multiply(double *a, double *b);
int conduction_evaluate(int target, int mode, double *in, double *out, double *sum,
			int *nexport, int *nsend_local);


void fof_get_group_center(double *cm, int gr);
void fof_get_group_velocity(double *cmvel, int gr);
int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local);
void fof_compute_group_properties(int gr, int start, int len);

void parallel_sort(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *));
int compare_IDs(const void *a, const void *b);
void test_id_uniqueness(void);

int io_compare_P_ID(const void *a, const void *b);
int io_compare_P_GrNr_SubNr(const void *a, const void *b);


void drift_particle(int i, int time1);
int ShouldWeDoDynamicUpdate(void);

void put_symbol(double t0, double t1, char c);
void write_cpu_log(void);
void *mymalloc_msg(size_t n, char *message);

int get_timestep_bin(int ti_step);

const char* svn_version(void);

void find_particles_and_save_them(int num);
void lineofsight_output(void);
void sum_over_processors_and_normalize(void);
void absorb_along_lines_of_sight(void);
void output_lines_of_sight(int num);
int find_next_lineofsighttime(int time0);
int find_next_gridoutputtime(int ti_curr);
void add_along_lines_of_sight(void);
double los_periodic(double x);
void do_the_kick(int i, int tstart, int tend, int tcurrent);


void x86_fix(void) ;

void *mymalloc(size_t n);
void myfree(void *p);
void *myrealloc(void *p, size_t n);
void mymalloc_init(void);

double get_shear_viscosity(int i);

void kinetic_feedback_mhm(void);
int kin_compare_key(const void *a, const void *b);
void kinetic_evaluate(int target, int mode);

void bubble(void);
void multi_bubbles(void);
void bh_bubble(double bh_dmass, MyFloat center[3], int BH_id);
void find_CM_of_biggest_group(void);
int compare_length_values(const void *a, const void *b);
double rho_dot(double z, void *params);
double bhgrowth(double z1, double z2);

void second_order_ics(void);


int smoothed_evaluate(int target, int mode, int *nexport, int *nsend_local);
void smoothed_values(void);

int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local);

double INLINE_FUNC hubble_function(double a);
#ifdef DARKENERGY
double DarkEnergy_a(double);
double DarkEnergy_t(double);
#ifdef TIMEDEPDE
void fwa_init(void);
double INLINE_FUNC fwa(double);
double INLINE_FUNC get_wa(double);
#ifdef TIMEDEPGRAV
double INLINE_FUNC dHfak(double a);
double INLINE_FUNC dGfak(double a);
#endif
#ifdef EXTERNALHUBBLE
double INLINE_FUNC hubble_function_external(double a);
#endif
#endif

#endif

void blackhole_accretion(void);
int blackhole_evaluate(int target, int mode, int *nexport, int *nsend_local);
int blackhole_evaluate_swallow(int target, int mode, int *nexport, int *nsend_local);

int  blackhole_compare_key(const void *a, const void *b);




void fof_fof(int num);
void fof_import_ghosts(void);
void fof_course_binning(void);
void fof_find_groups(void);
void fof_check_cell(int p, int i, int j, int k);
void fof_find_minids(void);
int fof_link_accross(void);
void fof_exchange_id_lists(void);
int fof_grid_compare(const void *a, const void *b);
void fof_compile_catalogue(void);
void fof_save_groups(int num);
void fof_save_local_catalogue(int num);
double fof_periodic(double x);
double fof_periodic_wrap(double x);
void fof_find_nearest_dmparticle(void);
int fof_find_nearest_dmparticle_evaluate(int target, int mode, int *nexport, int *nsend_local);

int fof_compare_key(const void *a, const void *b);
void fof_link_special(void);
void fof_link_specialpair(int p, int s);
void fof_make_black_holes(void);

int io_compare_P_GrNr_ID(const void *a, const void *b);

void write_file(char *fname, int readTask, int lastTask);

void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last);

int get_values_per_blockelement(enum iofields blocknr);

int get_datatype_in_block(enum iofields blocknr);
void get_dataset_name(enum iofields blocknr, char *buf);


int blockpresent(enum iofields blocknr);
void fill_write_buffer(enum iofields blocknr, int *pindex, int pc, int type);
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);

int get_particles_in_block(enum iofields blocknr, int *typelist);

int get_bytes_per_blockelement(enum iofields blocknr, int mode);

void read_file(char *fname, int readTask, int lastTask);

void get_Tab_IO_Label(enum iofields blocknr, char *label);


void long_range_init_regionsize(void);

int find_files(char *fname);

int metals_compare_key(const void *a, const void *b);
void enrichment_evaluate(int target, int mode);

int hydro_evaluate(int target, int mode, int *nexport, int *nsend_local);

void pm_init_nonperiodic_allocate(void);

void  pm_init_nonperiodic_free(void);

double get_random_number(unsigned int id);
void set_random_numbers(void);

int grav_tree_compare_key(const void *a, const void *b);
int dens_compare_key(const void *a, const void *b);
int hydro_compare_key(const void *a, const void *b);

int data_index_compare(const void *a, const void *b);
int peano_compare_key(const void *a, const void *b);

void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_domain(void *b, size_t n, size_t s);
void mysort_idlist(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_pmnonperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_peano(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));


int density_evaluate(int target, int mode, int *nexport, int *nsend_local);
int density_isactive(int n);

void GetMachNumberCR( struct sph_particle_data *Particle );
#ifndef CS_MODEL
void GetMachNumber( struct sph_particle_data* Particle );
#else
void GetMachNumber( struct sph_particle_data* Particle, struct particle_data* Particle_Hsml );
#endif
void GetShock_DtEnergy( struct sph_particle_data* Particle );

#ifdef MAGNETIC
#ifdef BFROMROTA
void rot_a(void);
void rot_a_evaluate(int i, int mode);
#endif
#endif
size_t sizemax(size_t a, size_t b);


void reconstruct_timebins(void);

void init_peano_map(void);
peanokey peano_hilbert_key(int x, int y, int z, int bits);
peanokey peano_and_morton_key(int x, int y, int z, int bits, peanokey *morton);
peanokey morton_key(int x, int y, int z, int bits);

void catch_abort(int sig);
void catch_fatal(int sig);
void terminate_processes(void);
void enable_core_dumps_and_fpu_exceptions(void);
void write_pid_file(void);

void pm_init_periodic_allocate(void);

void pm_init_periodic_free(void);

void move_particles(int time1);


void find_next_sync_point_and_drift(void);
void find_dt_displacement_constraint(double hfac);

void set_units_sfr(void);

void gravity_forcetest(void);

void allocate_commbuffers(void);
void allocate_memory(void);
void begrun(void);
void check_omega(void);
void close_outputfiles(void);
void compute_accelerations(int mode);
void compute_global_quantities_of_system(void);
void compute_potential(void);
void construct_timetree(void);
void cooling_and_starformation(void);
#ifdef GASRETURN
void gas_return(void);
#endif
void count_hot_phase(void);
void delete_node(int i);
void density(void);
void density_decouple(void);
void determine_interior(void);
int dissolvegas(void);
void do_box_wrapping(void);
void domain_Decomposition(void);
double enclosed_mass(double R);
void endrun(int);
void energy_statistics(void);
void ensure_neighbours(void);

void every_timestep_stuff(void);
void ewald_corr(double dx, double dy, double dz, double *fper);

void ewald_force(int ii, int jj, int kk, double x[3], double force[3]);
void ewald_force_ni(int iii, int jjj, int kkk, double x[3], double force[3]);

void ewald_init(void);
double ewald_psi(double x[3]);
double ewald_pot_corr(double dx, double dy, double dz);
int find_ancestor(int i);
int find_next_outputtime(int time);
void find_next_time(void);
int find_next_time_walk(int node);
void free_memory(void);
void advance_and_find_timesteps(void);
int get_timestep(int p, double *a, int flag);

void determine_PMinterior(void);

double get_starformation_rate(int i);
void gravity_tree(void);
void hydro_force(void);
void init(void);
#ifndef LT_STELLAREVOLUTION
void init_clouds(void);
void integrate_sfr(void);
#else
void init_clouds(int, double, double, double, double*, double*);
void integrate_sfr(double, double, double, double, double);
#endif
void insert_node(int i);
int mark_targets(void);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void open_outputfiles(void);
void write_outputfiles_header(void);
void peano_hilbert_order(void);
double pot_integrand(double xx);
void predict(double time);
void predict_collisionless_only(double time);
void predict_sph_particles(double time);
void prepare_decouple(void);
void read_ic(char *fname);
void read_ic_cluster(char *fname);
void read_ic_cluster_gas(char *fname);
void read_ic_cluster_wimp(char *fname);
int read_outputlist(char *fname);
void read_parameter_file(char *fname);
void rearrange_particle_sequence(void);
void reorder_gas(void);
void reorder_particles(void);
void restart(int mod);
void run(void);
void savepositions(int num);
void savepositions_ioformat1(int num);
double second(void);
void set_softenings(void);
void set_sph_kernel(void);
void set_units(void);
void setup_smoothinglengths(void);

void sumup_large_ints(int n, int *src, long long *res);
void sumup_longs(int n, long long *src, long long *res);

void statistics(void);
double timediff(double t0, double t1);
void veldisp(void);
void veldisp_ensure_neighbours(int mode);

void gravity_tree_shortrange(void);


double get_hydrokick_factor(int time0, int time1);
double get_gravkick_factor(int time0, int time1);
double drift_integ(double a, void *param);
double gravkick_integ(double a, void *param);
double growthfactor_integ(double a, void *param);
double hydrokick_integ(double a, void *param);
void init_drift_table(void);
double get_drift_factor(int time0, int time1);
#if defined(DISTORTIONTENSOR) || defined(DISTORTIONTENSORPS) 
double get_distortionkick_factor(int time0, int time1);
#endif
double measure_time(void);


/* on some DEC Alphas, the correct prototype for pow() is missing,
   even when math.h is included ! */

double pow(double, double);


void long_range_init(void);
void long_range_force(void);
void pm_init_periodic(void);
void pmforce_periodic(void);
void pm_init_regionsize(void);
void pm_init_nonperiodic(void);
int pmforce_nonperiodic(int grnr);

int pmpotential_nonperiodic(int grnr);
void pmpotential_periodic(void);

void readjust_timebase(double TimeMax_old, double TimeMax_new);

double enclosed_mass(double R);
void pm_setup_nonperiodic_kernel(void);

#ifdef CHEMISTRY
double dmax(double, double);
double dmin(double, double);
#endif

#ifdef LT_STELLAREVOLUTION
void fsolver_error_handler(const char *, const char *, int, int);
int get_Yset(int);

void init_SN(void);
int evolve_SN(int);
double INLINE_FUNC get_cost_SE(int);
int INLINE_FUNC get_chemstep(int, int, double*, double);
int INLINE_FUNC get_current_chem_bin(double, int);
int compare_steps(const void*, const void*);

void read_metals(void);

void init_clouds_cm(double*, double*, double, double, int, double*);

int write_eff_model(int);
int read_eff_model(int);

void read_SnIa_yields(void);
void read_SnII_yields(void);
void read_AGB_yields(void);

void read_metalcool_table(void);

void initialize_star_lifetimes(void);
double INLINE_FUNC get_age(double);

void recalculate_stellarevolution_vars(void *, int);
void recalc_eff_model();

void calculate_effective_yields(double, double, void*);
double calculate_FactorSN(double, double, void*);

double get_imf_params(double);
double INLINE_FUNC normalize_imf(double, double, void*);

double *get_meanZ(void);
void write_metallicity_stat(void);
void get_metals_sumcheck(int mode);
double INLINE_FUNC get_metallicity(int, int);
double INLINE_FUNC get_metallicity_solarunits(MyFloat);

double INLINE_FUNC get_entropy(int);

int INLINE_FUNC getindex(double*, int, int, double*, int*);

int INLINE_FUNC perturb(double*, double*);

#endif

#ifdef RADTRANSFER
/* Eddington tensor computation */
int eddington_treeevaluate(int target, int mode, int *nexport, int *nsend_local);
void eddington(void);

/* radiative transfer integration */
void radtransfer(void);
#ifndef CG
int radtransfer_evaluate(int target, int mode, int *nexport, int *nsend_local);
void radtransfer_begin(void);
#else
double radtransfer_vector_multiply(double *a, double *b);
double radtransfer_vector_sum(double *a);
void radtransfer_matrix_multiply(double *in, double *out, double *sum);
int radtransfer_evaluate(int target, int mode, double *in, double *out, double *sum, int *nexport, int *nsend_local);
#endif
void radtransfer_mean(void);
void radtransfer_set_simple_inits(void);
void radtransfer_simple_output(void);
void radtransfer_update_chemistry(void);
void radtransfer_get_active(void);

/* emission approximation */
int ngb_treefind_stars(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, 
		       int *nexport, int *nsend_local);
int star_density_evaluate(int target, int mode, int *nexport, int *nsend_local);
void star_density(void);
#endif


#ifdef AUTO_SWAP_ENDIAN_READIC
void swap_Nbyte(char *data, int n, int m);
void swap_header(void);
#endif

#ifdef DISTORTIONTENSORPS
void get_half_kick_distortion(int pindex, MyDouble half_kick_add[6][6]);
void analyse_phase_space(MyIDType pindex, 
                         MyDouble *s_1, MyDouble *s_2, MyDouble *s_3, 
                         MyDouble *smear_x, MyDouble *smear_y, MyDouble *smear_z,
                         MyDouble *D_xP_proj, MyDouble *second_deriv, MyDouble *sigma);
MyDouble get_analytic_annihilation(MyDouble s_1, MyDouble s_2, MyDouble s_3, 
                                   MyDouble s_1_prime, MyDouble s_2_prime, MyDouble s_3_prime,
                                   MyDouble second_deriv, MyDouble second_deriv_prime, 
                                   MyDouble dt, MyDouble sigma);
MyDouble get_max_caustic_density(MyDouble s_1, MyDouble s_2,
                                 MyDouble s_1_prime, MyDouble s_2_prime, 
                                 MyDouble second_deriv, MyDouble second_deriv_prime, 
                                 MyDouble sigma);
void get_current_ps_info(MyIDType pindex, MyDouble * flde, MyDouble * psde);
/* some math functions we need from phasespace_math.c */
void ludcmp(MyDouble ** a, int n, int *indx, MyDouble * d);
MyDouble **matrix(long nrl, long nrh, long ncl, long nch);
MyDouble *vector(long nl, long nh);
void free_matrix(MyDouble ** m, long nrl, long nrh, long ncl, long nch);
void free_vector(MyDouble * v, long nl, long nh);
void mult_matrix(MyDouble ** matrix_a, MyDouble ** matrix_b, int dimension, MyDouble ** matrix_result);
void mult_matrix_transpose_A(MyDouble ** matrix_a, MyDouble ** matrix_b, int dimension, MyDouble ** matrix_result);
void luinvert(MyDouble ** input_matrix, int n, MyDouble ** inverse_matrix);
void eigsrt(MyDouble d[], MyDouble ** v, int n);
void jacobi(MyDouble ** a, int n, MyDouble d[], MyDouble ** v, int *nrot);
#ifdef PMGRID
void pmtidaltensor_periodic(int component);
int pmtidaltensor_nonperiodic(int grnr, int component);
void check_tidaltensor(int particle_ID);
#endif
#endif



#ifdef EOS_DEGENERATE
extern void eos_init();
extern void eos_deinit();
extern int eos_calc_egiven( double rho, double *xnuc, double e, double *temp, double *p, double *dpdr );
extern int eos_calc_tgiven( double rho, double *xnuc, double temp, double *e, double *dedt );
extern int eos_calc_egiven_v( double rho, double *xnuc, double *dxnuc, double fac, double e, double *temp, double *p, double *dpdr );
extern int eos_calc_tgiven_v( double rho, double *xnuc, double *dxnuc, double fac, double temp, double *e, double *dedt );
#endif

#ifdef NUCLEAR_NETWORK
extern void network_init();
extern void network_integrate( double temp, double rho, double *x, double *dx, double dt, double *dedt );
#endif

#ifdef HEALPIX

void mk_xy2pix(int *x2pix, int *y2pix);
void vec2pix_nest( const long nside, double *vec, long *ipix); 
void healpix_halo(float *res);

#endif

