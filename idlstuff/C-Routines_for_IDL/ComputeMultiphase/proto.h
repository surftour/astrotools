
#ifndef ALLVARS_H
#include "allvars.h"
#endif
//#include "forcetree.h"
#include "cooling.h"


#ifdef HAVE_HDF5
#include <hdf5.h>
void write_header_attributes_in_hdf5(hid_t handle);
void read_header_attributes_in_hdf5(char *fname);
#endif

void hpm_pressureforce(void);
void hpm_find_eqs(void);
void hpm_find_eqs_simple(void);
void hpm_find_eqs_selfconsistent(void);

void kinetic_feedback_mhm(void);
int kin_compare_key(const void *a, const void *b);
void kinetic_evaluate(int target, int mode);

void bubble(void);
void find_CM_of_biggest_group(void);
int compare_length_values(const void *a, const void *b);
void fof_get_group_center(float *cm, int first, int len);


void conduction_smoothed_evaluate(int target, int mode);
void conduction_smoothed_temperature(void);

int fof_find_dmparticles_evaluate(int target, int mode);


#ifdef DARKENERGY

double DarkEnergy_a(double);
double DarkEnergy_t(double);

#ifdef TIMEDEPDE
void fwa_init(void);
double INLINE_FUNC fwa(double);
#endif

#endif

void blackhole_accretion(void);
void blackhole_evaluate(int target, int mode);
int  blackhole_compare_key(const void *a, const void *b);

int ngb_treefind_blackhole(FLOAT searchcenter[3], FLOAT hguess, int *startnode);


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
void fof_find_nearest_dmparticle_evaluate(int target, int mode);
int fof_compare_key(const void *a, const void *b);
void fof_link_special(void);
void fof_link_specialpair(int p, int s);
void fof_make_black_holes(void);

void ngb_treefind_flagexport(FLOAT searchcenter[3], FLOAT hguess);

void write_file(char *fname, int readTask, int lastTask);

void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last);

//int get_values_per_blockelement(enum iofields blocknr);

//int get_datatype_in_block(enum iofields blocknr);
//void get_dataset_name(enum iofields blocknr, char *buf);


//int blockpresent(enum iofields blocknr);
//void fill_write_buffer(enum iofields blocknr, int *pindex, int pc, int type);
//void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);

//int get_particles_in_block(enum iofields blocknr, int *typelist);

//int get_bytes_per_blockelement(enum iofields blocknr);

void read_file(char *fname, int readTask, int lastTask);
void fill_Tab_IO_Labels(void);

int ngb_treefind_darkmatter(FLOAT searchcenter[3], FLOAT hguess, int *startnode);

void long_range_init_regionsize(void);

int find_files(char *fname);

int metals_compare_key(const void *a, const void *b);
void enrichment_evaluate(int target, int mode);

void hydro_evaluate(int target, int mode);

void pm_init_nonperiodic_allocate(int dimprod);

void  pm_init_nonperiodic_free(void);

double get_random_number(int id);
void set_random_numbers(void);

int grav_tree_compare_key(const void *a, const void *b);
int dens_compare_key(const void *a, const void *b);
int hydro_compare_key(const void *a, const void *b);

void density_evaluate(int i, int mode);
#if defined(MAGNETIC) && defined(BSMOOTH)
void bsmooth_evaluate(int i, int mode);
#endif

void init_peano_map(void);
peanokey peano_hilbert_key(int x, int y, int z, int bits);

void catch_abort(int sig);
void catch_fatal(int sig);
void terminate_processes(void);
void enable_core_dumps_and_fpu_exceptions(void);
void write_pid_file(void);

void pm_init_periodic_allocate(int dimprod);

void pm_init_periodic_free(void);

void move_particles(int time0, int time1);

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
void count_hot_phase(void);
void delete_node(int i);
void density(void);
#if defined(MAGNETIC) && defined(BSMOOTH)
void bsmooth(void);
#endif
void density_decouple(void);
void determine_interior(void);
int dissolvegas(void);
void do_box_wrapping(void);
void DomainDecomposition(void);
double drand48();
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
void init_clouds(void);
void insert_node(int i);
void integrate_sfr(void);
int mark_targets(void);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void open_outputfiles(void);
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

void statistics(void);
double timediff(double t0, double t1);
void veldisp(void);
void veldisp_ensure_neighbours(int mode);

void gravity_tree_shortrange(void);


double get_hydrokick_factor(int time0, int time1);
double get_gravkick_factor(int time0, int time1);
double drift_integ(double a, void *param);
double gravkick_integ(double a, void *param);
double hydrokick_integ(double a, void *param);
void init_drift_table(void);
double get_drift_factor(int time0, int time1);

int ngb_clear_buf(FLOAT searchcenter[3], FLOAT hguess, int numngb);



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

void pmpotential_nonperiodic(int grnr);
void pmpotential_periodic(void);

void readjust_timebase(double TimeMax_old, double TimeMax_new);

double enclosed_mass(double R);
void pm_setup_nonperiodic_kernel(void);
