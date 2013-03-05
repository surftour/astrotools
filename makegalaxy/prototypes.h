
#include <stdio.h>

#include "forcetree.h"

int find_bin(double value, double* array, int size);
double interpolate_table(double f, int bin, double* array);
double interpolate_table_2d(double f1, double f2, int bin1, int bin2,
			    double** array);

void init_clouds(void);
double halo_mass(double r);

void integrate_gasdensity(int use_Dphi_z);
void integrate_and_adjust(int);

void read_parameterfile(char *fname);

void determine_cumulative_gasdensity(void);
void init_central_densities(void);
double surface_density_disk(double r);
double integrate_vertically(double rho0, double* Dphi_z, 
			    double* rho_arr, double* sigma_arr);

void dump_gas_density(char*);
void compute_vertical_force_field(void);
void compute_radial_force_field(void);
void compute_phi_field(void);

void compute_DF_lookuptable(void);
double comp_DF_halo(double E, double Lz);
double comp_DF_bulge(double E, double Lz);
double compute_ani_beta(double r);

void integrate_surfacedensity(void);

void dump_eqn_of_state(void);

double eqn_of_state(double rho);

void  set_dummy_particles(void);

void compute_vstream_gas(void);

void dump_veldisp_field(void);

double set_halo_positions(void);
double set_disk_positions(void);
double set_bulge_positions(void);
void set_gas_positions(void);

double set_halo_velocities(void);
double set_disk_velocities(void);
double set_bulge_velocities(void);
void set_gas_velocities(void);

void compute_velocity_dispersions_halo(void);
void compute_velocity_dispersions_disk(void);
void compute_velocity_dispersions_bulge(void);
void compute_local_escape_speed(void);
void allocate_memory(int NumPart);

void write_toomre(void);


void save_particles(void);

void init_units(void);
void init(void);

void structure(void);
void write_paramfile(void);

double second(void);

void prepare_cumlative_profile(void);
double mass_cumulative_disk(double);
double mass_cumulative_gas(double);
double mass_cumulative_bulge(double);
double disk_angmomentum(void);
double additional_mass_in_halo_cutoff(void);
double fc(double);
double gc(double);
void solve_mass_shells(void);
void setup_massprofile(void);


double halo_q_to_r(double q);
double disk_q_to_R(double q);
double disk_q_to_z(double q);
double bulge_q_to_r(double q);
double gas_q_to_R(double q);
double sample_gas_position(double* q, double *pos);


double comp_Dphi_z_halo(double R,double z);
double comp_Dphi_R_halo(double R,double z);
double comp_rho_halo(double R,double z);
double comp_phi_halo(double R, double z);


double comp_Dphi_z_bulge(double R,double z);
double comp_Dphi_R_bulge(double R,double z);
double comp_rho_bulge(double R,double z);
double comp_phi_bulge(double R,double z);


double comp_Dphi_R_disk_razorthin(double RR,double zz);
double comp_Dphi_z_disk(double R,double z);
double comp_Dphi_R_disk(double R,double z);
double comp_rho_disk(double R,double z);
double comp_Dphi_Z_disk_tree(double RR, double zz);
double comp_Dphi_R_disk_tree(double RR, double zz);
double comp_phi_disk_tree(double R, double z);

double comp_Dphi_z_gas(double R,double z);
double comp_Dphi_R_gas(double R,double z);
double comp_rho_gas(double R,double z);
double comp_Dphi_R_gas_razorthin(double RR,double zz);


double comp_Dphi_z(double R,double z);
double comp_Dphi_R(double R,double z);

double comp_phi(double R, double z);

double epicyclic_kappa2(double R);


double splint_zs_y_D2y(double z);
double splint_xl_yl_D2yl(double t);

int find_idx(double x, double *fx, int size);

double   drand48(void);

void write_cumulative_mass(void);



/* returns the maximum of two double
 */
double dmax(double x, double y);
double dmin(double x, double y);
int imax(int x, int y);
int imin(int x, int y);
