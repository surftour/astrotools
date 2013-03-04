
#define PI 3.1415926
#define G  43007.0                   /* gravitational constant */


extern char gal_output[100];

extern char gal_fname1[100];
extern double theta1;
extern double phi1;

extern char gal_fname2[100];
extern double theta2;
extern double phi2;

extern char combinetype[100];

extern double rperi;
extern double rstart;
extern double ecc;

extern double g1_x_start;
extern double g1_y_start;
extern double g1_z_start;
extern double g1_vx_start;
extern double g1_vy_start;
extern double g1_vz_start;

extern double g2_x_start;
extern double g2_y_start;
extern double g2_z_start;
extern double g2_vx_start;
extern double g2_vy_start;
extern double g2_vz_start;



extern struct galaxy_data {
  double Mdark, Scalefac;  
  int Ntot, Ngas, Ndm;
  double Mtot, Mdm;
  float **pos,**vel,*m,*u,*z;
  int   *id;
} gal[2];

/* Header for the standard file format.
 */
extern struct io_header
{
  int npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  int flag_stellarage;
  int flag_metals;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8 - 2*4];  /* fills to 256 Bytes */
} header[2];

void read_parameterfile(char *fname);

void manual_move_galaxies(struct galaxy_data *g1, struct galaxy_data *g2);

void kepler_move_galaxies(struct galaxy_data *g1, struct galaxy_data *g2);

void translate(struct galaxy_data *g, double vec[3]);

void vel_translate(struct galaxy_data *g, double vec[3]);

void load_particles(char *fname, struct galaxy_data *g, 
		    struct io_header *header);

void turn_galaxy(struct galaxy_data *g,double omega,double incl);

void save_combined(char *fname,struct galaxy_data *g1,struct galaxy_data *g2);




