

extern struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedbacktp;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
#ifdef TJ_VERSION
  int      flag_multiphase;
#endif
  int      flag_stellarage;
#ifdef TJ_VERSION
  int      flag_sfrhistogram;
  int      flag_stargens;
  int      flag_snaphaspot;
#endif
  int      flag_metals;
#ifdef TJ_VERSION
  int      flag_energydetail;
  int      flag_parentid;
  int      flag_starorig;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8 - 9*4];  /* fills to 256 Bytes */
#else
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8 - 2*4];  /* fills to 256 Bytes */
#endif
} header1;



extern int     NumPart, Ngas;

extern struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
  float  Nh, Mfs, Mclouds, hsml, sfr, meanz, TpU;
  float  Metallicity;
  float  totrad, totshock, totfb;
  int    ParentID;
  float  OrigMass, OrigHsml;
  float  Potential;
} *P;

extern int *Id;

extern double  Time, Redshift;
extern int     SnapInfo;
extern int     SnapFile;



int load_snapshot(char *fname, int files, int quiteon, int headeronly);
int allocate_memory(void);
int reordering(void);



