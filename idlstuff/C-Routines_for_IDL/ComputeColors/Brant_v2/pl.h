extern int nm;
extern int nz;
extern int na;
extern double *zz;
extern double *lage;
extern double *m;
void AllocateMagnitudes(void);
void FreeMagnitudes(void);
void GetMags(double lpa, double pz, double mags[]);
void LoadMagnitudeData(void);
int manual_gsl_interp_bsearch(double x_array[], double x, int index_lo, int index_hi);