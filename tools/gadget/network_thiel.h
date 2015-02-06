#ifndef NETWORK_THIEL_H
#define NETWORK_THIEL_H

#ifdef NUCLEAR_NETWORK

#define NETWORK_DIFFVAR 1e-6	/* variation for numerical derivatives */

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

void network_init();
void network_integrate( double temp, double rho, double *x, double *dx, double dt, double *dedt );
void network_getrhs( double temp, double rho, double *y, double *rhs );
void network_getjacob( double temp, double rho, double *y, double *rhs, double *jacob );

void network_part( double temp );
void network_getrates( double temp, double rho );
void network_checkswap( char* fname, int *swap );
double network_SwapDouble( double Val );
int network_SwapInt( int Val );

#pragma pack(push,1)

struct network_ratedata_1b {
	int n;
	int ni[3];
	int ic, ir, iv;
	double p[7];
	double q;
} *network_rates_1b;

struct network_ratedata_2b {
	int n;
	int ni[4];
	int ic, ir, iv;
	double p[7];
	double dummy;
} *network_rates_2b;

struct network_ratedata_3b {
	int n;
	int ni[3];
	int ic, ir, iv;
	double p[7];
	double dummy;
} *network_rates_3b;

#pragma pack(pop)

struct network_nucdata {
	int n;
	int *la;
	int *le;
	int nn, np;	/* neutron, proton number */
	double na, bi, exm, angm;
	char name[6];
	double dlg[24];
} *network_nucdata;

struct network_data {
	int nuc_count;
	double t9[24];
	double *a[3];
	int *mu[3];
	int mumax[3];
	double *gg;
	double *sig1b, *sig2b, *sig3b;
	double *x;
	double *y;
	double conv; /* unit conversion factor */
	double *na;
} network_data;

#endif

#endif /* NETWORK_THIEL_H */
