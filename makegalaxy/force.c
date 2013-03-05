#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_odeiv.h"

#include "prototypes.h"
#include "globvars.h"



void dump_gas_density(char* gsdfile)
{
  FILE *fd;
  int si, rbin;
  fd = fopen(gsdfile, "w");
  /*fd = fopen("gas_density.txt", "w");*/

  si = RSIZE + 1;
  fwrite(&si, sizeof(int), 1, fd);
  si = ZSIZE + 1;
  fwrite(&si, sizeof(int), 1, fd);

  fwrite(list_R, sizeof(double), RSIZE + 1, fd);
  fwrite(list_z, sizeof(double), ZSIZE + 1, fd);

  for(rbin=0; rbin<=RSIZE; ++rbin)
    fwrite(&RhoGas[rbin][0], sizeof(double), (ZSIZE + 1), fd);
  for(rbin=0; rbin<=RSIZE; ++rbin)
    fwrite(&CumulMassGas[rbin][0], sizeof(double), (ZSIZE + 1), fd);

  fclose(fd);
}


double surface_density_gasdisk(double r)
{
  /*
    if(r / H > NumGasScaleLengths)
    return 0;
    else
    // M_DISK is now just the stellar disk 
    return (M_GAS / (2 * M_PI * H * H)) * exp(-r / H)
    / (1 - (1 + NumGasScaleLengths) * exp(-NumGasScaleLengths));
  */

  double Rd=H;
  switch(GasDistribution) {
  case 1:
    Rd= H*GasExpAlpha;
    // fall through
  case 0:
    return M_GAS/(2*PI*Rd*Rd)*exp(-r/Rd);
    break;
  case 2:
    return r<PowerLawCutoff ? 
      (2-PowerLawGamma)*M_GAS/pow(PowerLawCutoff, 2-PowerLawGamma) : 0;
    break;
  }
}


/* Initializes the gas densities in the midplane in RhoGas[rbin][0]
   from the (known) surface density specified by
   surface_density_gasdisk, using the rhocentral[] table calculated in
   effmodel.c */
void init_central_densities(void)
{
  int tabbin, rbin;
  double u, sigma;

  if(N_GAS == 0)
    return;

  tabbin = 0;

  for(rbin = RSIZE; rbin >= 0; rbin--)
    {
      RhoGas[rbin][0] = 0;

      sigma = surface_density_gasdisk(list_R[rbin]);

      // find sigma bin
      tabbin = find_idx(sigma, sigmatab, TABSIZE-1);
      if(tabbin>=0) {

	u = (sigma - sigmatab[tabbin]) / (sigmatab[tabbin+1] - sigmatab[tabbin]);

	RhoGas[rbin][0] = interpolate_table(u, tabbin, rhocentral);
      } else {
	RhoGas[rbin][0] = rhocentral[0] * (sigma / sigmatab[0]);
      }

      printf("rbin= %d   r=%g  RhoGas[rbin][0] = %g  sigma= %g  sigmatab[%d]= %g   u= %g\n",rbin, list_R[rbin], RhoGas[rbin][0], sigma, tabbin, sigmatab[tabbin],u);
    }

  fflush(stdout);
}


/* Calculates the cumulative gas surface density in CumulMassGas[][]
   by summing the RhoGas[][] array. XXX This should NOT be used,
   instead we determine the surface density directly in
   integrate_vertically, because that's a better estimate than we get
   just integrating through the points afterward. */
void determine_cumulative_gasdensity(void)
{
  assert(0);
  int rbin, zbin;

  if(N_GAS == 0)
    return;

  for(rbin = 0; rbin <= RSIZE; rbin++)
    {
      // XXX really only need to reset zbin=0
      for(zbin = 0; zbin <= ZSIZE; zbin++)
	CumulMassGas[rbin][zbin] = 0;

      for(zbin = 1; zbin <= ZSIZE; zbin++) {
	assert(RhoGas[rbin][zbin]>=0);
	CumulMassGas[rbin][zbin] = CumulMassGas[rbin][zbin - 1] +
	  0.5 * (RhoGas[rbin][zbin] + RhoGas[rbin][zbin - 1]) * (list_z[zbin] - list_z[zbin - 1]);
      }

      for(zbin = 0; zbin <= ZSIZE; zbin++)
	CumulMassGas[rbin][zbin] *= 2;	/* to account for other side */

      /*
	if(surface_density_gasdisk(list_R[rbin]) > 0)
	printf("target=%g  found=%g\n", surface_density_gasdisk(list_R[rbin]), CumulMassGas[rbin][ZSIZE]);
      */
    }
}


// Figures out the vertical structure of the gas? 
void integrate_and_adjust(int use_Dphi_z)
{
  int rep, rbin, zbin;
  double minfactor=1, maxfactor=1;

  if(N_GAS == 0)
    return;


  printf("integrate and adjust surface density...\n");

  int n=0;
  while(1)
    {

      integrate_gasdensity(use_Dphi_z);

      ++n;

      minfactor= maxfactor =1;
      int rbinmax=0, rbinmin=0;
      for(rbin = 0; rbin <= RSIZE; rbin++)
	{
	  if(RhoGas[rbin][0] != 0)
	    {
	      double corrfactor=
		surface_density_gasdisk(list_R[rbin])/CumulMassGas[rbin][ZSIZE];

	      if(corrfactor<minfactor) {
		minfactor=corrfactor;
		rbinmin=rbin;
	      }
	      if(corrfactor>maxfactor) {
		maxfactor=corrfactor;		 
		rbinmax=rbin;
	      }
	      minfactor = dmin(minfactor, corrfactor);
	      maxfactor = dmax(maxfactor, corrfactor);

	      // to avoid overshooting the correction, we only apply a
	      // small part of it
	      for(zbin=0; zbin<ZSIZE; ++zbin)
		RhoGas[rbin][zbin] = RhoGas[rbin][zbin] * (fabs(corrfactor-1)<.001 ? corrfactor : pow(corrfactor, .25));

	    }
	}
      printf("Maximum surface density correction factors: %f,%f at bins %d, %d\n",
	     minfactor, maxfactor, rbinmin, rbinmax);

      if((minfactor>.98)&&(maxfactor<1.02)&&(n>=5))
	break;
      if(n>100) {
	printf("Convergence problem, stopping\n");
	break;
      }
    }
}




/* This calculates the vertical gas structure in RhoGas[][] by
   integrating up from the midplane, using RhoGas[rbin][0] as the
   initial value.

   For the first iteration, we use the cumulative surface density to
   estimate the gravity (use_Dphi_z=0), but later on when called from
   integrate_and_adjust, we use the actual computed force field.
*/
void integrate_gasdensity(int use_Dphi_z)
{
  double sigmacheck; //, rhoold;
  int rbin; //, zbin;

  if(N_GAS == 0)
    return;


  for(rbin = 0; rbin <= RSIZE; rbin++)
    {
      sigmacheck=integrate_vertically(RhoGas[rbin][0], 
				      use_Dphi_z ? Dphi_z[rbin] : NULL, 
				      RhoGas[rbin], CumulMassGas[rbin]);
      if(!use_Dphi_z && sigmacheck>0) {
	//printf("rbin=%d sigmacheck=%g\n", rbin, sigmacheck);
	//fflush(stdout);
      }
    }
}

typedef struct {
  int zbin;
  double u4;
  double* Dphi_z;
} drho_dz_pars;


/* integrand for vertical integration of the gas density. 
   x={rho, sigma}, f={drho_dz, dsigma_dz} */
int drho_dz (double z, const double x[], double f[],
	     void *params)
{
  drho_dz_pars* p = (drho_dz_pars*)params;
  double u, dphi_dz;

  if(p->Dphi_z) {
    u = (z-list_z[p->zbin])/(list_z[p->zbin+1]-list_z[p->zbin]);
    dphi_dz = interpolate_table(u, p->zbin, p->Dphi_z);
  }
  else {
    dphi_dz = 4 * PI * G * x[1];
  }

  double P, P1, P2, gam;

  if(x[0] > PhysDensThresh)
    {
      P = P1 = eqn_of_state(x[0]);
      
      P2 = eqn_of_state(1.1*x[0]);
      
      gam = log(P2 / P1) / log(1.1);
    }
  else
    {
      P = GAMMA_MINUS1 * x[0] * p->u4;
      gam = 1.0;
    }

  f[0] = -dphi_dz * x[0] * x[0] / (P*gam);
  f[1] = x[0];

  return GSL_SUCCESS;
}


  

/* Integrates the gas density upward from the midplane using the
   list_z grid, starting with the initial value rho0 for the midplane
   density. The total surface density is returned. If rho_arr is
   supplied, the gas density field is written into that.  If sigma_arr
   is supplied, the cumulative surface density is written into
   that. If Dphi_z is supplied, the force field in that array is used,
   otherwise the cumulative surface density is used to estimate the
   gravity. */
double integrate_vertically(double rho0, double* Dphi_z, double* rho_arr, double* sigma_arr)
{
  double rho, rho2, dz, gam, P1;
  double sigma, rhoold;
  double P, P2, drho, Dphi_z_here;
  double meanweight, u4, z;
  int zbin;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= UnitMass_in_g / UnitEnergy_in_cgs;

  /* GB */
#ifdef U4MOD
  //	u4 = u4/5.0;
  //	u4 = u4/2.0;
  u4 = u4/10.0;
  //	u4 = u4/50.0;
#endif

  zbin = 0;
  z = 0;
  rho = rho0;
  drho = 0;    
  sigma = 0;
 
  // gsl integrator and variables
  // we could probably optimize this and not realloc it every time
  gsl_odeiv_step * ode = gsl_odeiv_step_alloc (gsl_odeiv_step_rkf45, 2);
  gsl_odeiv_step_reset(ode);
  double ode_vars[2];
  ode_vars[0]=rho0;
  ode_vars[1]=0;
  double ode_err[2];
  double ode_dvar_dz[2];
  drho_dz_pars ode_pars;
  ode_pars.u4 = u4;
  ode_pars.Dphi_z = Dphi_z;
  gsl_odeiv_system ode_sys = {drho_dz, 0, 2, &ode_pars};
  gsl_odeiv_control * odestepper = gsl_odeiv_control_y_new (1e-7*rho0, 1e-7);
  gsl_odeiv_evolve * odeevolver = gsl_odeiv_evolve_alloc (2);
  gsl_odeiv_evolve_reset(odeevolver);
  double cur_z, cur_h;

  if(rho_arr) {
    for(zbin = 0; zbin <= ZSIZE; zbin++)
      rho_arr[zbin] = 0;
    rho_arr[0]=rho0;
  }
  if(sigma_arr)
    sigma_arr[0] = 0;
  
  if(rho0 > 0) 
    {
      zbin=0;
      while(zbin < ZSIZE ) 
	{

	  /* The function being integrated is:

	     drho/dz = - dphi/dz * rho^2 /(P*gamma)

	     Use a higer-order gsl integrator instead of this Euler
	     forward? 
	  */


	  dz = list_z[zbin + 1] - list_z[zbin];

	  ode_pars.zbin=zbin;
	  cur_z= list_z[zbin];
	  cur_h=dz;
	  gsl_odeiv_evolve_apply (odeevolver, odestepper, ode, &ode_sys, &cur_z, list_z[zbin+1], &cur_h, ode_vars);
	  rho=ode_vars[0];
	  sigma=ode_vars[1];
	  if(rho<0)
	    break;

	  zbin++;

	  if(rho_arr)
	    rho_arr[zbin] = rho;
	  if(sigma_arr)
	    sigma_arr[zbin] = sigma*2;
	}

      // fill out the sigma_arr
      if(sigma_arr)
	for(;zbin<ZSIZE;++zbin)
	  sigma_arr[zbin+1]=sigma_arr[zbin];

      sigma *= 2;
    }

  gsl_odeiv_step_free (ode);
  gsl_odeiv_control_free (odestepper);
  gsl_odeiv_evolve_free (odeevolver);

  return sigma;
}


void set_dummy_particles(void)
{
  int i, j, k, n, rbin, zbin;
  double q, R, z, phi, u;
  double qq[3], pos[3];
  /*
    FILE *fd;
  */


  printf("set dummy particles...\n");

  /* first, let's set the gas particles */

  n = 0;

  if(N_GAS > 0)
    {
      for(i = 0; i < RMASSBINS; i++) {
	
	qq[0] = (i + 0.5) / RMASSBINS;
	
	for(j = 0; j < ZMASSBINS; j++) {
	  // q[1] goes from -1 to 1
	  qq[1] = 2*(j+0.5)/ZMASSBINS-1;
	  
	  for(k = 0; k < PHIMASSBINS; k++) {
	    qq[2] = (k+0.5)/PHIMASSBINS;
	    sample_gas_position(qq, pos);
	    P[n].Pos[0]=pos[0];
	    P[n].Pos[1]=pos[1];
	    P[n].Pos[2]=pos[2];
	    P[n].Mass = M_GAS / (RMASSBINS * ZMASSBINS * PHIMASSBINS);
	    n++;
	  }
	}
      }  
    }

  /* now let's set the star particles (note: n continues to count higher) */

  if(N_DISK > 0)
    {
      for(i = 0; i < RMASSBINS; i++)
	{
	  q = (i + 0.5) / RMASSBINS;
	  R = H * disk_q_to_R(q);

	  for(j = 0; j < ZMASSBINS; j++)
	    {
	      q = (j + 0.5) / ZMASSBINS;
	      z = disk_q_to_z(q);

	      for(k = 0; k < PHIMASSBINS; k++)
		{
		  phi = (k + 0.5) / PHIMASSBINS * 2 * M_PI;

		  P[n].Pos[0] = R * cos(phi);
		  P[n].Pos[1] = R * sin(phi);
		  P[n].Pos[2] = z;

		  /*  M_DISK is now the stellar disk only 
		      P[n].Mass = (1 - GasFraction) * M_DISK / (RMASSBINS * ZMASSBINS * PHIMASSBINS); */
		  P[n].Mass = M_DISK / (RMASSBINS * ZMASSBINS * PHIMASSBINS);
		  n++;
		}
	    }
	}

    }


  force_treebuild();

  /* now, we read this from the parameter file
     ErrTolTheta = 0.15;
  */


  /*
    fd = fopen("particles.dat", "w");
    fwrite(&NumPart, sizeof(int), 1, fd);
    for(n = 0; n < NumPart; n++)
    fwrite(&P[n].Pos[0], sizeof(float), 3, fd);
    fclose(fd);
  */
}



void compute_phi_field(void)
{
  FILE *fd;
  char phifile[50]="";
  int i, j;

  printf("Start computing phi field  (GravSoftening= %g, ErrTolTheta= %g).\n",GravSoftening,ErrTolTheta); fflush(stdout);

  PhiMin= +1.0e+20;

  for(i = 0; i <= RSIZE; i++)
    {
      //printf("phi %d(%d)\r", i, RSIZE); fflush(stdout);
      printf("phi %d(%d)\n", i, RSIZE); fflush(stdout);

      for(j = 0; j <= ZSIZE; j++)
        {
          phi[i][j] = comp_phi(list_R[i], list_z[j]);

	  if(phi[i][j] < PhiMin)
	    PhiMin = phi[i][j];
        }
    }
  
  printf("\nPhiMin= %g\n",PhiMin);  fflush(stdout);
  printf("-G M / RH= %g\n", -G * M_HALO / RH);



  /* write this to a file so we have a record of it */
  if(strcmp(&OutputFile[strlen(OutputFile)-5],".hdf5") == 0)
    strncpy(phifile, OutputFile, strlen(OutputFile)-5);
  else
    strcpy(phifile, OutputFile);
  strcat(phifile, ".phi");
  if((fd = fopen(phifile,"w")))
    {
      fprintf(fd,"# \n");
      fprintf(fd,"# \n");
      fprintf(fd,"#      R             phi             phi             phi             phi             phi    \n");
      fprintf(fd,"#   (kpc/h)       @ %7.5g     @ %7.5g      @ %7.5g       @ %7.5g        @ %7.5g   \n",list_z[1],list_z[RSIZE/5],list_z[2*RSIZE/5],list_z[3*RSIZE/5],list_z[4*RSIZE/5]);
      fprintf(fd,"# \n");
      for(i = 0; i <= RSIZE; i++)
	fprintf(fd," %10.5f      %8.5e     %8.5e     %8.5e     %8.5e     %8.5e  \n",list_R[i],phi[1][i],phi[RSIZE/5][i],phi[2*RSIZE/5][i],phi[3*RSIZE/5][i],phi[4*RSIZE/5][i]);

      fclose(fd);
    }
}




void compute_vertical_force_field(void)
{
  int i, j;

  printf("Start computing vertical force field  (GravSoftening= %g, ErrTolTheta= %g).\n",GravSoftening,ErrTolTheta); fflush(stdout);


  for(i = 0; i <= RSIZE; i++)
    {
      //printf("\tcomputing radial bin %d of %d...\r", i, RSIZE);
      printf("\tcomputing radial bin %d of %d...\n", i, RSIZE);
      fflush(stdout);

      for(j = 0; j <= ZSIZE; j++)
	{
	  if(j == 0)
	    Dphi_z[i][j] = 0;
	  else
	    Dphi_z[i][j] = comp_Dphi_z(list_R[i], list_z[j]);
	}

    }
  printf("\n"); 

  /* write this to a file so we have a record of it */
  /*
    FILE *fd;
    char dffile[50]="";

    strncpy(dffile, OutputFile, strlen(OutputFile)-4);
    strcat(dffile, ".dphi_z");
    if((fd = fopen(dffile,"w")))
    {
    fprintf(fd,"# \n");
    fprintf(fd,"# \n");
    fprintf(fd,"#      z            Dphi_z          Dphi_z          Dphi_z          Dphi_z          Dphi_z  \n");
    fprintf(fd,"#   (kpc/h)       @ %7.5g     @ %7.5g      @ %7.5g       @ %7.5g        @ %7.5g   \n",list_R[1],list_R[100],list_R[200],list_R[300],list_R[400]);
    fprintf(fd,"# \n");
    for(i = 0; i <= ZSIZE; i++)
    fprintf(fd," %10.5f      %8.5e      %8.5e    %8.5e     %8.5e     %8.5e  \n",list_z[i],Dphi_z[i][1],Dphi_z[i][100],Dphi_z[i][200],Dphi_z[i][300],Dphi_z[i][400]);
    fclose(fd);
    }
  */

}



void compute_radial_force_field(void)
{
  int i, j;
  double k2, dphi_R_dr;
  FILE *fd;
  char dffile[100]="";

  printf("Start computing radial and vertical force field     (GravSoftening= %g, ErrTolTheta= %g).\n",GravSoftening,ErrTolTheta);  fflush(stdout);


  for(i = 0; i <= RSIZE; i++)
    {
      //printf("\tcomputing radial bin %d of %d...\r", i, RSIZE);
      printf("\tcomputing radial bin %d of %d...\n", i, RSIZE);
      fflush(stdout);

      for(j = 0; j <= ZSIZE; j++)
	{
	  if(j == 0)
	    {
	      Dphi_z[i][j] = 0;
	      Dphi_z_dR[i][j] = 0;
	    }
	  else
	    {
	      Dphi_z[i][j] = comp_Dphi_z(list_R[i], list_z[j]);
	      Dphi_z_dR[i][j] = comp_Dphi_z(list_RplusdR[i], list_z[j]);
	    }

	  if(i == 0)
	    Dphi_R[i][j] = 0;
	  else
	    Dphi_R[i][j] = comp_Dphi_R(list_R[i], list_z[j]);
	}
    }

  /* write this to a file so we have a record of it */
  strcpy(dffile, "");
  if(strcmp(&OutputFile[strlen(OutputFile)-5],".hdf5") == 0)
    strncpy(dffile, OutputFile, strlen(OutputFile)-5);
  else
    strcpy(dffile, OutputFile);
  strcat(dffile, ".dphi_R");
  printf("writing file: %s  .... \t", dffile);
  if((fd = fopen(dffile,"w")))
    {
      fprintf(fd,"# \n");
      fprintf(fd,"# \n");
      fprintf(fd,"#      R            Dphi_R          Dphi_R          Dphi_R          Dphi_R          Dphi_R  \n");
      fprintf(fd,"#   (kpc/h)       @ %7.5g     @ %7.5g      @ %7.5g       @ %7.5g        @ %7.5g   \n",list_z[1],list_z[RSIZE/5],list_z[2*RSIZE/5],list_z[3*RSIZE/5],list_z[4*RSIZE/5]);
      fprintf(fd,"# \n");
      for(i = 0; i <= RSIZE; i++)
	fprintf(fd," %10.5f      %8.5e     %8.5e     %8.5e     %8.5e     %8.5e  \n",list_R[i],Dphi_R[1][i],Dphi_R[RSIZE/5][i],Dphi_R[2*RSIZE/5][i],Dphi_R[3*RSIZE/5][i],Dphi_R[4*RSIZE/5][i]);
      fclose(fd);
    }
  printf("done.\n");



  /* write this to a file so we have a record of it */
  strcpy(dffile, "");
  if(strcmp(&OutputFile[strlen(OutputFile)-5],".hdf5") == 0)
    strncpy(dffile, OutputFile, strlen(OutputFile)-5);
  else
    strcpy(dffile, OutputFile);
  strcat(dffile, ".dphi_z");
  printf("writing file: %s  ....  \t", dffile);
  if((fd = fopen(dffile,"w")))
    {
      fprintf(fd,"# \n");
      fprintf(fd,"# \n");
      fprintf(fd,"#      z            Dphi_z          Dphi_z          Dphi_z          Dphi_z          Dphi_z  \n");
      fprintf(fd,"#   (kpc/h)       @ %7.5g     @ %7.5g      @ %7.5g       @ %7.5g        @ %7.5g   \n",list_R[1],list_R[RSIZE/5],list_R[2*RSIZE/5],list_R[3*RSIZE/5],list_R[4*RSIZE/5]);
      fprintf(fd,"# \n");
      for(i = 0; i <= ZSIZE; i++)
	fprintf(fd," %10.5f      %8.5e      %8.5e    %8.5e     %8.5e     %8.5e  \n",list_z[i],Dphi_z[i][1],Dphi_z[i][RSIZE/5],Dphi_z[i][2*RSIZE/5],Dphi_z[i][3*RSIZE/5],Dphi_z[i][3*RSIZE/5]);
      fclose(fd);
    }
  printf("done.\n");


  for(i = 1, epi_gamma2[0] = 1; i <= RSIZE; i++)
    {
      dphi_R_dr = comp_Dphi_R(list_RplusdR[i], 0);

      k2 = 3 / list_R[i] * Dphi_R[i][0] + (dphi_R_dr - Dphi_R[i][0]) / (list_RplusdR[i] - list_R[i]);

      epi_gamma2[i] = 4 / list_R[i] * Dphi_R[i][0] / k2;
      epi_kappa2[i] = k2;
    }

  epi_kappa2[0] = epi_kappa2[1];

  printf("Force field finished.\n");
}






double comp_Dphi_z(double R, double z)
{
  return comp_Dphi_Z_disk_tree(R, z) + comp_Dphi_z_halo(R, z) + comp_Dphi_z_bulge(R, z);
  /*
    return comp_Dphi_z_gas(R, z) + comp_Dphi_z_disk(R, z) + comp_Dphi_z_halo(R, z) + comp_Dphi_z_bulge(R, z);
  */
}

double comp_Dphi_R(double R, double z)
{

  return comp_Dphi_R_disk_tree(R, z) + comp_Dphi_R_halo(R, z) + comp_Dphi_R_bulge(R, z);

  /*
    return comp_Dphi_R_gas(R, z) + comp_Dphi_R_disk(R, z) + comp_Dphi_R_halo(R, z) + comp_Dphi_R_bulge(R, z);
  */
  /*
    double Dp1, Dp2, Dp3;

    Dp1= comp_Dphi_R_disk_tree(R, z);
    Dp2= comp_Dphi_R_halo(R, z);
    Dp3= comp_Dphi_R_bulge(R, z);
    printf("R= %g  z=%g      Dp1= %g  Dp2= %g  Dp3= %g\n",R,z,Dp1,Dp2,Dp3); fflush(stdout);
    return Dp1 + Dp2 + Dp3;
  */
}






double comp_phi(double R, double z)
{
  return comp_phi_disk_tree(R, z) + comp_phi_halo(R, z) + comp_phi_bulge(R, z);
}



