#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"


double comp_DF_halo_exact(double E, double L);
double comp_DF_bulge_exact(double E, double L);

double comp_DF_H_iso(double E);
double comp_DF_H_ani(double Q);

double comp_DF_halo_Eddington(double Q);
double comp_DF_bulge_Eddington(double Q);

void comp_DF_init(void);
void compute_d2rhodpsi2_halo(void);
void compute_d2rhodpsi2_bulge(void);

double eddington_integrand(double t, void *params);


#define WORKSIZE 100000
gsl_integration_workspace *Workspace;

double *xi, *yi, *eddint;
double *list_radius;
double *list_E;                            /* tabulated energies */
double *DistFunc_halo, *DistFunc_bulge;    /* distribution function(s) */
double *psi_R;
double *rho_R_halo,*rho_R_bulge;
double *drhodpsi_halo,*d2rhodpsi2_halo;
double *drhodpsi_bulge,*d2rhodpsi2_bulge;


/* size of DF look-up table */
int DFSIZE= RSIZE * 2;


void compute_DF_lookuptable(void)
{
  FILE *fd;
  char dffile[100]="";
  int i;

  printf("Start computing distribution function.\n"); fflush(stdout);

#ifdef MAXWELLIAN
  AnisotropyRadius= -1;
  printf("MAXWELLIAN is turned on\n");
#else
  printf("MAXWELLIAN is turned off\n");
#endif

#ifdef DF_H_MODEL
  printf("DF_H_MODEL is turned on\n");
#else
  printf("DF_H_MODEL is turned off\n");
#endif

#ifdef DF_EDDINGTON
  printf("DF_EDDINGTON is turned on\n");
#else
  printf("DF_EDDINGTON is turned off\n");
#endif

  printf("AnisotropyRadius= %g\n", AnisotropyRadius);
  printf("ra= %g\n", AnisotropyRadius * RH);

  comp_DF_init();

  for(i = 0; i <= DFSIZE; i++)
    {
	printf("Computing DF table bin %d/%d (e= %g, r=%g)\n", i, DFSIZE, -1.0*psi_R[i],list_radius[i]); fflush(stdout);

	list_E[i]= -1.0 * psi_R[i];

	DistFunc_halo[i] = comp_DF_halo_exact(list_E[i], 0.0);
	DistFunc_bulge[i] = comp_DF_bulge_exact(list_E[i], 0.0);
    }
  printf("\n");


  /* write this to a file so we have a record of it */
  if(strcmp(&OutputFile[strlen(OutputFile)-5],".hdf5") == 0)
    strncpy(dffile, OutputFile, strlen(OutputFile)-5);
  else
    strcpy(dffile, OutputFile);
  strcat(dffile, ".df");
  printf("DF look-up table saved in file: %s\n",dffile);
  if((fd = fopen(dffile,"w")))
    {
	fprintf(fd,"   Energy       f_halo(E)    f_bulge(E)\n");
	for(i= 0; i<= DFSIZE; i++)
	  fprintf(fd," %10.5e   %10.5e   %10.5e\n",list_E[i],DistFunc_halo[i],DistFunc_bulge[i]);
	fclose(fd);
    }

}










double comp_DF_halo_exact(double E, double L)
{
  double Q, f_e= 0.0;

  if((N_HALO == 0) || (M_HALO == 0))
    return 0;


  if (AnisotropyRadius > 0.0)
	Q= - E - (L*L/2./(AnisotropyRadius * RH)/(AnisotropyRadius * RH));
  else
	Q= - E;


#ifdef DF_H_MODEL
	/* calculate f_e exactly from the H90 analytic formula */
      if (AnisotropyRadius > 0.0)
        f_e= comp_DF_H_ani(Q);
      else
        f_e= comp_DF_H_iso(E);
#endif


#ifdef DF_EDDINGTON
	/* calculate f_e using Eddington's formula */
	f_e= comp_DF_halo_Eddington(Q);
	printf("f_e_edd(Q)/f_H_iso(E)= %g\n",f_e/comp_DF_H_iso(E));
	printf("f_e_edd(Q)/f_H_ani(Q)= %g\n",f_e/comp_DF_H_ani(Q)); fflush(stdout);
	printf("------\n");
#endif


  return f_e;
}





double comp_DF_bulge_exact(double E, double L)
{
  double Q, f_e= 0.0;

  if((N_BULGE == 0) || (M_BULGE == 0))
    return 0;


/* for the moment, we only assume isotropic DF and use Eddington's formula */

/*
  if (AnisotropyRadius > 0.0)
        Q= - E - (L*L/2./(AnisotropyRadius * RH)/(AnisotropyRadius * RH));
  else
*/
        Q= - E;

  f_e= comp_DF_bulge_Eddington(Q);

  return f_e;
}











/*
   interpolate the halo f_e from the pre-computed table
*/
double comp_DF_halo(double E, double L)
{
  double f_e= 0.0, Q, ee, ue;
  int ie;

  ee= E;

  if(AnisotropyRadius > 0.0)
    {
      Q= - E - (L*L/2./(AnisotropyRadius * RH)/(AnisotropyRadius * RH));

      /* in this case, table is stored in terms of -Q */
      ee= -1.0 * Q;
    }

  /* if e>0 or q<0 then dist. func. is zero */
  if(ee>0.0)
    return 0.0;

  ie= find_idx(ee, list_E, DFSIZE);

/*
  printf("\nE= %g  Q= %g ee= %g\n",E, Q, ee);
  printf("ie= %d  list_E[0]= %g   list_E[DFSIZE]= %g \n",ie, list_E[0], list_E[DFSIZE]);
  printf("list_E[%d]= %g\n", ie-1, list_E[ie-1]);
  printf("list_E[%d]= %g\n", ie, list_E[ie]);
  printf("list_E[%d]= %g\n", ie+1, list_E[ie+1]);
  fflush(stdout);
*/

  if(ie < 0 || ie > DFSIZE)
    ie= 0;

  ue = (ee - list_E[ie]) / (list_E[ie + 1] - list_E[ie]);

  if(ee>list_E[0])
    {
      /*
      printf("\nparticle has energy very close to zero\n");
      printf("ee=  %g, table min list_E[0]= %g\n",ee,list_E[0]);
      printf("return 0.\n");
      */
      ue= 0.0;
    }

  if(ie < 0 || ie >= DFSIZE || ue < 0 || ue > 1)
    {
      printf("\nin comp_DF_halo: potential problem with lookup table\n");
      printf("      ee= %g   ie= %d   ue=%d\n", ee, ie, ue);
      printf("      list_E[0]= %g  list_E[DFSIZE]= %g\n",list_E[0], list_E[DFSIZE]);
      printf("      try to do this manually .....\n");

      f_e= comp_DF_halo_exact(ee, 0.0);

      printf("      f_e= %g\n\n",f_e);
      if(f_e<0)
        exit(0);
    }
  else
    {
      f_e= DistFunc_halo[ie] * (1-ue) + DistFunc_halo[ie+1] * ue;
    }

  return f_e;
}






/*
   interpolate the bulge f_e from the pre-computed table
*/
double comp_DF_bulge(double E, double L)
{
  double f_e= 0.0, Q, ee, ue;
  int ie;
  
  ee= E;
  
  /* for the moment, we'll assume the bulge is isotropic */
/*
  if(AnisotropyRadius > 0.0)
    {
      Q= - E - (L*L/2./(AnisotropyRadius * RH)/(AnisotropyRadius * RH));
      ee= -1.0 * Q;
    }
*/

  /* if e>0 or q<0 then dist. func. is zero */
  if(ee>0.0)
    return 0.0;

  ie= find_idx(ee, list_E, DFSIZE);

  if(ie < 0 || ie > DFSIZE)
    ie= 0;

  ue = (ee - list_E[ie]) / (list_E[ie + 1] - list_E[ie]);

  if(ee>list_E[0])
    {
      /*
      printf("\nparticle has energy very close to zero\n");
      printf("ee=  %g, table min list_E[0]= %g\n",ee,list_E[0]);
      printf("return 0.\n");
      */
      ue= 0.0;
    }

  if(ie < 0 || ie >= DFSIZE || ue < 0 || ue > 1)
    {
      printf("\nin comp_DF_bulge: potential problem with lookup table\n");
      printf("      ee= %g   ie= %d   ue=%d\n", ee, ie, ue);
      printf("      list_E[0]= %g  list_E[DFSIZE]= %g\n",list_E[0], list_E[DFSIZE]);
      printf("      try to do this manually .....\n");

      f_e= comp_DF_bulge_exact(ee, 0.0);

      printf("      f_e= %g\n\n",f_e);
      if(f_e<0)
        exit(0);
    }
  else
    {
      f_e= DistFunc_bulge[ie] * (1-ue) + DistFunc_bulge[ie+1] * ue;
    }
    
  return f_e;
}








/* ----------------------------------------------------------------- */



/* the following are analytic forms of the distribution
   function, currently, it's only Hernquist 1990 */

/* Isotropic Hernquist Model */

double comp_DF_H_iso(double E)
{
  double vg, q, prefac, f_e;

  if (E >= 0.0)
	return 0.0;

  vg= sqrt(G * M_HALO / RH);
  q= sqrt(-1.0 * E * RH / G / M_HALO);

  if (q >= 1.0)
	return (1.0e+30);

  prefac= M_HALO / (8.0 * sqrt(2.0) * PI * PI * PI * RH * RH * RH * vg * vg * vg);

  f_e= pow((1-q*q), -5./2.) * ( 3.*asin(q) + q*pow((1.-q*q),0.5)*(1.-2.*q*q) * (8.*q*q*q*q - 8.*q*q - 3.));
/*
  printf("               E= %g  vg= %g  q= %g  prefac= %g f_e %g\n", E, vg, q, prefac, f_e); fflush(stdout);
*/

  return prefac * f_e;
}



/* Anisotropic Hernquist Model */

double comp_DF_H_ani(double Q)
{

  if (AnisotropyRadius <= 0.0)
	return 0;

  double vg, qbar, prefac, f_e, f_iso;
  double ra;

  ra= AnisotropyRadius * RH;

  if (Q <= 0.0)
	return 0;

  vg= sqrt(G * M_HALO / RH);
  qbar= sqrt(RH * Q / G / M_HALO);

  if (qbar >= 1.0)
	return 1.0e+30;

  prefac= M_HALO / (sqrt(2.0) * PI * PI * PI * RH * RH * RH * vg * vg * vg);

  f_e= (RH * RH / ra / ra) * qbar * (1. - 2.*qbar*qbar);

  f_iso= comp_DF_H_iso(-1.0*Q);

/*
printf("DF_ani  Q= %g   qbar= %g  f_iso= %g  prefac= %g  f_e= %g\n",Q,qbar,f_iso,prefac,f_e); fflush(stdout);
*/

  return f_iso + prefac * f_e;
}







/* ----------------------------------------------------------------- 

     Now, calculate DF using the Eddington formula.

   ----------------------------------------------------------------- */

gsl_spline * d2rhodpsi2_spline_halo;
gsl_interp_accel * d2rhodpsi2_spline_acc_halo;

gsl_spline * d2rhodpsi2_spline_bulge;
gsl_interp_accel * d2rhodpsi2_spline_acc_bulge;


// Parameter struct for the Eddington integration function.
typedef struct {
  const gsl_spline* s;
  gsl_interp_accel* a;
  double Q;
} spline_params;



/* The integrand for the Eddington integration. Now uses 
   the GSL spline functionality rather than the NR, which is evil. */
double eddington_integrand(double t, void *params)
{
  spline_params* p = params;

//printf("t= %g\t", t);
//printf("Q= %g\t",p->Q);
//printf("sqrt(p->Q-t)= %g   \t", sqrt(p->Q-t));
//printf("d2rhodpsi2= %g  \t", gsl_spline_eval(p->s, t, p->a));
  double v2 = gsl_spline_eval(p->s, t, p->a) / sqrt(p->Q - t);
//printf("v2= %g  \n",v2); fflush(stdout);
   
  return v2;
}







double comp_DF_halo_Eddington(double Q)
{
  double f_e;
  int i;
  double result, abserr;


  /* first, we'll generate the spline for d2rho/dpsi2 term in integrand */


  if(Q<=0)
    return 0.0;


  /* load spline for d2rho/dpsi2 term into function parameters */
  spline_params p;
  p.s = d2rhodpsi2_spline_halo;
  p.a = d2rhodpsi2_spline_acc_halo;
  p.Q = Q;



  /* error checking: write entire DF integrand and the spline */
/*
  FILE *fd;
  char tempfile[100]="";
  sprintf(tempfile, "df_integrands/%g.txt", Q);
  if((fd = fopen(tempfile,"w")))
    {
        fprintf(fd,"Q= %g \n",Q);
        fprintf(fd,"       xi            yi          d2rho_dpsi2            Q           Q-psiR     sqrt(Q-psiR) \n");
        for(i = 0; i <= DFSIZE; i++)
          fprintf(fd," %8.5e     %8.5e     %8.5e    %8.5e    %8.5e      %8.5e     \n",xi[i],yi[i],d2rho_dpsi2[i], Q, Q-psi_R[i], sqrt(Q-psi_R[i]));
        fclose(fd);
    }
*/



  gsl_function F;
  F.function = &eddington_integrand;
  F.params = &p;


  double epsabs=0;
  double epsrel=1e-4;
  int intstatus;

  /* do our own error handling, so that we can adjust the accuracy if needed (i.e., mainly at large r) */
  gsl_error_handler_t * old_handler = gsl_set_error_handler_off(); /* turn off the error-handler */

  do
  {
    intstatus= gsl_integration_qags(&F, 0.0, Q, epsabs, epsrel, WORKSIZE, Workspace, &result, &abserr);
    printf("epsrel = %g     result = %g +/- %g   (%g %%)\n", epsrel, result, abserr, 100.0*abserr/result);
    epsrel *= 2.0;
  } while(intstatus == GSL_EROUND);

  //printf("Workspace.intervals = %d\n", Workspace->size);
  //printf("Workspace.maximum_level = %d\n", Workspace->maximum_level);


  gsl_set_error_handler(old_handler);          /*  turn it on again */


  f_e = result / sqrt(8.0) / PI / PI;

  if(f_e<0)
    f_e= 0.0;


  return f_e;
}





double comp_DF_bulge_Eddington(double Q)
{
  double f_e;
  int i;
  double result, abserr;



  if(Q<=0)
    return 0.0;

  spline_params p;
  p.s = d2rhodpsi2_spline_bulge;
  p.a = d2rhodpsi2_spline_acc_bulge;
  p.Q = Q;


  gsl_function F;
  F.function = &eddington_integrand;
  F.params = &p;


  double epsabs=0;
  double epsrel=1e-4;
  int intstatus;

  gsl_error_handler_t * old_handler = gsl_set_error_handler_off(); /* turn off the error-handler */

  do
  {
    intstatus= gsl_integration_qags(&F, 0.0, Q, epsabs, epsrel, WORKSIZE, Workspace, &result, &abserr);
    //printf("epsrel = %g     result = %g +/- %g   (%g %)\n", epsrel, result, abserr, 100.0*abserr/result);
    //printf("error = %g\n", abserr);
    epsrel *= 2.0;
  } while(intstatus == GSL_EROUND);

  gsl_set_error_handler(old_handler);          /*  turn it on again */


  f_e = result / sqrt(8.0) / PI / PI;

  if(f_e<0)
    f_e= 0.0;


  return f_e;
}









/*
   Computes d2_rho / dpsi2, which is the primary term
   in the integral when using the Eddington formulation
   to calculate the distribution function.

   We do this for the halo and bulge components separately.

*/
void compute_d2rhodpsi2_halo()
{
  FILE *fd;
  char drhodpsifile[100]="";
  int i;

  double R, RplusdR, RminusdR, dR;
  double rho_RplusdR, rho_RminusdR;
  double psi_RplusdR, psi_RminusdR;
  double slope1, slope2, ra;

  // normalize just in case we get really small densities
  double rho_norm, psi_norm;


  for(i = 0; i <= DFSIZE; i++)
    {
	if(i < RSIZE)
	    /* continues same scaling as list_R 
            R= list_R[RSIZE] * exp((RSIZE-i) * (log(LL) - log(Baselen)) / (RSIZE - 1));   */
	    /* more modest scaling to zero energy (ie larger radii) */
            R= list_R[RSIZE] * pow(10., 5.0 * (RSIZE-i) / RSIZE);
            //R= list_R[RSIZE] * pow(10., 7.0 * (RSIZE-i) / RSIZE);
	else
            R= list_R[DFSIZE-i];

        if(i == DFSIZE)
            R= list_R[1] * 0.5;


	/* interval of integration */
	dR= R * exp(1.0e-2 * (log(LL) - log(Baselen)) / (RSIZE - 1)) - R;
	if(dR>1.5)
	  dR= 1.5;
	if((dR/R)<1.0e-4)
	  dR=1.0e-4*R;
        RplusdR= R + dR;
        RminusdR= R - dR;

	list_radius[i]= R;

        /* compute psi (=-phi, i.e., negative the potential) along the z axis, although
           the full DF assumes spherical symmetry */
        psi_R[i]= -1.0 * comp_phi(0, R);;
        psi_RplusdR= -1.0 * comp_phi(0, RplusdR);
        psi_RminusdR= -1.0 * comp_phi(0, RminusdR);

        /* -------------------------------
	    DF for halo component  */
        rho_R_halo[i]= comp_rho_halo(R,0);
        rho_RplusdR= comp_rho_halo(RplusdR,0);
        rho_RminusdR= comp_rho_halo(RminusdR,0);

	if(AnisotropyRadius>0)
	  {
	    ra= AnisotropyRadius * RH;
            rho_R_halo[i] *= (1.0 + (R * R / ra / ra ));
            rho_RplusdR *= (1.0 + (RplusdR * RplusdR / ra / ra));
            rho_RminusdR *= (1.0 + (RminusdR * RminusdR / ra / ra));
	  }

	// normalize
	rho_norm= rho_R_halo[i];
	rho_R_halo[i] /= rho_norm;
	rho_RplusdR /= rho_norm;
	rho_RminusdR /= rho_norm;
	psi_norm= psi_R[i];
	psi_R[i] /= psi_norm;
	psi_RplusdR /= psi_norm;
	psi_RminusdR /= psi_norm;

        slope1= (rho_RplusdR - rho_R_halo[i]) / (psi_RplusdR - psi_R[i]);
        slope2= (rho_R_halo[i] - rho_RminusdR) / (psi_R[i] - psi_RminusdR);

        drhodpsi_halo[i]= 0.5*(slope1+slope2);

        d2rhodpsi2_halo[i]= (rho_RplusdR+rho_RminusdR-2.*rho_R_halo[i])/(psi_RplusdR-psi_R[i])/(psi_R[i]-psi_RminusdR);

	// now put units back in
	psi_R[i] *= psi_norm;
	rho_R_halo[i] *= rho_norm;
	drhodpsi_halo[i] *= rho_norm / psi_norm;
	d2rhodpsi2_halo[i] *= rho_norm / psi_norm / psi_norm;

//printf("i= %d, R= %g,  second derivative denominators: %g, %g    diff= %g\n",i,R,(psi_RplusdR-psi_R[i]),(psi_R[i]-psi_RminusdR),(psi_RplusdR-psi_R[i])/(psi_R[i]-psi_RminusdR));
//printf("i= %d, R= %g,  rho= %g   drhodpsi_halo= %g    d2rhodpsi2_halo= %g     second derivative denominators: %g, %g    diff= %g\n",i,R,rho_R_halo[i],drhodpsi_halo[i],d2rhodpsi2_halo[i],(psi_RplusdR-psi_R[i]),(psi_R[i]-psi_RminusdR),(psi_RplusdR-psi_R[i])/(psi_R[i]-psi_RminusdR));
//printf("i= %d, R= %g,    rho_0= %g   ani_factor= %g  rho_R_halo= %g\n",i,R,comp_rho_halo(R,0),(1.0 + (R * R / ra / ra )),rho_R_halo[i]);


    }




  /* generate spline for the d2rho / dpsi2 term */
  for(i = 0; i <= DFSIZE; i++)
    {
	xi[i + 1] = psi_R[i];
	yi[i + 1] = d2rhodpsi2_halo[i];
    }

  gsl_spline_init(d2rhodpsi2_spline_halo, xi, yi, DFSIZE+1);




  /* write this to a file so we have a record of it */
  strcpy(drhodpsifile,"");
  if(strcmp(&OutputFile[strlen(OutputFile)-5],".hdf5") == 0)
    strncpy(drhodpsifile, OutputFile, strlen(OutputFile)-5);
  else
    strcpy(drhodpsifile, OutputFile);
  strcat(drhodpsifile, ".drhodpsi_halo");
  if((fd = fopen(drhodpsifile,"w")))
    {
        fprintf(fd,"# drhodpsi file, halo component\n");
        fprintf(fd,"# n= %d \n",DFSIZE);
        fprintf(fd,"# \n");
        fprintf(fd,"#   R (kpc)         psi               rho         drho/dpsi      d2rho/dpsi2   spline(psi)   \n");
        fprintf(fd,"# \n");
        for(i = 0; i <= DFSIZE; i++)
          fprintf(fd," %8.5e     %8.5e      %8.5e     %8.5e     %8.5e    %8.5e  \n",list_radius[i],psi_R[i],rho_R_halo[i],drhodpsi_halo[i],d2rhodpsi2_halo[i],gsl_spline_eval(d2rhodpsi2_spline_halo, 0.999*psi_R[i], d2rhodpsi2_spline_acc_halo));

        fclose(fd);
    }

}




void compute_d2rhodpsi2_bulge()
{
  FILE *fd;
  char drhodpsifile[100]="";
  int i;

  double R, RplusdR, RminusdR, dR;
  double rho_RplusdR, rho_RminusdR;
  double psi_RplusdR, psi_RminusdR;
  double slope1, slope2, ra;

  // normalize just in case we get really small densities
  double rho_norm, psi_norm;


  if(M_BULGE == 0)
    return;


  for(i = 0; i <= DFSIZE; i++)
    {
	if(i < RSIZE)
	    /* continues same scaling as list_R 
            R= list_R[RSIZE] * exp((RSIZE-i) * (log(LL) - log(Baselen)) / (RSIZE - 1));   */
	    /* more modest scaling to zero energy (ie larger radii) */
            R= list_R[RSIZE] * pow(10., 5.0 * (RSIZE-i) / RSIZE);
            //R= list_R[RSIZE] * pow(10., 7.0 * (RSIZE-i) / RSIZE);
	else
            R= list_R[DFSIZE-i];

        if(i == DFSIZE)
            R= list_R[1] * 0.5;


	/* interval of integration */
	dR= R * exp(1.0e-2 * (log(LL) - log(Baselen)) / (RSIZE - 1)) - R;
	if(dR>1.5)
	  dR= 1.5;
	if((dR/R)<1.0e-4)
	  dR=1.0e-4*R;
        RplusdR= R + dR;
        RminusdR= R - dR;

	list_radius[i]= R;

        /* compute psi (=-phi, i.e., negative the potential) along the z axis, although
           the full DF assumes spherical symmetry */
        psi_R[i]= -1.0 * comp_phi(0, R);;
        psi_RplusdR= -1.0 * comp_phi(0, RplusdR);
        psi_RminusdR= -1.0 * comp_phi(0, RminusdR);

        /* -------------------------------
	    DF for bulge component  */
        rho_R_bulge[i]= comp_rho_bulge(R,0);
        rho_RplusdR= comp_rho_bulge(RplusdR,0);
        rho_RminusdR= comp_rho_bulge(RminusdR,0);

	if(AnisotropyRadius>0)
	  {
	    ra= AnisotropyRadius * RH;
            rho_R_bulge[i] *= (1.0 + (R * R / ra / ra ));
            rho_RplusdR *= (1.0 + (RplusdR * RplusdR / ra / ra));
            rho_RminusdR *= (1.0 + (RminusdR * RminusdR / ra / ra));
	  }

	// normalize
	rho_norm= rho_R_bulge[i];
	rho_R_bulge[i] /= rho_norm;
	rho_RplusdR /= rho_norm;
	rho_RminusdR /= rho_norm;
	psi_norm= psi_R[i];
	psi_R[i] /= psi_norm;
	psi_RplusdR /= psi_norm;
	psi_RminusdR /= psi_norm;

        slope1= (rho_RplusdR - rho_R_bulge[i]) / (psi_RplusdR - psi_R[i]);
        slope2= (rho_R_bulge[i] - rho_RminusdR) / (psi_R[i] - psi_RminusdR);

        drhodpsi_bulge[i]= 0.5*(slope1+slope2);

        d2rhodpsi2_bulge[i]= (rho_RplusdR+rho_RminusdR-2.*rho_R_bulge[i])/(psi_RplusdR-psi_R[i])/(psi_R[i]-psi_RminusdR);

	// now put units back in
	psi_R[i] *= psi_norm;
	rho_R_bulge[i] *= rho_norm;
	drhodpsi_bulge[i] *= rho_norm / psi_norm;
	d2rhodpsi2_bulge[i] *= rho_norm / psi_norm / psi_norm;

//printf("i= %d, R= %g,  second derivative denominators: %g, %g    diff= %g\n",i,R,(psi_RplusdR-psi_R[i]),(psi_R[i]-psi_RminusdR),(psi_RplusdR-psi_R[i])/(psi_R[i]-psi_RminusdR));
//printf("i= %d, R= %g,  rho= %g   drhodpsi_bulge= %g    d2rhodpsi2_bulge= %g     second derivative denominators: %g, %g    diff= %g\n",i,R,rho_R_bulge[i],drhodpsi_bulge[i],d2rhodpsi2_bulge[i],(psi_RplusdR-psi_R[i]),(psi_R[i]-psi_RminusdR),(psi_RplusdR-psi_R[i])/(psi_R[i]-psi_RminusdR));
//printf("i= %d, R= %g,    rho_0= %g   ani_factor= %g  rho_R_bulge= %g\n",i,R,comp_rho_bulge(R,0),(1.0 + (R * R / ra / ra )),rho_R_bulge[i]);


    }




  /* generate spline for the d2rho / dpsi2 term */
  for(i = 0; i <= DFSIZE; i++)
    {
	xi[i + 1] = psi_R[i];
	yi[i + 1] = d2rhodpsi2_bulge[i];
    }

  gsl_spline_init(d2rhodpsi2_spline_bulge, xi, yi, DFSIZE+1);




  /* write this to a file so we have a record of it */
  strcpy(drhodpsifile,"");
  if(strstr(OutputFile,".hdf5"))
        strncpy(drhodpsifile, OutputFile, strlen(OutputFile)-5);
  if(strstr(OutputFile,".dat"))
        strncpy(drhodpsifile, OutputFile, strlen(OutputFile)-4);
  strcat(drhodpsifile, ".drhodpsi_bulge");
  if((fd = fopen(drhodpsifile,"w")))
    {
        fprintf(fd,"# drhodpsi file, bulge component\n");
        fprintf(fd,"# n= %d \n",DFSIZE);
        fprintf(fd,"# \n");
        fprintf(fd,"#   R (kpc)         psi               rho         drho/dpsi      d2rho/dpsi2   spline(psi)   \n");
        fprintf(fd,"# \n");
        for(i = 0; i <= DFSIZE; i++)
          fprintf(fd," %8.5e     %8.5e      %8.5e     %8.5e     %8.5e    %8.5e  \n",list_radius[i],psi_R[i],rho_R_bulge[i],drhodpsi_bulge[i],d2rhodpsi2_bulge[i],gsl_spline_eval(d2rhodpsi2_spline_bulge, 0.999*psi_R[i], d2rhodpsi2_spline_acc_bulge));

        fclose(fd);
    }

}









void comp_DF_init(void)
{
  /* stores spline on the Eddington integrand */
  xi= vector(1, DFSIZE+1);
  yi= vector(1, DFSIZE+1);
  eddint= vector(1, DFSIZE+1);

  /* arrays for energies, radii, d2rho/dpsi2 term in Eddington integrand, and pre-computed DF  */
  list_radius= vector(0, DFSIZE);
  list_E= vector(0, DFSIZE);
  psi_R= vector(0, DFSIZE);

  /* halo specific */
  rho_R_halo= vector(0, DFSIZE);
  drhodpsi_halo= vector(0, DFSIZE);
  d2rhodpsi2_halo= vector(0, DFSIZE);
  DistFunc_halo= vector(0, DFSIZE);

  /* bulge specific */
  rho_R_bulge= vector(0, DFSIZE);
  drhodpsi_bulge= vector(0, DFSIZE);
  d2rhodpsi2_bulge= vector(0, DFSIZE);
  DistFunc_bulge= vector(0, DFSIZE);


  /* setup the d2rho_dpsi2 array */

  /* for the halo */
  d2rhodpsi2_spline_halo = gsl_spline_alloc(gsl_interp_cspline, DFSIZE+1);
  d2rhodpsi2_spline_acc_halo = gsl_interp_accel_alloc();
  compute_d2rhodpsi2_halo();

  /* for the bulge */
  d2rhodpsi2_spline_bulge = gsl_spline_alloc(gsl_interp_cspline, DFSIZE+1);
  d2rhodpsi2_spline_acc_bulge = gsl_interp_accel_alloc();
  compute_d2rhodpsi2_bulge();

  Workspace = gsl_integration_workspace_alloc(WORKSIZE);
}



void comp_DF_Eddington_close(void)
{
  free_vector(xi, 1, DFSIZE+1);
  free_vector(yi, 1, DFSIZE+1);
  free_vector(eddint, 1, DFSIZE+1);
  free_vector(psi_R, 1, DFSIZE);
  free_vector(list_radius, 1, DFSIZE);
  free_vector(list_E, 1, DFSIZE);
  free_vector(rho_R_halo, 1, DFSIZE);
  free_vector(rho_R_bulge, 1, DFSIZE);
  free_vector(drhodpsi_halo, 1, DFSIZE);
  free_vector(drhodpsi_bulge, 1, DFSIZE);
  free_vector(d2rhodpsi2_halo, 1, DFSIZE);
  free_vector(d2rhodpsi2_bulge, 1, DFSIZE);
  free_vector(DistFunc_halo, 1, DFSIZE);
  free_vector(DistFunc_bulge, 1, DFSIZE);


  gsl_integration_workspace_free(Workspace);

  gsl_spline_free(d2rhodpsi2_spline_halo);
  gsl_interp_accel_free(d2rhodpsi2_spline_acc_halo);

  gsl_spline_free(d2rhodpsi2_spline_bulge);
  gsl_interp_accel_free(d2rhodpsi2_spline_acc_bulge);
}




/* ----------------------------------------------------------------- */



double compute_ani_beta(double r)
{
  double beta;
  
  /* isotropic model */
  beta = 0.0;
  
  if (AnisotropyRadius > 0)
    {
        double ra;

        ra= AnisotropyRadius * RH;

        /* Osipkov & Merritt anisotropy dependency */
        beta = r * r / (r * r + ra * ra);
    }
    
  return beta;
} 




