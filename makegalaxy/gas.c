#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"


static double R,z;

void compute_vstream_gas(void)
{
  int i, j;
  double P1, P2, vphipress, vphi;


  printf("comp vstream gas ... \t");

  for(i = 0; i < RSIZE; i++)
    {
      for(j = 0; j <= ZSIZE; j++)
	{
	  if(i == 0 || i == RSIZE)
	    VelStreamGas[i][j] = 0;
	  else
	    {
	      VelStreamGas[i][j] = 0;

	      vphi = Dphi_R[i][j];

	      P2 = eqn_of_state(RhoGas[i][j]);
	      P1 = eqn_of_state(1.1*RhoGas[i][j]);

	      if(RhoGas[i][j] > 0 && RhoGas[i+1][j]> 0 )
		{
		  vphipress = - log(P1 / P2) / log(1.1) * P2 / RhoGas[i][j] / H;

                  if(vphi+vphipress<0) 
                    {
                      /*
                      printf("r=%g vphi=%g vphipress=%g j=%d  %g %g %g %g  %g %g %g %g\n", list_R[i], sqrt(vphi), sqrt(-vphipress), j, 
                             P1,  P2, RhoGas[i][j], RhoGas[i+1][j],
                             log(P1 / P2) / log(1.1),
                             log(RhoGas[i+1][j] / RhoGas[i][j]), log(list_R[i+1]/list_R[i]),
                             list_R[i]);
                      */
                      vphi = 0;
                    }
                  else
                    vphi += vphipress;
		}

	      if(vphi > 0)
		{
		  VelStreamGas[i][j] = sqrt(vphi * list_R[i]);
		      
		  /*
		    if(j==0)
		    printf("R=%g  vphi=%g  %g P2=%g P1=%g rho2=%g rho1=%g  delta=%g\n", list_R[i], 
		    VelStreamGas[i][j], sqrt(Dphi_R[i][j] * list_R[i]),
		    P2, P1, RhoGas[i][j], RhoGas[i-1][j],
		    log(P2/P1)/log(list_R[i]/list_R[i-1]) * P2/ RhoGas[i][j]/list_R[i]);
		  */
		}
	    }
	}
    }

  for(j = 0; j <= ZSIZE; j++)
    {
      VelStreamGas[1][j] = VelStreamGas[2][j];
    }

  printf("done\n");
}






double mass_cumulative_gas(double R)
{
      switch(GasDistribution) {
          case 0:   /* Standard Exponential Distribution */
                return M_GAS*(1-(1+R/H)*exp(-R/H));
                break;
          case 1:  /* exponential with Rd x times the Disk one */
                return M_GAS*(1-(1+R/(GasExpAlpha*H))*exp(-R/(GasExpAlpha*H)));
                break;
          case 2:  /* Power-law distribution - gamma=1 is Mestel, or 1/R distribution */
                if(R>PowerLawCutoff)
                  return M_GAS;
                else
                  return M_GAS*pow((R/PowerLawCutoff), 2-PowerLawGamma);
                break;
        }

	return 0.0;
}






// XXX This function is actually never used!
double comp_rho_gas(double R,double z)
{
  double x=0;

  if(fabs(z)>6*Z0)
        return x;


  switch(GasDistribution) {
        case 0:
          x=(M_GAS)/(4*PI*H*H*Z0)*exp(-R/H)*pow(2/(exp(z/Z0)+exp(-z/Z0)),2);
          break;
        case 1:
          x=(M_GAS)/(4*PI*GasExpAlpha*H*GasExpAlpha*H*Z0)*exp(-R/(GasExpAlpha*H))*pow(2/(exp(z/Z0)+exp(-z/Z0)),2);
          break;
        case 2:
          x=(M_GAS)/(4*PI*Z0*PowerLawGamma*PowerLawGamma)*(2-PowerLawGamma)*pow((R/PowerLawCutoff),PowerLawGamma)*pow(2/(exp(z/Z0)+exp(-z/Z0)),2);
          break;
    }

  return x;
}




/* -------------------------------------
 *    Bessel function integrads for
 *    calculating potential of a 
 *    sheet of mass (as a function of 
 *    surface density).
 *
 *    binney & tremaine, pp. 74-79
 *      equation 2-167, with derivatives
 *      taken with respect to z and R
 *
 * ------------------------------------- */


double intz_g(double k);
double intz_g_abs(double);

double intR_g(double k);
double intR_g_abs(double k);





/* -------------------------------------
     Compute:  dphi/dz
   ------------------------------------- */
  
double comp_Dphi_z_gas(double RR,double zz)
{ 
  double comp_Dphi_z_gas_sph(double RR,double zz);
  double comp_Dphi_z_gas_exact(double RR,double zz);
  double Rd=0.0;
    
  switch(GasDistribution) {
        case 0:
          Rd= H;
          break;
        case 1:
          Rd= H*GasExpAlpha;
          break;
        case 2:
          break;
  }

  if(sqrt(RR*RR+zz*zz)>10*Rd)
    return comp_Dphi_z_gas_sph(RR,zz);
  else
    return comp_Dphi_z_gas_exact(RR,zz);

} 


/* assume spherical mass density */
double comp_Dphi_z_gas_sph(double RR,double zz)
{
  double m;
  double r;


  r=sqrt(RR*RR+zz*zz);

  m=mass_cumulative_gas(r);

  return G*zz/(r*r*r)*m;
}



/* This is only for the exponential distributions at the moment,
   we will prepare this for the power law one later  */
double comp_Dphi_z_gas_exact(double RR,double zz)
{
  int i;
  double dphiz;
  double Sigma0;
  double in1,in2,in3,bb;
  double deltaz,zpos;
  double Rd=0.0;


  if(N_GAS==0) return 0;


  switch(GasDistribution) {
        case 0:
          Rd= H;
          break;
        case 1:
          Rd= H*GasExpAlpha;
          break;
        case 2:
          return comp_Dphi_z_gas_sph(RR,zz);
          break;
    }

  if(fabs(zz)<4*Z0)
    {
      deltaz=(6.0*Z0)/NSHEETS;

      dphiz=0;

      for(i=0;i<NSHEETS;i++)
        {
          zpos=-3.0*Z0 + (i+0.5)*Z0*6.0/NSHEETS;

          R=RR;
          z=zz-zpos;

          Sigma0=(M_GAS)/(2*PI*Rd*Rd) * deltaz/(2*Z0) * pow(2/(exp(zpos/Z0)+exp(-zpos/Z0)),2);

          in1=qromb(intz_g,0,2/Rd);

          bb=2;
          do
            {
              in2=qromb(intz_g,bb/Rd,(bb+2)/Rd);
              in3=qromb(intz_g_abs,bb/Rd,(bb+2)/Rd);
              in1+=in2;
              bb+=2;
            }
          while(fabs(in3/in1)>1e-2);

          dphiz += 2*PI*G*Sigma0*Rd*Rd*( in1 );
        }
      return dphiz;
    }
  else
    {
      R=RR;
      z=zz;
      Sigma0=(M_GAS)/(2*PI*Rd*Rd);

      in1=qromb(intz_g,0,2/Rd);

      bb=2;
      do
        {
          in2=qromb(intz_g,bb/Rd,(bb+2)/Rd);
          in3=qromb(intz_g_abs,bb/Rd,(bb+2)/Rd);
          in1+=in2;
          bb+=2;
        }
      while(fabs(in3/in1)>1e-2);

      dphiz = 2*PI*G*Sigma0*Rd*Rd*( in1 );

      return dphiz;

    }
}




double intz_g(double k)
{
  double Rd=0.0;

  switch(GasDistribution) {
        case 0:
          Rd= H;
          break;
        case 1:
          Rd= H*GasExpAlpha;
          break;
        case 2:
          break;
    }

  if(z>0)
    return ( bessj0(k*R)*k*exp(-z*k)/pow(1+k*k*Rd*Rd,1.5));
  else
    return (-bessj0(k*R)*k*exp(z*k)/pow(1+k*k*Rd*Rd,1.5));
}


double intz_g_abs(double k)
{
  double Rd=0.0;

  switch(GasDistribution) {
        case 0:
          Rd= H;
          break;
        case 1:
          Rd= H*GasExpAlpha;
          break;
        case 2:
          break;
    }

  if(z>0)
    return fabs( bessj0(k*R)*k*exp(-z*k)/pow(1+k*k*Rd*Rd,1.5));
  else
    return fabs(-bessj0(k*R)*k*exp(z*k)/pow(1+k*k*Rd*Rd,1.5));
}




/* -------------------------------------
     Compute:  dphi/dR
   ------------------------------------- */

double comp_Dphi_R_gas(double RR,double zz)
{ 
  double comp_Dphi_R_gas_sph(double RR,double zz);
  double comp_Dphi_R_gas_exact(double RR,double zz);
  
  if(RR>0)
    {
      if(sqrt(RR*RR+zz*zz)>10*H)
        return comp_Dphi_R_gas_sph(RR,zz);
      else
        return comp_Dphi_R_gas_exact(RR,zz);
    }
  else
    return 0;
}




double comp_Dphi_R_gas_sph(double RR,double zz)
{ 
  double m;
  double r;
  
  r=sqrt(RR*RR+zz*zz);
  
  m=mass_cumulative_gas(r);
  
  return G*RR/(r*r*r)*m;
}


double comp_Dphi_R_gas_exact(double RR,double zz)
{ 
  int i; 
  double dphiR;
  double Sigma0;
  double in1,in2,in3,bb;
  double deltaz,zpos;
  double Rd=0.0;


  if(N_GAS==0) return 0;

  switch(GasDistribution) {
        case 0:
          Rd= H;
          break;
        case 1:
          Rd= H*GasExpAlpha;
          break;
        case 2:
          return comp_Dphi_R_gas_sph(RR,zz);
          break;
    }


  if(fabs(zz)<4*Z0)
    {
      deltaz=(6.0*Z0)/NSHEETS;
      dphiR=0;

      for(i=0;i<NSHEETS;i++)
        {
          zpos=-3.0*Z0 + (i+0.5)*Z0*6.0/NSHEETS;
          R=RR;
          z=zz-zpos;

          Sigma0=(M_GAS)/(2*PI*Rd*Rd) * deltaz/(2*Z0) * pow(2/(exp(zpos/Z0)+exp(-zpos/Z0)),2);

          in1=qromb(intR_g,0,2/Rd);
          bb=2;

          do
          {
            in2=qromb(intR_g,bb/Rd,(bb+2)/Rd);
            in3=qromb(intR_g_abs,bb/Rd,(bb+2)/Rd);
            in1+=in2;
            bb+=2;
          }
          while(fabs(in3/in1)>1e-2);

          dphiR += 2*PI*G*Sigma0*Rd*Rd*( in1 );
        }

      return dphiR;
    }
  else
    {
      R=RR;
      z=zz;

      Sigma0=(M_GAS)/(2*PI*Rd*Rd);

      in1=qromb(intR_g,0,2/Rd);
      bb=2;

      do
      {

        in2=qromb(intR_g,bb/Rd,(bb+2)/Rd);
        in3=qromb(intR_g_abs,bb/Rd,(bb+2)/Rd);
        in1+=in2;
        bb+=2;
       }
       while(fabs(in3/in1)>1e-2);

       dphiR = 2*PI*G*Sigma0*Rd*Rd*( in1 );

       return dphiR;
  }
}




double comp_Dphi_R_gas_razorthin(double RR,double zz)
{
  double Sigma0,y;
  double dphidR;
  double Rd=0;

  if(N_GAS==0) return 0;

  switch(GasDistribution) {
        case 0:
          Rd= H;
          break;
        case 1:
          Rd= H*GasExpAlpha;
          break;
        case 2:
          return comp_Dphi_R_gas_sph(RR,zz);
          break;
    }


  if(RR>0)
    {

      Sigma0=(M_GAS)/(2*PI*Rd*Rd);
      y=RR/(2*Rd);

      if(y>1e-4)
        dphidR = 2*PI*G*Sigma0*y*(bessi0(y)*bessk0(y)-bessi1(y)*bessk1(y));
      else
        dphidR =0;

      return dphidR;
    }
  else
    return 0;
}


double intR_g(double k)
{
  double Rd=0.0;

  switch(GasDistribution) {
        case 0:
          Rd= H;
          break;
        case 1:
          Rd= H*GasExpAlpha;
          break;
        case 2:
          break;
  }

  if(z>=0)
    return bessj1(k*R)*k*exp(-z*k)/pow(1+k*k*Rd*Rd,1.5);
  else
    return bessj1(k*R)*k*exp( z*k)/pow(1+k*k*Rd*Rd,1.5);
}

double intR_g_abs(double k)
{
  double Rd=0.0;

  switch(GasDistribution) {
        case 0:
          Rd= H;
          break;
        case 1:
          Rd= H*GasExpAlpha;
          break;
        case 2:
          break;
  }

  if(z>=0)
    return fabs(bessj1(k*R)*k*exp(-z*k)/pow(1+k*k*Rd*Rd,1.5));
  else
    return fabs(bessj1(k*R)*k*exp( z*k)/pow(1+k*k*Rd*Rd,1.5));
}


double gas_q_to_R(double q)
{
  int iter= 0;
  double R=1.0, Rold, Rd= 0.0, f, f_, pw;
  double Rnew;

  switch(GasDistribution) {
        case 0:
          Rd= H;
          break;
        case 1:
          Rd= H*GasExpAlpha;
          break;
        case 2:
          break;
  }

  // XXX Must include truncation
  switch(GasDistribution) {
          case 0:   /* Exponential Distribution (Mass is distributed
		       as R*exp(-R) due to surface element.) */
          case 1:

                do
                  {
                    f=(1+R)*exp(-R)+q-1;
                    f_=-R*exp(-R);

                    Rold=R;
                    R=R-f/f_;
                    if(iter > 200)
                      fprintf(stderr, "R=%g q=%g rold=%g iter=%d\n", R, q, Rold, iter);
                    if(iter > 210)
                      break;
                    iter++;
                  }
                while((fabs(R-Rold)/R> 1e-6) || R<0);
                R*=Rd;

                break;
          case 2:  /* Power-law distribution - gamma=1 is Mestel, or 1/R distribution */
                pw= 1/(2-PowerLawGamma);
                R=pow(q,pw);
                R*=PowerLawCutoff;                  /* nomalized by PowerLawCutoff */
                break;
        }

   return R;

}







