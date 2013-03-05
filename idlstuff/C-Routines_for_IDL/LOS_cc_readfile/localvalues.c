#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "overhead.h"
#include "proto.h"

/* This is the greatly-stripped version of this code used for the 
	quick tree calculations of column density */


/* Uses the two-phase model of SH03 to get the respective properties of the 
 *   two phases of a two-phase medium based on given density, specific energy,
 *   and electron fraction (ne_given)
 *   (the inputs should be in INTERNAL UNITS)
 */
void get_two_phase_breakdown(double u_tot, double rho_tot, double ne_given, LOCALVAL *local)
{
  int i, flag, stars_spawned, tot_spawned, stars_converted, tot_converted;
  unsigned int gen, bits;
  double dt, dtime, ascale = 1, hubble_a = 0, a3inv, ne = 1;
  double time_hubble_a, unew;
  double sum_sm, total_sm, sm, rate, sum_mass_stars, total_sum_mass_stars;
  double p, prob;
  double cloudmass;
  double factorEVP;
  double tsfr, trelax;
  double egyhot, egyeff, egycurrent, tcool, x, y, rate_in_msunperyear;
  double sfrrate, totsfrrate, dmax1, dmax2;

  double temperature_hot,egy_spec_hot,egy_spec_cold,density_hot,density_cold;
  double filling_factor_cold,entropy,intrinsic_density_hot,intrinsic_density_cold;
  double mass_fraction_cold,density_ratio,vff_factor,xc0,xh0,sfr_per_unit_mass;
  double egyhot_current,density_tot,u_cold,ne_temp,twophase_density_threshold;
  int f;
  ne_temp = ne_given;
  
  /* check whether we're in the region with thermal instability
   *  
   * f=1  normal cooling
   * f=0  star formation
   */
   if (rho_tot < All.PhysDensThresh) {

		f = 1;
   		temperature_hot         = convert_u_to_temp(u_tot, rho_tot, &ne_temp);
   		entropy                 = GAMMA_MINUS1*u_tot*(pow(rho_tot,GAMMA_MINUS1));
   		egy_spec_hot            = u_tot;
   		egy_spec_cold           = 0.0;
   		intrinsic_density_hot   = rho_tot;
   		intrinsic_density_cold  = 0.0;
   		mass_fraction_cold      = 0.0;
   		filling_factor_cold     = 0.0;
   		sfr_per_unit_mass		= 0.0;
   		x						= 0.0;
   		xh0 					= 1.0;
   		xc0						= 0.0;

   } else {
   
   		f = 0;

		ne = ne_given;
		density_tot = rho_tot;
		a3inv = 1.0;
		u_cold = All.EgySpecCold;   
   		entropy = GAMMA_MINUS1*u_tot*(pow(rho_tot,GAMMA_MINUS1));

		/* Calculate the star formation timescale, evaporation coefficient, etc 
		 * following SH03 -- that gives the equilibrium specific energy of the hot component
		 */
		tsfr = sqrt(All.PhysDensThresh / (density_tot * a3inv)) * All.MaxSfrTimescale;
		factorEVP = pow(density_tot * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

		egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold; 
		
		/* Use u_hot and cooling time to get fraction of gas in the cold phase (x) */  
		tcool = GetCoolingTime(egyhot, density_tot * a3inv, &ne);

		//printf(" %e %e %e %e %e %e %e %e \n",All.EgySpecCold,All.PhysDensThresh,All.MaxSfrTimescale,All.FactorEVP,All.EgySpecSN,a3inv,ne,ne);

		y =  tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN -
					   (1.0 - All.FactorSN) * All.EgySpecCold);

		x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

		egyeff = egyhot * (1 - x) + All.EgySpecCold * x;
		/* Recall that u is a mass-weighted average, so 
		 *   rho_hot*u_hot = rho*u_tot - rho_cold*u_cold
		 */
		egyhot_current  = (u_tot - x*u_cold) / (1.0 - x);


		/* Now convert this to the quantities above */
		egy_spec_hot   = egyhot_current;
		egy_spec_cold  = u_cold;
		mass_fraction_cold = x;
		
		/* Assume the two phases are in pressure equilibrium, so that
		 *   P_hot = P_cold  
		 *     -->  rho_hot_intrinsic * u_hot = rho_cold_intrinsic * u_cold
		 *     (we assume u_hot = u_hot_intrinsic, u_cold = u_cold_intrinsic)
		 */
		 density_ratio = egy_spec_cold / egy_spec_hot;
		 vff_factor    = (x / (1.0 - x)) * density_ratio;

		 filling_factor_cold = vff_factor / (1.0 + vff_factor);
		 xc0 = (1.0 + vff_factor) * x / vff_factor;		/* rho_cold_intrinsic/rho_tot  */
		 xh0 = density_ratio * xc0;						/* rho_hot_intrinsic/rho_tot   */

   		intrinsic_density_cold  = xc0 * rho_tot;
   		intrinsic_density_hot   = xh0 * rho_tot;

        temperature_hot = convert_u_to_temp(egy_spec_hot, intrinsic_density_hot, &ne_temp);

		sfr_per_unit_mass = (1 - All.FactorSN) * x / tsfr *
		  (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
		  /* (this multiplied by the total particle mass gives the total star 
		   *   formation rate in the particle)  */
		//printf("egy_spec_cold=%e x=%e egyhot_current=%e \n",egy_spec_cold,x,egyhot_current);
		//printf("density_ratio=%e y=%e egyeff=%e egyhot=%e \n",density_ratio,y,egyeff,egyeff);
		//printf("cold_massfrac=%e Vff_cold=%e X_COLD=%e X_HOT=%e \n",mass_fraction_cold,vff_factor,xc0,xh0);
   }

   /* Now dump everything we're interesting in keeping into a LOCALVALS struct */
   local->u_hot = egy_spec_hot;
   local->cloud_massfrac = x;
   local->cloud_fillingfactor = filling_factor_cold;
   local->rho_hot = xh0;
   local->rho_cold = xc0;

   return;
}

