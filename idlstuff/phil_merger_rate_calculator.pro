;;
;; An IDL function to return the galaxy (or halo) merger rate, 
;;   as a function of various parameters including e.g.
;;   redshift, mass (stellar, baryonic, or halo), mass ratio, 
;;   and gas-richness of the galaxies. All results are based on 
;;   the 'default' semi-empirical or HOD-based model presented in 
;;   Hopkins et al., ApJ, 2009 (arXiv:0906.5357). 
;; 
;;   The function can output either: 
;;      (1) the merger rate per galaxy (i.e. N_mergers per galaxy per Gyr)
;;		(2) the total number density of mergers (N_mergers per Mpc^3 per Gyr)
;; 		(3) the expected merger fraction in *pair-selected* samples
;;		(4) the expected merger fraction in *morphology-selected* samples
;;
;;	 Options exist to change the input stellar mass functions used in the 
;;     HOD-based approach (used to populate galaxies into halos), 
;;     or to change the input cosmological parameters
;;
;;
;;	We stress that the code here makes heavy use of various fitting functions 
;;    and simplifications in order to be self-contained and fast. As a result, it 
;;    is not exactly equivalent to the full calculations in Hopkins et al. 2009. 
;;    However, we have checked that these various fitting functions, or 
;;    eliminated steps in the calculation, all change the output by much less than 
;;    the uncertainties discussed in the paper. 
;;
function merger_rate_calculator, REDSHIFT, LOG_M_MIN, LOG_M_MAX, $
	MASSRATIO_MIN, MASSRATIO_MAX, GASFRACTION_MIN, GASFRACTION_MAX, $
	RETURN_N_MERGERS_PER_GALAXY=RETURN_N_MERGERS_PER_GALAXY, $
	RETURN_N_MERGERS_PER_VOLUME=RETURN_N_MERGERS_PER_VOLUME, $
	RETURN_MERGER_FRACTION_PAIRS=RETURN_MERGER_FRACTION_PAIRS, $
	RETURN_MERGER_FRACTION_MORPHOLOGY=RETURN_MERGER_FRACTION_MORPHOLOGY, $
	USE_STELLAR_MASSES=USE_STELLAR_MASSES, $
	USE_BARYONIC_MASSES=USE_BARYONIC_MASSES, $
	USE_HALO_MASSES=USE_HALO_MASSES, $
	USE_ALTERNATIVE_MASSFUNCTION=USE_ALTERNATIVE_MASSFUNCTION, $
	SET_SIGMA_8=SET_SIGMA_8, $
	QUIET=QUIET

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;		INPUTS 											  ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  REQUIRED :: 
;;
;;	 (1) REDSHIFT		
;;		The redshift at which to evaluate the merger rate
;;	 (2) LOG_M_MIN		
;;		Minimum galaxy (or halo, if the appropriate options are set) mass for 
;;		which to calculate the merger rate. Units are log(M/M_sun), so i.e. 
;;		a desired mass of 10^9 M_sun means LOG_M_MIN=9.0
;;	 (3) LOG_M_MAX
;;		Maximum mass, in the same manner. 
;;		The output merger rates will be the rates for galaxies in the mass 
;;		interval: LOG_M_MIN <= log(M_gal/M_sun) <= LOG_M_MAX
;;	 (4) MASSRATIO_MIN
;;		Minimum mass ratio to calculate the merger rate down to. 
;;		Mass ratio is defined as MASSRATIO = M_secondary/M_primary, 
;;      such that 0 < MASSRATIO < 1 (i.e. integrating down to 1:10 mergers 
;;		means setting MASSRATIO_MIN = 0.1). 'Major' mergers are usually 
;;   	defined down to a MASSRATIO_MIN ~ 0.2-0.3. 
;;	 (5) MASSRATIO_MAX
;;		Maximum mass ratio, in the same manner. If this is >1, 
;;		the value of 1 will be used. The code will output the merger rate 
;; 		for mass ratios in the interval: 
;;		MASSRATIO_MIN <= massratio <= MASSRATIO_MAX
;;		Generally, MASSRATIO_MAX=1 is the value of interest -- including 
;;		all mergers 'more major' than MASSRATIO_MIN (unless e.g. 
;;		minor mergers alone are desired). 
;;	 (6) GASFRACTION_MIN
;;		Minimum gas fraction f_gas of merging galaxies, for which 
;;		to calculate the merger rate. The gas fraction is 
;;		defined as f_gas = M_gas/(M_gas+M_stars), where M_gas is the cold 
;;		mass in a galactic disk, and M_stars is the stellar mass 
;;		(the gas fraction here is really the *net* gas fraction of the 
;;		merging pair). Clearly 0 <= f_gas <= 1. 
;;	 (7) GASFRACTION_MAX
;;		Maximum gas fraction to include. As for the MASSRATIO inputs, 
;;		the output includes all mergers where the merging pair 
;;		satisfies: GASFRACTION_MIN <= f_gas <= GASFRACTION_MAX
;;		To simply obtain all mergers irrespective of gas fraction, 
;;		set GASFRACTION_MIN=0 and GASFRACTION_MAX=1. 
;;		To isolate e.g. only gas-rich mergers, set some minimum 
;;		such as GASFRACTION_MIN=0.05, with GASFRACTION_MAX=1. 
;;
;;  OPTIONS :: 
;;
;;	RETURN_N_MERGERS_PER_GALAXY
;;		if set (RETURN_N_MERGERS_PER_GALAXY=1), then the output of the 
;;		code will be the number of mergers (in the appropriate mass, 
;;		mass ratio, and gas fraction range) at the given redshift, 
;;		*per galaxy*, per unit time. 
;;			units are Gyr^-1
;;		**this is the default code output, if no options here are set**
;;	RETURN_N_MERGERS_PER_VOLUME
;;		if set (=1), the code will output the total number of mergers 
;;		(in the given mass, mass ratio, and gas fraction range) 
;;		at the given redshift, *per unit volume* per unit time. 
;;		this is the number above, folded in with the abundance of 
;;		galaxies in the given range of mass and other properties.
;;			units are Mpc^-3 Gyr^-1 
;;
;;	** to turn either of these into an observable merger fraction or 
;;       number of mergers, the output must be multiplied by an 
;;		 appropriate 'observable merger time' **
;;  RETURN_MERGER_FRACTION_PAIRS
;;		if set(to any >0 value) returns the merger *fraction* 
;;		expected in a pair-selected sample, where the value
;;		x in RETURN_MERGER_FRACTION_PAIRS=x is the projected 
;;		maximum pair separation used, in h^-1 kpc
;;		(for a pair sample with separations <25 h^-1 kpc, for 
;;		 example, set RETURN_MERGER_FRACTION_PAIRS=25)
;;		this then returns the merger fraction for 
;;		galaxies in this mass range
;;  RETURN_MERGER_FRACTION_MORPHOLOGY
;;		as the above, set(=1) to return the pair fraction 
;;		for galaxies selected by morphological disturbances
;;
;;
;;	USE_STELLAR_MASSES
;;		if set(=1), all 'masses' referred to 
;;			(in particular the inputs LOG_M_MIN and LOG_M_MAX, 
;;			 and the masses that define merger mass ratios, i.e. 
;;			 MASSRATIO_MIN and MASSRATIO_MAX)
;;		will be galaxy *STELLAR* masses
;;	USE_BARYONIC_MASSES
;;		as the above, but if set(=1), all masses refer to BARYONIC masses
;;		** this is the default code option, if no 'USE_X_MASSES' is set**
;;	USE_HALO_MASSES
;;		as the above, but masses refer to HALO masses (i.e. halo-halo mergers)
;;
;;
;;	USE_ALTERNATIVE_MASSFUNCTION
;;		default code outputs use the redshift-dependent stellar 
;; 		mass functions from the compilation in 
;;		Perez-Gonzalez et al., ApJ, 2008, 675, 234 to determine 
;;		galaxy abundances and populate them in halos. 
;;		setting this option (=1) will switch to using the 
;;		observed mass functions in Marchesini et al. 2009 (arXiv:0811.1773)
;;		There are other possible choices of course, but we find 
;;		that these bracket the range of possibilities.
;;	SET_SIGMA_8
;;		sets the value of SIGMA_8 to whatever value this parameter is 
;;		set to. the default value is ~0.8. 
;;		This is the most important cosmological parameter affecting 
;;		the merger rates (other parameters make little or no difference). 
;;		But other parameters can be adjusted by manually altering their 
;;		values in the code below.
;;
;;  QUIET 
;;		set this (=1) to suppress both error messages and 
;;		the normal text block at the end of the script that 
;;		explains what is being output. do so with caution!
;;
;;    
;;   Example of code use: 
;;     Say I want to know the total rate of major (mu>1/3) mergers I should 
;;       observe at z=0.6, from galaxy stellar masses of 10^10 - 10^11 M_sun, 
;;       with any gas fraction. Then I input :: 
;;       IDL> n = merger_rate_calculator(0.6,10.0,11.0,0.33,1.00,0.0,1.0,/RETURN_N_MERGERS_PER_VOLUME,/USE_STELLAR_MASSES)
;;	   Now I want to do the same, but only for gas-rich mergers, which I 
;;		 estimate as f_gas > 0.2, say. Then just :: 
;;       IDL> n = merger_rate_calculator(0.6,10.0,11.0,0.33,1.00,0.2,1.0,/RETURN_N_MERGERS_PER_VOLUME,/USE_STELLAR_MASSES)
;;     Ok, now lets divide out the number of galaxies -- number of mergers 
;;		 *per galaxy* expected :: 
;;       IDL> n = merger_rate_calculator(0.6,10.0,11.0,0.33,1.00,0.2,1.0,/RETURN_N_MERGERS_PER_GALAXY,/USE_STELLAR_MASSES)
;;			or (since mergers per galaxy is the default, just)
;;       IDL> n = merger_rate_calculator(0.6,10.0,11.0,0.33,1.00,0.2,1.0,/USE_STELLAR_MASSES)
;;	   That's the rate... but say I've got a pair-selected sample, 
;;		 with separations <25 h^-1 kpc (projected).. lets output the expected 
;;		 pair fraction at these separations ::
;;       IDL> n = merger_rate_calculator(0.6,10.0,11.0,0.33,1.00,0.2,1.0,RETURN_MERGER_FRACTION_PAIRS=25.0,/USE_STELLAR_MASSES)
;;
;;     OK, now say that I'm looking at a morphological sample. First off, since I'm 
;;		measuring the perturbations to galaxies, the mass ratio that probably 
;;		matters more is the total baryonic (not just stellar) mass ratio... 
;;      so I can get my per-galaxy rate with :: 
;;       IDL> n = merger_rate_calculator(0.6,10.0,11.0,0.33,1.00,0.2,1.0,/RETURN_N_MERGERS_PER_GALAXY,/USE_BARYONIC_MASSES)
;;			or, since both the 'per galaxy' and 'baryonic' flags are defaults, just :: 
;;       IDL> n = merger_rate_calculator(0.6,10.0,11.0,0.33,1.00,0.2,1.0)
;;    and if I want the merger fraction in my sample, assuming a (rough) 
;;     1 Gyr observable timescale, just :: 
;;       IDL> n = merger_rate_calculator(0.6,10.0,11.0,0.33,1.00,0.2,1.0,/RETURN_MERGER_FRACTION_MORPHOLOGY)
;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
DEBUG_MODE = 0

	;; error block
	REDSHIFT=REDSHIFT(0)
	if keyword_set(QUIET) eq 0 then begin
	if keyword_set(USE_HALO_MASSES) then begin
		if ((LOG_M_MIN lt 9.) or (LOG_M_MIN ge 16.) or $
			(LOG_M_MAX le 9.) or (LOG_M_MAX gt 16.)) then begin
			print, ' For Halos, require : '
			print, '   9 <= LOG_M_MIN < 16, 9 < LOG_M_MAX <= 16 '
			return, 0			
		endif	
	endif else begin
		if ((LOG_M_MIN lt 8.) or (LOG_M_MIN ge 13.) or $
			(LOG_M_MAX le 8.) or (LOG_M_MAX gt 13.)) then begin
			print, ' For Galaxies, require : '
			print, '   8 <= LOG_M_MIN < 13, 8 < LOG_M_MAX <= 13 '
			return, 0			
		endif	
	endelse
	if (LOG_M_MIN ge LOG_M_MAX) then begin
		print, ' LOG_M_MIN must be less than LOG_M_MAX'
		return, 0
	endif
	if (REDSHIFT lt 0. or REDSHIFT gt 6.) then begin
		print, ' REDSHIFT must be in range 0 <= z <= 6'
		return, 0
	endif
	if (MASSRATIO_MIN lt 0.01) then begin
		print, ' MASSRATIO_MIN must be >= 0.01'
		return, 0
	endif
	if (MASSRATIO_MAX gt 1.0) then begin
		print, ' MASSRATIO_MAX must be <= 1'
		return, 0
	endif
	if (MASSRATIO_MIN ge MASSRATIO_MAX) then begin
		print, ' MASSRATIO_MIN must be less than MASSRATIO_MAX'
		return, 0
	endif
	if (GASFRACTION_MIN lt 0.0) then begin
		print, ' GASFRACTION_MIN must be >= 0'
		return, 0
	endif
	if (GASFRACTION_MAX gt 1.0) then begin
		print, ' GASFRACTION_MAX must be <= 1'
		return, 0
	endif
	if (GASFRACTION_MIN ge GASFRACTION_MAX) then begin
		print, ' GASFRACTION_MIN must be less than GASFRACTION_MAX'
		return, 0
	endif
	if keyword_set(SET_SIGMA_8) then begin
	if (SET_SIGMA_8 le 0.6 or SET_SIGMA_8 ge 1.0) then begin
		print, ' Values of sigma_8 in SET_SIGMA_8 must be in the range 0.6 < sigma_8 < 1.0'
		return, 0
	endif
	endif
	endif

EXCEPT_STORE=!EXCEPT
!EXCEPT=0


;; define cosmological parameters, default to WMAP5
OMEGA_LAMBDA=0.726
OMEGA_MATTER=0.274
HUBBLE=0.705
SPECTRAL_INDEX=0.960
OMEGA_BARYON=0.0456
SIGMA_8=0.812
if keyword_set(SET_SIGMA_8) then SIGMA_8=SET_SIGMA_8
z = REDSHIFT

	
;; evaluate quantities of interest at the given redshift
DMASS_MIN_0=0.5
DMASS_MIN=DMASS_MIN_0
DMASS_TMP=LOG_M_MAX-LOG_M_MIN
if DMASS_TMP lt DMASS_MIN then DMASS_MIN=DMASS_TMP
r_corr_dmass=DMASS_MIN/DMASS_MIN_0

;; sheth-tormen halo mass function
	dlog_mhalo_grid = 0.01
	dlog_mhalo_grid = 0.0125*r_corr_dmass
	log_m_halo_grid = 6.+dlog_mhalo_grid*findgen((17.-6.)/dlog_mhalo_grid+1.)
	;; inputs: cosmological mass density, growth factor, collapse threshold, 
	;;				fluctuation amplitudes, etc.
    rho_m = 2.76135d11*HUBBLE*HUBBLE*OMEGA_MATTER
	omega_matter_z = OMEGA_MATTER*((1.+z)^3.)/(OMEGA_MATTER*(1.+z)^3 + OMEGA_LAMBDA)
	omega_lambda_z = OMEGA_LAMBDA/(OMEGA_MATTER*(1.+z)^3 + OMEGA_LAMBDA)
		;; approximation from Carroll et al. 1992 for growth factor
		omz=omega_matter_z 
		olz=omega_lambda_z
		gz  = (5./2.)*omz/(omz^(4./7.) - olz + (1.+omz/2.)*(1.+olz/70.))
		omz=omega_matter
		olz=omega_lambda
		g0  = (5./2.)*omz/(omz^(4./7.) - olz + (1.+omz/2.)*(1.+olz/70.))
		Dz  = (gz/g0)/(1.+z)
	growth_factor = Dz
	d_c = 0.15*((12.*!PI)^(2./3.))*((omega_matter_z)^(0.0055))/growth_factor
      mi = 2. + 0.025*findgen((17.-2.)/0.025+1.)
	;; need fluctuation amplitude SIGMA_OF_MHALO[halo mass,z=0]
	;; fitting formulae from van den Bosch 2002 (mnras, 331, 98), 
	;;   modified with a tiny tweak to accommadate a non-unity spectral index, calibrated 
	;;   against both the WMAP1, WMAP3, WMAP5 and concordance cosmologies
	;;   -- gives s(M) to within a factor 0.005 error, for 10^6 < M < 10^16 
	log_m_halo=mi
	gamma = omega_matter*hubble*exp(-omega_baryon * (1.+sqrt(2.*hubble)/omega_matter))
	c = 3.804d-4 & x = (c*gamma*((10^(log_m_halo/3.))))/((omega_matter)^(1./3.))
	g1 = 64.087 * ((1. + 1.074*(x^(0.3)) - 1.581*(x^(0.4)) + 0.954*(x^(0.5)) - 0.185*(x^(0.6)))^(-10.))
	x = (32.*gamma)
	g2 = 64.087 * ((1. + 1.074*(x^(0.3)) - 1.581*(x^(0.4)) + 0.954*(x^(0.5)) - 0.185*(x^(0.6)))^(-10.))
	f = (g1*g1)/(g2*g2)
	s = SQRT(f * sigma_8 * sigma_8)
	s_of_m=s*(10^((log_m_halo-14.09)*(1.00-spectral_index)/9.2))/(1.+(1.00-spectral_index)/9.2)
	;; sheth-tormen fit parameters
		a = 0.707
		b = 0.3222
		p = 0.3
	  si = s_of_m
      sig = interpol(si,mi,log_m_halo_grid)
      ln_mi = mi * alog(10.)
      ln_si = alog(si)
      dlns_dlnM = -deriv(ln_mi,ln_si)
      dlns_dlnM = interpol(dlns_dlnM,mi,log_m_halo_grid)
	x = sqrt(a) * d_c/sig
	f_x = b*sqrt(2./!PI)*(1.+x^(-2.*p))*(x)*exp(-x*x/2.)
	M_dn_dM = (rho_m/(10^(log_m_halo_grid)/HUBBLE)) * (dlns_dlnM) * f_x	
	dn_dlogM_HALOS = M_dn_dM * alog(10.)

	if DEBUG_MODE then $
	plot,log_m_halo_grid,dn_dlogM_HALOS,thick=4.01,color=0,/ylog,$
		xstyle=1,ystyle=1,xrange=[8.,15.],yrange=[1.0d-10,1.]

;; galaxy stellar mass function
	;; Perez-Gonzalez compilation
 	z0=[  0.1,   0.3,   0.5,   0.7,   0.9,   1.1,  1.45,   1.8,  2.25,  2.75,  3.25,  3.75]
 	a0=[-1.18, -1.19, -1.22, -1.26, -1.23, -1.26, -1.29, -1.27, -1.26, -1.20, -1.14, -1.23]
 	m0=[11.16, 11.20, 11.26, 11.25, 11.27, 11.31, 11.34, 11.40, 11.46, 11.34, 11.33, 11.36]-0.15
 	p0=[-2.47, -2.65, -2.76, -2.82, -2.91, -3.06, -3.27, -3.49, -3.69, -3.64, -3.74, -3.94]
 	;; interpolate through their best-fit parameters as a function of redshift. 
 	;;		correct for the stellar IMF (Chabrier adopted here, as in manuscript)
	 phi_ast=10^interpol(p0,z0,z)
	 alpha_ast=interpol(a0,z0,z)
	 m_ast=interpol(m0,z0,z)
	if keyword_set(USE_ALTERNATIVE_MASSFUNCTION) then begin
	;; marchesini et al. quoted fits (use these instead of conroy+wechsler -- 
	;;   the former under-predicts the z>2 massive galaxy population significantly
	   z0=[  0.1,  1.65,  2.50,  3.50]
	   a0=[-1.18, -0.99, -1.01, -1.39]
	   m0=[10.96, 10.91, 10.96, 11.38] 
	   p0=alog10([30.87,10.17,3.95,0.53]*1.0d-4)
		 phi_ast=10^interpol(p0,z0,z)
		 alpha_ast=interpol(a0,z0,z)
		 m_ast=interpol(m0,z0,z)
	endif
	d_mgal_grid = 0.01
	d_mgal_grid = 0.025*r_corr_dmass
	log_m_gal_grid = 8.+d_mgal_grid*findgen((13.-8.)/d_mgal_grid+1.)
	x = 10^(log_m_gal_grid - m_ast)
  	phi = phi_ast * (x^alpha_ast) * exp(-x) *(x*alog(10.))
	dn_dlogM_GALAXY_STELLAR = phi

	if DEBUG_MODE then $
	plot,log_m_gal_grid,dn_dlogM_GALAXY_STELLAR,thick=4.01,color=0,/ylog,$
		xstyle=1,ystyle=1,xrange=[8.,13.],yrange=[1.0d-10,1.]

;; use abundance matching to link halos+subhalos to galaxies, determine Mgal(Mhalo)
	n_gal =reverse(total(reverse(d_mgal_grid*dn_dlogM_GALAXY_STELLAR),/cumulative))
  	n_halo=reverse(total(reverse(dlog_mhalo_grid*dn_dlogM_HALOS),/cumulative))
    mh_mg=interpol(log_m_halo_grid,alog10(n_halo),alog10(n_gal))
    mg_mh=interpol(log_m_gal_grid,alog10(n_gal),alog10(n_halo)) 
 	  mg_mh_STELLAR = mg_mh

;; assign median gas masses, use to determine mg_mh for galaxy baryonic mass
;; (not folding over distribution here b/c it doesnt change the statistics)
	 mg_tmp = mg_mh
	 ;; fit to the fractional lookback time to redshift z
	 p=[-1.495,1.292,-0.3913,-1.830]
	 frac_lookback_time = 1.-(1.+z)^(p(0))-exp(-p(1)*z^(p(2)))*(1.+z)^(p(3))
	 fg_z0_tmp = 1./(1. + 10^(0.4*(mg_tmp-9.15)))
	 fg_z_tmp=fg_z0_tmp*(1.-frac_lookback_time*(1.-fg_z0_tmp^1.5))^(-0.67)
		fg_z_tmp=fg_z_tmp*(fg_z_tmp gt 0.)*(fg_z_tmp lt 0.99) + 0.99*(fg_z_tmp ge 0.99)
		mgas_grid_tmp=(fg_z_tmp/(1.-fg_z_tmp)) * 10^(mg_tmp)
		mg_tmp = alog10( 10^mg_tmp + mgas_grid_tmp )
		mg_mh = mg_tmp
		  mg_mh_BARYONIC = mg_mh

	if DEBUG_MODE then begin
	plot,log_m_halo_grid,mg_mh_STELLAR,thick=4.01,color=0,$
		xstyle=1,ystyle=1,xrange=[9.,16.],yrange=[8.,13.]
	oplot,log_m_halo_grid,mg_mh_BARYONIC,thick=4.01,color=80
	endif

	mg_mh = mg_mh_BARYONIC ;; default
	if keyword_set(USE_STELLAR_MASSES) then mg_mh=mg_mh_STELLAR
	if keyword_set(USE_BARYONIC_MASSES) then mg_mh=mg_mh_BARYONIC
	if keyword_set(USE_HALO_MASSES) then mg_mh=log_m_halo_grid
		baddies=where((mg_mh lt 0.) or (mg_mh gt 16.) or (finite(mg_mh) eq 0) or $
					  (finite(mg_mh,/nan) eq 1),n_baddies)
		if (n_baddies gt 0) then begin
		  mg_mh[baddies] = -1.0
		  ok=where(mg_mh ne -1.0)
		  mg_mh_ok=mg_mh[ok]
		  mh_ok=log_m_halo_grid[ok]
		  mg_mh=interpol(mg_mh_ok,mh_ok,log_m_halo_grid)
		endif

;; now comes the real step -- integrating to determine the merger rates
	mh=log_m_halo_grid & dmh=dlog_mhalo_grid & n_mh=dn_dlogM_HALOS
	mgal_min=LOG_M_MIN
	mgal_max=LOG_M_MAX
	massratio_cut=MASSRATIO_MIN
	t_Hubble_z0 = 9.7813/HUBBLE
	dt_dz = t_Hubble_z0/((1.+z) * SQRT(OMEGA_MATTER*(1.+z)^3 + OMEGA_LAMBDA))
	mgr_z_norm=2.89d-2*(0.925+0.18*(1.-exp(-z^(0.75)/0.8)))
	mgr_t_norm=mgr_z_norm/dt_dz
	GAS_FRAC_SCATTER = 0.15 
		;; in dex -- looking at observations, at the masses of interest, this is quite a good 
		;;    approximation to the actual scatter. 
		;; assuming f_gas is lognormal with some median f0 and 
		;;   lognormal scatter s, then the fraction with f>fx 
		;;   is simply   0.5*( 1. + erf[-(fx-f0)/(sqrt[2]*s)] )
		i_min=round(interpol(findgen(n_elements(mh)),mg_mh,mgal_min))
		i_max=round(interpol(findgen(n_elements(mh)),mg_mh,mgal_max))
			if (mgal_min le min(mg_mh)) then i_min=0
			if (mgal_max ge max(mg_mh)) then i_max=n_elements(mh)-1
		if DEBUG_MODE then begin
			print, i_min, i_max, mgal_min, mgal_max
			plot,mg_mh,findgen(n_elements(mh))
			bad=where(finite(mg_mh) eq 0 or finite(mg_mh,/nan) eq 1,n_bad)
			print, 'BAD POINTS = ',n_bad
		endif
		 if (i_min ge n_elements(mh)-1 or i_min le 0 or $
		 	    finite(i_min) eq 0 or finite(i_min,/nan) eq 1) then i_min=0
		 if (i_max ge n_elements(mh)-1 or i_max le 0 or $
		  		finite(i_max) eq 0 or finite(i_max,/nan) eq 1) then i_max=n_elements(mh)-1
		n_gal_tot=0. & n_mgr_tot=0. & n_mgr_wfg_tot=0.
		for j=i_min,i_max do begin
 			xsi=10^(log_m_halo_grid - log_m_halo_grid[j])
 			  mg_tmp = mg_mh_STELLAR[j] 
				 fg_z0_tmp = 1./(1. + 10^(0.4*(mg_tmp-9.15)))
				 fg_z_tmp=fg_z0_tmp*(1.-frac_lookback_time*(1.-fg_z0_tmp^1.5))^(-0.67)
				 fg_z_tmp=fg_z_tmp*(fg_z_tmp gt 0.)*(fg_z_tmp lt 0.99) + 0.99*(fg_z_tmp ge 0.99)
 				 fg_exp_0=fg_z_tmp(0) 			
    		d_log_xsi=dmh & d_ln_xsi=d_log_xsi*alog(10.) & d_xsi=d_ln_xsi*xsi
    		;; now use fitting functions for halo merger rates (tweaked to more reflect model 
    		;;	after convolving with appropriate distributon functions)
    		xsi_bar=0.098 & beta=-2.01 & gamma=0.409
    		dp_dxsi=(xsi^(beta))*exp((xsi/xsi_bar)^gamma)*(xsi le 1.) 
     			bad=where((dp_dxsi le 0.) or (finite(dp_dxsi) eq 0) or (finite(dp_dxsi,/nan) eq 1),n_bad)
     			if (n_bad gt 0) then dp_dxsi[bad]=10^(-40.0) & dp=dp_dxsi*d_xsi
    		    p_xsi_int=reverse(total(reverse(dp),/cumulative))*mgr_t_norm
    		     f=max(p_xsi_int,imax) & g=min(p_xsi_int,imin)
    		    xsi_gal = 10^(mg_mh - mg_mh[j])
    		     p_above_massratio_min = $
    		    	10^interpol(alog10(p_xsi_int),alog10(xsi_gal),alog10(massratio_cut)) 
    		    	p_above_massratio_min=p_above_massratio_min(0)
    		    	p_merger=p_above_massratio_min
    		    if (MASSRATIO_MAX lt 1.) then begin
    		     p_above_massratio_max = $
    		    	10^interpol(alog10(p_xsi_int),alog10(xsi_gal),alog10(MASSRATIO_MAX)) 
    		    	p_above_massratio_max=p_above_massratio_max(0)
    		    	p_merger=p_above_massratio_min-p_above_massratio_max
    		    endif
			n_gal_j=dmh*n_mh[j]
			n_gal_tot = n_gal_tot + n_gal_j
			n_mgr_tot = n_mgr_tot + n_gal_j*p_merger
			
			;; now correct to fraction in the appropriate range of gas fractions
			if (GASFRACTION_MIN le 0.01) then begin
				p_fg_min=1.0
			endif else begin
				p_fg_min=0.5*(1.+erf(-(GASFRACTION_MIN-fg_exp_0)/(sqrt(2.)*GAS_FRAC_SCATTER)))
			endelse
			if (GASFRACTION_MAX ge 0.99) then begin
				p_fg_max=0.0
			endif else begin
				p_fg_max=0.5*(1.+erf(-(GASFRACTION_MAX-fg_exp_0)/(sqrt(2.)*GAS_FRAC_SCATTER)))
			endelse
			p_fg = p_fg_min(0) - p_fg_max(0) 
			if (p_fg le 0.) then p_fg=0.
			n_mgr_wfg_tot = n_mgr_wfg_tot + n_gal_j*p_merger*p_fg
		endfor
		n_mgr_total  = n_mgr_wfg_tot
		n_mgr_pergal = n_mgr_wfg_tot/n_gal_tot
		if DEBUG_MODE then begin
			print, n_gal_tot
			print, n_mgr_total
			print, n_mgr_pergal
		endif

		if keyword_set(QUIET) eq 0 then begin
		s1='Returning merger rate *per galaxy*, in Gyr^-1'
		if KEYWORD_SET(RETURN_N_MERGERS_PER_GALAXY) then $
			s1='Returning merger rate *per galaxy*, in Gyr^-1'
		if KEYWORD_SET(RETURN_N_MERGERS_PER_VOLUME) then $
			s1='Returning merger rate *per unit volume*, in Mpc^-3 Gyr^-1'
		if KEYWORD_SET(RETURN_MERGER_FRACTION_PAIRS) then $
			s1='Returning merger *fraction*, for pairs with projected separations <= '+$
				string(RETURN_MERGER_FRACTION_PAIRS)+' h^-1 kpc'
		if KEYWORD_SET(RETURN_MERGER_FRACTION_MORPHOLOGY) then $
			s1='Returning merger *fraction*, for morphology-selected systems (assumed t_obs~1 Gyr)'

		print, s1
		print, '  at redshift z = ',REDSHIFT
		print, '  for galaxies with masses in the range: '
		print, '     ',LOG_M_MIN,' <= mass <= ',LOG_M_MAX
		MASSKEY = 'BARYONIC mass'
			if keyword_set(USE_STELLAR_MASSES) then MASSKEY='STELLAR mass'
			if keyword_set(USE_BARYONIC_MASSES) then MASSKEY='BARYONIC mass'
			if keyword_set(USE_HALO_MASSES) then MASSKEY='HALO mass'
		print, '     in terms of their '+MASSKEY
		MASSKEY = 'BARYONIC-BARYONIC'
			if keyword_set(USE_STELLAR_MASSES) then MASSKEY='STELLAR-STELLAR'
			if keyword_set(USE_BARYONIC_MASSES) then MASSKEY='BARYONIC-BARYONIC'
			if keyword_set(USE_HALO_MASSES) then MASSKEY='HALO-HALO'
		print, '  where a merger is defined by a '+MASSKEY+' mass ratio in the range: '
		print, '     ',MASSRATIO_MIN,' <= M_secondary/M_primary <= ',MASSRATIO_MAX
		print, '  and galaxy gas fractions in the range: '
		print, '     ',GASFRACTION_MIN,' <= f_gas(merger) <= ',GASFRACTION_MAX
		if keyword_set(USE_ALTERNATIVE_MASSFUNCTION) then begin
		print, '  (using alternative observed stellar mass function from '
		print, '    Marchesini et al. 2009 to populate the HOD, instead of '
		print, '    default mass function from Perez-Gonzalez et al. 2008)'
		endif
		if keyword_set(SET_SIGMA_8) then begin
		print, '  (switching from default WMAP5 value to SIGMA_8 = ',SET_SIGMA_8,' )'
		endif
		endif

		;!EXCEPT=EXCEPT_STORE

		if KEYWORD_SET(RETURN_N_MERGERS_PER_GALAXY) then $
			return, n_mgr_pergal
		if KEYWORD_SET(RETURN_N_MERGERS_PER_VOLUME) then $
			return, n_mgr_total
		;; merger observable timescales calibrated from high-resolution 
		;;	simulations and mock imaging in Lotz et al., MNRAS, 2008, 391, 1137
		;;  (note that as in many places above, we're simply using mean 
		;;    values here to correct, statistically, for what is a broad 
		;;    distribution that in the full model is folded over)
		if KEYWORD_SET(RETURN_MERGER_FRACTION_PAIRS) then begin
			f = n_mgr_pergal * 0.35 * (RETURN_MERGER_FRACTION_PAIRS/30.)
			return, 1.-exp(-f) 
			;; this is necessary if f becomes large (multiple mergers per~t_Hubble), 
			;;	to correct for the fact that multiple mergers only appear as one -- 
			;;  i.e. merger fraction is, by definition, < 1
		endif
		if KEYWORD_SET(RETURN_MERGER_FRACTION_MORPHOLOGY) then begin
			f = n_mgr_pergal * 1.0  
			return, 1.-exp(-f)
		endif
		return, n_mgr_pergal  ;; default (same as RETURN_N_MERGERS_PER_GALAXY) 

	return, 0
end

