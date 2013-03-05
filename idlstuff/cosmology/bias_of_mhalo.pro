;; function to return the bias of a halo of mass Mhalo (in Msun/h) at redshift z
function bias_of_Mhalo, log_M_halo, REDSHIFT
	COMMON COSMOLOGICAL_PARAMETERS
	forward_function collapse_threshold, sigma_of_M

	;; return analytical Sheth-Tormen bias for a given halo mass at some redshift
;	a = 0.75
;	p = 0.3
;	delta_c = collapse_threshold(REDSHIFT)
;	sigma_M = sigma_of_Mhalo_z0(log_M_halo) * growth_factor(REDSHIFT)
;	nu = delta_c/sigma_M
;	b = 1. + (a*nu*nu - 1.)/delta_c + (2.*p)/(delta_c*(1. + (a*nu*nu)^p))

;; version from Croom et al. clustering paper -- Sheth, Mo, & Tormen (2001)
	a = 0.707
	b = 0.5
	c = 0.6
	delta_c = 1.69
	delta_c = 0.15*((12.*!PI)^(2./3.))*((OMEGA_MATTER_Z(REDSHIFT))^(0.0055))	
	;; slight correction from NFW '97
					;;	another factor of D(z) here would be double-counting, 
					;;    since it is included in sigma_M
	sigma_M = sigma_of_Mhalo_z0(log_M_halo) * growth_factor(REDSHIFT)


	nu = delta_c/sigma_M
	;nu = 1.
	denom_1 = SQRT(a) * delta_c
	numer_1 = a*SQRT(a) * nu*nu
	numer_2 = SQRT(a)*b*((a*nu*nu)^(1.-c))
	numer_3 = -(a*nu*nu)^c
	denom_3 = ((a*nu*nu)^c) + b*(1.-c)*(1.-c/2.)
		numer_4 = numer_3 / denom_3	
	bias = 1. + (numer_1 + numer_2 + numer_4)/denom_1	
	return, bias 
end
