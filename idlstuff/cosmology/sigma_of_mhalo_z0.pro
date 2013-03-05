COMMON sigma_of_Mhalo_z0_TO_PASS, R_of_M_to_pass

;; function to return the local sigma(M) as a function of Mhalo (in Msun/h)
;;   !! designed to be done in h_inverse coordinates !! (but can switch to h's below)
function sigma_of_Mhalo_z0, log_M_halo
	COMMON COSMOLOGICAL_PARAMETERS
	COMMON sigma_of_Mhalo_z0_TO_PASS
	forward_function sigma_of_Mhalo_z0_deriv
		int_floor = -3.
		int_roof  = 3.
	;; assumes the passed mass is in SOLAR MASSES
	rho_crit = (2.78d11 * OMEGA_MATTER );* Hubble*Hubble) 		;; rho_crit in M_sun/Mpc^3

	;; normalize to sigma8
	R_of_M_to_pass = 8. ;/Hubble
		k_char = alog(1./R_of_M_to_pass)
		sigma_R_squared = QROMO('sigma_of_Mhalo_z0_deriv',k_char-3.,k_char+3.)
		sigma_M_8 = SQRT(sigma_R_squared)
	R_of_M_to_pass = ((3./(4.*!PI))*(10^(log_M_halo)/rho_crit))^(1./3.)
		k_char = alog(1./R_of_M_to_pass)
			if ((k_char-3.) LT int_floor) then int_floor = k_char-3.
			if ((k_char+3.) GT int_roof)  then int_roof  = k_char+3.
		sigma_R_squared = QROMO('sigma_of_Mhalo_z0_deriv',int_floor,int_roof)
		sigma_M = SQRT(sigma_R_squared) * (SIGMA_8 / sigma_M_8)
	return, sigma_M
end
