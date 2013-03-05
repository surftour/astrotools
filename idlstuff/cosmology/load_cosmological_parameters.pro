COMMON COSMOLOGICAL_PARAMETERS, $
	Omega_LAMBDA, Omega_MATTER, Omega_BARYON, HUBBLE, theta_cmb, $
	spectral_index, sigma_8, z_matter_rad_eq, Omega_RADIATION, $
	tau_CMB
	
	
function load_cosmological_parameters, dummy, $
	SET_Omega_LAMBDA=SET_Omega_LAMBDA, SET_OMEGA_MATTER=SET_OMEGA_MATTER, $
	SET_OMEGA_BARYON=SET_OMEGA_BARYON, SET_HUBBLE=SET_HUBBLE, $
	SET_THETA_CMB=SET_THETA_CMB, SET_SPECTRAL_INDEX=SET_SPECTRAL_INDEX, $
	SET_SIGMA_8=SET_SIGMA_8,SET_OMEGA_RADIATION=SET_OMEGA_RADIATION, $
	SET_TAU_CMB=SET_TAU_CMB, SET_Z_MATTER_RAD_EQ=SET_Z_MATTER_RAD_EQ

	COMMON COSMOLOGICAL_PARAMETERS

	;; default to WMAP-ish values
	Omega_Lambda = 0.73
		if (keyword_set(SET_OMEGA_LAMBDA)) then Omega_Lambda=SET_OMEGA_LAMBDA
	Omega_Matter = 0.27
		if (keyword_set(SET_OMEGA_MATTER)) then Omega_Matter=SET_OMEGA_MATTER
	Omega_Baryon = 0.04
		if (keyword_set(SET_OMEGA_BARYON)) then Omega_Baryon=SET_OMEGA_BARYON
	Hubble = 0.71
		if (keyword_set(SET_HUBBLE)) then Hubble=SET_HUBBLE
	theta_cmb = 1.0
		if (keyword_set(SET_theta_cmb)) then theta_cmb=SET_THETA_CMB
	spectral_index = 1.0
		if (keyword_set(SET_SPECTRAL_INDEX)) then spectral_index=SET_SPECTRAL_INDEX
	sigma_8 = 0.84
		if (keyword_set(SET_SIGMA_8)) then sigma_8=SET_SIGMA_8
	z_matter_rad_eq = 2.50d4 * Omega_Matter * Hubble*Hubble
		if (keyword_set(SET_Z_MATTER_RAD_EQ)) then z_matter_rad_eq=SET_z_matter_rad_eq
	Omega_Radiation = Omega_Matter/z_matter_rad_eq
		if (keyword_set(SET_OMEGA_RADIATION)) then Omega_Radiation=SET_OMEGA_RADIATION
	tau_CMB = 0.09
		if (keyword_set(SET_tau_CMB)) then tau_CMB=SET_tau_CMB

	return, 1
end
