;; compute the critical overdensity for collapse vs. redshift z
function collapse_threshold, z
	COMMON COSMOLOGICAL_PARAMETERS
	forward_function growth_factor
	delta_c_0 = 0.15*((12.*!PI)^(2./3.))*((OMEGA_MATTER_Z(z))^(0.0055))	
	return, delta_c_0/growth_factor(z)
end

