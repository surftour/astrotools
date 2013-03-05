;; returns omega_radiation as a function of redshift
function omega_radiation_z, z
	COMMON COSMOLOGICAL_PARAMETERS
	forward_function hubble_E_z
	return, OMEGA_RADIATION * ((1.+z)^4) / ((hubble_E_z(z))^2)
end
