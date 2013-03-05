;; returns omega_matter as a function of redshift
function omega_matter_z, z
	COMMON COSMOLOGICAL_PARAMETERS
	forward_function hubble_E_z
	return, OMEGA_MATTER * ((1.+z)^3) / ((hubble_E_z(z))^2)
end
