;; returns the dimensionless Hubble scaling E(z)
function hubble_E_z, z
	COMMON COSMOLOGICAL_PARAMETERS
	return, SQRT(OMEGA_MATTER*(1.+z)^3 + OMEGA_LAMBDA + OMEGA_RADIATION*(1.+z)^4)
end