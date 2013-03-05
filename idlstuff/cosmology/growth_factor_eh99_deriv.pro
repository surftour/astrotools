function growth_factor_EH99_deriv, z
	COMMON COSMOLOGICAL_PARAMETERS
	g = SQRT(OMEGA_MATTER*(1.+z)^3 + (1.-OMEGA_MATTER-OMEGA_LAMBDA)*(1.+z)^2 + OMEGA_LAMDBDA)
	return, (1.+z)/g^3
end
