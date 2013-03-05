function growth_factor_deriv, a
	COMMON COSMOLOGICAL_PARAMETERS
	numerator = a^(3./2.)
	denominator = (Omega_LAMBDA*a*a*a + Omega_MATTER)^(3./2.)
	return, numerator/denominator
end
