;; growth factor, using the EH99 normalization

function growth_factor_EH99, z
	COMMON COSMOLOGICAL_PARAMETERS
	forward_function d_growth_factor_EH99
	Omega_M = OMEGA_MATTER
	Omega_L = OMEGA_LAMBDA
	h = HUBBLE

	g = SQRT(Omega_M*(1.+z)^3 + (1.-Omega_M-Omega_L)*(1.+z)^2 + Omega_L)
	zeq = 2.50d4 * Omega_M *h*h* (theta_cmb^(-4))
	Omega_z = Omega_M*((1.+z)^3)/(g*g)
	Omega_L_z = Omega_L/(g*g)
	D1 = ((1.+zeq)/(1.+z))*(5.*Omega_z/2.)/((Omega_z^(4./7.)) - Omega_L_z $
				+ (1.+Omega_z/2.)*(1.+Omega_L_z/70.))
	return, D1
end
