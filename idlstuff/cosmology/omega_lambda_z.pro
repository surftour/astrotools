;; returns omega_lambda as a function of redshift
function omega_lambda_z, z
	COMMON COSMOLOGICAL_PARAMETERS
	forward_function hubble_E_z
	return, OMEGA_LAMBDA / ((hubble_E_z(z))^2)
end
