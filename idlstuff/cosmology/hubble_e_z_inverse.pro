function hubble_E_z_inverse, z
	COMMON COSMOLOGICAL_PARAMETERS
	forward_function hubble_E_z
	return, 1./hubble_E_z(z)
end