;; return the dimensionless time element per redshift element
function dt_dz, z
	forward_function hubble_E_z
	return, 1./((1.+z) * hubble_E_z(z))
end
