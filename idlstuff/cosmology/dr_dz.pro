;; return the dimensionless comoving distance element per redshift element
function dr_dz, z
	forward_function hubble_E_z
	return, 1./hubble_E_z(z)
end
