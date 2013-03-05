;; returns the dimensionless comoving distance (need to multiply by hubble 
;;   distance to get units
;;
function comoving_distance, z
	forward_function hubble_E_z_inverse
	return, QROMB('hubble_E_z_inverse',0.,z)
end
