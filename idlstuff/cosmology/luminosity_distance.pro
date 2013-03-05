;; returns the dimensionless luminosity distance (need to multiply by hubble 
;;   distance to get units)
;;
function luminosity_distance, z
	forward_function comoving_distance
	return, (1.+z) * comoving_distance(z)
end
