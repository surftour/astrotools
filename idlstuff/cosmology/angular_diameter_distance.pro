;; returns the dimensionless angular diameter distance (need to multiply by hubble 
;;   distance to get units)
;;
function angular_diameter_distance, z
	forward_function comoving_distance
	return, comoving_distance(z) / (1. + z)
end
