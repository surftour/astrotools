;; returns the dimensionless comoving volume element (need to multiply by hubble 
;;   distance^3 to get units)
;;
;; comoving volume (dV/dz/d_Omega)
;;
function comoving_volume_element, z
	forward_function hubble_E_z, comoving_distance, angular_diameter_distance
	return, ((1.+z)^2) * (angular_diameter_distance(z)^2) / hubble_E_z(z)
end

