;; return the Hubble distance in Mpc (to give units to the various 
;;   distance measures, etc)
;; 
function hubble_distance_Mpc, dummy
	COMMON COSMOLOGICAL_PARAMETERS
	return, 3000./Hubble
end
