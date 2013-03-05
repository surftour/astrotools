;; return the Hubble time in yr (to give units to the various 
;;   time measures, etc)
;; 
function hubble_time_yr, dummy
	COMMON COSMOLOGICAL_PARAMETERS
	return, 9.7813d9 / HUBBLE
end
