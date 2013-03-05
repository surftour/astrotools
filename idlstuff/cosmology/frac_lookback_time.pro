;; returns the fractional lookback time (i.e. lookback time in units of 
;;   H0, multiply by 1/H0 to get the actual time back to a given point, 
;;	 or take 1/H0*(1-frac_lookback_time) to get the age of the universe at z
;;

function frac_lookback_time, z
	return, QROMO('dt_dz',0.0,z) / QROMO('dt_dz',0.0,1000.)
end
