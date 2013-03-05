;; return the dimensionless lookback time to some redshift
;; 
function lookback_time, z
	forward_function dt_dz
	return, QROMB('dt_dz',0.,z)	
end
