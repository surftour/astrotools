;; function to invert the "bias_of_Mhalo" function and give the halo mass 
;;   at a given bias at redshift z
function mhalo_of_bias, bias, REDSHIFT, NO_INTERPOL=NO_INTERPOL
	COMMON COSMOLOGICAL_PARAMETERS
	forward_function bias_of_Mhalo

	;; have to interpol across a grid; luckily only a narrow range in mass is relevant
	Mmin = 11.0
	Mmax = 15.0
	dM   = 0.25
	;; option to de-activate the more sophisticated re-interpolation algorithm, 
	;;   as it can run away in some cases
	if (keyword_set(NO_INTERPOL)) then begin
		Mmin = 10.0
		Mmax = 16.0
		dM   = 0.1
	endif
	M_grid = Mmin + dM*findgen((Mmax-Mmin)/dM+1.)
	b_grid = 0.0*M_grid
	for i = 0, n_elements(M_grid)-1 do b_grid[i] = bias_of_Mhalo(M_grid[i],REDSHIFT)
	log_Mhalo = INTERPOL(M_grid,b_grid,bias)
	
	if (keyword_set(NO_INTERPOL) EQ 0) then begin
	while ((log_Mhalo LT Mmin) OR (log_Mhalo GT Mmax)) do begin
		print, ' Halo mass outside of initial range, iterating...'
		Mmin = log_Mhalo - 2.
		Mmax = log_Mhalo + 2.
		M_grid = Mmin + dM*findgen((Mmax-Mmin)/dM+1.)
		b_grid = 0.0*M_grid
		for i = 0, n_elements(M_grid)-1 do b_grid[i] = bias_of_Mhalo(M_grid[i],REDSHIFT)
		log_Mhalo = INTERPOL(M_grid,b_grid,bias)
	endwhile
	endif
	
	return, log_Mhalo
end
