pro ir_band_lum
	forward_function attenuated_spectrum, cubic_fitfun
	L_bol_grid = 6.0  + 0.1*findgen(101)
	NH_grid    = 15.0 + 0.1*findgen(101)
	;;L_bol_grid = 10.0  + 1.0*findgen(8)
	
	
	;; 15 micron
	mu = 1.0d-6
	c  = 3.0d8
	nu = c/(15.*mu)

	l_bol_over_l_i = fltarr(n_elements(L_bol_grid),n_elements(NH_grid))
	for i=0, n_elements(L_bol_grid)-1 do begin
	print, 'i = ',i
	for j=0, n_elements(NH_grid)-1 do begin
		f = attenuated_spectrum(alog10(nu),L_bol_grid[i],NH_grid[j])
		l_bol_over_l_i[i,j] = 10^(L_bol_grid[i]-f)
	endfor
		epsilon = 0.01
		f = l_bol_over_l_i[i,*]
		plot,NH_grid,f
		df_dNH = DERIV(NH_grid,f)
		bad = where(df_dNH GE 0., n_bad)
		;;print, 'N_bad= ',n_bad
		if (n_bad GT 0) then begin
			bad_i = bad[0]-1
			f_bad_i = f[bad_i]
			f_gt_i  = f[bad_i+1:n_elements(f)-1]
			f_local_max = MAX(f_gt_i,max_i)
			f_gt_i  = f_gt_i[max_i:n_elements(f_gt_i)-1]
				bad_f = where(f_gt_i LT (f_bad_i-epsilon))
				bad_f = bad_i+1+max_i+bad_f[0]+1
				oplot,[NH_grid[bad_i],NH_grid[bad_i]],[-1.0d6,1.0d6]
				oplot,[NH_grid[bad_f],NH_grid[bad_f]],[-1.0d6,1.0d6]
			f_i = f[bad_i]
			f_f = f[bad_f]
			x_i = NH_grid[bad_i]
			x_f = NH_grid[bad_f]
			f_interpol = f_i + (f_f-f_i)/(x_f-x_i) * (NH_grid[bad_i+1:bad_f]-x_i)
			f[bad_i+1:bad_f] = f_interpol
			oplot,NH_grid,f,COLOR=250
			df_dNH = DERIV(NH_grid,f)
			bad = where(df_dNH GE 0., n_bad)
			print, 'N_bad (post-correction)= ',n_bad
			l_bol_over_l_i[i,*] = f
		endif
	endfor


	OPENW,1,'ir.15m.lum.dat'
	PRINTF,1,n_elements(L_bol_grid)
	PRINTF,1,L_bol_grid
	PRINTF,1,n_elements(NH_grid)
	PRINTF,1,NH_grid
	PRINTF,1,l_bol_over_l_i
	CLOSE,1

	plot,NH_grid,(l_bol_over_l_i[0,*]),/nodata,xstyle=1,xrange=[15.,26.]
		oplot,NH_grid,(l_bol_over_l_i[0,*]),COLOR=250,THICK=2.0

end
