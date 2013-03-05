pro bol_cor_bands
	forward_function attenuated_spectrum, cubic_fitfun
	L_bol_grid = 6.0  + 0.1*findgen(101)
		x0 = L_bol_grid-10.
	NH_grid    = 16.0 + 0.1*findgen(101)
	
	;; hard X-ray
	keV  = 2.418d17
	d_nu = 0.01
	nu   = 2.0 + d_nu*findgen((10.0-2.0)/d_nu+1.)
	nu   = nu*keV
	l_bol_over_l_B = fltarr(n_elements(L_bol_grid))
	for i=0, n_elements(L_bol_grid)-1 do begin
		f = attenuated_spectrum(alog10(nu),L_bol_grid[i],NH_grid[0])
		fnu = f - alog10(nu)
		f = alog10(TOTAL(10^(DOUBLE(fnu))) *d_nu*keV)
		l_bol_over_l_B[i] = 10^(L_bol_grid[i]-f)
	endfor

	plot,L_bol_grid,alog10(l_bol_over_l_B),/nodata
		oplot,L_bol_grid,alog10(l_bol_over_l_B),COLOR=250,THICK=2.0
	P_guess = [20.,10.,0.5,-0.01]
		P_fitted = MPFITFUN('cubic_fitfun',x0,alog10(l_bol_over_l_B),0.0*x0+1.,P_guess,/QUIET)
		print, P_fitted
		y0 = cubic_fitfun(x0,P_fitted)
		oplot,L_bol_grid,(y0),COLOR=250,THICK=2.0,LINESTYLE=2


	;; B band
	nu = [3.0d8/4.4d-7]	;; B-band at 4400 Angstroms
	l_bol_over_l_B = fltarr(n_elements(L_bol_grid))
	for i=0, n_elements(L_bol_grid)-1 do begin
		f = attenuated_spectrum(alog10(nu),L_bol_grid[i],NH_grid[0])
		l_bol_over_l_B[i] = 10^(L_bol_grid[i]-f)
	endfor

	oplot,L_bol_grid,alog10(l_bol_over_l_B),COLOR=80,LINESTYLE=0
	P_guess = [10.,10.,-0.5,-0.01]
		P_fitted = MPFITFUN('cubic_fitfun',x0,alog10(l_bol_over_l_B),0.0*x0+1.,P_guess,/QUIET)
		print, P_fitted
		y0 = cubic_fitfun(x0,P_fitted)
		oplot,L_bol_grid,(y0),COLOR=80,THICK=2.0,LINESTYLE=2



	;; soft X-ray
	keV  = 2.418d17
	d_nu = 0.01
	nu   = 0.5 + d_nu*findgen((2.0-0.5)/d_nu+1.)
	nu   = nu*keV
	l_bol_over_l_B = fltarr(n_elements(L_bol_grid))
	for i=0, n_elements(L_bol_grid)-1 do begin
		f = attenuated_spectrum(alog10(nu),L_bol_grid[i],NH_grid[0])
		fnu = f - alog10(nu)
		f = alog10(TOTAL(10^(DOUBLE(fnu))) *d_nu*keV)
		l_bol_over_l_B[i] = 10^(L_bol_grid[i]-f)
	endfor

	oplot,L_bol_grid,alog10(l_bol_over_l_B),COLOR=150,LINESTYLE=0
	P_guess = [20.,10.,0.5,-0.01]
		P_fitted = MPFITFUN('cubic_fitfun',x0,alog10(l_bol_over_l_B),0.0*x0+1.,P_guess,/QUIET)
		print, P_fitted
		y0 = cubic_fitfun(x0,P_fitted)
		oplot,L_bol_grid,(y0),COLOR=150,THICK=2.0,LINESTYLE=2

end



function cubic_fitfun, x, P
	;return, P[0]+P[1]*x+P[2]*x*x+P[3]*x*x*x
	return, alog10(P[0]*(10^(x*P[3])) + P[1]*(10^(x*P[2])))
end
