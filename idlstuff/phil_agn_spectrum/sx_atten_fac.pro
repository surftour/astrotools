function sx_atten_fac, log_NH
	keV=2.418d17
	d_nu = 0.01
	nu = (0.5 + 0.01*findgen(151))*keV

	f = attenuated_spectrum(alog10(nu),12.,log_NH)
		fnu = 10^(double(f)) / nu
	f0= attenuated_spectrum(alog10(nu),12.,0.)
		f0nu = 10^(double(f0)) / nu
	
	L = TOTAL(fnu)*d_nu
	L0= TOTAL(f0nu)*d_nu
	
	return, L/L0
end


pro compare_atten_fac
	forward_function fit_atten_fac
	
	NH = 16. + 0.1*findgen(101)
	f  = 0.0*NH
	for i=0,n_elements(NH)-1 do f[i] = sx_atten_fac(NH[i])
	plot,NH,10^alog10(f),/ylog ;,ystyle=1,yrange=[-10.1,0.1]
	sigma_eff = 5.0d-24
		oplot,NH,10^alog10(EXP(-(sigma_eff*10^(DOUBLE(NH)))^(0.8) )),COLOR=250

	P_guess = [-23.5, 0.8]
		P_fitted = MPFITFUN('fit_atten_fac',NH,alog10(f),0.0*NH+1.,P_guess,/QUIET)
		print, P_fitted
		y0 = fit_atten_fac(NH,P_fitted)
		oplot,NH,10^y0,COLOR=110,THICK=2.0,LINESTYLE=2
	
	bad = where((FINITE(f) EQ 0) OR (FINITE(f,/NAN) EQ 1), nbad)
	ok  = where((FINITE(f) NE 0) AND (FINITE(f,/NAN) NE 1), nok)
	if (nbad GT 0) then f[bad] = 10^(DOUBLE(INTERPOL(alog10(f[ok]),NH[ok],NH[bad])))
		
	OPENW,1,'sx_atten_fac.dat'
	PRINTF,1,n_elements(NH)
	PRINTF,1,NH
	PRINTF,1,alog10(f)
	CLOSE,1
	print, f

	oplot,NH,(atten_fac(NH,/SX)),LINESTYLE=2,THICK=2.,COLOR=150
	oplot,NH,(atten_fac(NH,/HX)),LINESTYLE=2,THICK=2.,COLOR=250
	oplot,NH,(atten_fac(NH,/BB)),LINESTYLE=2,THICK=2.,COLOR=80	
end

function fit_atten_fac, x, P
	return, alog10(EXP(-(10^(DOUBLE(P[0])) * (10^(DOUBLE(x))))^(P[1])))
end
