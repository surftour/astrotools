function intrinsic_spectrum, log_nu, log_l_bol
	;; from the full below calculation, use the fitting formula:
	;    (L_bol/L_uv) = P0 + P1*(L_uv/10^12 L_sun)^P2
	forward_function richards_return_spectrum
	
	;; recalc if change x_cut	-- currently using x_cut ~ 50 angstroms
	P0 = 4.5176299
	P1 = 0.63426829
	P2 = -0.36684631

	;; redone using the more recent Steffen et al. calibration of 
	;;   the luv-lx relation
	P0 = 4.4348811
	P1 = 0.69836704
	P2 = -0.27939027



		l_bol = 10^(DOUBLE(log_l_bol))
		l_uv_grid = 6. + 0.1*findgen(101)
		bolcorr_grid= P0 + P1*(10^((l_uv_grid-12.)*P2))
			bol_grid = alog10(bolcorr_grid) + (DOUBLE(l_uv_grid))
			bolcorr = 10^(INTERPOL(alog10(bolcorr_grid),bol_grid,log_l_bol))
		l_uv = l_bol / bolcorr
			nu_uv = 3.0d8/ (2500. * 1.0d-10)	;; 2500 AA in Hz
			l_uv_nu = l_uv/nu_uv
	alpha_OX = -0.137*(alog10(l_uv_nu)+alog10(3.9)+33.) + 2.638	
		;; Steffen et al. determination
	;;print, l_uv, alpha_OX

	;; now that that's out of the way, in a position to just load up the rest 
	;;   of the spectrum
	
	;; load intrinsic (blue optical -- closest to an un-reddened spectrum?)
	f_lt_c = richards_return_spectrum(log_nu,/BLUE)
		;; cutoff the high end at 500 angstroms
			nu_c = 3.0d8 / (500. * 1.0d-10)
			L_uv_0 = 10^(INTERPOL(alog10(f_lt_c),log_nu,alog10(nu_uv)))
			f_lt_c = f_lt_c * (l_uv/l_uv_0[0])
			L_c    = 10^(INTERPOL(alog10(f_lt_c),log_nu,alog10(nu_c)))
			f_lt_c = f_lt_c * (log_nu LE alog10(nu_c))

	;; load generic X-ray spectrum
	f_gt_c = load_xr_spec(log_nu)
		keV = 2.418d17
		nu_cx = (3.0d8 / 50.0d-10)
		nu_2keV = 2.0*keV
		L_2keV_0 = 10^(INTERPOL(alog10(f_gt_c),log_nu,alog10(2.0*keV)))

	fac_temp = alog10(nu_uv/(2.0*keV))
	;; instead of this, use the lx-luv direct scaling (bisector calibrated -- 
	;;    end point is that it's slightly steeper, giving a slightly more 
	;;    well-behaved L_2kev at high luminosities
	log_l_2keV_grid_nu = (0.72) * (alog10(l_uv_nu)+alog10(3.9)+33.) + 4.53
		;;;log_l_2keV = -alpha_OX*fac_temp + alog10(l_uv_nu)
	l_2keV     = 10^(DOUBLE(log_l_2keV_grid_nu + alog10(2.0*keV) - alog10(3.9)-33.))
	f_gt_c = f_gt_c * (l_2keV/l_2keV_0[0])
		ok = where(log_nu GE alog10(nu_cx), n_ok)
		if (n_ok GT 0) then $
		L_cx   = 10^(INTERPOL(alog10(f_gt_c[ok]),log_nu[ok],alog10(nu_cx)))
		;print, L_cx
		;print, L_c

	beta = alog10((L_cx/nu_cx)/(L_c/nu_c)) / alog10(nu_cx/nu_c)
	f_gap = L_c * ((10^(DOUBLE(log_nu)) / nu_c)^(beta+1.))
	
	f = f_lt_c*(log_nu LT alog10(nu_c)) + $
		f_gap*((log_nu GE alog10(nu_c)) AND (log_nu LE (alog10(nu_cx)))) + $
		f_gt_c*(log_nu GT alog10(nu_cx))

	;plot,log_nu,alog10(f),xstyle=1,xrange=[12.,21.],ystyle=1, $
	;	yrange=[alog10(l_uv)-3.,alog10(l_uv)+1.]
	;	
	;	nu0 = 3.0d8/1.5d-6
	;	oplot,alog10([nu0,nu0]),[-20.,20.]
	;
	;	nu0 = 3.0d8/2.5d-6
	;	oplot,alog10([nu0,nu0]),[-20.,20.]

	return, alog10(f)
end



function intrinsic_spectrum_calc_luv_to_lbol_full
	;; the detailed calculation for a range of spectra to do alpha_ox & 
	;;  use it to do the inversion problem and solve for an equation of 
	;;  l_bol(l_uv) ultimately in a truly self-consistent manner

	d_log_nu = 0.001
	nu_min = 11.0
	nu_max = 22.0
	log_nu = nu_min + d_log_nu*findgen((nu_max-nu_min)/d_log_nu+1.)

	;; load intrinsic (blue optical -- closest to an un-reddened spectrum?)
	f = richards_return_spectrum(log_nu,/BLUE)
	;; cutoff the high end at 500 angstroms
	nu_c = 3.0d8 / (500. * 1.0d-10)
	nu_uv = 3.0d8/ (2500. * 1.0d-10)
	L_uv  = 10^(INTERPOL(alog10(f),log_nu,[alog10(nu_uv)]))
	L_c   = 10^(INTERPOL(alog10(f),log_nu,[alog10(nu_c)]))
	alpha_c = (L_c/nu_c)/(L_uv/nu_uv)
	
	f = f * (log_nu LE alog10(nu_c))
	L_lt_c = TOTAL(DOUBLE(f)*d_log_nu*alog(10.))
	print, L_lt_c
	plot,log_nu,alog10(f)

	print, l_uv[0], L_lt_c/l_uv
	


	;; load generic X-ray spectrum
	f_x = load_xr_spec(log_nu)
	keV = 2.418d17
	;; cutoff the low end at 50 angstroms
	nu_cx = (3.0d8 / 50.0d-10)
	nu_2keV = 2.0*keV
	L_gt_cx = TOTAL(DOUBLE(f_x)*d_log_nu*alog(10.))
	L_2keV = 10^(INTERPOL(alog10(f_x),log_nu,[alog10(2.0*keV)]))
		ok = where(log_nu GE alog10(nu_cx))
	L_cx   = 10^(INTERPOL(alog10(f_x[ok]),log_nu[ok],[alog10(nu_cx)]))
	alpha_cx = (L_cx/nu_cx)/(L_2keV/nu_2keV)
	print, l_2keV[0], L_gt_cx/L_2keV


	l_uv_grid = 6.0 + 0.1*findgen(101)
	l_uv_grid_nu = 10^(DOUBLE(l_uv_grid+alog10(3.9)+33.))/nu_uv


	alpha_OX = -0.137*alog10(l_uv_grid_nu) + 2.638	;; Steffen et al. determination
	fac_temp = alog10(nu_uv/(2.0*keV))
	log_l_2keV_grid_nu = -alpha_OX*fac_temp + alog10(l_uv_grid_nu)

	;; instead of this, use the lx-luv direct scaling (bisector calibrated -- 
	;;    end point is that it's slightly steeper, giving a slightly more 
	;;    well-behaved L_2kev at high luminosities
	log_l_2keV_grid_nu = (0.72) * alog10(l_uv_grid_nu) + 4.53


	l_2keV_grid = log_l_2keV_grid_nu + alog10(2.0*keV) - (alog10(3.9) + 33.)
	l_rat = 10^(l_2keV_grid - l_uv_grid)
	l_rat_0 = 10^INTERPOL(alog10(l_rat),l_uv_grid,11.0)
	print, l_rat_0, l_uv, l_2keV


	renorm_fac = (L_uv[0] * l_rat_0[0] / l_2keV[0])
	f_x = f_x * renorm_fac
	L_cx   = L_cx * renorm_fac
	beta   = alog10((L_cx/nu_cx)/(L_c[0]/nu_c))/alog10(nu_cx/nu_c)


	plot,log_nu,alog10(f),ystyle=1,yrange=[-2.5,-0.5],xstyle=1,xrange=[12.,21.]
	oplot,log_nu,alog10(f_x)
	;;oplot,log_nu,alog10(richards_return_spectrum(log_nu,/BLUE)),COLOR=255


	f_gap = l_c[0] * ((10^(DOUBLE(log_nu))/nu_c)^(beta[0]))
	f_gap = f_gap * (10^(log_nu)/nu_c)
	

	f = f*(log_nu LT alog10(nu_c)) + $
		f_gap*((log_nu GE alog10(nu_c)) AND (log_nu LE (alog10(nu_cx)))) + $
		f_x*(log_nu GT alog10(nu_cx))

	oplot,log_nu,alog10(f),COLOR=255

	renorm_fac_uv = 10^(l_uv_grid) / l_uv[0]
	L_lt_c = L_lt_c * renorm_fac_uv
	L_uv = l_uv[0] * renorm_fac_uv
	L_c  = L_c[0] * renorm_fac_uv

	renorm_fac_x = (L_uv * l_rat / l_2keV[0])
	L_cx   = L_cx[0] * renorm_fac_x
	beta   = alog10((L_cx/nu_cx)/(L_c/nu_c))/alog10(nu_cx/nu_c)
	L_gt_cx = L_gt_cx * renorm_fac_x

	L_gap = (1./(1.+beta)) * (nu_c/nu_uv) * alpha_c[0] * l_uv * ((nu_cx/nu_uv)^(beta+1.)-1.)
	L_tot = L_lt_c + L_gap + L_gt_cx	;; now have L_bol as a function of l_2500

	print, L_tot/L_uv
	bol_uv_corr = L_tot/L_uv
	plot,l_uv_grid,bol_uv_corr

	P_guess = [4.,1.,-3.];,-1.]
	x0 = l_uv_grid-12.
	P_fitted = MPFITFUN('cubic_fitfun',x0,bol_uv_corr,0.0*x0+1.,P_guess,/QUIET)
	print, P_fitted
	y0 = cubic_fitfun(x0,P_fitted)
	oplot,l_uv_grid,y0,COLOR=255

end


function cubic_fitfun, x, P
	;return, P[0]+P[1]*x+P[2]*x*x+P[3]*x*x*x
	return, P[0] + P[1]*(10^(x*P[2]))
end

