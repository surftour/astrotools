function atten_fac, log_NH, BB=BB, SX=SX, HX=HX
	;; gives f in L_band_obs = f * L_band_intrinsic for a given NH & band
	forward_function return_idl_routines_homedir

	if (keyword_set(BB)) then begin
		NH = 16. + 0.1*findgen(101)
		;sigma = cross_section(3.0d8/4.4d-7)
		sigma = 7.4763487d-22	;; answer from above (just to save time)
		f = -sigma * 10^(DOUBLE(NH)) * alog10(exp(1.))
	endif
	if (keyword_set(SX)) then begin
		OPENR,1,return_idl_routines_homedir(0)+'/agn_spectrum/sx_atten_fac.dat'
		N_NH = 0
		READF,1,N_NH
			NH = fltarr(N_NH)
			f  = fltarr(N_NH)
		READF,1,NH
		READF,1,f
		CLOSE,1
	endif
	if (keyword_set(HX)) then begin
		OPENR,1,return_idl_routines_homedir(0)+'/agn_spectrum/hx_atten_fac.dat'
		N_NH = 0
		READF,1,N_NH
			NH = fltarr(N_NH)
			f  = fltarr(N_NH)
		READF,1,NH
		READF,1,f
		CLOSE,1
	endif

	y = 10^(DOUBLE(INTERPOL(f,NH,log_NH)))
	return, y	
end
