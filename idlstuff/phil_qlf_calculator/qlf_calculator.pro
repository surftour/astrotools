;;
;; returns the QLF from HRH 2006, at some list of luminosities and some redshift, 
;;    for arbitrary band, and any of the fitting formulae given therein
;; 
;; see qlf_calculator_callable.c for a detailed description of what all these 
;;    keywords actually mean
;;
;;
function qlf_calculator, log_l_band_in_solar, redshift, $
	L_IS_L_BOL=L_IS_L_BOL, VERBOSE=VERBOSE, $
	BB=BB, IR=IR, SX=SX, HX=HX, BOL=BOL, NU=NU, LAMBDA=LAMBDA, $
	FULL=FULL, PLE=PLE, FAINT=FAINT, BRIGHT=BRIGHT, SCATTER=SCATTER, $
		SCHECHTER=SCHECHTER, LDDE=LDDE, PDE=PDE, FIT_NUM=FIT_NUM 

	exec_call=return_idl_routines_homedir(0)+'/luminosity_functions/qlf_calculator_callable.so'

	;; default to the bolometric QLF
	nu_in_Hz = 0.0
		if (keyword_set(BOL)) then nu_in_Hz = 0.0
		if (keyword_set(BB)) then nu_in_Hz = -1.0
		if (keyword_set(IR)) then nu_in_Hz = -2.0
		if (keyword_set(SX)) then nu_in_Hz = -3.0
		if (keyword_set(HX)) then nu_in_Hz = -4.0
		if (keyword_set(NU)) then nu_in_Hz = NU
		if (keyword_set(LAMBDA)) then nu_in_Hz = 2.998d8/LAMDA

	;; default to the best-fit total model
	fit_model = 0
		if (keyword_set(FULL))  	then fit_model = 0
		if (keyword_set(PLE))   	then fit_model = 1
		if (keyword_set(FAINT)) 	then fit_model = 2
		if (keyword_set(BRIGHT))	then fit_model = 3
		if (keyword_set(SCATTER))  	then fit_model = 4
		if (keyword_set(SCHECHTER))	then fit_model = 5
		if (keyword_set(LDDE))  	then fit_model = 6
		if (keyword_set(PDE))  		then fit_model = 7
		if (keyword_set(FIT_NUM))	then fit_model = FIT_NUM

	redshift     = DOUBLE(redshift)
	nu_in_Hz     = DOUBLE(nu_in_Hz)
	fit_model    = LONG(fit_model)
	
	n_grid = 100L
	l_bol_pass = dblarr(n_grid)
	l_bnd_pass = dblarr(n_grid)
	p_bnd_pass = dblarr(n_grid)

	vocal = 0L
	if (keyword_set(VERBOSE)) then vocal = 1L

	S = CALL_EXTERNAL(exec_call, $
           'main',nu_in_Hz,redshift,fit_model,l_bol_pass,l_bnd_pass,p_bnd_pass,vocal)
	
	l_grid = l_bnd_pass
	if (keyword_set(L_IS_L_BOL)) then l_grid = l_bol_pass
	l_grid = l_grid - alog10(3.9) - 33.0	;; converts to solar
	
	return, INTERPOL(p_bnd_pass,l_grid,log_l_band_in_solar)
end
