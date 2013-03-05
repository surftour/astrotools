function fitted_band_lum, log_L_bol, BB=BB, SX=SX, HX=HX, IR=IR, JACOBIAN=JACOBIAN
	;; returns the band luminosity from the 
	;;	 fitted bolometric corrections from "bol_cor_bands.pro" & 
	;;   everything herein (log_L_bol in solar luminosities)
	
	if (keyword_set(BB) OR keyword_set(IR)) then $
		  P = [8.99833   ,   6.24800  ,  -0.370587  , -0.0115970]
	if (keyword_set(SX)) then $
		  P = [10.0287   ,   17.8653  ,   0.276804  , -0.0199558]
	if (keyword_set(HX)) then $
		  P = [6.08087   ,   10.8331  ,   0.276802  , -0.0199597]


	x = log_L_bol - 10.
	L_bol_over_L_band = P[0]*(10^(x*P[3])) + P[1]*(10^(x*P[2]))
	L_band = 10^(DOUBLE(log_L_bol)) / DOUBLE(L_bol_over_L_band)
	y = alog10(L_band)
	if (keyword_set(IR)) then y = y + (-0.073656070)	;; B-band to 15 micron conversion

	;; JACOBIAN -- if set, return instead the jacobian (d_log_L / d_log_L_band)
	;;					which is straightforward to calculate from the equations above
	if (keyword_set(JACOBIAN)) then begin
		D1 = P[0]*(1.+P[3])*(10^(x*P[3])) + P[1]*(1.+P[2])*(10^(x*P[2]))
		D2 = P[0]*(10^(x*P[3])) + P[1]*(10^(x*P[2]))
		dlogL_over_dlogLband = D1/D2
		y = dlogL_over_dlogLband
	endif

	return, y	
end


;; (if use alpha-ox instead of the luv-lx bisector)
	;; HX : 
		;; P = [6.26204   ,   6.77995  ,   0.354465  , -0.0210125]
	;; SX : 
      	;; P = [10.3274   ,   11.1812   ,  0.354466  , -0.0210084]
	;; BB : 
      	;; P = [9.79932   ,   11.0396  ,  -0.530405  , -0.0198206]	

