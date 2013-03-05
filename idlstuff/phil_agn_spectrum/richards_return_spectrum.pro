function richards_return_spectrum, log_nu, ALL=ALL,BLUE=BLUE,RED=RED,OPTLUM=OPTLUM,OPTDIM=OPTDIM, $
									IRLUM=IRLUM,IRDIM=IRDIM,BB=BB,SX=SX,HX=HX
	OPENR,1,return_idl_routines_homedir(0)+'/agn_spectrum/richards_bol_corr_tables.dat'
	x = '1 2'
	y = STRSPLIT(x,/EXTRACT)
	i_line = 1
	i_line_max = 227
	nu = fltarr(i_line_max-1)
	f1 = nu
	f2 = nu
	f3 = nu
	f4 = nu
	f5 = nu
	f6 = nu
	f7 = nu
	while(i_line LT i_line_max) do begin
		READF,1,x
		y = STRSPLIT(x,/EXTRACT)
		nu[i_line-1] = FLOAT(y[0])
		f1[i_line-1] = FLOAT(y[1])
		f2[i_line-1] = FLOAT(y[2])
		f3[i_line-1] = FLOAT(y[3])
		f4[i_line-1] = FLOAT(y[4])
		f5[i_line-1] = FLOAT(y[5])
		f6[i_line-1] = FLOAT(y[6])
		f7[i_line-1] = FLOAT(y[7])
		i_line = i_line + 1
	endwhile
	CLOSE,1
	nuLnu = f1
		if (KEYWORD_SET(ALL)) then 		nuLnu = f1
		if (KEYWORD_SET(BLUE)) then 	nuLnu = f2
		if (KEYWORD_SET(RED)) then 		nuLnu = f3
		if (KEYWORD_SET(OPTLUM)) then 	nuLnu = f4
		if (KEYWORD_SET(OPTDIM)) then 	nuLnu = f5
		if (KEYWORD_SET(IRLUM)) then 	nuLnu = f6
		if (KEYWORD_SET(IRDIM)) then 	nuLnu = f7
	new_nu = 11.0 + 0.01*findgen(801)
	nuLnu = INTERPOL(nuLnu,nu,new_nu)
	nu = new_nu
		d_log_nu = 0.01
		l_bol = alog(10.) * TOTAL(10^(DOUBLE(nuLnu))) * d_log_nu	

	nuLnu_of_log_nu = INTERPOL(nuLnu,nu,log_nu) - alog10(l_bol)
	
	if (KEYWORD_SET(BB)) then begin
		log_nu = [alog10(6.818d14)]
		nuLnu_of_log_nu = INTERPOL(nuLnu,nu,log_nu) - alog10(l_bol)
	endif 
	if (KEYWORD_SET(SX)) then begin
		keV = 2.418d17
		l_band = alog(10.) * TOTAL(10^(DOUBLE(nuLnu*((10^nu GE 0.5*keV) AND (10^nu LE 2.0*keV))))) * d_log_nu
		nuLnu_of_log_nu = alog10(l_band) - alog10(l_bol)
	endif 
	if (KEYWORD_SET(HX)) then begin
		keV = 2.418d17
		l_band = alog(10.) * TOTAL(10^(DOUBLE(nuLnu*((10^nu GE 2.0*keV) AND (10^nu LE 10.0*keV))))) * d_log_nu
		nuLnu_of_log_nu = alog10(l_band) - alog10(l_bol)
	endif 
	return, 10^(DOUBLE(nuLnu_of_log_nu))
end
