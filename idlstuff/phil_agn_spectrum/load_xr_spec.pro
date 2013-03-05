function load_xr_spec, log_nu
	;; loads XR spectrum for nu > 1 keV

;; unnecessary 20,000 point binning
	homedir = return_idl_routines_homedir(0)+'/agn_spectrum/'

	OPENR,7,homedir+'xr.spec.dat'
	READF,7,n_to_read
	nu=fltarr(n_to_read)
	lnu=fltarr(n_to_read)
	READF,7,nu
	READF,7,lnu
	CLOSE,7

	f = INTERPOL(lnu,nu,log_nu)
	keV = 2.418d17	;; keV in Hz
	f = 10^(DOUBLE(f))
	f = f*(log_nu GE MIN(nu))

	if (2 EQ 0) then begin
	nu0 = MIN(nu) + 0.05*findgen((MAX(nu)-MIN(nu))/0.05+1.)
	lnu0= INTERPOL(lnu,nu,nu0)
	OPENW,8,'xr.spec.downsampled.dat'
	PRINTF,8,n_elements(nu0)
	PRINTF,8,nu0
	PRINTF,8,lnu0
	CLOSE,8
	endif

	return, f	
end