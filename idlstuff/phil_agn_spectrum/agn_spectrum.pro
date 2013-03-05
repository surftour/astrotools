;;
;; function to return the AGN spectrum
;;   given nu (in Hz), and log(L_bol/L_sun), returns 
;;   log(nuL_nu/L_sun)
;;
;;	keywords to shortcut for various bands:
;;	  BB = B-band, IR = mid-IR (15 microns), 
;;    SX = soft X-rays (0.5-2 keV), 
;;    HX = hard X-rays (2-10 keV)
;;
;;  keywords for different spectra:
;;    MARCONI = Marconi et al. 2004-ish spectrum -- it's not *really* 
;;      the same (updated w. more recent estimates of the X-ray 
;;      input parameters, i.e. Gamma & reflection angle, from Tozzi et al. 2006 
;;		as well as the more recent Steffen et al. 2006 calibration of the 
;;      alpha_ox relation (specifically their l_uv-l_2keV bisector measurement, 
;;      which is a subtle difference but is very important for the bolometric 
;;      corrections of the most luminous X-ray sources
;;    HRH = Hopkins, Richards, & Hernquist 2006 template spectrum
;;	    compiled from a number of observations therein. It's a little closer to 
;;      a typical observed spectrum, with a more detailed modeling of a number of 
;;      parts of the spectrum (i.e. more continuum shape). the major difference is 
;;      that it includes a hot dust component, which is important to account for 
;;      since much of the energy comes out there. however, if you want to model 
;;      that generation self-consistently, you should use the "input" marconi et al. 
;;      spectrum. still, the galaxy-scale obscuration doesn't necessarily produce this 
;;      and the alpha-ox calibrations are for objects with this feature, so it is 
;;      generally more representative
;;    RICHARDS = Richards et al. 2006 all-quasar mean SED: if you want a 
;;      non-luminosity-dependent spectrum for whatever reason, this is an improved 
;;      version of Elvis et al. 1994 (which is too X-ray bright)
;;
;;    SDSS can be added as a keyword with any of the model spectra, and it will 
;;      overlay the vanden Berk et al. 2001 median SDSS SED over the relevant 
;;      portion of the spectrum :: basically does so with a sliding continuum determination
;;      that then overlays this spectrum, such that the integrated luminosity in the 
;;      entire bandpass and the continuum spectrum are conserved (i.e. you can put this 
;;      over an arbitrary continuum
;;
;;	unless you're using the Richards et al. 2006 spectra, if you want a reference for this, 
;;    it should be Hopkins, Richards, & Hernquist 2006. Even the "MARCONI" key spectrum 
;;    is substantially modified as described. If you want the additional relevant 
;;    observational compilations on which it's all based, the list is : 
;;    Richards et al. 2006, Steffen et al. 2006, Hatziminaoglou et al. 2005, 
;;    Tozzi et al. 2006, Strateva et al. 2005, Telfer et al. 2002, Vanden Berk et al. 2001, 
;;    George et al. 1998, Elvis et al. 1994, Marconi et al. 2004 (not really an observational 
;;    paper but the methodology does follow them), with the X-ray reflection component 
;;    following Ueda et al. 2003 in the PEXRAV code, Magdziarz & Zdziarski 1995
;;
;;
;;
function agn_spectrum, nu_in_Hz, log_l_bol, $
	BB=BB, IR=IR, SX=SX, HX=HX, $
	HRH=HRH, MARCONI=MARCONI, RICHARDS=RICHARDS, $
	SDSS=SDSS

	exec_call=return_idl_routines_homedir(0)+'/agn_spectrum/agn_spectrum.so'

	nu_in_Hz     = DOUBLE(nu_in_Hz)
		if (keyword_set(BB)) then nu_in_Hz = -1.0d0
		if (keyword_set(IR)) then nu_in_Hz = -2.0d0
		if (keyword_set(SX)) then nu_in_Hz = -3.0d0
		if (keyword_set(HX)) then nu_in_Hz = -4.0d0	
		N_nu  = LONG(n_elements(nu_in_Hz))

	log_l_bol    = DOUBLE(log_l_bol)
		N_lum = LONG(n_elements(log_l_bol))


	spectrum_key = 0L
		if (keyword_set(HRH)) 		then spectrum_key = 0L
		if (keyword_set(MARCONI)) 	then spectrum_key = 1L
		if (keyword_set(RICHARDS)) 	then spectrum_key = 2L	
	sloan_key    = 0L
		if (keyword_set(SDSS))		then sloan_key = 1L

	l_band_all = DOUBLE(fltarr(N_nu,N_lum))

	for i=0, N_lum-1 do begin
	l_band_vec = DOUBLE(fltarr(N_nu))
	log_l_bol_pass = DOUBLE(log_l_bol[i])
	S = CALL_EXTERNAL(exec_call, $
           'main',N_nu,nu_in_Hz,log_l_bol_pass,spectrum_key,sloan_key,l_band_vec)
	l_band_all[*,i] = l_band_vec
	endfor


;;	for i=0,n_lum-1 do begin
;;		plot,alog10(nu_in_Hz),l_band_all[*,i]-log_l_bol[i],ystyle=1,yrange=[-4.,0.]
;;	endfor


	return, l_band_all
end

