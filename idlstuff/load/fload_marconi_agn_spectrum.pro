function fload_marconi_agn_spectrum, l_bol, redshift, n_spectrum, l_nu_min, l_nu_max
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;	mas.pro
;
;	Returns marconi spectrum for a given
;	bolometric AGN luminosity
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
log_nu_min = l_nu_min
log_nu_max = l_nu_max


c_light = 3.0e8
nu_uv =  alog10(c_light/(2.5e-7))
;print,"nu_uv",nu_uv

kev_lambda = c_light/(2.418e17)
kev_nu  = 2.418e+17 ;kev in freq
kev_erg = 1.602e-9 ;kev in ergs
h_planck = 6.625e-27

nu_two_kev = (2.0*2.418e17)


;if (log_nu_min gt nu_uv-1.0) then begin
;	log_nu_min = nu_uv - 1.0
;	;print,nu_uv,log_nu_min,log_nu_max
;endif
print,"In marconi_spectrum "
print,"L_bol=      ",l_bol
print,"Redshift=   ",redshift
print,"n_spectrum= ",n_spectrum
print,"log_nu_min= ",log_nu_min
print,"log_nu_max= ",log_nu_max





; define what we're looking for
; -------------------------------
nu = fltarr(n_spectrum)
nu = (log_nu_max-log_nu_min)*findgen(n_spectrum)/float(n_spectrum) + log_nu_min
;print, "nu=        ",nu

; wavelength

lambda = fltarr(n_spectrum)
lambda(*) = c_light/(10.0^(nu))
;print, "lambda=    ",lambda






;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make complete Marconi spectrum
;;;;;;;;;;;;;;;;;;;;;;;;;;;


; 1) make non-normalized spectrum

sL     = l_bol-12.0
l_B    = l_bol - (0.80 - 0.067*sL + 0.017*sL*sL - 0.0023*sL*sL*sL)
nu_B   = alog10(c_light/(4.4e-7))
l_B    = alog10((10.0^(l_B))/10.0^(nu_B))
;print,'L_B ',l_B
l_nu   = fltarr(n_spectrum)



;  > 1 micron
; --------------
; alpha = 2.0 for lambda > 1 micron
rj_index = where(lambda ge 1.0e-6, nrj)
if nrj gt 0 then l_nu(rj_index) = 2.0*( nu(rj_index) - nu(rj_index(nrj-1)))


;  1 um - 1300 Angstroms
; --------------------------
; alpha = -0.44 for  1 micron > lambda > 1300 Ang
ouv_index = where((lambda lt 1.0e-6) and (lambda ge 1.3e-7), nouv)
if nouv gt 0 then l_nu(ouv_index) = -0.44*( nu(ouv_index)-nu(ouv_index(0)))



; do luminosity correction
lind = where( abs(nu-nu_B) lt 0.001, nl)
if nl gt 0 then begin
	l_corr = l_nu(lind(0))-l_B 
	l_nu(*) = l_nu(*)-l_corr
endif


;  UV and X-ray corrections
; --------------------------
nu_uv =  alog10(c_light/(2.5e-7))
linuv = where( abs(nu-nu_uv) lt 0.001, nluv)
if nluv gt 0 then begin
	L_uv = l_nu(linuv(0))
	;print,"Lum at 2500 Ang = ",L_uv

	alpha_OX = 0.11*l_uv + 1.85
	print,"alpha_OX = ",alpha_OX
endif


;print, "L_B*nu_B = ",alog10((10.0^l_B)*(10.0^nu_B))
;print, "L_uv*nu_uv = ",alog10((10.0^l_uv)*(10.0^nu_uv))


;   1300 - 1200  Angstroms
;  --------------------------
; flat at the top
f_index = where((lambda lt 1.3e-7) and (lambda ge 1.2e-7), nf)
if nf gt 0 then l_nu(f_index) = l_nu(ouv_index(nouv-1))



;   1200 - 500  Angstroms
; ---------------------------
; alpha = -1.76 for  1200 > lambda > 500 Ang
fuv_index = where((lambda lt 1.2e-7) and (lambda ge 5.0e-8), nfuv)
if nfuv gt 0 then l_nu(fuv_index) = -1.76*( nu(fuv_index)-nu(fuv_index(0))) + l_nu(f_index(nf-1))




;    X-ray Emission
; ---------------------
;do_xray_emission= 0
do_xray_emission= 1
if log_nu_max lt kev_nu then do_xray_emission= 0
if do_xray_emission eq 1 then begin
	;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; Get X-ray spectrum
	; from pexrav model
	; ref: Magdziarz and Zdziarski, MNRAS, 273, 837, 1995
	;;;;;;;;;;;;;;;;;;;;;;;;;;;
		print, "Calc X-ray emission"

		Nen  = 500 ; number of energies
		XEar  = fltarr(Nen+1)
		XSpec = fltarr(Nen+1)
		emax = 2.0e5 ;keV
		emin = 1.0    ;keV
		XEar = 10.0^((alog10(emax)-alog10(emin))*findgen(Nen+1)/float(Nen-1) + alog10(emin))

	        S = CALL_EXTERNAL('/home/brant/code/idl/marconi_agn_spectrum/marconi_agn_spectrum', $
	        ;S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/AGN_Spectrum/marconi_agn_spectrum', $
	                'marconi_agn_spectrum', $
			Nen, $
			XSpec, $
			XEar, $
	                redshift,/F_VALUE)
		xspec = xspec/xear
	;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;print,"Reflection spectrum calculated"
	;;;;;;;;;;;;;;;;;;;;;;;;;;;

	; use pexrav for E > 1keV
	;kev_lambda = c_light/(2.418e17)
	;kev_nu  = 2.418e17 ;kev in freq
	;kev_erg = 1.602e-9 ;kev in ergs
	;h_planck = 6.625e-27

	;find x-ray lum at 2 keV 
	;nu_two_kev = (2.0*2.418e17)
	f_nu = nu_uv - alog10(nu_two_kev)
	L_two_kev = alpha_OX*f_nu +L_uv

	xray_index = where((lambda lt kev_lambda) and (lambda ge kev_lambda/emax), nx)
	if nx gt 0 then begin
		e_nu = fltarr(nx)
		e_nu = h_planck*(10.0^nu(xray_index))/kev_erg
		l_nu(xray_index) = alog10(interpol(xspec,xear,e_nu))

		two_kev_index = where( abs(nu(xray_index)-alog10(nu_two_kev)) lt 0.001,ntkv)
		l_corr = l_nu(xray_index(two_kev_index(0))) - L_two_kev
		l_nu(xray_index) = l_nu(xray_index) - l_corr
	endif
endif


;    UV Emission
;  ----------------
; connect fuv to xray with power-law
c_index = where( (lambda lt 5.0e-8) and (lambda ge kev_lambda), nc)
if nc gt 0 then l_nu(c_index) = (l_nu(xray_index(0))-l_nu(fuv_index(nfuv-1)))*(nu(c_index)-nu(fuv_index(nfuv-1)))/(nu(xray_index(0))-nu(fuv_index(nfuv-1)))    + l_nu(fuv_index(nfuv-1))




;   Plot Spectrum
; ------------------
;window,3,xsize=400,ysize=400
;;plot,Xear,Xspec,/xlog,/ylog,yrange=[1.0e-3,0.2],ystyle=1
;plot,nu,(10^(l_nu))*(10.0^(nu)),xrange=[12.0,21.0],xstyle=1,/ylog,yrange=[10.0^(9.5),10.0^(12.0)],ystyle=1

;set_plot,'ps'
;device,filename="marconi_spectrum.eps",/encapsulated
;device,xsize=5.0,ysize=5.0,/inches
;plot,nu,(10^(l_nu))*(10.0^(nu)),xrange=[12.0,21.0],xstyle=1,/ylog,yrange=[10.0^(9.5),10.0^(12.0)],ystyle=1,xtitle='log !Mn [Hz]',ytitle='log !Mn L!L!Mn!N [L!LSun!N]',font=1,charsize=1.2
;xyouts,0.66,0.87,'L!Lbol!N = 1.0 x 10!U12!N L!Lsun!N',font=1,charsize=1.2,/normal
;device,/close
;set_plot,'x'




;  Done
; -------
return,[nu,l_nu]


end



