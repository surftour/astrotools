function oplot_coldtemp_cutoff_log, dummy, rhocut=rhocut

; ----------------------------------
; these numbers are for h=0.75
; but we want to allow for any h or
; taking h out all together
; ----------------------------------

if not keyword_set(rhocut) then begin
        ;hubble= fload_cosmology('h')
        hubble= 1.0
        rhocut= 0.0171*hubble*hubble     ; units of Msolar/pc3
endif

tempcut = 150       ; in GADGET units

if dummy eq 1 then begin
	x=1e-5*findgen(10) + 1e-15
	idx=where(x gt rhocut)
	if idx(0) gt 0 then x(idx)= rhocut else x(9)= rhocut
	y= tempcut/0.012381322  + 0*x        ; convert to K

	; convert to log
	x= alog10(x)
	y= alog10(y)

	oplot, x, y, linestyle=0, color=0
endif

print, "oplot_coldtemp_cutoff is using rhocut=",rhocut, "  and tempcut=",tempcut

return, 1

end



