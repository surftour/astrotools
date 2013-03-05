function oplot_rho_cutoff_y_log, dummy, rhocut=rhocut

; ----------------------------------
; these numbers are for h=0.75
; but we want to allow for any h or
; taking h out all together
; ----------------------------------

;hubble= fload_cosmology('h')
hubble= 1.0
if not keyword_set(rhocut) then rhocut= 0.0171*hubble*hubble     ; units of Msolar/pc3



if dummy eq 1 then begin
	y=1e8*findgen(10)+0.1
	x= rhocut  + 0*y        ; convert to Msolar/pc3

	; convert to log
	x= alog10(x)
	y= alog10(y)

	oplot, x, y, linestyle=0, color=0
endif

;print, "oplot_rho_cutoff_y is using rhocut=",rhocut

return, 1

end



