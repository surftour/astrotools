function oplot_smoothing_length, smoothinglength, linear=linear, devac=devac, plog=plog

if not keyword_set(smoothinglength) then smoothinglength=0.1

; ----------------------------------
; ----------------------------------

;smoothinglength = 0.1     ; in GADGET units  = kpc/h

xmax=1e12

if smoothinglength GT 0 then begin
	y=alog10(xmax)*findgen(20)+0.0001 - 3.0
	x= smoothinglength + 0.*y        ; convert to Msolar/pc3
;print, "x= ", x
;print, "y= ", y
	if keyword_set(linear) then x= x             ; default
	if keyword_set(devac) then x= x^(0.25)       ; de vacouleurs
	if keyword_set(plog) then x= alog10(x)       ; log plot
	oplot, x, y, linestyle=1, color= 0
endif

;print, "oplot_smoothing_length is using smoothinglength=",smoothinglength

return, 1

end



