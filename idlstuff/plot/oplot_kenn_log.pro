function oplot_kenn_log, sendto, linear=linear

; ----------------------------------
; kennicutt uses h=0.75
; but we want to allow for any h or
; taking h out all together
; ----------------------------------

h = 0.70     ; 1.0 means in units of h, h-2 etc.
; h = 1.0

x=100000*findgen(10)+0.1
y=(2.5e-4)*((h/0.75)^0.6)*(x^(1.4))
;y1=(3.2e-4)*((h/0.75)^0.6)*(x^(1.4))
;y2=(1.8e-4)*((h/0.75)^0.6)*(x^(1.4))

if keyword_set(linear) then begin
	x= alog10(x)
	y= alog10(y)
endif

if sendto eq 'ps' then begin
	oplot, x, y, psym=-3, linestyle=1, color= 0, thick= 2.0
endif else begin
	oplot, x, y, psym=-3, linestyle=0, color= getcolor('black')
endelse



; -------------------------------------------------
; now plot the rough outline of the datapoints
;    this is for the linear keyword begin on
;y1= y+alog10(10)
;oplot, x, y1, psym=-3, linestyle=1, color= 0, thick= 2.0
;y2= y-alog10(10)
;oplot, x, y2, psym=-3, linestyle=1, color= 0, thick= 2.0


; -------------------------------------------------
; do 1 sigma if we don't want to do the outline
;    this is for the linear keyword begin on
;
;    2.5+/-0.7, so essentially x1.28
;
;y1= y+alog10(1.28)
;oplot, x, y1, psym=-3, linestyle=1, color= 0, thick= 2.0
;y2= y-alog10(1.28)
;oplot, x, y2, psym=-3, linestyle=1, color= 0, thick= 2.0

; -------------------------------------------------

;print, "oplot_kenn_log is using h=",h

return, 1

end



