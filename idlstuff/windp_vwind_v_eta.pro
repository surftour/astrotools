pro allowed, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "allowed, junk"
   return
endif




initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename='windallowed.eps', colortable=3



;--------------------------------------
;  Print the Shit
;--------------------------------------


; wind velocity
; 
ymax = 1000.0
ymin = 0.0
yaxistitle="!6Wind Velocity, !8v!6!DW!N (km s!E-1!n)"

; eta
; (wind mass - sfr conversion factor)
xmax = 10.0
xmin = 0.03
;xmin = 0.003
xaxistitle="!6Wind Efficiency, !7g!6"


;---------------------------

!p.position= [0.19, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
	ystyle=1, $
	;ystyle=8, $                     ; suppress right y axis
	;/ylog, $
	/xlog, $
        xcharsize=1.5, ycharsize=1.5, $
        charthick=2.0, $
        xthick=4.0, ythick=4.0, $
	;ytickformat='noexp_label', $
	;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata






	; ----------------------
        ; designate regions
	; ----------------------

	fload_newcolortable, 1

	; no burst region
	; ---------------  
	;etareg= [8.0,0.2,xmax,xmax,8.0]
	;vwindreg= [ymin,ymax,ymax,ymin,ymin]
	;polyfill, etareg, vwindreg, /data, color= 50, /line_fill, linestyle=0, $
	                               thick=3.0

	; white out area, and label it
	;xs= [0.6,0.6,1.8,1.8,0.6]
	;ys= [690,770,770,690,690]
	;polyfill, xs, ys, /data, color= 1, /fill
	;xyouts, 0.8, 740, 'no', /data, size=1.3, charthick=2.0, color= 0
	;xyouts, 0.7, 700, 'burst', /data, size=1.3, charthick=2.0, color= 0


	; nothing happens region
	; ----------------------
	; old
	;xs= [0.2, 0.013]
	;ys= [0, 1200]
	;oplot, xs, ys, linestyle= 1, color= 0, thick= 3.0, psym=-3 
	; new
	;etareg= [xmin,xmin,0.013,0.2,xmin]
        ;vwindreg= [ymin,ymax,ymax,ymin,ymin]
        ;polyfill, etareg, vwindreg, /data, color= 100, /line_fill, linestyle=0, $
        ;                       thick=3.0, orientation=45.0


	; white out and label
	;xs= [0.01,0.01,0.07,0.07,0.01]
        ;ys= [340,400,400,340,340]
        ;polyfill, xs, ys, /data, color= 1, /fill
	;xyouts, 0.012, 360, 'no wind', /data, size=1.3, charthick=2.0, color= 0



	; black hole is too large
        ; ------------------------
	;ets= (alog10(xmax)-alog10(xmin))*findgen(30)/29.0 + alog10(xmin)
	;ets= 10^(ets)
	;ets= [ets,xmax]

	;vws= [140, 150, 160, 170, 180, $
;		200, 230, 260, 300, 340, $
	;	380, 420, 460, 485, 500, $
	;	485, 465, 445, 425, 400, $
	;	375, 350, 325, 300, 275, $
	;	260, 250, 240, 230, 220, 210]

	;etareg= [xmin, ets, xmax, xmin]
        ;vwindreg= [ymin, vws, ymin,ymin]
        ;polyfill, etareg, vwindreg, /data, color= 200, /fill, linestyle=0
;                               thick=3.0, orientation=90.0

	;xyouts, 0.05, 40, 'BH is too massive', /data, size=1.3, charthick=2.0, color= 0



        ; no burst region
        ; ---------------  
        ;etareg= [2.6,0.8,xmax,xmax,2.6]
        ;vwindreg= [ymin,ymax,ymax,ymin,ymin]
        ;polyfill, etareg, vwindreg, /data, color= 50, /line_fill, linestyle=0, $
        ;                               thick=3.0
	;
	y=findgen(1001)
	wigzero= 1.3
	wignorm= 0.1
	period= 250
	x= wignorm * sin(2*!PI*y/period)
        etareg= [wigzero*10^(x),xmax,xmax,x(0)]
        vwindreg= [y,ymax,ymin,ymin]
        polyfill, etareg, vwindreg, /data, color= 50, /line_fill, linestyle=0, $
                                       thick=3.0

        ; white out area, and label it
        xs= [1.2,1.2,3.4,3.4,1.2]
        ys= [540,620,620,540,540]
        polyfill, xs, ys, /data, color= 1, /fill
	oplot, xs, ys, psym=-3, color= 0, thick=1.0
        xyouts, 1.6, 590, 'no', /data, size=1.3, charthick=2.0, color= 0
        xyouts, 1.4, 550, 'burst', /data, size=1.3, charthick=2.0, color= 0






        ; parameter sets
	; ----------------------

	fload_newcolortable, 3

	; nothing happened
        ; ---
	; sb 1, 2, 3, 4, 5, 6
        eta= [0.005,0.005,0.005,0.005,0.005,0.005]
        vwind= [837,1675,418,209,104,5298]
        usersym,cos(findgen(49)/49*2*!pi),sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        oplot, eta, vwind, psym=8, color= 0, symsize=1.1, thick= 6.0

	; sb 7, 9, 19, 18
        eta= [0.05, 0.05, 0.05, 0.05]
        vwind= [132, 209, 418, 837]
        usersym,cos(findgen(49)/49*2*!pi),sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        oplot, eta, vwind, psym=8, color= 0, symsize=1.1, thick= 6.0


        ; strange oscillations
        ; ---
	; sb 13
        eta= [2.0]
        vwind= [209]
        oplot, eta, vwind, psym=7, color= 150, symsize=1.5, thick= 6.0

	; sb 14, 21, 22
        eta= [2.0, 5.0, 5.0]
        vwind= [105, 209, 105]
        usersym,cos(findgen(49)/49*2*!pi),sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        oplot, eta, vwind, psym=8, color= 0, symsize=1.1, thick= 6.0


        ; large blowout
        ; ---
	; sb 8
        eta= [2.0]
        vwind= [837]
        oplot, eta, vwind, psym=7, color= 150, symsize=1.5, thick= 6.0

	; sb 11, 20
        eta= [5.0, 5.0]
        vwind= [837, 418]
        usersym,cos(findgen(49)/49*2*!pi),sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        oplot, eta, vwind, psym=8, color= 0, symsize=1.1, thick= 6.0


	; moderate blowout
	; ---
	; sb 10
        eta= [0.5]
        vwind= [837]
        oplot, eta, vwind, psym=7, color= 150, symsize=1.5, thick= 6.0


	; small blowout
	; ---
	; sb 12, 15, 16, 17
        eta= [2.0, 0.5, 0.5, 0.5]
        vwind= [418, 418, 209, 105]
        ;oplot, eta, vwind, psym=7, color= 100, symsize=1.5, thick= 6.0
        usersym,cos(findgen(49)/49*2*!pi),sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        oplot, eta, vwind, psym=8, color= 0, symsize=1.1, thick= 6.0



;--------------------------------------
; 
;alpha= 3.51
;M_100= 10^(9.79-(2*alpha))

;x=[xmin,xmax]
;y= M_100 * (x^alpha)
;oplot, x, y, psym=-3, thick=3.0, linestyle=3, color= 0


;--------------------------------------
; 
;  Black out region with Chi > 1
;
;x=[xmin,xmax]
x=2.0 * findgen(100)/99.0 - 1.0
x=10^x
chi= 1.0
y= 1185.0 * sqrt(chi/x)
;oplot, x, y, psym=-3, thick=8.0, linestyle=0, color= 0

idx=where(y le ymax)
x= x(idx)
y= y(idx)

idx= where(x le xmax)
x= x(idx)
y= y(idx)

xat= (1000.0/1185.0)^(-2.0)

xs= [xmax, xat, x, xmax]
ys= [ymax, ymax, y, ymax]

; now fill it in black
polyfill, xs, ys, /data, color= 0, /fill

; white out area, and label it
xs= [3.0,3.0,8.0,8.0,3.0]
ys= [810,890,890,810,810]
polyfill, xs, ys, /data, color= 1, /fill
xyouts, 3.3, 840, '!7v!6 > 1', /data, size=1.3, charthick=2.0, color= 0


	; redraw this symbol so that it's on top of black region
	; sb 8
        eta= [2.0]
        vwind= [837]
        oplot, eta, vwind, psym=7, color= 150, symsize=1.5, thick= 6.0


;
;  Chi= 0.1 line
;x=2.0 * findgen(100)/99.0 - 1.0
x=2.0 * findgen(10)/9.0 - 1.0
x=10^x
chi= 0.1
y= 1185.0 * sqrt(chi/x)
oplot, x, y, psym=-3, thick=4.0, linestyle=1, color= 0
; white out area, and label it
xs= [0.1,0.1,0.3,0.3,0.1]
ys= [850,930,930,850,850]
polyfill, xs, ys, /data, color= 1, /fill
oplot, xs, ys, color= 0, psym=-3, linestyle= 0
xyouts, 0.115, 880, '!7v!6=0.1', /data, size=1.5, charthick=2.0, color= 0

;
;  Chi= 0.01 line
;x=3.0 * findgen(100)/99.0 - 2.0
x=3.0 * findgen(10)/9.0 - 2.0
x=10^x
chi= 0.01
y= 1185.0 * sqrt(chi/x)
oplot, x, y, psym=-3, thick=4.0, linestyle=2, color= 0
; white out area, and label it
xs= [0.035,0.035,0.127,0.127,0.035]
ys= [490,570,570,490,490]
polyfill, xs, ys, /data, color= 1, /fill
oplot, xs, ys, color= 0, psym=-3, linestyle= 0
xyouts, 0.04, 520, '!7v!6=0.01', /data, size=1.5, charthick=2.0, color= 0




;--------------------------------------

!p.position= [0.19, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $                     ; suppress right y axis
        ;/ylog, $
        /xlog, $
        xcharsize=1.5, ycharsize=1.5, $
        charthick=2.0, $
        xthick=4.0, ythick=4.0, $
        ;ytickformat='noexp_label', $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata, /noerase


;--------------------------------------

device, /close








end







