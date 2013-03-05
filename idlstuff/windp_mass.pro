;------------------------------------------------------------------------
;
;    Random Procedures related to Metals
;
;
;
;
;------------------------------------------------------------------------



pro windm, junk, filename=filename


if not keyword_set(junk) then begin
   print, "  "
   print, " windm, junk, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='windm.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

ymax= 10.9
ymin= 7.4
yaxistitle='!6Log Wind Mass (M!D!9n!6!N)'

xmax= 950.0
xmin= 0.0
xaxistitle='!6Wind Velocity, !8v!6!DW!N (km s!E-1!n)'



; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
	;/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; -----------------------------------------------


; universal wind speed 
vwind= [105, 209, 418, 837]

xpt1=  60.0
xpt2= 110.0
xpt3= 150.0



; eta= 0.005
; ------------------
;sfrmax= [167.7, 121.8, 134.4, 103.6]

;ylbl= ymin+20.0
;pointt= 1
;select_thispoint, pointt, thispsym, thiscolor
;oplot, vwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl, '!7g!6= 0.005', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 0.05  (sb7, sb19, sb18, sb9)
; ------------------
mwind= [0.00365, 0.00467, 0.00554, 0.07575]    ; fake data
mwind= alog10(mwind * 1.0d+10 / 0.7)

ylbl= 10.3  ;ymin+15.0
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, mwind, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.03, '!7g!6= 0.05', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 0.5   (sb17, sb16, sb15, sb10)
; ------------------
mwind= [0.00469, 0.00498, 0.01699, 0.79305]
mwind= alog10(mwind * 1.0d+10 / 0.7)

ylbl= 10.0  ; ymin+12.0
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, mwind, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.03, '!7g!6= 0.5', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 2.0   (sb14, sb13, sb12, sb8)
; ------------------
mwind= [0.00551, 0.00898, 0.04481, 1.80527]
mwind= alog10(mwind * 1.0d+10 / 0.7)

ylbl= 9.7  ;ymin+4.0
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, mwind, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.03, '!7g!6= 2.0', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 5.0    (sb22, sb21, sb20, sb11)
; ------------------
mwind= [0.00515, 0.01665, 0.07407, 2.51611]
mwind= alog10(mwind * 1.0d+10 / 0.7)

ylbl= 9.4  ;ymin+2.0
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, mwind, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.03, '!7g!6= 5.0', size=1.2, color=thiscolor, /data, charthick=4.0




; print extras
; -------------
x= [xmin,xmax]

mwind= [0.004, 0.004]
mwind= alog10(mwind * 1.0d+10 / 0.7)
oplot, x, mwind, psym=-3, color= 0, linestyle= 1, thick= 3.0
;xyouts, 0.82, 0.90, 'no wind', size=1.1, color=0, /normal
;xyouts, 0.82, 0.87, 'max SFR', size=1.1, color=0, /normal


xpt4= 330.0 
ylbl= alog10(3.12387 * 1.0d+10 / 0.7)
xyouts, xpt3-30.0, ylbl-0.03, 'initial gas mass', size=1.2, color=0, /data, charthick=4.0
arrow, xpt4+200.0, ylbl, xpt4+340.0, ylbl, $
                        COLOR=zcolor, THICK=3.0, hthick=3.0, /data


;x= [685.0, 900.0]
;y= [ylbl, ylbl]
;oplot, x, y, psym=-3, thick=3.0, color= 0, linestyle= 2



; done
; -----
device, /close


end
   




;--------------------------------------------------------------------------







pro windm2, junk, filename=filename


if not keyword_set(junk) then begin
   print, "  "
   print, " windm2, junk, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='windm2.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

;ymax= 10.9
;ymin= 7.4
;yaxistitle='!6Log Wind Mass (M!D!9n!6!N)'
ymax= 1.0
ymin= 0.0
yaxistitle='!6Wind Mass Fraction'

xmax= 2.0
;xmin= 0.0001
xmin= 0.0
xaxistitle='Wind Energy Fraction, !7v!6'



; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
	;/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
	;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; -----------------------------------------------



xpt1=  60.0
xpt2= 110.0
xpt3= 150.0



; eta= 0.005
; ------------------
;sfrmax= [167.7, 121.8, 134.4, 103.6]

;ylbl= ymin+20.0
;pointt= 1
;select_thispoint, pointt, thispsym, thiscolor
;oplot, vwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl, '!7g!6= 0.005', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 0.05  (sb7, sb19, sb18, sb9)
; ------------------
vwind= [0.000625, 0.0015625, 0.00625, 0.025]
mwind= [0.00365, 0.00467, 0.00554, 0.07575]    ; fake data
;mwind= alog10(mwind * 1.0d+10 / 0.7)
mwind= mwind/3.12387

ylbl= 10.3  ;ymin+15.0
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, mwind, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.03, '!7g!6= 0.05', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 0.5   (sb17, sb16, sb15, sb10)
; ------------------
vwind= [0.00390625, 0.015625, 0.0625, 0.25]
mwind= [0.00469, 0.00498, 0.01699, 0.79305]
;mwind= alog10(mwind * 1.0d+10 / 0.7)
mwind= mwind/3.12387

ylbl= 10.0  ; ymin+12.0
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, mwind, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.03, '!7g!6= 0.5', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 2.0   (sb14, sb13, sb12, sb8)
; ------------------
vwind= [0.015625, 0.0625, 0.25, 1.0]
mwind= [0.00551, 0.00898, 0.04481, 1.80527]
;mwind= alog10(mwind * 1.0d+10 / 0.7)
mwind= mwind/3.12387

ylbl= 9.7  ;ymin+4.0
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, mwind, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.03, '!7g!6= 2.0', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 5.0    (sb22, sb21, sb20, sb11)
; ------------------
vwind= [0.03125, 0.15625, 0.625, 2.5]
mwind= [0.00515, 0.01665, 0.07407, 2.51611]
;mwind= alog10(mwind * 1.0d+10 / 0.7)
mwind= mwind/3.12387

ylbl= 9.4  ;ymin+2.0
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, mwind, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.03, '!7g!6= 5.0', size=1.2, color=thiscolor, /data, charthick=4.0




; print extras
; -------------
x= [xmin,xmax]

mwind= [0.004, 0.004]
;mwind= alog10(mwind * 1.0d+10 / 0.7)
mwind= mwind/3.12387
oplot, x, mwind, psym=-3, color= 0, linestyle= 1, thick= 3.0
;xyouts, 0.82, 0.90, 'no wind', size=1.1, color=0, /normal
;xyouts, 0.82, 0.87, 'max SFR', size=1.1, color=0, /normal


xpt4= 330.0 
ylbl= alog10(3.12387 * 1.0d+10 / 0.7)
xyouts, xpt3-30.0, ylbl-0.03, 'initial gas mass', size=1.2, color=0, /data, charthick=4.0
arrow, xpt4+200.0, ylbl, xpt4+340.0, ylbl, $
                        COLOR=zcolor, THICK=3.0, hthick=3.0, /data


;x= [685.0, 900.0]
;y= [ylbl, ylbl]
;oplot, x, y, psym=-3, thick=3.0, color= 0, linestyle= 2



; done
; -----
device, /close


end
   




;--------------------------------------------------------------------------






