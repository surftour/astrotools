;------------------------------------------------------------------------
;
;    Random Procedures related to Metals
;
;
;
;
;------------------------------------------------------------------------




;--------------------------------------------------------------------------



; routine to select points
; ---------------------------
pro select_thispoint, pointselection, thispsym, thiscolor


;  open (or filled) square
; --------------------------
if pointselection eq 1 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 150
endif

;  x's
; -----
if pointselection eq 2 then begin
        symsel= 7
        symcolor= 100
endif

;  open circle
; --------------
if pointselection eq 3 then begin
        symsize= 1.0
        ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        symsel= 8
        symcolor= 50
endif

;  open triangle
; ----------------
if pointselection eq 4 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
        usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
        symsel= 8
        symcolor= 0
endif

;  asterisk
; ----------
if pointselection eq 5 then begin
        symsel= 2
        symcolor= 200
endif

;  filled square
; --------------------------
if pointselection eq 6 then begin
        symsize= 1.0
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 20
endif

;  diamond
; ----------
if pointselection eq 7 then begin
        symsel= 4
        symcolor= 120
endif

;  triangle (open)
; -----------------
if pointselection eq 8 then begin
        symsel= 5
        symcolor= 180
endif

;  plus sign
; ----------
if pointselection eq 9 then begin
        symsel= 1
        symcolor= 170
endif


thispsym=symsel
thiscolor=symcolor


end












;=====================================================================================





pro plot_sfrmax, junk, filename=filename


if not keyword_set(junk) then begin
   print, "  "
   print, " plot_sfrmax, junk, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfrmax.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
;setup_plot_stuff, 'ps', filename=filename, colortable=2
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; set up constants
; -----------------------------

xmax= 1000.0
xmin= 0.0
xaxistitle='!6Wind Velocity, !8v!6!DW!N (km s!E-1!n)'

ymax = 300.0
ymin = 10.0
yaxistitle= '!6Maximum SFR (M!D!9n!6!N Yr!E-1!N)'



; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
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


; eta= 0.005
; ------------------
;sfrmax= [167.7, 121.8, 134.4, 103.6]

xpt1= 700.0
xpt2= 760.0
xpt3= 785.0

;ylbl= ymin+20.0
;pointt= 1
;select_thispoint, pointt, thispsym, thiscolor
;oplot, vwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl, '!7g!6= 0.005', size=1.2, color=thiscolor, /data, charthick=4.0



; eta= 0.05
; ------------------
;sfrmax= [169.2, 147.7, 255.8, 151.9]
sfrmax= [169.2, 175.7, 255.8, 151.9]

ylbl= 75.0  ;ymin+15.0
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-3.0, '!7g!6= 0.05', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 0.5
; ------------------
sfrmax= [89.6, 113.8, 174.6, 126.5]

ylbl= 60.0  ; ymin+12.0
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-2.0, '!7g!6= 0.5', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 2.0
; ------------------
sfrmax= [56.6, 108.3, 43.6, 20.0]

ylbl= 48.0  ;ymin+4.0
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-2.0, '!7g!6= 2.0', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 5.0
; ------------------
sfrmax= [23.5, 35.9, 20.0, 20.0]

ylbl= 37.0  ;ymin+2.0
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-2.0, '!7g!6= 5.0', size=1.2, color=thiscolor, /data, charthick=4.0







; print extras
; -------------
x= [xmin,xmax]

y= [172.0, 172.0]
oplot, x, y, psym=-3, color= 0, linestyle= 1
xyouts, 0.82, 0.90, 'no wind', size=1.1, color=0, /normal
xyouts, 0.82, 0.87, 'max SFR', size=1.1, color=0, /normal

y= [100.0, 100.0]
oplot, x, y, psym=-3, color= 0, linestyle= 0
xyouts, 0.52, 0.73, 'burst criteria', size=1.1, color=0, /normal

y= [20.0, 20.0]
oplot, x, y, psym=-3, color= 0, linestyle= 2
xyouts, 0.25, 0.26, '2x quiescent', size=1.1, color=0, /normal



; done
; -----
device, /close


end
   




;--------------------------------------------------------------------------






pro sfrmax_v_windm, junk, filename=filename


if not keyword_set(junk) then begin
   print, "  "
   print, " sfrmax_v_windm, junk, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfrmax.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
;setup_plot_stuff, 'ps', filename=filename, colortable=2
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; set up constants
; -----------------------------

xmax= 10.3
xmin= 6.0
xaxistitle='!6Log Wind Mass (M!D!9n!6!N)'

ymax = 300.0
ymin = 0.0
yaxistitle= '!6Maximum SFR (M!D!9n!6!N Yr!E-1!N)'



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



; eta= 0.005
; ------------------
;sfrmax= [167.7, 121.8, 134.4, 103.6]

xpt1= 700.0
xpt2= 760.0
xpt3= 785.0

;ylbl= ymin+20.0
;pointt= 1
;select_thispoint, pointt, thispsym, thiscolor
;oplot, vwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl, '!7g!6= 0.005', size=1.2, color=thiscolor, /data, charthick=4.0



; eta= 0.05  (sb7, sb19, sb18, sb9)
; ------------------
;sfrmax= [169.2, 147.7, 255.8, 151.9]
sfrmax= [169.2, 175.7, 255.8, 151.9]
;mwind= [0.00365, xx, xx, 0.07575]
mwind= [0.00365, 0.008, 0.040, 0.07575]    ; fake data
mwind= alog10(mwind * 1.0d+10 / 0.7)

ylbl= 75.0  ;ymin+15.0
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, mwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-3.0, '!7g!6= 0.05', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 0.5   (sb17, sb16, sb15, sb10)
; ------------------
sfrmax= [89.6, 113.8, 174.6, 126.5]
;mwind= [xx, xx, xx, 0.79305]
mwind= [0.001, 0.02, 0.1, 0.79305]
mwind= alog10(mwind * 1.0d+10 / 0.7)

ylbl= 60.0  ; ymin+12.0
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, mwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-2.0, '!7g!6= 0.5', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 2.0   (sb14, sb13, sb12, sb8)
; ------------------
sfrmax= [56.6, 108.3, 43.6, 20.0]
mwind= [0.00551, 0.00898, 0.04481, 1.80527]
mwind= alog10(mwind * 1.0d+10 / 0.7)

ylbl= 48.0  ;ymin+4.0
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, mwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-2.0, '!7g!6= 2.0', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 5.0    (sb22, sb21, sb20, sb11)
; ------------------
sfrmax= [23.5, 35.9, 20.0, 20.0]
;mwind= [xx, xx, xx, 2.51611]
mwind= [0.001, 0.01, 0.1, 2.51611]
mwind= alog10(mwind * 1.0d+10 / 0.7)

ylbl= 37.0  ;ymin+2.0
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, mwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-2.0, '!7g!6= 5.0', size=1.2, color=thiscolor, /data, charthick=4.0







; print extras
; -------------
;x= [xmin,xmax]

;y= [172.0, 172.0]
;oplot, x, y, psym=-3, color= 0, linestyle= 1
;xyouts, 0.82, 0.90, 'no wind', size=1.1, color=0, /normal
;xyouts, 0.82, 0.87, 'max SFR', size=1.1, color=0, /normal

;y= [100.0, 100.0]
;oplot, x, y, psym=-3, color= 0, linestyle= 0
;xyouts, 0.52, 0.73, 'burst criteria', size=1.1, color=0, /normal

;y= [20.0, 20.0]
;oplot, x, y, psym=-3, color= 0, linestyle= 2
;xyouts, 0.25, 0.26, '2x quiescent', size=1.1, color=0, /normal



; done
; -----
device, /close


end
   




;--------------------------------------------------------------------------






