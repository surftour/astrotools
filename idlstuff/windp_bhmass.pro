;------------------------------------------------------------------------
;
;    BH mass vs. V_wind
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





pro plot_bhm, junk, filename=filename


if not keyword_set(junk) then begin
   print, "  "
   print, " plot_bhm, junk, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='bhm.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
;setup_plot_stuff, 'ps', filename=filename, colortable=2
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; set up constants
; -----------------------------

xmax= 925.0
xmin= 0.0
xaxistitle='!6Wind Velocity, !8v!6!DW!N (km s!E-1!n)'

ymax = 9.8
ymin = 6.2
yaxistitle= '!6Log Black Hole Mass (M!D!9n!6!N)'



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

xpt1= 590.0
xpt2= 640.0
xpt3= 675.0

;ylbl= ymin+20.0
;pointt= 1
;select_thispoint, pointt, thispsym, thiscolor
;oplot, vwind, sfrmax, psym=-thispsym, color=thiscolor, thick=6.0
;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl, '!7g!6= 0.005', size=1.2, color=thiscolor, /data, charthick=4.0



; eta= 0.05  (sb7, sb19, sb18, sb9)
; ------------------
bhmass= [169.2, 175.7, 255.8, 151.9]
bhmass= alog10(bhmass * 1.0d+10 / 0.7)

ylbl= 9.6  ;ymin+15.0
pointt= 2
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, bhmass, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.05, '!7g!6= 0.05', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 0.5   (sb17, sb16, sb15, sb10)
; ------------------
bhmass= [89.6, 113.8, 174.6, 126.5]
bhmass= alog10(bhmass * 1.0d+10 / 0.7)

ylbl= 9.3  ; ymin+12.0
pointt= 3
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, bhmass, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.05, '!7g!6= 0.5', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 2.0   (sb14, sb13, sb12, sb8)
; ------------------
bhmass= [56.6, 108.3, 43.6, 20.0]
bhmass= alog10(bhmass * 1.0d+10 / 0.7)

ylbl= 9.0  ;ymin+4.0
pointt= 4
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, bhmass, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.05, '!7g!6= 2.0', size=1.2, color=thiscolor, /data, charthick=4.0




; eta= 5.0    (sb22, sb21, sb20, sb11)
; ------------------
bhmass= [23.5, 35.9, 20.0, 20.0]
bhmass= alog10(bhmass * 1.0d+10 / 0.7)

ylbl= 8.7  ;ymin+2.0
pointt= 5
select_thispoint, pointt, thispsym, thiscolor
oplot, vwind, bhmass, psym=-thispsym, color=thiscolor, thick=6.0, symsize=1.3
oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl-0.05, '!7g!6= 5.0', size=1.2, color=thiscolor, /data, charthick=4.0





; done
; -----
device, /close


end
   




;--------------------------------------------------------------------------









