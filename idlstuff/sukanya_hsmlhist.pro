;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Histogram of Gas Smoothing Lengths
;     ------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro hsmlhist, frun, snapnum, filename=filename

if not keyword_set(frun) then begin
   print, "  "
   print, "hsmlhist, frun, snapnum, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='hsmlhist.eps'
if not keyword_set(snapnum) then snapnum=25


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; smoothing length
xaxistitle= "Log !6Hsml (!8h!6!E-1!N pc)"
xmax = 3.0
xmin = -2.0

; number (histogram)
;yaxistitle= "!6V!Dr!N (km s!E-1!N)"
;yaxistitle= ' '
;ymax = 1.2
yaxistitle= ' Number '
ymax = 800.0
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.95

!p.position= [x0, y0, x1, y1]

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
	;ytickformat='(a1)', $
	ytitle=yaxistitle, $
	/nodata, /noerase


; -----------------------------------------------


; ------------------
; compute histogram
; ------------------


ok=fload_snapshot_bh(frun,snapnum)

gas_hsml= fload_gas_hsml(1)
print, "gas hsml max/min= ", max(gas_hsml), min(gas_hsml)

gas_hsml= alog10(gas_hsml)

temp= process_histogram(gas_hsml, xmax=xmax, xmin=xmin, levels=100, oplotit=50, /nonorm)
print, min(temp), max(temp)

; -----------------------------------------------

avg_gashsml= mean(gas_hsml)
print, "avg= ", avg_gashsml



;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.2, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0
xyouts, 0.25, 0.90, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0
snaplbl= strcompress(string(snapnum),/remove_all)
xyouts, 0.25, 0.85, 'snap='+snaplbl, /normal, size= 1.2, charthick=3.0, color= 0

;xyouts, 0.10, 0.85, "!6V!D200!N = 50 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 500 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 160 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0



; print extras
; -------------


; done
; -----
device, /close


end






;===========================================================================================





pro hsml_v_r, frun, snapnum, filename=filename

if not keyword_set(frun) then begin
   print, "  "
   print, "hsml_v_r, frun, snapnum, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='hsml_v_r.eps'
if not keyword_set(snapnum) then snapnum=25


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; smoothing length
xaxistitle= "Log !6Hsml (!8h!6!E-1!N pc)"
xmax = 3.0
xmin = -2.0

; number (histogram)
;yaxistitle= "!6V!Dr!N (km s!E-1!N)"
;yaxistitle= ' '
;ymax = 1.2
yaxistitle= 'Radius (!8h!6!E-1!N kpc)'
ymax = 800.0
ymin = 0.01



; ------------------
; Plot this up
; ------------------
x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.95

!p.position= [x0, y0, x1, y1]

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
	;ytickformat='(a1)', $
	ytitle=yaxistitle, $
	/nodata, /noerase


; -----------------------------------------------


; ------------------
; compute histogram
; ------------------


ok=fload_snapshot_bh(frun,snapnum)

gas_hsml= fload_gas_hsml(1)
print, "gas hsml max/min= ", max(gas_hsml), min(gas_hsml)
gas_hsml= alog10(gas_hsml)


r= fload_gas_xyz('r')
;r= alog10(r)


oplot, gas_hsml, r, psym=3, color= 150

; -----------------------------------------------


;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.2, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0
xyouts, 0.25, 0.90, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0
snaplbl= strcompress(string(snapnum),/remove_all)
xyouts, 0.25, 0.85, 'snap='+snaplbl, /normal, size= 1.2, charthick=3.0, color= 0

;xyouts, 0.10, 0.85, "!6V!D200!N = 50 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 500 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 160 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0



; print extras
; -------------


; done
; -----
device, /close


end






