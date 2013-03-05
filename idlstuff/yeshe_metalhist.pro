;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Histogram of Stellar Particle Metallicities
;     ------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro metal_hist, frun, snapnum, filename=filename

if not keyword_set(frun) then begin
   print, "  "
   print, "metal_hist, frun, snapnum, filename=filename"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='metalhist.eps'
if not keyword_set(snapnum) then snapnum=25


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius, linear scale
xaxistitle= "Log !6Z (Z!D!9n!6!N)"
xmax = 1.0
xmin = -3.0

; number (histogram)
;yaxistitle= "!6V!Dr!N (km s!E-1!N)"
yaxistitle= ' '
ymax = 1.2
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.05
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
	ytickformat='(a1)', $
	;ytitle=yaxistitle, $
	/nodata, /noerase


; -----------------------------------------------


; ------------------
; compute histogram
; ------------------


ok=fload_snapshot_bh(frun,snapnum)

newstar_metals= fload_allstars_z(1)
print, "new star metallicity max/min= ", max(newstar_metals), min(newstar_metals)
newstar_metals= fload_allstars_z(1) / 0.02
print, "    (in solar)                ", max(newstar_metals), min(newstar_metals)

newstar_metals= alog10(newstar_metals)

temp= process_histogram(newstar_metals, xmax=xmax, xmin=xmin, levels=100, oplotit=50)
print, min(temp), max(temp)

; -----------------------------------------------

avg_metallicity= mean(newstar_metals)
print, "avg= ", avg_metallicity
avglbl= strcompress(string(avg_metallicity),/remove_all)
xyouts, 0.10, 0.75, "<Z(Z!D!9n!6!N)> = "+avglbl, size=1.2, color=0, /normal, charthick=3.0



;xyouts, 0.75, 0.90, fload_timelbl(1,2), size=1.2, color=0, /normal, charthick=4.0
;xyouts, 0.55, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0
;xyouts, 0.75, 0.94, fload_fid(1), /normal, size= 1.2, charthick=3.0, color= 0

;xyouts, 0.10, 0.85, "!6V!D200!N = 50 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
;xyouts, 0.10, 0.85, "!6V!D200!N = 500 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0
xyouts, 0.10, 0.85, "!6V!D200!N = 160 km s!E-1!N", size=1.5, color=0, /normal, charthick=3.0



; print extras
; -------------


; done
; -----
device, /close


end






