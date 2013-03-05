;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Histogram of Gas Densities
;     ------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------



pro rhohist, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "rhohist, junk"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='rhohist.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius, log scale
xaxistitle= "!6Log Gas Desntiy (M!D!9n!6!N pc!E-3!N)"
;xaxistitle= "!6Log Gas Desntiy (M!D!9n!6!N pc!E-3!N)"
xmax = 5.0
xmin = -7.0

; number (histogram)
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

plot, [0], [0], psym=-3, linestyle=0, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
	xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0,$
	xtitle=xaxistitle, ytickformat='(a1)', /nodata, /noerase


; -----------------------------------------------

;snapnum= 6
snapnum= 8

do_one_rhohist, 'yuexing', snapnum, xmin=xmin, xmax=xmax, thisc= 50
xyouts, 0.10, 0.85, "yuexing", size=1.5, color=50, /normal, charthick=3.0

; -----------------------------------------------

do_one_rhohist, 'yuexing2', snapnum, xmin=xmin, xmax=xmax, thisc= 150
xyouts, 0.10, 0.80, "yuexing2", size=1.5, color=150, /normal, charthick=3.0

; -----------------------------------------------


; print extras
; -------------


; done
; -----
device, /close


end





;=======================================================================






pro do_one_rhohist, frun, snapnum, xmin=xmin, xmax=xmax, thisc=thisc

ok=fload_snapshot_bh(frun, snapnum)

;rho= fload_gas_rho(1)

; convert this to Msolar pc-3
rho= 10.0*fload_gas_rho(1)

; convert to cm-3
;UnitDensity_in_cgs = 6.76991d-22
;ProtonMass = 1.6726d-24
;rho = fload_gas_rho(1) * UnitDensity_in_cgs / ProtonMass

print, "gas densities max/min= ", max(rho), min(rho)

rho= alog10(rho)

temp= process_histogram(rho, xmax=xmax, xmin=xmin, levels=100, oplotit=thisc, mannorm= 500.0)
print, min(temp), max(temp)


end



;=======================================================================

