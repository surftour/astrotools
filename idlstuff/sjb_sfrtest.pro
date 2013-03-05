pro sfrtest, junk, filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "  "
   print, "sfrtest, junk"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename = 'sjb_sfrtest.eps'

initialize_plotinfo, 1

setup_plot_stuff, 'ps', filename=filename, colortable=4
; should I try and send a different size for this?



;-------------
;  Plot axes 
;-------------

xaxistitle= "Gas Density (Gadget Units)"
;yaxistitle= "SFR (M!D!9n!3!N Yr!E-1!N)"
yaxistitle= "SFR/V (M!D!9n!3!N Yr!E-1!N / V)"
xmax = 900.0
;xmax = 800.0
xmin = 1.0e-5
ymax = 0.1
ymin = 1.0e-5



!p.position= [0.20, 0.15, 0.98, 0.98]

plot, [10.0], [10.0], psym=-3, linestyle=0, $
        color= 0, /xlog, /ylog, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




;-----------------------
;  Put data on graph
;-----------------------2

frun= "data/ds/vc3vc3e_2"
snapnum= 48


ok= fload_snapshot_bh(frun, snapnum, /nopot_in_snap)

rho= fload_gas_rho(1)
sfr= fload_gas_sfr(1)
hsml= fload_gas_hsml(1)
vol= hsml * hsml * hsml

;oplot, rho, sfr, psym=2, color= 150
oplot, rho, sfr/vol, psym=2, color= 150

;stop

x= [xmin,xmax]
y= 1.0e-4 * x^(1.5)

oplot, x, y, psym=-3, color= 50






;--------------------------------------
device, /close


end




