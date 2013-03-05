pro tdwarfs_find, frun, snapnum

if not keyword_set(frun) then begin
   print, "  "
   print, "  "
   print, "tdwarfs_find, frun, snapnum"
   print, "  "
   print, "  "
   return
endif


filename='findtdwarfs.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4




;---------------
;  Print it Up
;---------------

xaxistitle="!6R (kpc)"
xmax = 200.0
xmin = 0.0

;yaxistitle="!6Log Gas Density (GU) (10!E10!N M!D!9n!6!N)"
yaxistitle="!6Log Gas Density (M!D!9n!6!N pc!E-1!N)"
ymax = 1.5
ymin = -8.0



!p.position= [0.17, 0.14, 0.98, 0.98]

plot, [10.0], [10.0], psym=-3, linestyle=0, $
        color= 0, $
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
;-----------------------


ok=fload_snapshot(frun,snapnum)

radius= fload_gas_xyz('r')
gasdensity= alog10(10.0*fload_gas_rho(1))
u= fload_gas_u(1)

idx= where(u gt 500.0)
oplot, radius(idx), gasdensity(idx), psym=3, color= 150

idx= where(u le 500.0)
oplot, radius(idx), gasdensity(idx), psym=3, color= 50


xyouts, 0.22, 0.90, frun, color= 0, charthick=2.0, size= 1.2, /normal
xyouts, 0.22, 0.85, fload_timelbl(1,2), color= 0, charthick=2.0, size= 1.0, /normal


;-------------
device, /close


end




