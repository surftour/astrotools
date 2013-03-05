

;--------------------------------------------------------------------------


function grab_bh_growth, frun


;-------------------------------------
;   Get BH Info from fid.bh file
;-------------------------------------
open_blackholes_file, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass
lasti= n_elements(bh_mass)-1

; final/initial
return, bh_mass[lasti]/bh_mass[0]

end




;--------------------------------------------------------------------------



pro overplot_onerun, frun, mratio, thisptype=thisptype

bhgrowth= grab_bh_growth(frun)

; grab point information
select_thispoint, thisptype, thispsym, thiscolor

oplot, [mratio], [bhgrowth], psym=thispsym, color= thiscolor, symsize= 1.5, thick=4.0

end




;--------------------------------------------------------------------------


pro overplot_onelabel, xx0, yy0, thislabel, thisptype=thisptype


; grab point information
select_thispoint, thisptype, thispsym, thiscolor

xyouts, xx0, yy0, thislabel, /normal, color= thiscolor, size= 1.2, charthick=4.0

end




;--------------------------------------------------------------------------







pro mj_bhgrowth, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "mj_bhgrowth, junk"
   return
endif




;============================================
;============================================
;
;
;
;============================================
;============================================


xaxistitle= "Mass Ratio (M!Dsat!N/M!Dprimary!N)"
xmax= 2.0
;xmax= 1.05
xmin= 0.01
;xmin= 0.0

yaxistitle= "BH Growth (M!Dfinal!N/M!Dseed!N)"
ymax= 1300.0
ymin= 0.2

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename='bhgrowth.eps', colortable=4

!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, ystyle=1, $
        xcharsize=1.5, ycharsize=1.5, charthick=3.0, xthick=4.0, ythick=4.0, color= 0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /xlog, /ylog, ytickformat='exp_label'



;----------------------------------------


; minors (josh's)
overplot_onerun, "/raid4/tcox/minor/min_0", 1/8.0, thisptype= 2
overplot_onerun, "/raid4/tcox/minor/min_30", 1/8.0, thisptype= 2
overplot_onerun, "/raid4/tcox/minor/min_90", 1/8.0, thisptype= 2
overplot_onerun, "/raid4/tcox/minor/min_150", 1/8.0, thisptype= 2
overplot_onerun, "/raid4/tcox/minor/min_180", 1/8.0, thisptype= 2

overplot_onerun, "/raid4/tcox/minor/Sbfg0.4Sdfg0.4_000", 1/4.1, thisptype= 2
overplot_onerun, "/raid4/tcox/minor/Sbfg0.4Sdfg0.4_030", 1/4.1, thisptype= 2
overplot_onerun, "/raid4/tcox/minor/Sbfg0.4Sdfg0.4_090", 1/4.1, thisptype= 2

overplot_onelabel, 0.23, 0.92, 'Josh', thisptype= 2


; minors (mine)
;overplot_onerun, "/raid4/tcox/vcs/minors/vc2vc1", 1/2.0, thisptype= 5
overplot_onerun, "/raid4/tcox/vcs/minors/vc3bvc1", 1/8.0, thisptype= 5
overplot_onerun, "/raid4/tcox/vcs/minors/vc3vc1", 1/8.0, thisptype= 5
overplot_onerun, "/raid4/tcox/vcs/minors/vc3vc2", 1/3.5, thisptype= 5
overplot_onerun, "/raid4/tcox/vcs/minors/vc3bvc2", 1/3.5, thisptype= 5

overplot_onelabel, 0.23, 0.88, 'TJ (old)', thisptype= 5


; majors
overplot_onerun, "/raid4/tcox/vc3vc3e", 1/1.0, thisptype= 3
overplot_onerun, "/raid4/tcox/vc3vc3h", 1/1.0, thisptype= 3
overplot_onerun, "/raid4/tcox/vc3vc3j", 1/1.0, thisptype= 3
overplot_onerun, "/raid4/tcox/vc3vc3i", 1/1.0, thisptype= 3
overplot_onerun, "/raid4/tcox/vc3vc3k", 1/1.0, thisptype= 3

overplot_onelabel, 0.23, 0.84, 'TJ (majors)', thisptype= 3


; isolated
overplot_onerun, "/raid4/tcox/isolated/d3", 1/50.0, thisptype= 4
overplot_onerun, "/raid4/tcox/isolated/b3", 1/50.0, thisptype= 4
overplot_onerun, "/raid4/tcox/isolated/vc3c", 1/50.0, thisptype= 4
overplot_onerun, "/raid4/tcox/isolated/e3_1", 1/50.0, thisptype= 4
overplot_onerun, "/raid4/tcox/isolated/e3_2", 1/50.0, thisptype= 4
overplot_onerun, "/raid4/tcox/isolated/e3_bulge1", 1/50.0, thisptype= 4
overplot_onerun, "/raid4/tcox/isolated/e3_bulge2", 1/50.0, thisptype= 4

overplot_onelabel, 0.23, 0.80, 'isolated', thisptype= 4


;----------------------------------------


;x0= 0.30
;y0= 0.82
;xyouts, x0+0.05, y0+0.05, msg0, /normal, charthick=3.0, size=1.3, color= 200
;xyouts, x0+0.05, y0, msg1, /normal, charthick=3.0, size=1.3, color= 150
;xyouts, x0+0.05, y0-0.05, msg2, /normal, charthick=3.0, size=1.3, color= 100
;xyouts, x0+0.05, y0-0.10, msg3, /normal, charthick=3.0, size=1.3, color= 50
;xyouts, x0+0.05, y0-0.15, msg4, /normal, charthick=3.0, size=1.3, color= 0
;oplot, [0.018], [0.73], psym=4, thick=3.0, symsize=1.5, color=200
;oplot, [0.018], [0.68], psym=2, thick=3.0, symsize=1.5, color=150
;oplot, [0.018], [0.63], psym=5, thick=3.0, symsize=1.5, color=100
;oplot, [0.018], [0.58], psym=7, thick=3.0, symsize=1.5, color=50
;oplot, [0.018], [0.53], psym=6, thick=3.0, symsize=1.5, color= 0



oplot, [xmin,xmax], [1.0, 1.0], psym=-3, color=0, thick= 3.0, linestyle= 1



;============================================
device, /close




end




