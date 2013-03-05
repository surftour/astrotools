

; ---------------------------------------------------

function get_bh_mass, frun
	open_blackholes_file, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass
	return, [max(bh_mass)]
end

; ---------------------------------------------------

function get_sigma, frun
	read_sigma_file, '/raid4/tcox/'+frun, time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
                                Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, $
                                gsigavg, gsigerr
	slstidx= n_elements(time)-1
	return, [Asigavg[slstidx]]
end

; ---------------------------------------------------







;==================================================================
;==================================================================
;
;      Wind Model M_bh v. sigma
;
;
;==================================================================
;==================================================================







;-------------------------------------------------
;
;-------------------------------------------------
pro msig, junk


if not keyword_set(junk) then begin
	print, " "
	print, " msig, junk"
	print, " "
	print, " needs:  .run bh_multi"
	print, "         .run time_m_sigma_re"
	print, " "
	print, " "
	return
endif

filename='jymsig.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4




xaxistitle = "!7r!6!N(km s!E-1!N)"
yaxistitle = "!6M!DBH!N (M!D!9n!6!N)"
xmax = 600
xmin = 30
ymax = 8.0e+9
ymin = 2.0e+5


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
	/xlog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; -----------------------------------
; -----------------------------------


;xpt1= 37.0
;xpt2= 43.0
xpt2= 41.0
xpt3= 48.0





; no winds
; ----------
ylbl= 4.0d+9
pointt= 9
select_thispoint, pointt, thispsym, thiscolor

oplot, [xpt2], [ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl*0.85, 'major', size=1.2, color=thiscolor, /data, charthick=4.0

frun= 'minor/Sbfg0.4Sbfg0.4_000'
oplot, get_sigma(frun), 2.0 * get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Sbfg0.4_030'
oplot, get_sigma(frun), 2.0 * get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Sbfg0.4_090'
oplot, get_sigma(frun), 2.0 * get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Sbfg0.4_150'
oplot, get_sigma(frun), 2.0 * get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Sbfg0.4_180'
oplot, get_sigma(frun), 2.0 * get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor



ylbl= 2.0d+9
pointt= 2
select_thispoint, pointt, thispsym, thiscolor

;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl*0.85, '!7g!6= 0.05', size=1.2, color=thiscolor, /data, charthick=4.0
oplot, [xpt2], [ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl*0.85, '1:2', size=1.2, color=thiscolor, /data, charthick=4.0

frun= 'minor/Sbfg0.4Scfg0.4_000'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Scfg0.4_030'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Scfg0.4_090'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Scfg0.4_150'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Scfg0.4_180'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor





; eta= 0.5   (sb17, sb16, sb15, sb10)
; ------------------
ylbl= 1.1d+9
pointt= 3
select_thispoint, pointt, thispsym, thiscolor

;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl*0.85, '!7g!6= 0.5', size=1.2, color=thiscolor, /data, charthick=4.0
oplot, [xpt2], [ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl*0.85, '1:4', size=1.2, color=thiscolor, /data, charthick=4.0


frun= 'minor/Sbfg0.4Sdfg0.4_000'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Sdfg0.4_030'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Sdfg0.4_090'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Sdfg0.4_150'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Sdfg0.4_180'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor








; eta= 2.0   (sb14, sb13, sb12, sb8)
; ------------------
ylbl= 6.0d+8
pointt= 4
select_thispoint, pointt, thispsym, thiscolor

;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl*0.85, '!7g!6= 2.0', size=1.2, color=thiscolor, /data, charthick=4.0
oplot, [xpt2], [ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl*0.85, '1:8', size=1.2, color=thiscolor, /data, charthick=4.0


frun= 'minor/Sbfg0.4Imfg0.4_000'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Imfg0.4_030'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Imfg0.4_090'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Imfg0.4_150'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'minor/Sbfg0.4Imfg0.4_180'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor






; -----------------------------------

mblackhole_sigma_relation, 1

;xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0
;if bhmsg NE '' then xyouts, 0.25, 0.80, bhmsg, /normal, charthick=3.0, size=1.33, color=0




; -----------------------------------

device, /close


end




; -----------------------------------




;   m-sigma relation
; -----------------------------------
pro mblackhole_sigma_relation, junk

; original Gebhardt et al. (2000)
sigm=[0.1, 1.0, 10.0, 100.0, 1000.0]
mbh= (1.2e+8)*((sigm/200.0)^(3.75))
oplot, sigm, mbh, psym=-3, linestyle=1, thick=5.0, color=0
oplot, [100.,130.], [1.6e+6,1.6e+6], psym=-3, linestyle=1, thick=5.0, color=0
xyouts, 0.60, 0.3, 'Gebhardt et al. (2000)', color= 0, charsize=0.9, /normal

; original Ferrarese & Merritt (2000)
sigm=[0.1, 1.0, 10.0, 100.0, 1000.0]
mbh= (1.2e+8)*((sigm/200.0)^(4.8))
oplot, sigm, mbh, psym=-3, linestyle=2, thick=2.0, color=0
oplot, [100.,130.], [8.0e+5,8.0e+5], psym=-3, linestyle=2, thick=2.0, color=0
xyouts, 0.60, 0.25, 'Ferrarese & Merritt (2000)', color= 0, charsize=0.9, /normal

; Tremaine et al. (2002)
sigm=[0.1, 1.0, 10.0, 100.0, 1000.0]
mbh= (10^8.13)*((sigm/200.0)^(4.02))
oplot, sigm, mbh, psym=-3, linestyle=0, thick=3.0, color=0
oplot, [100.,130.], [4.0e+5,4.0e+5], psym=-3, linestyle=0, thick=3.0, color=0
xyouts, 0.60, 0.2, 'Tremaine et al. (2002)', color= 0, charsize=0.9, /normal

end



;===================================================================================











