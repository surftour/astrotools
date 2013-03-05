

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
pro wmsig, junk


if not keyword_set(junk) then begin
	print, " "
	print, " wmsig, junk"
	print, " "
	print, " needs:  .run bh_multi"
	print, "         .run time_m_sigma_re"
	print, " "
	print, " "
	return
endif

filename='msigma.eps'

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
xyouts, xpt3, ylbl*0.85, 'no SB wind', size=1.2, color=thiscolor, /data, charthick=4.0

frun= 'ds/d0e2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d1e2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d2e2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d3e7'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'vc3vc3e_2'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d4e2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d5e2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d6e2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor

frun= 'bs/b0e'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b1e'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b2e'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b3e'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b4e'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b4e'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b6e'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor



; -----------------------------------



; eta= 0.05  (sb7, sb19, sb18, sb9)
; ------------------
ylbl= 2.0d+9
pointt= 2
select_thispoint, pointt, thispsym, thiscolor

;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl*0.85, '!7g!6= 0.05', size=1.2, color=thiscolor, /data, charthick=4.0
oplot, [xpt2], [ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl*0.85, '!7g!6= 0.05', size=1.2, color=thiscolor, /data, charthick=4.0

; sb7BH
;sig= [157.3]
;bh= [0.0433465]
;oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=0.05, !8v!6=132", /data, charthick=3.0, size=1, color= 50

frun= 'sbw/sb19BH'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=0.05, !8v!6=209", /data, charthick=3.0, size=1, color= 50

frun= 'sbw/sb18BH'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=0.05, !8v!6=418", /data, charthick=3.0, size=1, color= 50

frun= 'sbw/sb9BH'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=0.05, !8v!6=837", /data, charthick=3.0, size=1, color= 50



; kind of cheating, but will look good - and it is correct
frun= 'ds/d0h2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d1h2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d2h2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d3h7'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
;frun= 'vc3vc3h_2'
frun= 'vc3vc3h'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d4h2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d5h2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'ds/d6h2_q'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor

frun= 'bs/b0h'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b1h'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b2h'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b3h'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b4h'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b4h'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor
frun= 'bs/b6h'
oplot, get_sigma(frun), get_bh_mass(frun) * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color=thiscolor







; -----------------------------------



; eta= 0.5   (sb17, sb16, sb15, sb10)
; ------------------
ylbl= 1.1d+9
pointt= 3
select_thispoint, pointt, thispsym, thiscolor

;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl*0.85, '!7g!6= 0.5', size=1.2, color=thiscolor, /data, charthick=4.0
oplot, [xpt2], [ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl*0.85, '!7g!6= 0.5', size=1.2, color=thiscolor, /data, charthick=4.0

; sb17BH
; is also too big
;bh= [0.159028]

; sb16BH
;sig= [165.2]
;bh= [0.0880842]
;oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=0.5, !8v!6=209", /data, charthick=3.0, size=1, color= 50

; sb15BH
;sig= [148.8]
;bh= [0.0217242]
;oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=0.5, !8v!6=418", /data, charthick=3.0, size=1, color= 50

frun= 'sbw/sb10BH'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
;xyouts, sig+10, bh * 1e+10 / 0.7, "sb10", /data, charthick=3.0, size=1, color= 150
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=0.5, !8v!6=837", /data, charthick=3.0, size=1, color= 50

frun= 'sbw/sb10BHtr1'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

frun= 'sbw/sb10BHtr2'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

frun= 'sbw/sb10BHtr2I'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

frun= 'sbw/sb10BHtr3'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor



;different masses
frun= 'sb10_mass/b0e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/b1e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/b2e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/b3e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/b4e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/b5e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/b6e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor



;different masses
frun= 'sb10_mass/d0e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/d1e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/d2e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/d4e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/d5e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb10_mass/d6e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor






; -----------------------------------






; eta= 2.0   (sb14, sb13, sb12, sb8)
; ------------------
ylbl= 6.0d+8
pointt= 4
select_thispoint, pointt, thispsym, thiscolor

;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl*0.85, '!7g!6= 2.0', size=1.2, color=thiscolor, /data, charthick=4.0
oplot, [xpt2], [ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl*0.85, '!7g!6= 2.0', size=1.2, color=thiscolor, /data, charthick=4.0



; sb14BH
;sig= [127.8]
;bh= [0.137036]
;oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
;xyouts, sig+10, bh * 1e+10 / 0.7, "sb13", /data, charthick=3.0, size=1, color= 150
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=2.0, !8v!6=105", /data, charthick=3.0, size=1, color= 50

; sb13BH
;sig= [158.8]
;bh= [0.158933]
;oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=2.0, !8v!6=209", /data, charthick=3.0, size=1, color= 50

frun= 'sbw/sb13BHtr1'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

frun= 'sbw/sb13BHtr2'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

frun= 'sbw/sb13BHtr3'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

; different masses
frun= 'sb13_mass/b0e'




frun= 'sbw/sb12BH'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=2.0, !8v!6=418", /data, charthick=3.0, size=1, color= 50

frun= 'sbw/sb8BH'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
;xyouts, sig+10, bh * 1e+10 / 0.7, "sb8", /data, charthick=3.0, size=1, color= 150
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=2.0, !8v!6=837", /data, charthick=3.0, size=1, color= 50

frun= 'sbw/sb8BHtr0'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

frun= 'sbw/sb8BHtr1'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

frun= 'sbw/sb8BHtr2'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

frun= 'sbw/sb8BHtr3'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor



; different masses
frun= 'sb8_mass/b0e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb8_mass/b1e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb8_mass/b2e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb8_mass/b4e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb8_mass/b5e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb8_mass/b6e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor


; different masses
frun= 'sb8_mass/d0e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb8_mass/d1e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb8_mass/d2e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb8_mass/d4e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
;frun= 'sb8_mass/d5e'
;bh= get_bh_mass(frun)
;sig= get_sigma(frun)
;oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
frun= 'sb8_mass/d6e'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor







; -----------------------------------




; eta= 5.0    (sb22, sb21, sb20, sb11)
; ------------------
ylbl= 3.2d+8
pointt= 5
select_thispoint, pointt, thispsym, thiscolor

;oplot, [xpt1,xpt2], [ylbl,ylbl], psym=-thispsym, color=thiscolor, thick=3.0
;xyouts, xpt3, ylbl*0.85, '!7g!6= 5.0', size=1.2, color=thiscolor, /data, charthick=4.0
oplot, [xpt2], [ylbl], psym=-thispsym, color=thiscolor, thick=3.0
xyouts, xpt3, ylbl*0.85, '!7g!6= 5.0', size=1.2, color=thiscolor, /data, charthick=4.0



; sb22BH
;sig= [121.3]
;bh= [0.294031]
;oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

; sb22BHtr1
;sig= [152.1]
;bh= [0.0354629]
;oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

frun= 'sbw/sb22BHtr2'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

frun= 'sbw/sb22BHtr3'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor

; sb21BH
;sig= [130.0]
;bh= [0.0817862]
;oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=5.0, !8v!6=209", /data, charthick=3.0, size=1, color= 50

frun= 'sbw/sb20BH'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=5.0, !8v!6=418", /data, charthick=3.0, size=1, color= 50

frun= 'sbw/sb11BH'
bh= get_bh_mass(frun)
sig= get_sigma(frun)
oplot, sig, bh * 1e+10 / 0.7, thick=3.0, psym=-thispsym, color= thiscolor
;xyouts, sig+10, 0.9 * bh * 1e+10 / 0.7, "!7g!6=5.0, !8v!6=837", /data, charthick=3.0, size=1, color= 50





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











