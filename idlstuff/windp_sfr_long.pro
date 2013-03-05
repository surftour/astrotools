

; --------------------------------
;  Read Bolometric Lum File
; ----------------------------------
pro read_bololum_file, frun, time, sbololum, bhbololum, $
                bolofile=bolofile

if not keyword_set(bolofile) then bolofile= frun+'/lum_bolo.txt'

spawn, "wc "+bolofile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then re_data= fltarr(3,lines)

openr, 1, bolofile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, re_data
close, 1


time= re_data[0,*]
sbololum= re_data[1,*]
bhbololum= re_data[2,*]


end




;===================================================================






pro wsfr, junk, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		h=h



if not keyword_set(junk) then begin
   print, "  "
   print, "wsfr, junk, filename=filename, /h"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='wsfr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize= 23.0, newysize= 5.5
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
if keyword_set(cumulative) then yaxistitle="!6Mass (10!e10!n h!E-1!NM!D!9n!6!N)"
if keyword_set(gasmass) then yaxistitle="!6Gas Mass (10!e10!n h!E-1!NM!D!9n!6!N)"

;xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
xaxistitle = "!6Time (Gyr)"  & h= 0.7

;xmax = 14.0
;xmax = 12.0
;xmax = 10.0
;xmax = 7.5
;xmax = 6.0
;xmax = 5.0
;xmax = 4.25
;xmax = 4.0
xmax = 3.1
;xmax = 3.0
;xmax = 2.8
;xmax = 2.0
;xmax = 1.5
;xmax = 1.3
;xmax = 1.0
;xmax = 0.5
xmin = 0

;ymax = 8000
;ymax = 2500
;ymax = 1750
;ymax = 1500
;ymax = 1000
;ymax = 400
;ymax = 300
;ymax = 250.0
ymax = 220
;ymax = 180
;ymax = 120
;ymax = 100
;ymax = 90
;ymax = 80
;ymax = 60
;ymax = 50.0
;ymax = 20.0
;ymax = 16.0
;ymax = 15
;ymax = 13
;ymax = 12
;ymax = 10.0
;ymax = 8.0
;ymax = 6.0
;ymax = 5.0
;ymax = 3.5
;ymax = 2.0
;ymax = 1.5
;ymax = 1.3
;ymax = 1.2
;ymax = 0.75
;ymax = 0.6
;ymax = 0.5
;ymin = 0 
;ymin = 0.1
;ymin= 0.01
ymin= 0.006
;ymin= 0.001
;ymin= 0.0001
;ymin= 0.00001


; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;---------------------------

;!p.position= [0.09, 0.30, 0.99, 0.98]
!p.position= [0.09, 0.30, 0.91, 0.98]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        ;xstyle=1, ystyle=1, $
        xstyle=1, ystyle=9, $
	xticklen= 0.08, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=6.0, $
	ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata





;-------------------------------------------
;   Load the runs to display
;-------------------------------------------


;xyouts, 0.70, 0.92, 'Log(M!Dtot!N/M!D!9n!6!N)', /normal, charthick=3, size=1.33, color=0
process_one_sfr, "sbw/sb10BH", lcolor=150, lthick= 4.0, msg=' ', x0= 0.80, y0= 0.63, h=h, lpsym=-3




; ---------------------------------------------------



        ;-------------------------------------------
        ;   Get SFR rate from txt - for each file
        ;-------------------------------------------
        open_sfr_file, "sbw/sb10BH", sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass


	sfrtime = sfrtime / h

	ulirg_time0= 1.155 / h
	ulirg_time1= 1.25 / h
	;ulirg_time1= 1.35 / h

	idx= where((sfrtime ge ulirg_time0) and (sfrtime le ulirg_time1))

	ulirg_sfrtime= sfrtime(idx)
	ulirg_sfrsfr= sfrsfr(idx)
	nidx= n_elements(idx)-1

	xrange= [ulirg_sfrtime(0),ulirg_sfrtime,ulirg_sfrtime(nidx),ulirg_sfrtime(0)]
	;yrange= [ulirg_sfrsfr(0),ulirg_sfrsfr,ulirg_sfrsfr(nidx),ulirg_sfrsfr(0)]
	;yrange= [1.0e-5,ulirg_sfrsfr,1.0e-5,1.0e-5]
	yrange= [ymin,ulirg_sfrsfr,ymin,ymin]

	polyfill, xrange, yrange, /data, color= 150, /fill



; ---------------------------------------------------



; std
std= 0
;std= 1
if std eq 1 then begin
    x0= 0.70
    ;x0= 0.65
    ;x0= 0.55
    ;y0= 0.40
    ;x0= 0.35
    y0= 0.92
    zcolor= zcolor_orig
    for i=0, n_elements(fruns)-1 do begin
        y0= y0-0.04
        if n_elements(lbls) gt 0 then begin
                xyouts, x0, y0, lbls[i], /normal, charthick=3, size=1.33, color=zcolor
                    zcolor= zcolor+deltacolor
        endif else begin
                xyouts, x0, y0, fruns[i], /normal, charthick=3, size=1.33, color=zcolor
                    zcolor= zcolor+deltacolor
        endelse
    endfor
endif





; ---------------------------------------------------


add_right_axis= 1
if add_right_axis eq 1 then begin

	read_bololum_file, '/raid4/tcox/sbw/sb10BH', time, sbololum, bhbololum, bolofile=bolofile
	;S_bololum= sbololum + 33.591
	;BH_bololum= bhbololum + 33.591
	S_bololum= 10.0^(sbololum) * 3.99d+33
	BH_bololum= 10.0^(bhbololum) * 3.99d+33
	;bololum= alog10(10^(sbololum) + 10^(bhbololum))

	time= time/0.7


	;---------------------------
	;yaxistitle="!6Log L!DBolo!N (ergs s!E-1!N)"
	yaxistitle="!6L!DBolo!N (ergs s!E-1!N)"
	;ymax= max(bololum)
	;ymin= min(bololum)
	;ymax= 45.5
	;ymin= 41.5
	ymax= 3.1d+45
	ymin= 3.1d+41

	!p.position= [0.09, 0.30, 0.91, 0.98]

	plot, [1.0],[1.0], psym=-3, $
	        xrange=[xmin,xmax], yrange=[ymin,ymax], $
	        color= 0, $
	        /ylog, $
	        xstyle=5, ystyle=5, $
	        xticklen= 0.05, $
	        xcharsize=1.50, ycharsize=1.50, $
	        xthick=4.0, ythick=4.0, $
	        charthick=6.0, $
	        ytickformat='exp_label', $
	        ;xtitle=xaxistitle, $
	        ;ytitle=yaxistitle, $
	        /nodata, /noerase

	axis, yaxis=1, xstyle=1,ystyle=1,color=0,charsize=1.5, xthick=4.0,ythick=4.0, $
		ytitle= yaxistitle, charthick= 6.0, /ylog, ytickformat= 'exp_label'

        oplot, time, S_bololum, thick=4.0, color= 0, psym=-3, linestyle=1
        oplot, time, BH_bololum, thick=4.0, color= 0, psym=-3, linestyle=2



endif




; ---------------------------------------------------








;--------------------------------------
;--------------------------------------

device, /close


end

    
















