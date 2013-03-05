

; -------------------------------------------------------------------------------------------------




; ----------------------------
;  Read hotgas.txt file
; ----------------------------
pro read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, $
				gas_tot, gas_hot, gas_cold, gas_sf

hgasfile= frun+'/hotgas.txt'

spawn, "wc "+hgasfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then hgas_data= fltarr(9,lines)

openr, 1, hgasfile

junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, hgas_data
close, 1


time= hgas_data[0,*]
temp_keV_X= hgas_data[1,*]
temp_keV= hgas_data[2,*]
temp_K= hgas_data[3,*]
entropy= hgas_data[4,*]
gas_tot= hgas_data[5,*]
gas_hot= hgas_data[6,*]
gas_cold= hgas_data[7,*]
gas_sf= hgas_data[8,*]


end








; -------------------------------------------------------------------------------------------






;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro xraylum_vs_time, junk


if not keyword_set(junk) then begin
	print, " "
	print, " xraylum_vs_time, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='xraylum.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4
;setup_plot_stuff, 'ps', filename=filename, colortable=0



;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "!6Log L!DX!N (ergs s!E-1!N)"
xaxistitle = "!6Time (Gyr)"
;xmax = max(time)
xmax = 4.25
;xmax = 3.0
xmin = 0
ymax = 42.0
;ymax = 41.0
ymin = 37.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; Galaxy 1
;----------
;frun="/raid4/tcox/As/A1"
;frun="/raid4/tcox/vc1vc1"
;frun="pool/vc3vc3"
;frun="/raid4/tcox/vc3bvc3b"
;frun="/raid4/tcox/vc3avc3a_3"
;frun="/raid4/tcox/vc3ba_2_3"
;frun="/raid4/tcox/vc3vc3e_2"
frun="/raid4/tcox/vc3vc3h_2"
;frun="/raid4/tcox/gfs/vc3vc3x2_h"
read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

time= time/0.7
symsize= 1.2
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
;oplot, time, xray, thick=3.0, psym=-8, color= 50
oplot, time, xray, thick=5.0, psym=-3, color= 50
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
oplot, time, xray_rs_s, thick=3.0, psym=-8, color= 50, linestyle= 2
;oplot, time, xray_rs_s, thick=3.0, psym=-3, color= 50, linestyle= 1

; upper right
;oplot, [2.4, 2.6], [41.25,41.25], thick=3.0, psym=-3, color= 50
;xyouts, 0.69, 0.82, 'BH', /normal, charthick=3.0, size=1.33, color=50
; lower right
oplot, [2.4, 2.6], [38.30,38.30], thick=3.0, psym=-3, color= 50
xyouts, 0.69, 0.35, 'BH', /normal, charthick=3.0, size=1.33, color=50

usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
oplot, [2.4, 2.6], [38.65,38.65], thick=3.0, psym=-8, color= 50, linestyle=2
xyouts, 0.69, 0.40, 'BH (RS)', /normal, charthick=3.0, size=1.33, color=50


; Galaxy 2
;----------
;frun="/raid4/tcox/As/A2"
;frun="/raid4/tcox/vc2vc2"
;frun="pool/vc3vc3_wBH"
;frun="/raid4/tcox/vc3bvc3b_no"
;frun="/raid4/tcox/vc3avc3a_2"
;frun="/raid4/tcox/vc3ba_2_2"
;frun="/raid4/tcox/vc3vc3e_no"
frun="/raid4/tcox/vc3vc3h_no"
read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

time= time/0.7
symsize= 1.0
;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
;oplot, time, xray, thick=3.0, psym=-8, color= 150
oplot, time, xray, thick=3.0, psym=-2, color= 150, linestyle= 1
;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
oplot, time, xray_rs_s, thick=3.0, psym=-8, color= 150, linestyle= 0
;oplot, time, xray_rs_s, thick=3.0, psym=-3, color= 150, linestyle= 1

oplot, [2.4, 2.6], [38.00,38.00], thick=3.0, psym=-8, color= 150
xyouts, 0.69, 0.30, 'std (RS)', /normal, charthick=3.0, size=1.33, color=150

oplot, [2.4, 2.6], [37.68,37.68], thick=3.0, psym=-2, color= 150, linestyle= 1
xyouts, 0.69, 0.25, 'std', /normal, charthick=3.0, size=1.33, color=150


; Galaxy 3
;----------
;frun="/raid4/tcox/As/A3"
;frun="/raid4/tcox/vc3vc3"
;frun="/raid4/tcox/vc3avc3a_1"
;frun="/raid4/tcox/vc3ba_2_1"
;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
;;read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
;oplot, time/0.7, xray, thick=3.0, psym=-6, color= 220
;oplot, [2.2, 2.4], [41.08,41.08], thick=3.0, psym=-6, color= 150
;xyouts, 0.82, 0.77, '1e5', /normal, charthick=3.0, size=1.33, color=220


; Galaxy 4
;----------
;frun="/raid4/tcox/As/A4"
;frun="/raid4/tcox/vc4vc4a"
;frun="/raid4/tcox/vc3avc3a_4"
;frun="/raid4/tcox/vc3ba_2_4"
;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
;;read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
;oplot, time/0.7, xray, thick=3.0, psym=-6, color= 100
;oplot, [2.2, 2.4], [41.08,41.08], thick=3.0, psym=-6, color= 150
;xyouts, 0.82, 0.72, '1e6', /normal, charthick=3.0, size=1.33, color=100


; Galaxy 5
;----------
;frun="/raid4/tcox/As/A5"
;frun="/raid4/tcox/vc5vc5a"
;frun="/raid4/tcox/vc3avc3a_6"
;frun="/raid4/tcox/vc3ba_2_5"
;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
;;read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
;oplot, time/0.7, xray, thick=3.0, psym=-6, color= 0
;oplot, [2.2, 2.4], [41.08,41.08], thick=3.0, psym=-6, color= 150
;xyouts, 0.82, 0.67, '1e7', /normal, charthick=3.0, size=1.33, color=0


; Galaxy 6
;----------
;frun="/raid4/tcox/vc6vc6a"
;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
;;read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
;oplot, time/0.7, xray, thick=3.0, psym=-7, color= 180
;oplot, [2.2, 2.4], [41.08,41.08], thick=3.0, psym=-6, color= 150
;xyouts, 0.82, 0.67, '1e7', /normal, charthick=3.0, size=1.33, color=0




; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0

; draw arrow for merger time
timemerge=1.05/0.7
arrow, timemerge, 37.2, timemerge, 37.9, COLOR=0, THICK=3.0, hthick=3.0, /data

device, /close




end




;================================================================================





;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro xraylum_many_vs_time, junk


if not keyword_set(junk) then begin
	print, " "
	print, " xraylum_vs_time, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='xraylum.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=4
setup_plot_stuff, 'ps', filename=filename, colortable=1
;setup_plot_stuff, 'ps', filename=filename, colortable=0



;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "!6Log L!DX!N (ergs s!E-1!N)"
xaxistitle = "!6Time (Gyr)"
;xmax = max(time)
;xmax = 7.25
xmax = 4.25
;xmax = 3.0
xmin = 0
ymax = 44.0
;ymax = 41.0
ymin = 36.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; Gas Fractions
;---------------
do_gfs= 0
if do_gfs eq 1 then begin
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3u_e', linecolor=20
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3v_e', linecolor=60
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3w_e', linecolor=100
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3x2_e', linecolor=140
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3y_e', linecolor=180
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3z_e', linecolor=220
	; ------
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3u_h', linecolor=20
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3v_h', linecolor=60
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3w_h', linecolor=100
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3x2_h', linecolor=140
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3y_h', linecolor=180
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3z_h', linecolor=220
	; ------
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3u_k', linecolor=220, lthick= 12.0
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3v_k', linecolor=180, lthick= 10.0
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3w_k', linecolor=140, lthick= 8.0
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3x2_k', linecolor=100, lthick= 6.0
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3y_k', linecolor=60, lthick= 4.0
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3z_k', linecolor=20, lthick= 3.0

	;process_one_sim, '/raid4/tcox/zs/z3h', linecolor=20
	;process_one_sim, '/raid4/tcox/es/e3h', linecolor=80
	;process_one_sim, '/raid4/tcox/ds/d3h7', linecolor=120
	;process_one_sim, '/raid4/tcox/bs/b3h', linecolor=160
	;process_one_sim, '/raid4/tcox/As/A3', linecolor=200
	process_one_sim, '/raid4/tcox/zs/z3e', linecolor=50, lthick= 3.0
	process_one_sim, '/raid4/tcox/es/e3e', linecolor=100, lthick= 5.0
	process_one_sim, '/raid4/tcox/ds/d3e7', linecolor=145, lthick= 7.0
	process_one_sim, '/raid4/tcox/vc3vc3e_2', linecolor=150, lthick= 9.0
	process_one_sim, '/raid4/tcox/bs/b3e', linecolor=200, lthick= 12.0

	xyouts, 0.24, 0.86, "Gas Fraction", /normal, charthick=4, size=1.53, color=0
endif


; Masses
; --------
do_ms= 0
if do_ms eq 1 then begin
	;process_one_sim, '/raid4/tcox/ds/d0k2_q', linecolor=20, lthick= 3.0
	;process_one_sim, '/raid4/tcox/ds/d1k2_q', linecolor=50, lthick= 5.0
	;process_one_sim, '/raid4/tcox/ds/d2k2_q', linecolor=80, lthick= 6.0
	;process_one_sim, '/raid4/tcox/ds/d3k7',   linecolor=110, lthick= 7.0
	;process_one_sim, '/raid4/tcox/ds/d4k2_q', linecolor=140, lthick= 8.0
	;process_one_sim, '/raid4/tcox/ds/d5k2_q', linecolor=170, lthick= 9.0
	;process_one_sim, '/raid4/tcox/ds/d6k2_q', linecolor=200, lthick= 10.0

	;process_one_sim, '/raid4/tcox/zs/z0e', linecolor=20
	;process_one_sim, '/raid4/tcox/zs/z1e', linecolor=40
	;process_one_sim, '/raid4/tcox/zs/z2e', linecolor=60
	;process_one_sim, '/raid4/tcox/zs/z3e', linecolor=80
	;process_one_sim, '/raid4/tcox/zs/z4e', linecolor=100
	;process_one_sim, '/raid4/tcox/zs/z5e', linecolor=120
	;process_one_sim, '/raid4/tcox/zs/z6e', linecolor=150

	process_one_sim, '/raid4/tcox/bs/b0e', linecolor=20, lthick= 3.0
	process_one_sim, '/raid4/tcox/bs/b1e', linecolor=50, lthick= 5.0
	process_one_sim, '/raid4/tcox/bs/b2e', linecolor=80, lthick= 6.0
	process_one_sim, '/raid4/tcox/bs/b3e', linecolor=110, lthick= 7.0
	process_one_sim, '/raid4/tcox/bs/b4e', linecolor=140, lthick= 8.0
	process_one_sim, '/raid4/tcox/bs/b5e', linecolor=170, lthick= 9.0
	process_one_sim, '/raid4/tcox/bs/b6e', linecolor=200, lthick= 10.0

	xyouts, 0.54, 0.89, "Progenitor Mass", /normal, charthick=4, size=1.53, color=0
endif



; Orientations
; ------------
do_ori= 1
if do_ori eq 1 then begin
	process_one_sim, '/raid4/tcox/ds/d3e7', linecolor=50, lthick= 4.0
	process_one_sim, '/raid4/tcox/ds/d3h7', linecolor=100, lthick= 4.0
	process_one_sim, '/raid4/tcox/ds/d3k7', linecolor=150, lthick= 4.0
	process_one_sim, '/raid4/tcox/ds/d3f7', linecolor=200, lthick= 4.0

	;process_one_sim, '/raid4/tcox/gfs/vc3vc3x2_e', linecolor=50, lstyle= 1, lthick= 6.0
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3x2_h', linecolor=100, lstyle= 1, lthick= 6.0
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3x2_k', linecolor=150, lstyle= 1, lthick= 6.0
	;process_one_sim, '/raid4/tcox/gfs/vc3vc3_f', linecolor=200, lstyle= 1, lthick= 6.0

        xyouts, 0.24, 0.86, "Disk Orientations", /normal, charthick=4, size=1.53, color=0
endif



; Orbits
; ------------
do_orb= 0
if do_orb eq 1 then begin
	process_one_sim, '/raid4/tcox/ds/d3k1', linecolor=20, lthick= 1.0
	process_one_sim, '/raid4/tcox/ds/d3k7', linecolor=50, lthick= 3.0
	process_one_sim, '/raid4/tcox/ds/d3k2', linecolor=80, lthick= 5.0
	process_one_sim, '/raid4/tcox/ds/d3k3', linecolor=110, lthick= 7.0
	process_one_sim, '/raid4/tcox/ds/d3k4', linecolor=140, lthick= 9.0
	process_one_sim, '/raid4/tcox/ds/d3k6', linecolor=170, lthick= 11.0
	process_one_sim, '/raid4/tcox/ds/d3k5', linecolor=200, lthick= 13.0

	xyouts, 0.24, 0.86, "Orbital Angular Momentum", /normal, charthick=4, size=1.53, color=0
	;oplot, [1.65,2.2],[13.0,13.0], psym=-3, color= 0, linestyle= 0
	;xyouts, 0.7, 0.80, '4', /normal, charthick=1, size=1.33, color=50
	;xyouts, 0.7, 0.76, '2', /normal, charthick=1, size=1.33, color=75
	;xyouts, 0.7, 0.72, '1', /normal, charthick=1, size=1.33, color=100
	;xyouts, 0.7, 0.68, '1/2', /normal, charthick=1, size=1.33, color=125
	;xyouts, 0.7, 0.64, '1/4', /normal, charthick=1, size=1.33, color=150
	;xyouts, 0.7, 0.60, '1/8', /normal, charthick=1, size=1.33, color=175

endif






; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0

; draw arrow for merger time
;timemerge=1.05/0.7
;arrow, timemerge, 37.2, timemerge, 37.9, COLOR=0, THICK=3.0, hthick=3.0, /data

device, /close




end









pro process_one_sim, frun, linecolor=linecolor, $
			lthick=lthick, $
			lstyle=lstyle

if not keyword_set(lthick) then lthick= 5.0
if not keyword_set(lstyle) then lstyle= 0

read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
time= time/0.7

;xrays_to_plot= 30.0 + alog10(0.5*(10^(xray_rs_s-30.0) + 10^(xray-30.0)))
xrays_to_plot= xray_rs_s
;xrays_to_plot= xray

oplot, time, xrays_to_plot, thick=lthick, psym=-3, color= linecolor, linestyle= lstyle


end







;================================================================================





;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro xraydivB_vs_time, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_xraydivB_vs_time, junk"
	print, " "
	print, " "
	return
endif

filename='xrayB.eps'
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


yaxistitle = "!6Log L!DX!N / L!DB!N (erg s!e-1!n L!DB,!9n!3!N!E-1!N)"
xaxistitle = "!6Time (Gyr)"

;xmax = max(time)
xmax = 4.25
xmin = 0
ymax = 31.0
ymin = 26.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; Galaxy 1
;----------
;frun="pool/vc3vc3"
;frun="pool/vc3bvc3b"
;frun="/raid4/tcox/vc3bvc3b"
frun="/raid4/tcox/vc3vc3h_2"
read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

read_colors_file, frun, time,  mags
Mb= mags[2,*]
Lb= -0.4*(Mb-5.51)

time= time/0.7

; xray and Lb are both in log here
;xrayB= xray - Lb
xrayB= xray_rs_s - Lb

; make it log
;idx= where(xrayB le 0)
;if idx(0) eq -1 then begin
;	xrayB= alog10(xrayB)
;endif else begin
;	xrayB(idx)= 1e10
;	xrayB= alog10(xrayB)
;	xrayB(idx)= 0.0
;endelse

symsize= 1.2
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
oplot, time, xrayB, thick=3.0, psym=-8, color= 50, linestyle= 2

;;xyouts, 0.65, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=50
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
;oplot, [2.2, 2.4], [41.25,41.25], thick=3.0, psym=-8, color= 50
;;oplot, [2.2, 2.4], [40.65,40.65], thick=3.0, psym=-8, color= 50
;xyouts, 0.65, 0.82, 'black hole', /normal, charthick=3.0, size=1.33, color=50
;xyouts, 0.82, 0.87, '1e3', /normal, charthick=3.0, size=1.33, color=50

usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
oplot, [2.4, 2.6], [27.00,27.00], thick=3.0, psym=-8, color= 50, linestyle=2
xyouts, 0.69, 0.30, 'BH (RS)', /normal, charthick=3.0, size=1.33, color=50



; Galaxy 2
;----------
;frun="pool/vc3vc3_wBH"
;frun="pool/vc3bvc3b_wBH"
;frun="/raid4/tcox/vc3bvc3b_no"
frun="/raid4/tcox/vc3vc3h_no"
read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

read_colors_file, frun, time,  mags
Mb= mags[2,*]
Lb= -0.4*(Mb-5.51)

time= time/0.7

;xrayB= xray-Lb
xrayB= xray_rs_s - Lb


symsize= 1.0
;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
oplot, time, xrayB, thick=3.0, psym=-8, color= 150
;oplot, time, xrayB, thick=3.0, psym=-2, color= 150

oplot, [2.4, 2.6], [26.70,26.70], thick=3.0, psym=-8, color= 150
xyouts, 0.69, 0.25, 'std (RS)', /normal, charthick=3.0, size=1.33, color=150


;;xyouts, 0.7, 0.85, fload_getid(frun), /normal, charthick=1, size=1.33, color=150



; --------------------------
; draw arrow for merger time
timemerge=1.05/0.7
arrow, timemerge, 26.4, timemerge, 27.1, COLOR=0, THICK=3.0, hthick=3.0, /data


; --------------------------
; OFP '01 slope
x= [2.5,3.0,3.5]
y= 28.4 + 0.063*x
oplot, x, y, psym=-3, thick=5.0, color= 0
xyouts, 0.65, 0.49, 'OFP01', /normal, charthick=3.0, size=1.33, color=0



device, /close




end









;----------------------------------------------------------------------------------------


