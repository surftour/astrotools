pro sfr_multi, junk, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		ynotlog=ynotlog, $
		h=h


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_multi, junk, filename=filename, /h, "
   print, "           /cumulative, /gasmass, /ynotlog"
   print, "  "
   print, "  "
   print, "  WARNING: cumulative and gasmass do not currently work!  "
   print, "  "
   print, "  default filename: sfr.eps"
   print, "  "
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'
if keyword_set(gasmass) then filename='gmass.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
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
;xmax = 3.5
;xmax = 3.0
;xmax = 2.8
;xmax = 2.5
;xmax = 2.4
xmax = 2.2
;xmax = 2.0
;xmax = 1.5
;xmax = 1.4
;xmax = 1.3
;xmax = 1.0
;xmax = 0.5
xmin = 0

;ymax= 8e+5
;ymax= 8e+4
;ymax= 2e+4
;ymax = 8000
;ymax = 4000
;ymax = 2500
;ymax = 1750
;ymax = 1500
ymax = 1000
;ymax = 800
;ymax = 600
;ymax = 400
;ymax = 300
;ymax = 250.0
;ymax = 202
;ymax = 180
;ymax = 120
;ymax = 100
;ymax = 90
;ymax = 80
;ymax = 60
;ymax = 50.0
;ymax = 40.0
;ymax = 30.0
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
;ymin= 8.0
;ymin= 1.0
;ymin = 0  & ynotlog= 1
;ymin = 0.8
;ymin = 0.4
;ymin = 0.1
;ymin= 0.07
ymin= 0.01
;ymin= 0.001
;ymin= 0.0001
;ymin= 0.00001

if keyword_set(ynotlog) then ymin = 0 

; physical units
if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	;h = fload_cosmology('h')
	;xmax = xmax / h
	;xmin = xmin / h
endif


;---------------------------

!p.position= [0.16, 0.14, 0.98, 0.98]

if keyword_set(ynotlog) then begin
    plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	;ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata
endif else begin
    plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata
endelse




;-------------------------------------------
;   Load the runs to display
;-------------------------------------------

; -------------------------------------
;    VC Galaxies
; -------------------------------------

; v1
; -----
;fruns= ["vc1vc1","vc1"]
;       tmerg=[1.1,1.1]
;       lbls= ['vc1 Major Merger','vc1 Isolated']
;       msg= ' '

;fruns= ["vc1vc1","vc1vc1_no", "vc1vc1_wBH"]
;       tmerg=[1.1,1.1,1.1]
;fruns= ["vc1vc1", "vc1vc1_wBH"]
;       tmerg=[1.1,1.1]
       ;lbls= ['vc1 Major Merger','vc1 Isolated']
       ;lbls= ['vc1vc1 Major Merger','vc1 Isolated']
;       msg= ' '

; v2
; -----
;fruns= ["vc2vc2","vc2vc2_no", "vc2vc2_wBH"]
;       tmerg=[1.1,1.1,1.1]
;fruns= ["vc2vc2","vc2vc2_wBH"]
;       tmerg=[1.1,1.1]
;       msg= ' '


; v3
; -----
;fruns= ["vc3vc3", "vc3vc3a"]             ; h_dm = 0.4 and 0.2, respectively
;fruns= ["vc3vc3e", "vc3vc3e_2"]
;fruns= ["vc3vc3e", "vc3vc3h"]
;fruns= ["vc3vc3e", "vc3vc3e_2","vc3vc3e_no"]
;fruns= ["vc3vc3e_2","vc3vc3e_no"] & lbls= ["black hole","no black hole"]
;fruns= ["vc3vc3", "vc3vc3_wBH"]
;fruns= ["vc3bvc3b", "vc3bvc3b_wBH"]
;fruns= ["vc3vc3_wBH", "vc3vc3_sauron"]
;	tmreg=[1.1,1.1,1.1]
;	lbls=['vc3vc3','vc3vc3_wBH']
;	lbls=['mako','sauron']
	msg=''

;fruns= ["vc3bvc3b", "vc3bvc3b_wBH"]
;       tmreg=[1.1,1.1,1.1]
;       lbls=['vc3bvc3b','vc3bvc3b_w']
;       msg=''



; vc3 - orbits
; --------------
;fruns= ["vc3vc3b", "vc3vc3c","vc3vc3d","vc3vc3e","vc3vc3f","vc3vc3g","vc3vc3h"]
;fruns= ["vc3vc3i", "vc3vc3j","vc3vc3k","vc3vc3l","vc3vc3m","vc3vc3n","vc3vc3o","vc3vc3p"]
;	tmreg=[1.1,1.1,1.1,1.1,1.1,1.1,1.1]
;	msg=''


;fruns= ["vc3vc3d","vc3vc3i","vc3vc3l","vc3vc3n"] ; embedded disks
;fruns= ["vc3vc3k","vc3vc3m","vc3vc3p"]  ; multiple components
;fruns= ["vc3vc3h","vc3vc3b", "vc3vc3c"]  ; S0's
;fruns= ["vc3vc3e","vc3vc3j"]  ;non-rotating
;fruns= ["vc3vc3f","vc3vc3g","vc3vc3o"]  ; rotating E's

;--------------------------------------------------------------------------------

;process_and_plot_sfr_history, "ds/vc3vc3b", lcolor=40,  lthick= 2.0, msg='b', x0= 0.85, y0= 0.78
;process_and_plot_sfr_history, "ds/vc3vc3c", lcolor=80,  lthick= 2.0, msg='c', x0= 0.85, y0= 0.75
;process_and_plot_sfr_history, "ds/vc3vc3d", lcolor=120, lthick= 2.0, msg='d', x0= 0.85, y0= 0.72
;process_and_plot_sfr_history, "ds/vc3vc3e", lcolor=160, lthick= 2.0, msg='e', x0= 0.85, y0= 0.69
;process_and_plot_sfr_history, "ds/vc3vc3f", lcolor=200, lthick= 2.0, msg='f', x0= 0.85, y0= 0.66

;process_and_plot_sfr_history, "ds/vc3vc3g", lcolor=40,  lthick= 2.0, msg='g', x0= 0.85, y0= 0.78
;process_and_plot_sfr_history, "ds/vc3vc3h", lcolor=80,  lthick= 2.0, msg='h', x0= 0.85, y0= 0.75
;process_and_plot_sfr_history, "ds/vc3vc3i", lcolor=120, lthick= 2.0, msg='i', x0= 0.85, y0= 0.72
;process_and_plot_sfr_history, "ds/vc3vc3j", lcolor=160, lthick= 2.0, msg='j', x0= 0.85, y0= 0.69
;process_and_plot_sfr_history, "ds/vc3vc3k", lcolor=200, lthick= 2.0, msg='k', x0= 0.85, y0= 0.66

;process_and_plot_sfr_history, "ds/vc3vc3l", lcolor=40,  lthick= 2.0, msg='l', x0= 0.85, y0= 0.78
;process_and_plot_sfr_history, "ds/vc3vc3m", lcolor=80,  lthick= 2.0, msg='m', x0= 0.85, y0= 0.75
;process_and_plot_sfr_history, "ds/vc3vc3n", lcolor=120, lthick= 2.0, msg='n', x0= 0.85, y0= 0.72
;process_and_plot_sfr_history, "ds/vc3vc3o", lcolor=160, lthick= 2.0, msg='o', x0= 0.85, y0= 0.69
;process_and_plot_sfr_history, "ds/vc3vc3p", lcolor=200, lthick= 2.0, msg='p', x0= 0.85, y0= 0.66

;--------------------------------------------------------------------------------
;process_and_plot_sfr_history, "ds/vc3vc3e", lcolor=150, lthick= 2.0, msg='e', x0= 0.85, y0= 0.85
;process_and_plot_sfr_history, "ds/vc3vc3h", lcolor=0, lthick= 2.0, msg='h', x0= 0.85, y0= 0.80
;--------------------------------------------------------------------------------

;process_and_plot_sfr_history, "vc3vc3l", lcolor=50, lthick= 2.0, msg='l', x0= 0.25, y0= 0.32
;process_and_plot_sfr_history, "vc3vc3l_bg", lcolor=150, lthick= 2.0, msg='l_bg', x0= 0.25, y0= 0.28

;fruns= ["vc3vc3e"]
;lbls=[" "]
;tmerg= 1.1
;msg= ' '

;fruns= ["vc3vc3e","vc3vc3h"]
;msg= ' '


; vc3 - gas fractions
; --------------------
;fruns=["vc3vc3","vc3vc3h"]
;fruns= ["vc3vc3e","gfs/vc3vc3x2_e","ds/d3e7"]
;fruns= ["vc3vc3h","gfs/vc3vc3x2_h","ds/d3h7"]
;fruns= ["vc3vc3k","gfs/vc3vc3x2_k","ds/d3k7"]

;fruns= ["gfs/vc3vc3u_e","gfs/vc3vc3v_e","gfs/vc3vc3w_e","gfs/vc3vc3x2_k","gfs/vc3vc3y_e","gfs/vc3vc3z_e"]
;fruns= ["gfs/vc3vc3u_h","gfs/vc3vc3v_h","gfs/vc3vc3w_h","gfs/vc3vc3x2_k","gfs/vc3vc3y_h","gfs/vc3vc3z_h"]
;fruns= ["gfs/vc3vc3u_k","gfs/vc3vc3v_k","gfs/vc3vc3w_k","gfs/vc3vc3x2_k","gfs/vc3vc3y_k","gfs/vc3vc3z_k"]
;lbls=['100%','80%','60%','40%','20%','5%']

;fruns= ["gfs/vc3vc3w_e","gfs/vc3vc3x2_e","gfs/vc3vc3z_e"]


;fruns=["vc3vc3","vc3vc3h"]
;fruns=["vc3bvc3b","As/A3duprr1"]
;fruns=["vc3vc3e","vc3vc3e_gf6","vc3vc3e_gf8","vc3vc3e_gf8a"]

;process_and_plot_sfr_history, "bs/b3e", lcolor=50, lthick= 2.0, msg='b3e', x0= 0.25, y0= 0.32
;process_and_plot_sfr_history, "vc3vc3e_gf8", lcolor=150, lthick= 2.0, msg='vc3vc3e_gf8', x0= 0.25, y0= 0.29
;process_and_plot_sfr_history, "vc3vc3e_gf8a", lcolor=100, lthick= 2.0, msg='vc3vc3e_gf8a', x0= 0.25, y0= 0.26

;msg=' '

;fload_newcolortable, 1
;process_and_plot_sfr_history, "zs/z3e", lcolor=20, lthick= 2.0, msg='5%', x0= 0.75, y0= 0.82
;process_and_plot_sfr_history, "es/e3e", lcolor=70, lthick= 2.0, msg='20%', x0= 0.75, y0= 0.78
;process_and_plot_sfr_history, "ds/vc3vc3e", lcolor=120, lthick= 2.0, msg='40%', x0= 0.75, y0= 0.74
;process_and_plot_sfr_history, "ds/d3e7", lcolor=120, lthick= 2.0, msg='40%', x0= 0.75, y0= 0.70
;process_and_plot_sfr_history, "bs/b3e", lcolor=170, lthick= 2.0, msg='80%', x0= 0.75, y0= 0.66



; vc3 - halo concentration
; -------------------------
;fruns=["vc3vc3e_c5","vc3vc3e_c7","vc3vc3e", "vc3vc3e_c11", "vc3vc3e_c13"]
;msg=' '




; ------------------



;fruns= ["ds/d0e", "ds/d0h"]
;fruns= ["ds/d1e", "ds/d1h"]
;fruns= ["ds/d2e", "ds/d2h"]


; ------------------



; larger vcs
; -----------
;fruns= ["vc4", "vc4a"]
;msg=' '
;fruns= ["vc5", "vc5a"]
;msg=' '
;fruns= ["vc6", "vc6a"]
;msg=' '


;fruns= ["vc4avc4a", "vc5avc5a","vc6avc6a"]
;	tmreg=[1.1,1.1,1.1]
;	msg=' '


;xyouts, 0.20, 0.36, 'Orbit', /normal, charthick=3, size=1.33, color=0
;process_and_plot_sfr_historprocess_and_plot_sfr_history, "bs/b2e", lcolor=200, lthick= 2.0, msg='e', x0= 0.25, y0= 0.32
;process_and_plot_sfr_history, "bs/b2h", lcolor=150, lthick= 2.0, msg='h', x0= 0.25, y0= 0.28
;process_and_plot_sfr_history, "bs/b2f", lcolor=100, lthick=2.0, msg='f', x0= 0.25, y0= 0.24
;process_and_plot_sfr_history, "bs/b2k", lcolor=50, lthick=2.0, msg='k', x0= 0.25, y0= 0.20

; reposition on pot min - does it matter?
;fruns= ["vc4vc4a","vc4vc4b"]
;msg=' '	



;fruns= ["vc4vc4e","vc4vc4a","vc4avc4a","vc4vc4_no"]
;msg= ' '
;fruns= ["vc5vc5e","vc5vc5a","vc5avc5a","vc5vc5_no"]
;msg= ' '
;fruns= ["vc6vc6e","vc6vc6a"]
;msg= ' '



; a3 - z comparison
; -----------------
;xyouts, 0.70, 0.92, 'Log(M!Dtot!N/M!D!9n!6!N)', /normal, charthick=3, size=1.33, color=0
;process_and_plot_sfr_history, "As/A3", lcolor=150, lthick= 2.0, msg='z= 0.0', x0= 0.80, y0= 0.63
;process_and_plot_sfr_history, "z3/a3", lcolor=50, lthick=2.0, msg='z= 3.0', x0= 0.80, y0= 0.67
;process_and_plot_sfr_history, "z3/a4.5", lcolor=50, lthick=2.0, msg='z= 3.0', x0= 0.80, y0= 0.67
;process_and_plot_sfr_history, "bs/b2e", lcolor=140, ctab= 1, lthick=2.0, msg='11.7', x0= 0.80, y0= 0.71



; z3 galaxies (brant v. mine)
; ------------------------------
;xyouts, 0.70, 0.92, 'Log(M!Dtot!N/M!D!9n!6!N)', /normal, charthick=3, size=1.33, color=0
;process_and_plot_sfr_history, "brantsruns/full_model/z3/v3g2q1z3_v3g2q1z3_p", lcolor=0, lthick= 2.0, msg='v3g2q1z3', x0= 0.60, y0= 0.85
;process_and_plot_sfr_history, "z3/b3h", lcolor=150, lthick= 2.0, msg='z3/b3h', x0= 0.60, y0= 0.80
;process_and_plot_sfr_history, "z3/b3h2", lcolor=100, lthick= 2.0, msg='z3/b3h2', x0= 0.60, y0= 0.75
;process_and_plot_sfr_history, "z3/b3h3", lcolor=200, lthick= 2.0, msg='z3/b3h3', x0= 0.60, y0= 0.70
;process_and_plot_sfr_history, "z3/b3e", lcolor=50, lthick=2.0, msg='z3/b3e', x0= 0.60, y0= 0.87

;process_and_plot_sfr_history, "z3/b3h2", lcolor=100, lthick= 2.0, msg='z3/b3h2 - raw', x0= 0.60, y0= 0.85, /raw
;process_and_plot_sfr_history, "z3/b3h2", lcolor=150, lthick= 2.0, msg='z3/b3h2', x0= 0.60, y0= 0.75

; 80% gas
; set xmax= 202 above
;process_and_plot_sfr_history, "brantsruns/full_model/z3/v3g2q1z3_v3g2q1z3_p", lcolor=0, lthick= 2.0, msg='v3g2q1z3', x0= 0.60, y0= 0.85
;process_and_plot_sfr_history, "z3/b3h3", lcolor=200, lthick= 2.0, msg='z3/b3h3', x0= 0.60, y0= 0.80
;process_and_plot_sfr_history, "z3/b3e2", lcolor=150, lthick= 2.0, msg='z3/b3e2', x0= 0.60, y0= 0.75

; set xmax= 600 above
;process_and_plot_sfr_history, "brantsruns/full_model/z3/v4g2q1z3_v4g2q1z3_p", lcolor=0, lthick= 2.0, msg='v4g2q1z3', x0= 0.60, y0= 0.85
;process_and_plot_sfr_history, "z3/b4h", lcolor=200, lthick= 2.0, msg='z3/b4h', x0= 0.60, y0= 0.80
;process_and_plot_sfr_history, "z3/b4e", lcolor=150, lthick= 2.0, msg='z3/b4e', x0= 0.60, y0= 0.75

; set xmax= 2000 above
;process_and_plot_sfr_history, "brantsruns/full_model/z3/v5g2q1z3_v5g2q1z3_p", lcolor=0, lthick= 2.0, msg='v5g2q1z3', x0= 0.60, y0= 0.85
;process_and_plot_sfr_history, "z3/b5h", lcolor=200, lthick= 2.0, msg='z3/b5h', x0= 0.60, y0= 0.80
;process_and_plot_sfr_history, "z3/b5e", lcolor=150, lthick= 2.0, msg='z3/b5e', x0= 0.60, y0= 0.75

;process_and_plot_sfr_history, "brantsruns/full_model/z3/v6g2q1z3_v6g2q1z3_p", lcolor=0, lthick= 2.0, msg='v6g2q1z3', x0= 0.60, y0= 0.85
;process_and_plot_sfr_history, "z3/b6h", lcolor=200, lthick= 2.0, msg='z3/b6h', x0= 0.60, y0= 0.80
;process_and_plot_sfr_history, "z3/b6e", lcolor=150, lthick= 2.0, msg='z3/b6e', x0= 0.60, y0= 0.85
;process_and_plot_sfr_history, "z3/b6f", lcolor=100, lthick= 2.0, msg='z3/b6f', x0= 0.60, y0= 0.80
;process_and_plot_sfr_history, "z3/b6h", lcolor=50, lthick= 2.0, msg='z3/b6h', x0= 0.60, y0= 0.75
;process_and_plot_sfr_history, "z3/b6k", lcolor=200, lthick= 2.0, msg='z3/b6k', x0= 0.60, y0= 0.70

;process_and_plot_sfr_history, "z3/b6k", lcolor=200, lthick= 2.0, msg='z3/b6k', x0= 0.30, y0= 0.45
;process_and_plot_sfr_history, "z3/b6k1", lcolor=175, lthick= 2.0, msg='z3/b6k1', x0= 0.30, y0= 0.42
;process_and_plot_sfr_history, "z3/b6k2", lcolor=150, lthick= 2.0, msg='z3/b6k2', x0= 0.30, y0= 0.39
;process_and_plot_sfr_history, "z3/b6k3", lcolor=125, lthick= 2.0, msg='z3/b6k3', x0= 0.30, y0= 0.36
;process_and_plot_sfr_history, "z3/b6k4", lcolor=100, lthick= 2.0, msg='z3/b6k4', x0= 0.30, y0= 0.33
;process_and_plot_sfr_history, "z3/b6k5", lcolor=75, lthick= 2.0, msg='z3/b6k5', x0= 0.30, y0= 0.30
;process_and_plot_sfr_history, "z3/b6k6", lcolor=50, lthick= 2.0, msg='z3/b6k6', x0= 0.30, y0= 0.27
;process_and_plot_sfr_history, "z3/b6k7", lcolor=25, lthick= 2.0, msg='z3/b6k7', x0= 0.30, y0= 0.24

;
; new z=3 SMG models
;
;process_and_plot_sfr_history, "z3/b6e", lcolor=200, lthick= 2.0, msg='z3b6e', x0= 0.30, y0= 0.45
;process_and_plot_sfr_history, "z3/b6v1e", lcolor=150, lthick= 2.0, msg='z3b6v1e', x0= 0.30, y0= 0.40
;process_and_plot_sfr_history, "z3/b6v2e", lcolor=100, lthick= 2.0, msg='z3b6v2e', x0= 0.30, y0= 0.35
;process_and_plot_sfr_history, "z3/b6v3e_a", lcolor=50, lthick= 2.0, msg='z3b6v3e', x0= 0.30, y0= 0.30

;process_and_plot_sfr_history, "z3/b5e", lcolor=0, lthick= 5.0, msg='z3/b5e', x0= 0.70, y0= 0.91
;process_and_plot_sfr_history, "z3/b5e_2", lcolor=200, lthick= 2.0, msg='z3/b5e_2', x0= 0.70, y0= 0.88
;process_and_plot_sfr_history, "z3/b5e_3", lcolor=170, lthick= 2.0, msg='z3/b5e_3', x0= 0.70, y0= 0.83
;process_and_plot_sfr_history, "z3/b5e_4", lcolor=140, lthick= 2.0, msg='z3/b5e_4', x0= 0.70, y0= 0.80
;process_and_plot_sfr_history, "z3/b5e_5", lcolor=110, lthick= 2.0, msg='z3/b5e_5', x0= 0.70, y0= 0.77
;process_and_plot_sfr_history, "z3/b5e_6", lcolor=80, lthick= 2.0, msg='z3/b5e_6', x0= 0.70, y0= 0.74

;process_and_plot_sfr_history, "z3/b5e_7", lcolor=50, lthick= 2.0, msg='z3/b5e_7', x0= 0.70, y0= 0.83
;process_and_plot_sfr_history, "z3/b5e_7", lcolor=50, lthick= 2.0, msg='z3/b5e_7', x0= 0.70, y0= 0.74
;process_and_plot_sfr_history, "z3/b5e_8", lcolor=150, lthick= 2.0, msg='z3/b5e_8', x0= 0.70, y0= 0.71

;process_and_plot_sfr_history, "z3/b5e_9", lcolor=50, lthick= 2.0, msg='z3/b5e_9', x0= 0.70, y0= 0.74
;process_and_plot_sfr_history, "z3/b5e_10", lcolor=150, lthick= 2.0, msg='z3/b5e_10', x0= 0.70, y0= 0.71

;process_and_plot_sfr_history, "z3/b5e_11", lcolor=50, lthick= 2.0, msg='z3/b5e_11', x0= 0.70, y0= 0.74
;process_and_plot_sfr_history, "z3/b5e_12", lcolor=150, lthick= 2.0, msg='z3/b5e_12', x0= 0.70, y0= 0.71

;process_and_plot_sfr_history, "z3/b5e_13", lcolor=100, lthick= 2.0, msg='z3/b5e_13', x0= 0.70, y0= 0.80
;process_and_plot_sfr_history, "z3/b5e_13", lcolor=50, lthick= 2.0, msg='z3/b5e_13', x0= 0.70, y0= 0.74
;process_and_plot_sfr_history, "z3/b5e_14", lcolor=150, lthick= 2.0, msg='z3/b5e_14', x0= 0.70, y0= 0.71

;process_and_plot_sfr_history, "z3/b5e_15", lcolor=130, lthick= 2.0, msg='z3/b5e_15', x0= 0.70, y0= 0.77
;process_and_plot_sfr_history, "z3/b5e_16", lcolor=160, lthick= 2.0, msg='z3/b5e_16', x0= 0.70, y0= 0.74
;process_and_plot_sfr_history, "z3/b5e_17", lcolor=240, lthick= 2.0, msg='z3/b5e_17', x0= 0.70, y0= 0.71
;--------------

;process_and_plot_sfr_history, "z3/b4e", lcolor=0, lthick= 5.0, msg='z3/b4e', x0= 0.70, y0= 0.93
;process_and_plot_sfr_history, "z3/b4v1e", lcolor=220, lthick= 2.0, msg='z3/b4v1e', x0= 0.70, y0= 0.88
;process_and_plot_sfr_history, "z3/b4v2e", lcolor=200, lthick= 2.0, msg='z3/b4v2e', x0= 0.70, y0= 0.85
;process_and_plot_sfr_history, "z3/b4v3e", lcolor=180, lthick= 2.0, msg='z3/b4v3e', x0= 0.70, y0= 0.82
;process_and_plot_sfr_history, "z3/b4v4e", lcolor=160, lthick= 2.0, msg='z3/b4v4e', x0= 0.70, y0= 0.79
;process_and_plot_sfr_history, "z3/b4v5e", lcolor=140, lthick= 2.0, msg='z3/b4v5e', x0= 0.70, y0= 0.76
;process_and_plot_sfr_history, "z3/b4v6e", lcolor=120, lthick= 2.0, msg='z3/b4v6e', x0= 0.70, y0= 0.73
;process_and_plot_sfr_history, "z3/b4v7e", lcolor=100, lthick= 2.0, msg='z3/b4v7e', x0= 0.70, y0= 0.70


;--------------
;process_and_plot_sfr_history, "z3/b5e", lcolor=0, lthick= 5.0, msg='z3/b5e', x0= 0.70, y0= 0.91, gasmass=gasmass
;process_and_plot_sfr_history, "z3/b5e_2", lcolor=200, lthick= 2.0, msg='z3/b5e_2', x0= 0.70, y0= 0.88, gasmass=gasmass
;process_and_plot_sfr_history, "z3/b5e_7", lcolor=50, lthick= 2.0, msg='z3/b5e_7', x0= 0.70, y0= 0.83, gasmass=gasmass
;process_and_plot_sfr_history, "z3/b5e_13", lcolor=100, lthick= 2.0, msg='z3/b5e_13', x0= 0.70, y0= 0.80, gasmass=gasmass
;process_and_plot_sfr_history, "z3/b5e_15", lcolor=130, lthick= 2.0, msg='z3/b5e_15', x0= 0.70, y0= 0.77, gasmass=gasmass
;process_and_plot_sfr_history, "z3/b5e_16", lcolor=160, lthick= 2.0, msg='z3/b5e_16', x0= 0.70, y0= 0.74, gasmass=gasmass
;process_and_plot_sfr_history, "z3/b5e_17", lcolor=240, lthick= 2.0, msg='z3/b5e_17', x0= 0.70, y0= 0.71, gasmass=gasmass
;--------------

;process_and_plot_sfr_history, "z3/b5v1e", lcolor=150, lthick= 2.0, msg='z3/b5v1e', x0= 0.70, y0= 0.86
;process_and_plot_sfr_history, "z3/b5v2e", lcolor=250, lthick= 2.0, msg='z3/b5v2e', x0= 0.70, y0= 0.83
;process_and_plot_sfr_history, "z3/b5v3e", lcolor=220, lthick= 2.0, msg='z3/b5v3e', x0= 0.70, y0= 0.81
;process_and_plot_sfr_history, "z3/b5v4e", lcolor=190, lthick= 2.0, msg='z3/b5v4e', x0= 0.70, y0= 0.78
;process_and_plot_sfr_history, "z3/b5v4ea", lcolor=170, lthick= 2.0, msg='z3/b5v4ea', x0= 0.70, y0= 0.75
;process_and_plot_sfr_history, "z3/b5v4ea", lcolor=120, lthick= 2.0, msg='z3/b5v4ea', x0= 0.70, y0= 0.75
;process_and_plot_sfr_history, "z3/b5v9e", lcolor=150, lthick= 2.0, msg='z3/b5v9e', x0= 0.70, y0= 0.72
;process_and_plot_sfr_history, "z3/b5v10e", lcolor=130, lthick= 2.0, msg='z3/b5v10e', x0= 0.70, y0= 0.69
;process_and_plot_sfr_history, "z3/b5v11e", lcolor=110, lthick= 2.0, msg='z3/b5v11e', x0= 0.70, y0= 0.66
;process_and_plot_sfr_history, "z3/b5v11ea", lcolor=90, lthick= 2.0, msg='z3/b5v11ea', x0= 0.70, y0= 0.63
;process_and_plot_sfr_history, "z3/b5v11ea", lcolor=50, lthick= 2.0, msg='z3/b5v11ea', x0= 0.70, y0= 0.63
;process_and_plot_sfr_history, "z3/b5v12e", lcolor=70, lthick= 2.0, msg='z3/b5v12e', x0= 0.70, y0= 0.60
;process_and_plot_sfr_history, "z3/b5v13e", lcolor=50, lthick= 2.0, msg='z3/b5v13e', x0= 0.70, y0= 0.57
;process_and_plot_sfr_history, "z3/b5v13e", lcolor=50, lthick= 2.0, msg='z3/b5v13e', x0= 0.70, y0= 0.80
;process_and_plot_sfr_history, "z3/b5v13e_7", lcolor=80, lthick= 2.0, msg='z3/b5v13e_7', x0= 0.70, y0= 0.77
;process_and_plot_sfr_history, "z3/b5v14e", lcolor=30, lthick= 2.0, msg='z3/b5v14e', x0= 0.70, y0= 0.54
;process_and_plot_sfr_history, "z3/b5v14e", lcolor=120, lthick= 2.0, msg='z3/b5v14e', x0= 0.70, y0= 0.70
;process_and_plot_sfr_history, "z3/b5v14e_7", lcolor=150, lthick= 2.0, msg='z3/b5v14e_7', x0= 0.70, y0= 0.67

;process_and_plot_sfr_history, "brantsruns/full_model/z3/v6g2q1z3_v6g2q1z3_p", lcolor=50, lthick= 2.0, msg='v6g2q1z3', x0= 0.60, y0= 0.85, lstyle= 1
;process_and_plot_sfr_history, "brantsruns/full_model/z3/v6g2q0z3_v6g2q0z3_p", lcolor=50, lthick= 2.0, msg='v6g2q0z3', x0= 0.60, y0= 0.81
;process_and_plot_sfr_history, "brantsruns/full_model/z3/v6g1q1z3_v6g1q1z3_p", lcolor=150, lthick= 2.0, msg='v6g1q1z3', x0= 0.60, y0= 0.77, lstyle= 1
;process_and_plot_sfr_history, "brantsruns/full_model/z3/v6g1q0z3_v6g1q0z3_p", lcolor=150, lthick= 2.0, msg='v6g1q0z3', x0= 0.60, y0= 0.73

;process_and_plot_sfr_history, "brantsruns/full_model/z3/v5g2q1z3_v5g2q1z3_p", lcolor=50, lthick= 2.0, msg='v5g2q1z3', x0= 0.60, y0= 0.85, lstyle= 1
;process_one_sfr, "brantsruns/full_model/z3/v5g2q0z3_v5g2q0z3_p", lcolor=50, lthick= 2.0, msg='v5g2q0z3', x0= 0.60, y0= 0.81
;process_one_sfr, "brantsruns/full_model/z3/v5g1q1z3_v5g1q1z3_p", lcolor=150, lthick= 2.0, msg='v5g1q1z3', x0= 0.60, y0= 0.77, lstyle= 1
;process_one_sfr, "brantsruns/full_model/z3/v5g1q0z3_v5g1q0z3_p", lcolor=150, lthick= 2.0, msg='v5g1q0z3', x0= 0.60, y0= 0.73

;process_one_sfr, "brantsruns/full_model/z3/v4g2q1z3_v4g2q1z3_p", lcolor=50, lthick= 2.0, msg='v4g2q1z3', x0= 0.60, y0= 0.85, lstyle= 1
;process_one_sfr, "brantsruns/full_model/z3/v4g2q0z3_v4g2q0z3_p", lcolor=50, lthick= 2.0, msg='v4g2q0z3', x0= 0.60, y0= 0.81
;process_one_sfr, "brantsruns/full_model/z3/v4g1q1z3_v4g1q1z3_p", lcolor=150, lthick= 2.0, msg='v4g1q1z3', x0= 0.60, y0= 0.77, lstyle= 1
;process_one_sfr, "brantsruns/full_model/z3/v4g1q0z3_v4g1q0z3_p", lcolor=150, lthick= 2.0, msg='v4g1q0z3', x0= 0.60, y0= 0.73

;process_one_sfr, "z3/iso_b5", lcolor=50, lthick= 2.0, msg='z3/iso_b5', x0= 0.70, y0= 0.89 ; , multiby=2.0
;;process_one_sfr, "z3/iso_b5", lcolor=60, lthick= 2.0, msg='z3/iso_b5 (x2)', x0= 0.70, y0= 0.86, multiby=2.0
;process_one_sfr, "z3/iso_b5v1", lcolor=80, lthick= 2.0, msg='z3/iso_b5v1', x0= 0.70, y0= 0.83 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v2", lcolor=100, lthick= 2.0, msg='z3/iso_b5v2', x0= 0.70, y0= 0.80 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v3", lcolor=120, lthick= 2.0, msg='z3/iso_b5v3', x0= 0.70, y0= 0.77 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v4", lcolor=140, lthick= 2.0, msg='z3/iso_b5v4', x0= 0.70, y0= 0.74 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v5", lcolor=160, lthick= 2.0, msg='z3/iso_b5v5', x0= 0.70, y0= 0.71 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v6", lcolor=180, lthick= 2.0, msg='z3/iso_b5v6', x0= 0.70, y0= 0.68 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v7", lcolor=200, lthick= 2.0, msg='z3/iso_b5v7', x0= 0.70, y0= 0.65 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v8", lcolor=220, lthick= 2.0, msg='z3/iso_b5v8', x0= 0.70, y0= 0.62 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v9", lcolor=240, lthick= 2.0, msg='z3/iso_b5v9', x0= 0.70, y0= 0.59 ; , multiby=2.0


;process_one_sfr, "z3/iso_b5", lcolor=50, lthick= 2.0, msg='z3/iso_b5', x0= 0.70, y0= 0.89 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v1", lcolor=80, lthick= 2.0, msg='z3/iso_b5v1', x0= 0.70, y0= 0.86 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v4", lcolor=100, lthick= 2.0, msg='z3/iso_b5v4', x0= 0.70, y0= 0.83 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v7", lcolor=120, lthick= 2.0, msg='z3/iso_b5v7', x0= 0.70, y0= 0.80 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v9", lcolor=140, lthick= 2.0, msg='z3/iso_b5v9', x0= 0.70, y0= 0.77 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v10", lcolor=160, lthick= 2.0, msg='z3/iso_b5v10', x0= 0.70, y0= 0.74 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v11", lcolor=180, lthick= 2.0, msg='z3/iso_b5v11', x0= 0.70, y0= 0.71 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v12", lcolor=200, lthick= 2.0, msg='z3/iso_b5v12', x0= 0.70, y0= 0.68 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v13", lcolor=220, lthick= 2.0, msg='z3/iso_b5v13', x0= 0.70, y0= 0.65 ; , multiby=2.0
;process_one_sfr, "z3/iso_b5v14", lcolor=240, lthick= 2.0, msg='z3/iso_b5v14', x0= 0.70, y0= 0.62 ; , multiby=2.0



;process_one_sfr, "z3/b6b5e", lcolor=240, lthick= 2.0, msg='b6b5e', x0= 0.70, y0= 0.88 ; , multiby=2.0
;process_one_sfr, "z3/b6b5k", lcolor=200, lthick= 2.0, msg='b6b5k', x0= 0.70, y0= 0.84 ; , multiby=2.0
;process_one_sfr, "z3/b6b5f", lcolor=160, lthick= 2.0, msg='b6b5f', x0= 0.70, y0= 0.80 ; , multiby=2.0

;process_one_sfr, "z3/b6b4e", lcolor=120, lthick= 2.0, msg='b6b4e', x0= 0.70, y0= 0.76 ; , multiby=2.0
;process_one_sfr, "z3/b6b4k", lcolor=80, lthick= 2.0, msg='b6b4k', x0= 0.70, y0= 0.72 ; , multiby=2.0
;process_one_sfr, "z3/b6b4f", lcolor=40, lthick= 2.0, msg='b6b4f', x0= 0.70, y0= 0.68 ; , multiby=2.0

;process_one_sfr, "z3/b6e", lcolor=00, lthick= 2.0, msg='b6e', x0= 0.70, y0= 0.57 ; , multiby=2.0
;process_one_sfr, "z3/b6f", lcolor=50, lthick= 2.0, msg='b6f', x0= 0.70, y0= 0.53 ; , multiby=2.0
;process_one_sfr, "z3/b6h", lcolor=80, lthick= 2.0, msg='b6h', x0= 0.70, y0= 0.49 ; , multiby=2.0
;process_one_sfr, "z3/b6k", lcolor=120, lthick= 2.0, msg='b6k', x0= 0.70, y0= 0.45 ; , multiby=2.0


; 40% gas
;process_one_sfr, "brantsruns/full_model/z3/v3g1q1z3_v3g1q1z3_p", lcolor=50, lthick= 2.0, msg='v3g1q1z3', x0= 0.60, y0= 0.85
;process_one_sfr, "z3/d3h", lcolor=200, lthick= 2.0, msg='z3/d3h', x0= 0.60, y0= 0.80
;process_one_sfr, "z3/d3e", lcolor=150, lthick= 2.0, msg='z3/d3e', x0= 0.60, y0= 0.75

;process_one_sfr, "brantsruns/full_model/z3/v4g1q1z3_v4g1q1z3_p", lcolor=50, lthick= 2.0, msg='v4g1q1z3', x0= 0.60, y0= 0.85
;process_one_sfr, "z3/d4h", lcolor=200, lthick= 2.0, msg='z3/d4h', x0= 0.60, y0= 0.80
;process_one_sfr, "z3/d4e", lcolor=150, lthick= 2.0, msg='z3/d4e', x0= 0.60, y0= 0.75

;process_one_sfr, "brantsruns/full_model/z3/v5g1q1z3_v5g1q1z3_p", lcolor=50, lthick= 2.0, msg='v5g1q1z3', x0= 0.60, y0= 0.85
;process_one_sfr, "z3/d5h", lcolor=200, lthick= 2.0, msg='z3/d5h', x0= 0.60, y0= 0.80
;process_one_sfr, "z3/d5e", lcolor=150, lthick= 2.0, msg='z3/d5e', x0= 0.60, y0= 0.75

;process_one_sfr, "brantsruns/full_model/z3/v6g1q1z3_v6g1q1z3_p", lcolor=50, lthick= 2.0, msg='v6g1q1z3', x0= 0.60, y0= 0.85
;process_one_sfr, "z3/d6h", lcolor=200, lthick= 2.0, msg='z3/d6h', x0= 0.60, y0= 0.80
;process_one_sfr, "z3/d6e", lcolor=150, lthick= 2.0, msg='z3/d6e', x0= 0.60, y0= 0.75




;Large A's
;------------
;process_one_sfr, "As/A7", lcolor=150, lthick= 2.0, msg='A7', x0= 0.80, y0= 0.83
;process_one_sfr, "As/A7_q", lcolor=50, lthick= 2.0, msg='A7_q', x0= 0.80, y0= 0.80

;process_one_sfr, "As/A4", lcolor=150, lthick= 2.0, msg='A4', x0= 0.26, y0= 0.35
;process_one_sfr, "As/A5", lcolor=50, lthick= 2.0, msg='A5', x0= 0.26, y0= 0.30
;process_one_sfr, "As/A6", lcolor=100, lthick= 2.0, msg='A6', x0= 0.26, y0= 0.25


;yz= 0.5
;process_one_sfr, "As/A3", lcolor=0, lthick= 2.0, msg='A3', x0= 0.22, y0= yz+0.35, cumulative=cumulative
;process_one_sfr, "As/A3sb8", lcolor=50, lthick= 2.0, msg='A3sb8', x0= 0.22, y0= yz+0.31, cumulative=cumulative
;process_one_sfr, "As/A3sb9", lcolor=100, lthick= 2.0, msg='A3sb9', x0= 0.22, y0= yz+0.27, cumulative=cumulative
;process_one_sfr, "As/A3sb10", lcolor=150, lthick= 2.0, msg='A3sb10', x0= 0.22, y0= yz+0.23, cumulative=cumulative
;process_one_sfr, "As/A3sb13", lcolor=200, lthick= 2.0, msg='A3sb13', x0= 0.22, y0= yz+0.19, cumulative=cumulative


; size comparison
; -----------------
;xyouts, 0.70, 0.92, 'Log(M!Dtot!N/M!D!9n!6!N)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "bs/b0e", lcolor=200, ctab= 1, lthick= 2.0, msg='10.8', x0= 0.80, y0= 0.63
;process_one_sfr, "bs/b1e", lcolor=170, ctab= 1, lthick=2.0, msg='11.2', x0= 0.80, y0= 0.67
;process_one_sfr, "bs/b2e", lcolor=140, ctab= 1, lthick=2.0, msg='11.7', x0= 0.80, y0= 0.71
;process_one_sfr, "bs/b3e", lcolor=110, ctab= 1, lthick=2.0, msg='12.1', x0= 0.80, y0= 0.75
;process_one_sfr, "bs/b4e", lcolor=80, ctab= 1, lthick=2.0, msg='12.5', x0= 0.80, y0= 0.79
;process_one_sfr, "bs/b5e", lcolor=50, ctab= 1, lthick=2.0, msg='13.0', x0= 0.80, y0= 0.83
;process_one_sfr, "bs/b6e", lcolor=20, ctab= 1, lthick=2.0, msg='13.6', x0= 0.80, y0= 0.87
;xyouts, 0.20, 0.52, 'Log(M!Dtot!N/M!D!9n!6!N)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "bs/b0e", lcolor=200, ctab= 1, lthick= 2.0, msg='10.8', x0= 0.30, y0= 0.23, /normalize
;process_one_sfr, "bs/b1e", lcolor=170, ctab= 1, lthick=2.0, msg='11.2', x0= 0.30, y0= 0.27, /normalize
;process_one_sfr, "bs/b2e", lcolor=140, ctab= 1, lthick=2.0, msg='11.7', x0= 0.30, y0= 0.31, /normalize
;process_one_sfr, "bs/b3e", lcolor=110, ctab= 1, lthick=2.0, msg='12.1', x0= 0.30, y0= 0.35, /normalize
;process_one_sfr, "bs/b4e", lcolor=80, ctab= 1, lthick=2.0, msg='12.5', x0= 0.30, y0= 0.39, /normalize
;process_one_sfr, "bs/b5e", lcolor=50, ctab= 1, lthick=2.0, msg='13.0', x0= 0.30, y0= 0.43, /normalize
;process_one_sfr, "bs/b6e", lcolor=20, ctab= 1, lthick=2.0, msg='13.6', x0= 0.30, y0= 0.47, /normalize





;xyouts, 0.20, 0.25, 'Log(M!Dtot!N/M!D!9n!6!N)= 10.8', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.20, 0.85, 'Log(M!Dtot!N/M!D!9n!6!N)= 10.8', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "bs/b0e", lcolor=200, lthick= 2.0, msg=' ', x0= 0.30, y0= 0.23
;process_one_sfr, "bs/b0h", lcolor=150, lthick=2.0, msg=' ', x0= 0.30, y0= 0.27
;process_one_sfr, "bs/b0f", lcolor=100, lthick=2.0, msg=' ', x0= 0.30, y0= 0.31
;process_one_sfr, "bs/b0k", lcolor=50, lthick=2.0, msg=' ', x0= 0.30, y0= 0.35
;process_one_sfr, "ds/d0e2_q", lcolor=200, lthick= 2.0, msg=' ', x0= 0.30, y0= 0.23
;process_one_sfr, "ds/d0h2_q", lcolor=150, lthick=2.0, msg=' ', x0= 0.30, y0= 0.27
;process_one_sfr, "ds/d0f2_q", lcolor=100, lthick=2.0, msg=' ', x0= 0.30, y0= 0.31
;process_one_sfr, "ds/d0k2_q", lcolor=50, lthick=2.0, msg=' ', x0= 0.30, y0= 0.35

;xyouts, 0.50, 0.85, 'Log(M!Dtot!N/M!D!9n!6!N)= 12.1', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "bs/b3h", lcolor=150, lthick=2.0, msg=' ', x0= 0.30, y0= 0.27

;xyouts, 0.50, 0.85, 'Log(M!Dtot!N/M!D!9n!6!N)= 13.0', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "bs/b5h", lcolor=150, lthick=2.0, msg=' ', x0= 0.30, y0= 0.27

;process_one_sfr, "vcs/vc1vc1", lcolor=200, ctab= 1, lthick=2.0, msg='with BH', x0= 0.30, y0= 0.43
;process_one_sfr, "vcs/vc1vc1_no", lcolor=20, ctab= 1, lthick=2.0, msg='without BH', x0= 0.30, y0= 0.47



; with halo gas
; ---------------
;fruns=["v3c","v3h","v3hh","v3ha","v3hha"]
;fruns=["vc3HGe","vc3HGh"]
;fruns=["vc3HGe","vc3vc3e"] & lbls= ["with halo gas", "std"]
;fruns=["vc3HGh","vc3vc3h"] & lbls= ["with halo gas", "std"]


;process_one_sfr, "bs/b5e", lcolor=200, ctab= 1, lthick=2.0, msg='std', x0= 0.30, y0= 0.43
;process_one_sfr, "bs/b5e_igm1", lcolor=20, ctab= 1, lthick=2.0, msg='with IGM', x0= 0.30, y0= 0.47


;process_one_sfr, "bs/b2e", lcolor=0, lthick=2.0, msg='b2e', x0= 0.20, y0= 0.43
;process_one_sfr, "bs/b2e_igm1", lcolor=50, lthick=2.0, msg='b2e_igm1', x0= 0.20, y0= 0.39
;process_one_sfr, "bs/b2e_igm2", lcolor=100, lthick=2.0, msg='b2e_igm2', x0= 0.20, y0= 0.35
;process_one_sfr, "bs/b2e_igm3", lcolor=150, lthick=2.0, msg='b2e_igm3', x0= 0.20, y0= 0.31



;-------------------
; nosf's (for phil)

;process_one_sfr, "nosf/nod3e", lcolor=200, lthick=2.0, msg='nod3e (!8t!6!D*!N/10!E4!N, q=0.25)', x0= 0.25, y0= 0.92
;process_one_sfr, "nosf/nod3e2", lcolor=150, lthick=2.0, msg='nod3e2 (!8t!6!D*!N/10!E4!N, q=1.0)', x0= 0.25, y0= 0.88
;process_one_sfr, "nosf/nod3e3", lcolor=100, lthick=2.0, msg='nod3e3 (nosf, q=0.25)', x0= 0.25, y0= 0.84
;process_one_sfr, "nosf/nod3e4", lcolor=50, lthick=2.0, msg='nod3e4 (nosf, q=1.0)', x0= 0.25, y0= 0.80
;xyouts, 0.25, 0.76, '(note: no actual sf in nod3e3,4)', /normal, charthick=3, size=0.75, color=0

; essentially the same as nod3e
;process_one_sfr, "nosf/nod3h", lcolor=0, lthick=2.0, msg='nod3h (tau_sf)', x0= 0.25, y0= 0.70



;-------------------
; q_eos comparisons

;process_one_sfr, "bs/b5e", lcolor=200, lthick=2.0, msg='b5e (q=1.0)', x0= 0.65, y0= 0.92
;process_one_sfr, "ds/d5e", lcolor=200, lthick=2.0, msg='d5e (q=1.0)', x0= 0.65, y0= 0.92
;process_one_sfr, "ds/d5e2_q", lcolor=150, lthick=2.0, msg='d5e2_q (q=0.25)', x0= 0.65, y0= 0.88
;process_one_sfr, "ds/d5e2", lcolor=100, lthick=2.0, msg='d5e2 (q=1.0)', x0= 0.65, y0= 0.84
;process_one_sfr, "ds/d5e2_2", lcolor=50, lthick=2.0, msg='d5e2_2 (q=1.0)', x0= 0.65, y0= 0.80

;process_one_sfr, "ds/d4e", lcolor=200, lthick=2.0, msg='d4e (q=1.0)', x0= 0.65, y0= 0.92
;process_one_sfr, "ds/d4e2_q", lcolor=150, lthick=2.0, msg='d4e2_q (q=0.25)', x0= 0.65, y0= 0.88
;process_one_sfr, "ds/d4e2", lcolor=100, lthick=2.0, msg='d4e2 (q=1.0)', x0= 0.65, y0= 0.84

;process_one_sfr, "ds/d6e", lcolor=200, lthick=2.0, msg='d6e (q=1.0)', x0= 0.65, y0= 0.92
;process_one_sfr, "ds/d6e2_q", lcolor=150, lthick=2.0, msg='d6e2_q (q=0.25)', x0= 0.65, y0= 0.88
;process_one_sfr, "ds/d6e2", lcolor=100, lthick=2.0, msg='d6e2 (q=1.0)', x0= 0.65, y0= 0.84

;process_one_sfr, "ds/d4e2_q", lcolor=150, lthick=2.0, msg='low feedback', x0= 0.65, y0= 0.88
;process_one_sfr, "ds/d4e2", lcolor=50, lthick=2.0, msg='high feedback', x0= 0.65, y0= 0.84
;process_one_sfr, "ds/d5e2_q", lcolor=150, lthick=2.0, msg='low feedback', x0= 0.65, y0= 0.84
;process_one_sfr, "ds/d5e2", lcolor=50, lthick=2.0, msg='high feedback', x0= 0.65, y0= 0.88

;process_one_sfr, "bs/b3e", lcolor=200, lthick=2.0, msg='b3e (q=1.0)', x0= 0.62, y0= 0.92
;process_one_sfr, "ds/vc3vc3e_2", lcolor=150, lthick=2.0, msg='vc3vc3e_2 (q=0.25)', x0= 0.26, y0= 0.28
;process_one_sfr, "ds/d3e7", lcolor=100, lthick=2.0, msg='d3e7 (q=0.25)', x0= 0.62, y0= 0.84


;process_one_sfr, "ds/d4e2_q", lcolor=50, lthick=2.0, msg='d4e2_q', x0= 0.65, y0= 0.88
;process_one_sfr, "ds/d5e2_q", lcolor=100, lthick=2.0, msg='d5e2_q', x0= 0.65, y0= 0.84
;process_one_sfr, "ds/d6e2_q", lcolor=150, lthick=2.0, msg='d5e2_q', x0= 0.65, y0= 0.80


; with bulge
; ------------
;fruns=["vc3vc3e","vc3vc3e_bulge","vc3vc3f", "vc3vc3f_bulge","vc3vc3k","vc3vc3k_bulge"]
;fruns=["vc3vc3e","vc3vc3e_bulge"]
;fruns=["vc3vc3h","vc3vc3h_bulge"]




; Minors
; --------
;fruns= ["vc3bvc3b","vc3bvc3","vc3bvc2","vc3bvc1"]
;       tmreg=[1.1,1.1,1.1,1.1]
       ;lbls=['vc3bvc3b','vc3bvc3b_w']
;       msg=''

;fruns= ["vc3vc3","vc3vc2","vc3vc1"]
;       tmreg=[1.1,1.1,1.1]
       ;lbls=['vc3bvc3b','vc3bvc3b_w']
;       msg=''

;fruns= ["vc2vc2","vc2vc1"]
;       tmreg=[1.1,1.1]
       ;lbls=['vc3bvc3b','vc3bvc3b_w']
;       msg=''


;fruns= ["vc2vc2_wBH","vc2vc1_wBH"]
;       tmreg=[1.1,1.1]
       ;lbls=['vc3bvc3b','vc3bvc3b_w']
;       msg=''

;fruns= ["vc3vc3_wBH","vc3vc2_wBH","vc3vc1_wBH"]
;       tmreg=[1.1,1.1,1.1]
       ;lbls=['vc3bvc3b','vc3bvc3b_w']
;       msg=''

;fruns= ["vc3bvc3b_wBH","vc3bvc3_wBH","vc3bvc2_wBH","vc3bvc1_wBH"]
;       tmreg=[1.1,1.1,1.1,1.1]
       ;lbls=['vc3bvc3b','vc3bvc3b_w']
;       msg=''



; vary black hole mass
; ----------------------------
;fruns= ["vc3avc3a_1", "vc3avc3a_2", "vc3avc3a_3"]
;fruns= ["vc3avc3a_3","vc3avc3a_2","vc3avc3a_1","vc3avc3a_4", "vc3avc3a_6"]
;fruns= ["vc3ba_2_1", "vc3ba_2_2", "vc3ba_2_3"]
;fruns= ["vc3ba_2_3","vc3ba_2_2","vc3ba_2_1","vc3ba_2_4","vc3ba_2_5"]
        ;tmerg=[1.1,1.1,1.1,1.1,1.1]
        ;lbls=['1e3','1e4','1e5','1e6','1e7']
	;lbls=['1e3 (89%)','1e4 (77%)','1e5 (65%)','1e6 (59%)','1e7 (33%)']    ; vc3avc3a_* with gas consumption
	;lbls=['1e3 (94%)','1e4 (93%)','1e5 (78%)','1e6 (69%)','1e7 (42%)']    ; vc3ba_2_* with gas consumption
	;msg=''

;process_one_sfr, "vc3avc3a_3", lcolor= 20, ctab= 1, lthick=2.0, msg='1e3 (89%)', x0= 0.20, y0= 0.36
;process_one_sfr, "vc3avc3a_2", lcolor= 65, ctab= 1, lthick=2.0, msg='1e4 (77%)', x0= 0.20, y0= 0.32
;process_one_sfr, "vc3avc3a_1", lcolor= 110, ctab= 1, lthick=2.0, msg='1e5 (65%)', x0= 0.20, y0= 0.28
;process_one_sfr, "vc3avc3a_4", lcolor= 155, ctab= 1, lthick=2.0, msg='1e6 (59%)', x0= 0.20, y0= 0.24
;process_one_sfr, "vc3avc3a_6", lcolor= 200, ctab= 1, lthick=2.0, msg='1e7 (33%)', x0= 0.20, y0= 0.20



; Mihos Hernquist feedback models
; --------------------------------
;fruns= ["vc3vc3h", "vc3vc3h_mh_no"]
;lbls=['1e3','1e4','1e5','1e6','1e7']


;process_one_sfr, "ds/vc3vc3h", lcolor=0, lthick= 2.0, msg='std', x0= 0.20, y0= 0.34
;process_one_sfr, "ds/vc3vc3h_2", lcolor=150, lthick= 2.0, msg='std', x0= 0.20, y0= 0.34
;process_one_sfr, "ds/vc3vc3h_no", lcolor=100, lthick= 2.0, msg='no BH', x0= 0.20, y0= 0.30
;process_one_sfr, "vc3vc3h_mh_no", lcolor=150, lthick= 2.0, msg='no, MH', x0= 0.20, y0= 0.26


;-------------------------------------------------------------------


; Springel & Hernquist Winds model
; -----------------------------------

; std's
;process_one_sfr, "ds/vc3vc3e_2", lcolor= 0, lthick=2.0, msg='std', y0= 0.88
;process_one_sfr, "ds/vc3vc3e_no", lcolor= 100, lthick=2.0, msg='no BH', y0=0.84

; all windspeed= 837.6 km/sec (but different efficiency and energy fraction)
;xyouts, 0.18, 0.42, 'Wind Efficiency/', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.20, 0.38, 'Energy Fraction', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "vc3vc3e_sb1", lcolor=200, ctab= 1, lthick= 2.0, msg='0.005/0.0025', x0= 0.20, y0= 0.34
;process_one_sfr, "vc3vc3e_sb9", lcolor=150, ctab= 1, lthick=2.0, msg='0.05/0.025', x0= 0.20, y0= 0.30
;process_one_sfr, "vc3vc3e_sb10", lcolor=110, ctab= 1, lthick=2.0, msg='0.5/0.25', x0= 0.20, y0= 0.26
;process_one_sfr, "vc3vc3e_sb8", lcolor=60, ctab= 1, lthick=2.0, msg='2.0/1.0', x0= 0.20, y0= 0.22
;process_one_sfr, "vc3vc3e_sb11", lcolor=20, ctab= 1, lthick=2.0, msg='5.0/2.5', x0= 0.20, y0= 0.18



; all windefficiency=2.0 and various energyfractions (i.e., windspeeds)
;xyouts, 0.18, 0.40, 'Wind Speed', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.20, 0.36, '(km/sec)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "vc3vc3e_sb8", lcolor=200, ctab= 1, lthick= 2.0, msg='837', x0= 0.20, y0= 0.32
;process_one_sfr, "vc3vc3e_sb12", lcolor=140, ctab= 1, lthick=2.0, msg='418', x0= 0.20, y0= 0.28
;process_one_sfr, "vc3vc3e_sb13", lcolor=80, ctab= 1, lthick=2.0, msg='209', x0= 0.20, y0= 0.24
;process_one_sfr, "vc3vc3e_sb14", lcolor=20, ctab= 1, lthick=2.0, msg='105', x0= 0.20, y0= 0.20


; all windefficiency=0.005 and various energyfractions (i.e., windspeeds)
;xyouts, 0.18, 0.48, 'Wind Speed', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.20, 0.44, '(km/sec)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "vc3vc3e_sb6", lcolor=220, ctab= 1, lthick= 2.0, msg='5297', x0= 0.20, y0= 0.40
;process_one_sfr, "vc3vc3e_sb2", lcolor=180, ctab= 1, lthick= 2.0, msg='1675', x0= 0.20, y0= 0.36
;process_one_sfr, "vc3vc3e_sb1", lcolor=140, ctab= 1, lthick=2.0, msg='837', x0= 0.20, y0= 0.32
;process_one_sfr, "vc3vc3e_sb3", lcolor=100, ctab= 1, lthick=2.0, msg='418', x0= 0.20, y0= 0.28
;process_one_sfr, "vc3vc3e_sb4", lcolor=60, ctab= 1, lthick=2.0, msg='209', x0= 0.20, y0= 0.24
;process_one_sfr, "vc3vc3e_sb5", lcolor=20, ctab= 1, lthick=2.0, msg='105', x0= 0.20, y0= 0.20


;xyouts, 0.20, 0.37, 'WindSpeed', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.22, 0.33, '(km/s)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "sbw/sb14", lcolor=200, ctab= 1, lthick= 2.0, msg='105', x0= 0.20, y0= 0.28
;process_one_sfr, "sbw/sb13", lcolor=200, ctab= 1, lthick= 2.0, msg='210', x0= 0.20, y0= 0.28
;process_one_sfr, "vc3vc3e_sb8", lcolor=200, ctab= 1, lthick= 2.0, msg='837', x0= 0.20, y0= 0.28
;process_one_sfr, "vc3vc3e_sb8BH", lcolor=20, ctab= 1, lthick= 2.0, msg='837 (BH)', x0= 0.20, y0= 0.24
;xyouts, 0.20, 0.18, 'WindEfficiency=2.0', /normal, charthick=3, size=1.33, color=0


;xyouts, 0.20, 0.37, 'WindSpeed', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.22, 0.33, '(km/s)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "sbw/sb22", lcolor=200, ctab= 1, lthick= 2.0, msg='105', x0= 0.20, y0= 0.28
;process_one_sfr, "sbw/sb21", lcolor=200, ctab= 1, lthick= 2.0, msg='210', x0= 0.20, y0= 0.28
;xyouts, 0.20, 0.18, 'WindEfficiency=5.0', /normal, charthick=3, size=1.33, color=0


;xyouts, 0.20, 0.39, 'WindTravL', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.22, 0.35, '(kpc/h)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "sbw/sb10", lcolor=200, ctab= 1, lthick= 2.0, msg='20', x0= 0.20, y0= 0.30
;process_one_sfr, "sbw/sb10tr1", lcolor=140, ctab= 1, lthick= 2.0, msg=' 2', x0= 0.20, y0= 0.26
;process_one_sfr, "sbw/sb10tr2", lcolor=80, ctab= 1, lthick= 2.0, msg=' 0.2', x0= 0.20, y0= 0.22
;process_one_sfr, "sbw/sb10tr3", lcolor=20, ctab= 1, lthick= 2.0, msg=' 0.02', x0= 0.20, y0= 0.18
;process_one_sfr, "sbw/sb10BH", lcolor=200, ctab= 1, lthick= 2.0, msg='20 (BH)', x0= 0.20, y0= 0.30
;process_one_sfr, "sbw/sb10BHtr1", lcolor=140, ctab= 1, lthick= 2.0, msg=' 2 (BH)', x0= 0.20, y0= 0.26
;process_one_sfr, "sbw/sb10BHtr2", lcolor=80, ctab= 1, lthick= 2.0, msg=' 0.2 (BH)', x0= 0.20, y0= 0.22
;process_one_sfr, "sbw/sb10BHtr3", lcolor=20, ctab= 1, lthick= 2.0, msg=' 0.02 (BH)', x0= 0.20, y0= 0.18
;process_one_sfr, "sbw/sb13BH", lcolor=200, ctab= 1, lthick= 2.0, msg='20 (BH)', x0= 0.20, y0= 0.30
;process_one_sfr, "sbw/sb13BHtr1", lcolor=140, ctab= 1, lthick= 2.0, msg=' 2 (BH)', x0= 0.20, y0= 0.26
;process_one_sfr, "sbw/sb13BHtr2", lcolor=80, ctab= 1, lthick= 2.0, msg=' 0.2 (BH)', x0= 0.20, y0= 0.22
;process_one_sfr, "sbw/sb13BHtr3", lcolor=20, ctab= 1, lthick= 2.0, msg=' 0.02 (BH)', x0= 0.20, y0= 0.18
;process_one_sfr, "sbw/sb22BH", lcolor=200, ctab= 1, lthick= 2.0, msg='20 (BH)', x0= 0.20, y0= 0.30
;process_one_sfr, "sbw/sb22BHtr1", lcolor=140, ctab= 1, lthick= 2.0, msg=' 2 (BH)', x0= 0.20, y0= 0.26
;process_one_sfr, "sbw/sb22BHtr2", lcolor=80, ctab= 1, lthick= 2.0, msg=' 0.2 (BH)', x0= 0.20, y0= 0.22
;process_one_sfr, "sbw/sb22BHtr3", lcolor=20, ctab= 1, lthick= 2.0, msg=' 0.02 (BH)', x0= 0.20, y0= 0.18
;process_one_sfr, "sbw/sb8BH", lcolor=200, ctab= 1, lthick= 2.0, msg='20 (BH)', x0= 0.20, y0= 0.30
;process_one_sfr, "sbw/sb8BHtr1", lcolor=140, ctab= 1, lthick= 2.0, msg=' 2 (BH)', x0= 0.20, y0= 0.26
;process_one_sfr, "sbw/sb8BHtr2", lcolor=80, ctab= 1, lthick= 2.0, msg=' 0.2 (BH)', x0= 0.20, y0= 0.22
;process_one_sfr, "sbw/sb8BHtr3", lcolor=20, ctab= 1, lthick= 2.0, msg=' 0.02 (BH)', x0= 0.20, y0= 0.18


;xyouts, 0.20, 0.37, 'WindSpeed', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.22, 0.33, '(km/s)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "sbw/sb10", lcolor=200, ctab= 1, lthick= 2.0, msg='837', x0= 0.20, y0= 0.28
;process_one_sfr, "sbw/sb10BH", lcolor=20, ctab= 1, lthick= 2.0, msg='837 (BH)', x0= 0.20, y0= 0.24
;xyouts, 0.20, 0.18, 'WindEfficiency=0.5', /normal, charthick=3, size=1.33, color=0


;xyouts, 0.20, 0.33, 'Wind Kicks', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "sbw/sb10", lcolor=120, ctab= 1, lthick= 2.0, msg='Axial', x0= 0.20, y0= 0.28
;process_one_sfr, "sbw/sb10I", lcolor=20, ctab= 1, lthick= 2.0, msg='Isotropic', x0= 0.20, y0= 0.24
;process_one_sfr, "sbw/sb10tr2", lcolor=120, ctab= 1, lthick= 2.0, msg='Axial', x0= 0.20, y0= 0.28
;process_one_sfr, "sbw/sb10tr2I", lcolor=20, ctab= 1, lthick= 2.0, msg='Isotropic', x0= 0.20, y0= 0.24
;process_one_sfr, "sbw/sb10BHtr2", lcolor=120, ctab= 1, lthick= 2.0, msg='Axial', x0= 0.20, y0= 0.28
;process_one_sfr, "sbw/sb10BHtr2I", lcolor=20, ctab= 1, lthick= 2.0, msg='Isotropic', x0= 0.20, y0= 0.24




; no bh feedback models
;-----------------------
;process_one_sfr, "sbw/sb10", lcolor=200, ctab= 1, lthick= 2.0, msg='no BH', x0= 0.75, y0= 0.88
;process_one_sfr, "sbw/sb10BH", lcolor=20, ctab= 1, lthick= 2.0, msg='BH', x0= 0.75, y0= 0.84
;process_one_sfr, "sbw/sb10BHnofb", lcolor=100, ctab= 1, lthick= 2.0, msg='BH, no fb', x0= 0.75, y0= 0.80
;process_one_sfr, "sbw/sb8", lcolor=200, ctab= 1, lthick= 2.0, msg='no BH', x0= 0.75, y0= 0.88
;process_one_sfr, "sbw/sb8BH", lcolor=20, ctab= 1, lthick= 2.0, msg='BH', x0= 0.75, y0= 0.84
;process_one_sfr, "sbw/sb8BHnofb", lcolor=100, ctab= 1, lthick= 2.0, msg='BH, no fb', x0= 0.75, y0= 0.80
;process_one_sfr, "sbw/sb13", lcolor=200, ctab= 1, lthick= 2.0, msg='no BH', x0= 0.75, y0= 0.88
;process_one_sfr, "sbw/sb13BH", lcolor=20, ctab= 1, lthick= 2.0, msg='BH', x0= 0.75, y0= 0.84
;process_one_sfr, "sbw/sb13BHnofb", lcolor=100, ctab= 1, lthick= 2.0, msg='BH, no fb', x0= 0.75, y0= 0.80


;xyouts, 0.20, 0.37, 'WindSpeed', /normal, charthick=3, size=1.33, color=0
;xyouts, 0.22, 0.33, '(km/s)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "sbw/sb13", lcolor=200, ctab= 1, lthick= 2.0, msg='209', x0= 0.20, y0= 0.28
;process_one_sfr, "sbw/sb13BH", lcolor=20, ctab= 1, lthick= 2.0, msg='209 (BH)', x0= 0.20, y0= 0.24
;xyouts, 0.20, 0.18, 'WindEfficiency=2.0', /normal, charthick=3, size=1.33, color=0

;fruns= ["vc3vc3e_2", "vc3vc3e_no", "vc3vc3e_sb8","vc3vc3e_sb12","vc3vc3e_sb13","vc3vc3e_sb14"]
;fruns= ["vc3vc3e_2", "vc3vc3e_no", "vc3vc3e_sb1","vc3vc3e_sb2_test"]
;lbls=['std','no BH','Winds (837)','Winds (1675)']




; sb winds with various mass progenitors
;------------------------------------------
;process_one_sfr, "ds/d0e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88
;process_one_sfr, "sb8_mass/d0e", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84
;process_one_sfr, "sb10_mass/d0e", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.80
;process_one_sfr, "sb13_mass/d0e", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.76

;process_one_sfr, "ds/d1e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88
;process_one_sfr, "sb8_mass/d1e", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84
;process_one_sfr, "sb10_mass/d1e", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.80
;process_one_sfr, "sb13_mass/d1e", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.76

;process_one_sfr, "ds/d2e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88
;process_one_sfr, "sb8_mass/d2e", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84
;process_one_sfr, "sb10_mass/d2e", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.80
;process_one_sfr, "sb13_mass/d2e", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.76

;process_one_sfr, "ds/d4e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88
;process_one_sfr, "sb8_mass/d4e", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84
;process_one_sfr, "sb10_mass/d4e", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.80
;process_one_sfr, "sb13_mass/d4e", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.76

;process_one_sfr, "ds/d5e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88
;process_one_sfr, "sb8_mass/d5e", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84
;process_one_sfr, "sb10_mass/d5e", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.80
;process_one_sfr, "sb13_mass/d5e", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.76

;process_one_sfr, "ds/d6e2_q", lcolor=0, lthick= 2.0, msg='nosb', x0= 0.75, y0= 0.88
;process_one_sfr, "sb8_mass/d6e", lcolor=200, lthick= 2.0, msg='sb8', x0= 0.75, y0= 0.84
;process_one_sfr, "sb10_mass/d6e", lcolor=50, lthick= 2.0, msg='sb10', x0= 0.75, y0= 0.80
;process_one_sfr, "sb13_mass/d6e", lcolor=100, lthick= 2.0, msg='sb13', x0= 0.75, y0= 0.76


;process_one_sfr, "ds/d6e", lcolor=100, lthick= 2.0, msg='d6e', x0= 0.75, y0= 0.88
;process_one_sfr, "ds/d6e2", lcolor=50, lthick= 2.0, msg='d6e2', x0= 0.75, y0= 0.84
;process_one_sfr, "ds/d6e2_q", lcolor=150, lthick= 2.0, msg='d6e2_q', x0= 0.75, y0= 0.80

;process_one_sfr, "bs/b5e", lcolor=0, lthick= 2.0, msg='b5e', x0= 0.75, y0= 0.92
;process_one_sfr, "ds/d5e", lcolor=100, lthick= 2.0, msg='d5e', x0= 0.75, y0= 0.88
;process_one_sfr, "ds/d5e2", lcolor=50, lthick= 2.0, msg='d5e2', x0= 0.75, y0= 0.84
;process_one_sfr, "ds/d5e2_q", lcolor=150, lthick= 2.0, msg='d5e2_q', x0= 0.75, y0= 0.80

;process_one_sfr, "bs/b4e", lcolor=0, lthick= 2.0, msg='b4e', x0= 0.75, y0= 0.92
;process_one_sfr, "ds/d4e", lcolor=100, lthick= 2.0, msg='d4e', x0= 0.75, y0= 0.88
;process_one_sfr, "ds/d4e2", lcolor=50, lthick= 2.0, msg='d4e2', x0= 0.75, y0= 0.84
;process_one_sfr, "ds/d4e2_q", lcolor=150, lthick= 2.0, msg='d4e2_q', x0= 0.75, y0= 0.80





;--- xxxxx  55555  xxxxx ---
;process_one_sfr, "bs/b5e", lcolor=0, lthick= 2.0, msg='b5e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb10_mass/b5e", lcolor=100, lthick= 2.0, msg='sb10/b5e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb10_mass/b5e_no", lcolor=150, lthick= 2.0, msg='sb10/b5e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d5e2_q", lcolor=0, lthick= 2.0, msg='d5e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb10_mass/d5e", lcolor=100, lthick= 2.0, msg='sb10/d5e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb10_mass/d5e_no", lcolor=150, lthick= 2.0, msg='sb10/d5e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b5e", lcolor=0, lthick= 2.0, msg='b5e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb13_mass/b5e", lcolor=100, lthick= 2.0, msg='sb13/b5e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb13_mass/b5e_no", lcolor=150, lthick= 2.0, msg='sb13/b5e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d5e2_q", lcolor=0, lthick= 2.0, msg='d5e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb13_mass/d5e", lcolor=100, lthick= 2.0, msg='sb13/d5e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb13_mass/d5e_no", lcolor=150, lthick= 2.0, msg='sb13/d5e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b5e", lcolor=0, lthick= 2.0, msg='b5e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb8_mass/b5e", lcolor=100, lthick= 2.0, msg='sb8/b5e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb8_mass/b5e_no", lcolor=150, lthick= 2.0, msg='sb8/b5e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d5e2_q", lcolor=0, lthick= 2.0, msg='d5e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb8_mass/d5e", lcolor=100, lthick= 2.0, msg='sb8/d5e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb8_mass/d5e_no", lcolor=150, lthick= 2.0, msg='sb8/d5e_no', x0= 0.70, y0= 0.80


;process_one_sfr, "bs/b5e", lcolor=0, lthick= 2.0, msg='b5e', x0= 0.70, y0= 0.92
;process_one_sfr, "bs/b5e_no", lcolor=50, lthick= 2.0, msg='b5e_no', x0= 0.70, y0= 0.88
;process_one_sfr, "bs/b5e_no2", lcolor=80, lthick= 2.0, msg='b5e_no2', x0= 0.70, y0= 0.84
;process_one_sfr, "bs/b5e_no3", lcolor=110, lthick= 2.0, msg='b5e_no3', x0= 0.70, y0= 0.80
;process_one_sfr, "bs/b5e_no4", lcolor=140, lthick= 2.0, msg='b5e_no4', x0= 0.70, y0= 0.76
;process_one_sfr, "bs/b5e_no5", lcolor=170, lthick= 2.0, msg='b5e_no5', x0= 0.70, y0= 0.72
;process_one_sfr, "bs/b5e_no6", lcolor=200, lthick= 2.0, msg='b5e_no6', x0= 0.70, y0= 0.68




;--- xxxxx  33333  xxxxx ---
;process_one_sfr, "bs/b3e", lcolor=0, lthick= 2.0, msg='b3e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb10_mass/b3e", lcolor=100, lthick= 2.0, msg='sb10/b3e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb10_mass/b3e_no", lcolor=150, lthick= 2.0, msg='sb10/b3e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d3e7", lcolor=0, lthick= 2.0, msg='d3e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb10_mass/d3e", lcolor=100, lthick= 2.0, msg='sb10/d3e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb10_mass/d3e_no", lcolor=150, lthick= 2.0, msg='sb10/d3e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b3e", lcolor=0, lthick= 2.0, msg='b3e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb13_mass/b3e", lcolor=100, lthick= 2.0, msg='sb13/b3e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb13_mass/b3e_no", lcolor=150, lthick= 2.0, msg='sb13/b3e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d3e7", lcolor=0, lthick= 2.0, msg='d3e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb13_mass/d3e", lcolor=100, lthick= 2.0, msg='sb13/d3e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb13_mass/d3e_no", lcolor=150, lthick= 2.0, msg='sb13/d3e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b3e", lcolor=0, lthick= 2.0, msg='b3e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb8_mass/b3e", lcolor=100, lthick= 2.0, msg='sb8/b3e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb8_mass/b3e_no", lcolor=150, lthick= 2.0, msg='sb8/b3e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d3e7", lcolor=0, lthick= 2.0, msg='d3e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb8_mass/d3e", lcolor=100, lthick= 2.0, msg='sb8/d3e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb8_mass/d3e_no", lcolor=150, lthick= 2.0, msg='sb8/d3e_no', x0= 0.70, y0= 0.80

;
;process_one_sfr, "bs/b3e_2", lcolor=0, lthick= 2.0, msg='std', x0= 0.75, y0= 0.88
;process_one_sfr, "bs/b3e_no", lcolor=100, lthick= 2.0, msg='no BH', x0= 0.75, y0= 0.84

;
;process_one_sfr, "bs/b3e_no", lcolor=50, lthick= 2.0, msg='b3e_no', x0= 0.75, y0= 0.92
;process_one_sfr, "bs/b3e_2", lcolor=150, lthick= 2.0, msg='b3e_2', x0= 0.75, y0= 0.88
;process_one_sfr, "bs/b3e", lcolor=0, lthick= 2.0, msg='b3e', x0= 0.75, y0= 0.84
;process_one_sfr, "bs/b3e_q", lcolor=200, lthick= 2.0, msg='b3e_q', x0= 0.75, y0= 0.80
;process_one_sfr, "bs/b3e_q2", lcolor=100, lthick= 2.0, msg='b3e_q2', x0= 0.75, y0= 0.76

;------------------------------------

;process_one_sfr, "bs/b3e", lcolor=0, lthick= 2.0, msg=' ', x0= 0.75, y0= 0.84
;process_one_sfr, "bs/b3e_q2", lcolor=150, lthick= 2.0, msg=' ', x0= 0.75, y0= 0.76
;process_one_sfr, "ds/d5e2_q", lcolor=150, lthick=2.0, msg=' ', x0= 0.65, y0= 0.88
;process_one_sfr, "ds/d5e2_2", lcolor=0, lthick=2.0, msg=' ', x0= 0.65, y0= 0.80

;------------------------------------

;xyouts, 0.60, 0.92, 'ErrTolIntAccuracy', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "bs/b3e_i2", lcolor=  0, lthick= 2.0, msg='1.0', x0= 0.75, y0= 0.86
;process_one_sfr, "bs/b3e_i1", lcolor= 50, lthick= 2.0, msg='0.1', x0= 0.75, y0= 0.82
;process_one_sfr, "bs/b3e", lcolor=150, lthick= 2.0, msg='0.025', x0= 0.75, y0= 0.78
;process_one_sfr, "bs/b3e_i3", lcolor=100, lthick= 2.0, msg='0.0025', x0= 0.75, y0= 0.74


;process_one_sfr, "bs/b3e", lcolor=150, lthick= 2.0, msg='0.025', x0= 0.75, y0= 0.78
;process_one_sfr, "bs/b3e_i3", lcolor=100, lthick= 2.0, msg='0.0025', x0= 0.75, y0= 0.74
;process_one_sfr, "bs/b3e_i4", lcolor=50, lthick= 2.0, msg='P-Gadget3.2', x0= 0.75, y0= 0.70

;--- xxxxx  11111  xxxxx ---
;process_one_sfr, "bs/b1e", lcolor=0, lthick= 2.0, msg='b1e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb10_mass/b1e", lcolor=100, lthick= 2.0, msg='sb10/b1e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb10_mass/b1e_no", lcolor=150, lthick= 2.0, msg='sb10/b1e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d1e2_q", lcolor=0, lthick= 2.0, msg='d1e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb10_mass/d1e", lcolor=100, lthick= 2.0, msg='sb10/d1e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb10_mass/d1e_no", lcolor=150, lthick= 2.0, msg='sb10/d1e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b1e", lcolor=0, lthick= 2.0, msg='b1e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb13_mass/b1e", lcolor=100, lthick= 2.0, msg='sb13/b1e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb13_mass/b1e_no", lcolor=150, lthick= 2.0, msg='sb13/b1e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d1e2_q", lcolor=0, lthick= 2.0, msg='d1e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb13_mass/d1e", lcolor=100, lthick= 2.0, msg='sb13/d1e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb13_mass/d1e_no", lcolor=150, lthick= 2.0, msg='sb13/d1e_no', x0= 0.70, y0= 0.80
;---
;process_one_sfr, "bs/b1e", lcolor=0, lthick= 2.0, msg='b1e', x0= 0.70, y0= 0.88
;process_one_sfr, "sb8_mass/b1e", lcolor=100, lthick= 2.0, msg='sb8/b1e', x0= 0.70, y0= 0.84
;process_one_sfr, "sb8_mass/b1e_no", lcolor=150, lthick= 2.0, msg='sb8/b1e_no', x0= 0.70, y0= 0.80
; process_one_sfr, "ds/d1e2_q", lcolor=0, lthick= 2.0, msg='d1e', x0= 0.70, y0= 0.88
; process_one_sfr, "sb8_mass/d1e", lcolor=100, lthick= 2.0, msg='sb8/d1e', x0= 0.70, y0= 0.84
; process_one_sfr, "sb8_mass/d1e_no", lcolor=150, lthick= 2.0, msg='sb8/d1e_no', x0= 0.70, y0= 0.80




; -------------------------------------
;  Volker's Runs

;fruns= ["aA1","aA1w"] & lbls=['coll_1/vc1/nbh/','coll_1/vc1/bh']
;fruns= ["aA2","aA2w"] & lbls=['coll_1/vc2/nbh/','coll_1/vc2/bh']
;fruns= ["aA3","aA3w"] & lbls=['coll_1/vc3/nbh/','coll_1/vc3/bh']
;fruns= ["aA4","aA4w"] & lbls=['coll_1/vc4/nbh/','coll_1/vc4/bh']

;fruns= ["bA1","bA1w"] & lbls=['coll_2/vc1/nbh/','coll_2/vc1/bh']
;fruns= ["bA2","bA2w"] & lbls=['coll_2/vc2/nbh/','coll_2/vc2/bh']
;fruns= ["bA3","bA3w"] & lbls=['coll_2/vc3/nbh/','coll_2/vc3/bh']
;fruns= ["bA4","bA4w"] & lbls=['coll_2/vc4/nbh/','coll_2/vc4/bh']

;msg=''


;fruns= ["As/A3", "As/A3_10x"]
;fruns= ["As/A3_1x", "As/A3_10x", "As/A3_10xa"]
;fruns= ["As/A3", "As/A3_1x", "As/A3_10x"]
;fruns= ["As/A3", "As/A3nopot"]
;fruns= ["As/A3", "As/A3duprr1"]
msg= ' '




; -------------------------------------
;  minor merger

;process_one_sfr, "minor/circ2", lcolor=100, lthick= 2.0, msg='circ2', x0= 0.70, y0= 0.84
;process_one_sfr, "minor/circ3", lcolor=150, lthick= 2.0, msg='circ3', x0= 0.70, y0= 0.80
;process_one_sfr, "minor/circ4", lcolor=0, lthick= 2.0, msg='circ4', x0= 0.70, y0= 0.88


;process_one_sfr, "isolated/G3", lcolor=0, lthick= 2.0, msg='G3', x0= 0.70, y0= 0.92

;process_one_sfr, "isolated/G3bl", lcolor=50, lthick= 2.0, msg='G3bl', x0= 0.70, y0= 0.88
;process_one_sfr, "isolated/G3blv6", lcolor=100, lthick= 2.0, msg='G3blv6', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/G3blv7", lcolor=150, lthick= 2.0, msg='G3blv7', x0= 0.70, y0= 0.80

;process_one_sfr, "isolated/G3Rd1", lcolor=150, lthick= 2.0, msg='G3Rd1', x0= 0.70, y0= 0.88
;process_one_sfr, "isolated/G3Rd1a", lcolor=200, lthick= 2.0, msg='G3Rd1a', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/G3Rd1b", lcolor=50, lthick= 2.0, msg='G3Rd1b', x0= 0.70, y0= 0.80

;process_one_sfr, "isolated/G3Rd2", lcolor=150, lthick= 2.0, msg='G3Rd2', x0= 0.70, y0= 0.88
;process_one_sfr, "isolated/G3Rd2a", lcolor=200, lthick= 2.0, msg='G3Rd2a', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/G3Rd2b", lcolor=50, lthick= 2.0, msg='G3Rd2b', x0= 0.70, y0= 0.80

;process_one_sfr, "isolated/G3Rd3", lcolor=150, lthick= 2.0, msg='G3Rd3', x0= 0.70, y0= 0.88
;process_one_sfr, "isolated/G3Rd3a", lcolor=200, lthick= 2.0, msg='G3Rd3a', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/G3Rd3b", lcolor=50, lthick= 2.0, msg='G3Rd3b', x0= 0.70, y0= 0.80

;process_one_sfr, "isolated/G3Rd4", lcolor=150, lthick= 2.0, msg='G3Rd4', x0= 0.70, y0= 0.88
;process_one_sfr, "isolated/G3Rd4a", lcolor=200, lthick= 2.0, msg='G3Rd4a', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/G3Rd4b", lcolor=50, lthick= 2.0, msg='G3Rd4b', x0= 0.70, y0= 0.80


; these look pretty good!
;process_one_sfr, "isolated/G3Rd4d", lcolor=50, lthick= 2.0, msg='G3Rd4d', x0= 0.70, y0= 0.88
;process_one_sfr, "isolated/G3Rd4e", lcolor=150, lthick= 2.0, msg='G3Rd4e', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/G3Rd4e_n0", lcolor=100, lthick= 2.0, msg='G3Rd4e_n0', x0= 0.70, y0= 0.80

;process_one_sfr, "isolated/G3Rd4e", lcolor=50, lthick= 2.0, msg='G3Rd4e', x0= 0.70, y0= 0.88
;process_one_sfr, "isolated/G3Rd4f", lcolor=150, lthick= 2.0, msg='G3Rd4f', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/G3Rd4g", lcolor=100, lthick= 2.0, msg='G3Rd4g', x0= 0.70, y0= 0.80
;process_one_sfr, "isolated/G3Rd4h", lcolor=200, lthick= 2.0, msg='G3Rd4h', x0= 0.70, y0= 0.76
;
;process_one_sfr, "isolated/G3Rd4e_n0", lcolor=50, lthick= 2.0, msg='G3Rd4e_n0', x0= 0.70, y0= 0.88
;process_one_sfr, "isolated/G3Rd4f_n0", lcolor=150, lthick= 2.0, msg='G3Rd4f_n0', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/G3Rd4g_n0", lcolor=100, lthick= 2.0, msg='G3Rd4g_n0', x0= 0.70, y0= 0.80
;process_one_sfr, "isolated/G3Rd4h_n0", lcolor=200, lthick= 2.0, msg='G3Rd4h_n0', x0= 0.70, y0= 0.76

;process_one_sfr, "isolated/G3Rd4e", lcolor=0, lthick= 2.0, msg='G3Rd4e', x0= 0.70, y0= 0.92
;process_one_sfr, "isolated/G3Rd4e_bd1", lcolor=150, lthick= 2.0, msg='G3Rd4e_bd1', x0= 0.70, y0= 0.88
;process_one_sfr, "isolated/G3Rd4e_bd2", lcolor=100, lthick= 2.0, msg='G3Rd4e_bd2', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/G3Rd4e_bd3", lcolor=200, lthick= 2.0, msg='G3Rd4e_bd3', x0= 0.70, y0= 0.80
;process_one_sfr, "isolated/G3Rd4e_bd4", lcolor=50, lthick= 2.0, msg='G3Rd4e_bd4', x0= 0.70, y0= 0.76


; mergers
;process_one_sfr, "G3/G3Rd4G3Rd4", lcolor=50, lthick= 2.0, msg='G3Rd4G3Rd4', x0= 0.70, y0= 0.88
;process_one_sfr, "G3/G3Rd4G2", lcolor=150, lthick= 2.0, msg='G3Rd4G2', x0= 0.70, y0= 0.84
;process_one_sfr, "G3/G3Rd4G1", lcolor=100, lthick= 2.0, msg='G3Rd4G1', x0= 0.70, y0= 0.80
;process_one_sfr, "G3/G3Rd4G0", lcolor=200, lthick= 2.0, msg='G3Rd4G0', x0= 0.70, y0= 0.76

;process_one_sfr, "G3/G3Rd4eG3", lcolor=50, lthick= 2.0, msg='G3Rd4eG3', x0= 0.70, y0= 0.88
;process_one_sfr, "G3/G3Rd4eG2", lcolor=150, lthick= 2.0, msg='G3Rd4eG2', x0= 0.70, y0= 0.84
;process_one_sfr, "G3/G3Rd4eG1", lcolor=100, lthick= 2.0, msg='G3Rd4eG1', x0= 0.70, y0= 0.80
;process_one_sfr, "G3/G3Rd4eG0", lcolor=200, lthick= 2.0, msg='G3Rd4eG0', x0= 0.70, y0= 0.76

;process_one_sfr, "G3/G3Rd4eG3_n0", lcolor=50, lthick= 2.0, msg='G3Rd4eG3_n0', x0= 0.70, y0= 0.88
;process_one_sfr, "G3/G3Rd4eG2_n0", lcolor=150, lthick= 2.0, msg='G3Rd4eG2_n0', x0= 0.70, y0= 0.84
;process_one_sfr, "G3/G3Rd4eG1_n0", lcolor=100, lthick= 2.0, msg='G3Rd4eG1_n0', x0= 0.70, y0= 0.80
;process_one_sfr, "G3/G3Rd4eG0_n0", lcolor=200, lthick= 2.0, msg='G3Rd4eG0_n0', x0= 0.70, y0= 0.76

;process_one_sfr, "G3/G3Rd4eG2", lcolor=50, lthick= 2.0, msg='G3Rd4eG2', x0= 0.70, y0= 0.92
;process_one_sfr, "G3/G3Rd4fG2", lcolor=100, lthick= 2.0, msg='G3Rd4fG2', x0= 0.70, y0= 0.88
;process_one_sfr, "G3/G3Rd4gG2", lcolor=150, lthick= 2.0, msg='G3Rd4gG2', x0= 0.70, y0= 0.84
;process_one_sfr, "G3/G3Rd4hG2", lcolor=200, lthick= 2.0, msg='G3Rd4hG2', x0= 0.70, y0= 0.80

;process_one_sfr, "G3/G3bd3G3", lcolor=50, lthick= 2.0, msg='G3bd3G3', x0= 0.70, y0= 0.88
;process_one_sfr, "G3/G3bd3G2", lcolor=150, lthick= 2.0, msg='G3bd3G2', x0= 0.70, y0= 0.84
;process_one_sfr, "G3/G3bd3G1", lcolor=100, lthick= 2.0, msg='G3bd3G1', x0= 0.70, y0= 0.80
;process_one_sfr, "G3/G3bd3G0", lcolor=200, lthick= 2.0, msg='G3bd3G0', x0= 0.70, y0= 0.76

;process_one_sfr, "G3/G3bd3G3_n0", lcolor=50, lthick= 2.0, msg='G3bd3G3_n0', x0= 0.70, y0= 0.88
;process_one_sfr, "G3/G3bd3G2_n0", lcolor=150, lthick= 2.0, msg='G3bd3G2_n0', x0= 0.70, y0= 0.84
;process_one_sfr, "G3/G3bd3G1_n0", lcolor=100, lthick= 2.0, msg='G3bd3G1_n0', x0= 0.70, y0= 0.80
;process_one_sfr, "G3/G3bd3G0_n0", lcolor=200, lthick= 2.0, msg='G3bd3G0_n0', x0= 0.70, y0= 0.76


;process_one_sfr, "G3/G3bd3G1", lcolor=50, lthick= 2.0, msg='G3bd3G1', x0= 0.70, y0= 0.92
;process_one_sfr, "G3/G3bd3G1a", lcolor=100, lthick= 2.0, msg='G3bd3G1a', x0= 0.70, y0= 0.88
;process_one_sfr, "G3/G3bd3G1_n0", lcolor=150, lthick= 2.0, msg='G3bd3G1_n0', x0= 0.70, y0= 0.84
;process_one_sfr, "G3/G3bd3G1a_n0", lcolor=200, lthick= 2.0, msg='G3bd3G1a_n0', x0= 0.70, y0= 0.80


;process_one_sfr, "G3/G3Rd4eG1", lcolor=50, lthick= 2.0, msg='G3Rd4eG1', x0= 0.70, y0= 0.92
;process_one_sfr, "G3/G3Rd4eG1a", lcolor=100, lthick= 2.0, msg='G3Rd4eG1a', x0= 0.70, y0= 0.88
;process_one_sfr, "G3/G3Rd4eG1_n0", lcolor=150, lthick= 2.0, msg='G3Rd4eG1_n0', x0= 0.70, y0= 0.84
;process_one_sfr, "G3/G3Rd4eG1a_n0", lcolor=200, lthick= 2.0, msg='G3Rd4eG1a_n0', x0= 0.70, y0= 0.80


;process_one_sfr, "Gs/G3gf2G3gf2", lcolor=50, lthick= 2.0, msg='G3gf2G3gf2', x0= 0.70, y0= 0.92
;process_one_sfr, "Gs/G3gf2G3gf2_a", lcolor=150, lthick= 2.0, msg='G3gf2G3gf2_a', x0= 0.70, y0= 0.88
;process_one_sfr, "Gs/G3G3", lcolor=100, lthick= 2.0, msg='G3G3', x0= 0.70, y0= 0.84      ; b.s. orbit
;process_one_sfr, "Gs/G3G3b-u2", lcolor=200, lthick= 2.0, msg='G3G3b-u2', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/G3gf2_a", lcolor=0, lthick= 2.0, msg='G3gf2_a (iso)', x0= 0.70, y0= 0.80



; -------------------------------------


;process_one_sfr, "minor/min_0", lcolor=50, lthick= 2.0, msg='min_0', x0= 0.70, y0= 0.92, h=h
;process_one_sfr, "minor/min_30", lcolor=80, lthick= 2.0, msg='min_30', x0= 0.70, y0= 0.88, h=h
;process_one_sfr, "minor/min_90", lcolor=110, lthick= 2.0, msg='min_90', x0= 0.70, y0= 0.84, h=h
;process_one_sfr, "minor/min_150", lcolor=140, lthick= 2.0, msg='min_150', x0= 0.70, y0= 0.80, h=h
;process_one_sfr, "minor/min_180", lcolor=170, lthick= 2.0, msg='min_180', x0= 0.70, y0= 0.76, h=h

;process_one_sfr, "minor/iso_Sb", lcolor=150, lthick= 2.0, msg='iso_Sb', x0= 0.70, y0= 0.92, h=h

;process_one_sfr, "minor/lr_min_0", lcolor=100, lthick= 2.0, msg='lr_min_0', x0= 0.70, y0= 0.88, h=h
;process_one_sfr, "minor/bd1_lr_min_0", lcolor=150, lthick= 2.0, msg='bd1_lr_min_0', x0= 0.70, y0= 0.84, h=h
;process_one_sfr, "minor/bd2_lr_min_0", lcolor=200, lthick= 2.0, msg='bd2_lr_min_0', x0= 0.70, y0= 0.80, h=h

;process_one_sfr, "minor/min_30", lcolor=150, lthick= 2.0, msg=' ', x0= 0.70, y0= 0.88, h=h

;process_one_sfr, "minor/min_0",        lcolor=10,  lthick= 1.0, msg=' ', x0= 0.70, y0= 0.92, h=h, ctab=1
;process_one_sfr, "minor/min_30",       lcolor=20,  lthick= 1.0, msg=' ', x0= 0.70, y0= 0.88, h=h, ctab=1
;process_one_sfr, "minor/min_90",       lcolor=30,  lthick= 1.0, msg=' ', x0= 0.70, y0= 0.84, h=h, ctab=1
;process_one_sfr, "minor/min_150",      lcolor=40,  lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/min_180",      lcolor=50,  lthick= 1.0, msg=' ', x0= 0.70, y0= 0.76, h=h, ctab=1
;process_one_sfr, "minor/iso_Sb",       lcolor=60,  lthick= 1.0, msg=' ', x0= 0.70, y0= 0.92, h=h, ctab=1
;process_one_sfr, "minor/lr_min_0",     lcolor=70,  lthick= 1.0, msg=' ', x0= 0.70, y0= 0.88, h=h, ctab=1
;process_one_sfr, "minor/bd1_lr_min_0", lcolor=80,  lthick= 1.0, msg=' ', x0= 0.70, y0= 0.84, h=h, ctab=1
;process_one_sfr, "minor/bd2_lr_min_0", lcolor=90,  lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.4Sdfg0.4_000", lcolor=100, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.4Sdfg0.4_030", lcolor=110, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.4Sdfg0.4_090", lcolor=120, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.4Sdfg0.4_150", lcolor=130, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.4Sdfg0.4_180", lcolor=140, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.4Imfg0.4_000", lcolor=150, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.4Imfg0.4_030", lcolor=160, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.4Imfg0.4_090", lcolor=170, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.4Imfg0.4_150", lcolor=180, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.4Imfg0.4_180", lcolor=190, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.8Sdfg0.8_000", lcolor=200, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.8Sdfg0.8_030", lcolor=205, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.8Sdfg0.8_090", lcolor=210, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.8Sdfg0.8_150", lcolor=215, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.8Sdfg0.8_180", lcolor=220, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.8Imfg0.8_000", lcolor=225, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.8Imfg0.8_030", lcolor=230, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.8Imfg0.8_090", lcolor=235, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.8Imfg0.8_150", lcolor=240, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1
;process_one_sfr, "minor/Sbfg0.8Imfg0.8_180", lcolor=245, lthick= 1.0, msg=' ', x0= 0.70, y0= 0.80, h=h, ctab=1



; these have orbits similar to G3G1, not min_0
;process_one_sfr, "G3/SbIm_o1", lcolor=150, lthick= 2.0, msg='SbIm_o1', x0= 0.70, y0= 0.84
;process_one_sfr, "G3/SbIm_o1_n0", lcolor=200, lthick= 2.0, msg='SbIm_o1_n0', x0= 0.70, y0= 0.80
;process_one_sfr, "G3/SbIm_o2", lcolor=50, lthick= 2.0, msg='SbIm_o2', x0= 0.70, y0= 0.76
;process_one_sfr, "G3/SbIm_o2_n0", lcolor=10, lthick= 2.0, msg='SbIm_o2_n0', x0= 0.70, y0= 0.72

;process_one_sfr, "G3/SbSb_o2_n0", lcolor=150, lthick= 2.0, msg='SbSb_o2_n0', x0= 0.70, y0= 0.84
;process_one_sfr, "G3/SbG2_o2_n0", lcolor=200, lthick= 2.0, msg='SbG2_o2_n0', x0= 0.70, y0= 0.80
;process_one_sfr, "G3/SbG1_o2_n0", lcolor=50, lthick= 2.0, msg='SbG1_o2_n0', x0= 0.70, y0= 0.76
;process_one_sfr, "G3/SbG0_o2_n0", lcolor=10, lthick= 2.0, msg='SbG0_o2_n0', x0= 0.70, y0= 0.72


; -------------------------------------

; substr

;process_one_sfr, "substr/vc3c_1", lcolor=0, lthick= 2.0, msg='vc3c_1', x0= 0.70, y0= 0.84
;process_one_sfr, "substr/t_23_0", lcolor=100, lthick= 2.0, msg='t_23_0', x0= 0.70, y0= 0.80
;process_one_sfr, "substr/t_23_1", lcolor=150, lthick= 2.0, msg='t_23_1', x0= 0.70, y0= 0.76

; -------------------------------------

; compare to enzo

;process_one_sfr, "enzocomp/lr_twogals", lcolor=40, lthick= 1.0, lstyle= 2, msg='lr_twogals', x0= 0.70, y0= 0.88, h=h
;process_one_sfr, "enzocomp/lr_2gal_igm1", lcolor=60, lthick= 3.0, msg='lr_2gal_igm1', x0= 0.70, y0= 0.84, h=h
;process_one_sfr, "enzocomp/hr_twogals", lcolor=140, lthick= 1.0, lstyle= 2, msg='hr_twogals', x0= 0.70, y0= 0.80, h=h
;process_one_sfr, "enzocomp/hr_2gal_igm1", lcolor=160, lthick= 3.0, msg='hr_2gal_igm1', x0= 0.70, y0= 0.76, h=h


; -------------------------------------

; shock-sf models

;process_one_sfr, "isolated/vc3c", lcolor=150, lthick= 2.0, msg='vc3c', x0= 0.70, y0= 0.84
;process_one_sfr, "isolated/d3", lcolor=200, lthick= 2.0, msg='d3', x0= 0.70, y0= 0.80
;process_one_sfr, "ssft/iso", lcolor=50, lthick= 2.0, msg='iso', x0= 0.70, y0= 0.76
;process_one_sfr, "ssft/iso1a", lcolor=10, lthick= 2.0, msg='iso1a', x0= 0.70, y0= 0.72

;process_one_sfr, "ssft/iso", lcolor=200, lthick= 2.0, msg='iso', x0= 0.70, y0= 0.88
;process_one_sfr, "ssft/iso1a", lcolor=150, lthick= 2.0, msg='iso1a', x0= 0.70, y0= 0.84
;process_one_sfr, "ssft/iso1b", lcolor=100, lthick= 2.0, msg='iso1b', x0= 0.70, y0= 0.80

;process_one_sfr, "ssft/iso2a", lcolor=200, lthick= 2.0, msg='iso2a', x0= 0.70, y0= 0.88
;process_one_sfr, "ssft/iso2b", lcolor=150, lthick= 2.0, msg='iso2b', x0= 0.70, y0= 0.84
;process_one_sfr, "ssft/iso2c", lcolor=100, lthick= 2.0, msg='iso2c', x0= 0.70, y0= 0.80

;process_one_sfr, "ssft/iso3a", lcolor=200, lthick= 2.0, msg='iso3a', x0= 0.70, y0= 0.88
;process_one_sfr, "ssft/iso3b", lcolor=150, lthick= 2.0, msg='iso3b', x0= 0.70, y0= 0.84
;process_one_sfr, "ssft/iso3c", lcolor=100, lthick= 2.0, msg='iso3c', x0= 0.70, y0= 0.80


; -------------------------------------

; remnant-remnant mergers

; vc3 remnants
;fruns= ["vc3vc3_wBH", "vc3rem_vc3","vc3rem_vc3rem"]
;	tmreg=[1.1,1.1,1.1]
;	msg=''


;fruns= ["vc3vc3", "vc3rem_vc3","vc3rem_vc2"]
;fruns= ["vc3vc3", "vc3rem_vc3","vc3rem_vc2","vc3rem_vc1","vc3rem_vc1a","vc3rem_vc1b","vc3rem_vc1c"]
;	tmreg=[1.1,1.1,1.1]
;	msg=''



;process_one_sfr, "ds/vc3vc3f", lcolor=50, lthick= 2.0, msg='f', x0= 0.70, y0= 0.92
;process_one_sfr, "remergers/vc3fvc3f", lcolor=200, lthick= 2.0, msg='fxf', x0= 0.70, y0= 0.88
;process_one_sfr, "remergers/vc3fvc3f_polar", lcolor=150, lthick= 2.0, msg='fxf_polar', x0= 0.70, y0= 0.84
;process_one_sfr, "remergers/vc3mvc3m", lcolor=100, lthick= 2.0, msg='m', x0= 0.70, y0= 0.84
;process_one_sfr, "remergers/vc3fvc3f_polar", lcolor=50, lthick= 2.0, msg='m_polar', x0= 0.70, y0= 0.80

;process_one_sfr, "remergers/b2eb2e", lcolor=200, lthick= 2.0, msg='b2eb2e', x0= 0.70, y0= 0.88
;process_one_sfr, "remergers/b2eb2h", lcolor=180, lthick= 2.0, msg='b2eb2h', x0= 0.70, y0= 0.85
;process_one_sfr, "remergers/b2eb2k", lcolor=160, lthick= 2.0, msg='b2eb2k', x0= 0.70, y0= 0.82

;process_one_sfr, "remergers/b3eb3e", lcolor=120, lthick= 2.0, msg='b3eb3e', x0= 0.70, y0= 0.78
;process_one_sfr, "remergers/b3eb3f", lcolor=100, lthick= 2.0, msg='b3eb3f', x0= 0.70, y0= 0.75
;process_one_sfr, "remergers/b3eb3h", lcolor=80, lthick= 2.0, msg='b3eb3h', x0= 0.70, y0= 0.72

;process_one_sfr, "remergers/b4eb3e", lcolor=60, lthick= 2.0, msg='b4eb3e', x0= 0.70, y0= 0.68
;process_one_sfr, "remergers/b4eb4e", lcolor=50, lthick= 2.0, msg='b4eb4e', x0= 0.70, y0= 0.65
;process_one_sfr, "remergers/b4eb4k", lcolor=40, lthick= 2.0, msg='b4eb4k', x0= 0.70, y0= 0.62
;process_one_sfr, "remergers/b4fb4h", lcolor=30, lthick= 2.0, msg='b4fb4h', x0= 0.70, y0= 0.59

;process_one_sfr, "remergers/iso_b3e", lcolor=0, lthick= 2.0, msg='iso_b3e', x0= 0.70, y0= 0.88, h=h
;process_one_sfr, "remergers/iso_b3e_2", lcolor=0, lthick= 2.0, msg='iso_b3e_2', x0= 0.70, y0= 0.84, h=h
;process_one_sfr, "remergers/iso_sph", lcolor=50, lthick= 2.0, msg='iso_sph', x0= 0.70, y0= 0.80, h=h
;process_one_sfr, "remergers/iso_v4x3", lcolor=150, lthick= 2.0, msg='iso_v4x3', x0= 0.70, y0= 0.76, h=h
;process_one_sfr, "remergers/iso_v3x2", lcolor=100, lthick= 2.0, msg='iso_v3x2', x0= 0.70, y0= 0.72, h=h
;process_one_sfr, "remergers/iso_v2v2", lcolor=200, lthick= 2.0, msg='iso_v2v2', x0= 0.70, y0= 0.68, h=h

;process_one_sfr, "remergers/b3eb3e2", lcolor=0, lthick= 2.0, msg='b3eb3e2', x0= 0.70, y0= 0.78, h=h
;process_one_sfr, "remergers/b3esph", lcolor=50, lthick= 2.0, msg='b3esph', x0= 0.70, y0= 0.74, h=h
;process_one_sfr, "remergers/b3ev4", lcolor=150, lthick= 2.0, msg='b3ev4', x0= 0.70, y0= 0.70, h=h
;process_one_sfr, "remergers/b3ev3", lcolor=100, lthick= 2.0, msg='b3ev3', x0= 0.70, y0= 0.66, h=h
;process_one_sfr, "remergers/b3ev2", lcolor=200, lthick= 2.0, msg='b3ev2', x0= 0.70, y0= 0.62, h=h


; -------------------------------------


;process_one_sfr, "priya/v3b1v3b1_e", lcolor=20, lthick= 2.0, msg='v3b1v3b1_e', x0= 0.70, y0= 0.94
;process_one_sfr, "priya/v3b1v3b1_f", lcolor=20, lthick= 2.0, msg='v3b1v3b1_f', x0= 0.70, y0= 0.94
;process_one_sfr, "priya/v3b1v3b1_k", lcolor=20, lthick= 2.0, msg='v3b1v3b1_k', x0= 0.70, y0= 0.94

;process_one_sfr, "priya/v3b1v3b4_e", lcolor=40, lthick= 2.0, msg='v3b1v3b4_e', x0= 0.70, y0= 0.91
;process_one_sfr, "priya/v3b1v3b4_f", lcolor=40, lthick= 2.0, msg='v3b1v3b4_f', x0= 0.70, y0= 0.91
;process_one_sfr, "priya/v3b1v3b4_k", lcolor=40, lthick= 2.0, msg='v3b1v3b4_k', x0= 0.70, y0= 0.91

;process_one_sfr, "priya/v3b2v3b1_e", lcolor=60, lthick= 2.0, msg='v3b2v3b1_e', x0= 0.70, y0= 0.88
;process_one_sfr, "priya/v3b2v3b1_f", lcolor=60, lthick= 2.0, msg='v3b2v3b1_f', x0= 0.70, y0= 0.88
;process_one_sfr, "priya/v3b2v3b1_k", lcolor=60, lthick= 2.0, msg='v3b2v3b1_k', x0= 0.70, y0= 0.88

;process_one_sfr, "priya/v3b2v3b2_e", lcolor=80, lthick= 2.0, msg='v3b2v3b2_e', x0= 0.70, y0= 0.85
;process_one_sfr, "priya/v3b2v3b2_f", lcolor=80, lthick= 2.0, msg='v3b2v3b2_f', x0= 0.70, y0= 0.85
;process_one_sfr, "priya/v3b2v3b2_k", lcolor=80, lthick= 2.0, msg='v3b2v3b2_k', x0= 0.70, y0= 0.85

;process_one_sfr, "priya/v3b2v3b4_e", lcolor=100, lthick= 2.0, msg='v3b2v3b4_e', x0= 0.70, y0= 0.82
;process_one_sfr, "priya/v3b2v3b4_f", lcolor=100, lthick= 2.0, msg='v3b2v3b4_f', x0= 0.70, y0= 0.82
;process_one_sfr, "priya/v3b2v3b4_k", lcolor=100, lthick= 2.0, msg='v3b2v3b4_k', x0= 0.70, y0= 0.82

;process_one_sfr, "priya/v3b3v3b1_e", lcolor=120, lthick= 2.0, msg='v3b3v3b1_e', x0= 0.70, y0= 0.79
;process_one_sfr, "priya/v3b3v3b1_f", lcolor=120, lthick= 2.0, msg='v3b3v3b1_f', x0= 0.70, y0= 0.79
;process_one_sfr, "priya/v3b3v3b1_k", lcolor=120, lthick= 2.0, msg='v3b3v3b1_k', x0= 0.70, y0= 0.79

;process_one_sfr, "priya/v3b3v3b2_e", lcolor=140, lthick= 2.0, msg='v3b3v3b2_e', x0= 0.70, y0= 0.76
;process_one_sfr, "priya/v3b3v3b2_f", lcolor=140, lthick= 2.0, msg='v3b3v3b2_f', x0= 0.70, y0= 0.76
;process_one_sfr, "priya/v3b3v3b2_k", lcolor=140, lthick= 2.0, msg='v3b3v3b2_k', x0= 0.70, y0= 0.76

;process_one_sfr, "priya/v3b3v3b3_e", lcolor=160, lthick= 2.0, msg='v3b3v3b3_e', x0= 0.70, y0= 0.73
;process_one_sfr, "priya/v3b3v3b3_f", lcolor=160, lthick= 2.0, msg='v3b3v3b3_f', x0= 0.70, y0= 0.73
;process_one_sfr, "priya/v3b3v3b3_k", lcolor=160, lthick= 2.0, msg='v3b3v3b3_k', x0= 0.70, y0= 0.73

;process_one_sfr, "priya/v3b3v3b4_e", lcolor=180, lthick= 2.0, msg='v3b3v3b4_e', x0= 0.70, y0= 0.70
;process_one_sfr, "priya/v3b3v3b4_f", lcolor=180, lthick= 2.0, msg='v3b3v3b4_f', x0= 0.70, y0= 0.70
;process_one_sfr, "priya/v3b3v3b4_k", lcolor=180, lthick= 2.0, msg='v3b3v3b4_k', x0= 0.70, y0= 0.70


; -------------------------------------

;process_one_sfr, "tst/iso_1", lcolor=0, lthick= 2.0, msg='iso_1', x0= 0.58, y0= 0.93
;process_one_sfr, 'tst/iso_2', lcolor=50, lthick= 2.0, msg='iso_2 (no BH)', x0= 0.58, y0= 0.90
;process_one_sfr, "tst/iso_1", lcolor=0, lthick= 2.0, msg='Arepo', x0= 0.58, y0= 0.93
;process_one_sfr, 'tst/iso_2', lcolor=50, lthick= 2.0, msg='Arepo (no BH)', x0= 0.58, y0= 0.90

;process_one_sfr, "priya/iso_v3b1",   lcolor= 150, lthick= 2.0, msg='iso_v3b1', x0= 0.58, y0= 0.85
;process_one_sfr, 'priya/iso_v3b1_2', lcolor=100, lthick= 2.0, msg='iso_v3b1_2 (no BH)', x0= 0.58, y0= 0.82
;process_one_sfr, "priya/iso_v3b1",   lcolor= 150, lthick= 2.0, msg='Gadget', x0= 0.58, y0= 0.85
;process_one_sfr, 'priya/iso_v3b1_2', lcolor=100, lthick= 2.0, msg='Gadget (no BH)', x0= 0.58, y0= 0.82
;xyouts, 0.3, 0.3, "Isolated Disk", color= 0, charthick= 3.0, size= 1.2
;process_one_sfr, "priya/iso_v3b2",   lcolor= 150, lthick= 2.0, msg='iso_v3b2', x0= 0.70, y0= 0.79
;process_one_sfr, "priya/iso_v3b3",   lcolor= 200, lthick= 2.0, msg='iso_v3b3', x0= 0.70, y0= 0.76
;process_one_sfr, "priya/iso_v3b4",   lcolor= 180, lthick= 2.0, msg='iso_v3b4', x0= 0.70, y0= 0.73

;process_one_sfr, "priya/v3b1v3b1_e", lcolor=0, lthick= 2.0, msg='v3b1v3b1_e', x0= 0.70, y0= 0.94
;process_one_sfr, "priya/v3b1v3b1_e1", lcolor=0, lthick= 2.0, msg='v3b1v3b1_e1', x0= 0.70, y0= 0.91
;process_one_sfr, "priya/v3b1v3b1_e2", lcolor=40, lthick= 2.0, msg='v3b1v3b1_e2', x0= 0.70, y0= 0.88
;process_one_sfr, "tst/tst1", lcolor=100, lthick= 2.0, msg='tst1', x0= 0.70, y0= 0.85
;process_one_sfr, "tst/tst2", lcolor=50,  lthick= 2.0, msg='tst2', x0= 0.70, y0= 0.82
;process_one_sfr, "tst/tst3", lcolor=200, lthick= 2.0, msg='tst3', x0= 0.70, y0= 0.79
;process_one_sfr, "tst/tst3_hdf", lcolor=100, lthick= 2.0, msg='tst3_hdf', x0= 0.70, y0= 0.76
;process_one_sfr, "tst/tst3_2", lcolor=150, lthick= 2.0, msg='tst3_2', x0= 0.70, y0= 0.76
;process_one_sfr, "tst/tst3_3", lcolor=200, lthick= 2.0, msg='tst3_3', x0= 0.70, y0= 0.73
;process_one_sfr, "tst/tst3_4", lcolor=100, lthick= 2.0, msg='tst3_4', x0= 0.70, y0= 0.70
;process_one_sfr, "tst/tst3_5", lcolor=50, lthick= 2.0, msg='tst3_5', x0= 0.70, y0= 0.67
;process_one_sfr, "tst/tst3_6", lcolor=10, lthick= 2.0, msg='tst3_6', x0= 0.70, y0= 0.64
;process_one_sfr, "tst/tst4", lcolor=150, lthick= 2.0, msg='tst4', x0= 0.70, y0= 0.73
;process_one_sfr, "tst/tst5", lcolor=50, lthick= 2.0, msg='tst5', x0= 0.70, y0= 0.70

;process_one_sfr, "priya/v3b1v3b1_e1", lcolor=50, lthick= 2.0, msg='Gadget', x0= 0.30, y0= 0.40
;process_one_sfr, "tst/tst3", lcolor=150, lthick= 2.0, msg='Arepo', x0= 0.30, y0= 0.35


; -------------------------------------


; -------------------------------------

; gurtina's lmc

;process_one_sfr, "gurtina/lmcq0", lcolor=50, lthick= 2.0, msg='q=0', x0= 0.70, y0= 0.88
;process_one_sfr, "gurtina/lmcq005", lcolor=200, lthick= 2.0, msg='q=0.05', x0= 0.70, y0= 0.84
;process_one_sfr, "gurtina/lmcq01", lcolor=100, lthick= 2.0, msg='q=0.1', x0= 0.70, y0= 0.80
;process_one_sfr, "1yygurtina/lmcq03", lcolor=150, lthick= 2.0, msg='q=0.3', x0= 0.70, y0= 0.76




; -------------------------------------

; keck timing tests

;do_ody_timing= 1
do_ody_timing= 0
if do_ody_timing eq 1 then begin
	process_one_sfr, "tests/SbSbhs_e_8", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_8e", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_8h", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88

	process_one_sfr, "tests/SbSbhs_e_16", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_16a", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_16b", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_16c", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_16d", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_16e", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_16f", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_16g", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_16h", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_16i", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	;process_one_sfr, "tests/SbSbhs_e_16j", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88   ; didn't finish
	process_one_sfr, "tests/SbSbhs_e_16k", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	;process_one_sfr, "tests/SbSbhs_e_16l", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88   ; didn't finish

	process_one_sfr, "tests/SbSbhs_e_32", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_32e", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	;process_one_sfr, "tests/SbSbhs_e_32g", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88    ; didn't finish
	process_one_sfr, "tests/SbSbhs_e_32h", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_32h2", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88

	process_one_sfr, "tests/SbSbhs_e_48", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_48", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_48e", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_48h", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88

	process_one_sfr, "tests/SbSbhs_e_64", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_64", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_64e", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88
	process_one_sfr, "tests/SbSbhs_e_64h", lcolor=100, lthick= 1.0, msg='1x', x0= 0.25, y0= 0.88

	process_one_sfr, "tests/Sb2x1Sb2x1_e_16", lcolor=50, lthick= 1.0, msg='2x', x0= 0.25, y0= 0.84
	process_one_sfr, "tests/Sb2x1Sb2x1_e_32", lcolor=50, lthick= 1.0, msg='2x', x0= 0.25, y0= 0.84

	process_one_sfr, "tests/Sb3x1Sb3x1_e_16", lcolor=150, lthick= 1.0, msg='3x', x0= 0.25, y0= 0.80
	process_one_sfr, "tests/Sb3x1Sb3x1_e_32", lcolor=150, lthick= 1.0, msg='3x', x0= 0.25, y0= 0.80

	process_one_sfr, "tests/Sb5x1Sb5x1_e_32", lcolor=200, lthick= 1.0, msg='5x', x0= 0.25, y0= 0.76
	process_one_sfr, "tests/Sb5x1Sb5x1_e_48", lcolor=200, lthick= 1.0, msg='5x', x0= 0.25, y0= 0.76
	process_one_sfr, "tests/Sb5x1Sb5x1_e_64", lcolor=200, lthick= 1.0, msg='5x', x0= 0.25, y0= 0.76
	process_one_sfr, "tests/Sb5x1Sb5x1_e_96", lcolor=200, lthick= 1.0, msg='5x', x0= 0.25, y0= 0.76

	process_one_sfr, "tests/Sb10xSb10xhs_e_8", lcolor=0, lthick= 1.0, msg='10x', x0= 0.25, y0= 0.72
	process_one_sfr, "tests/Sb10xSb10xhs_e_16", lcolor=0, lthick= 1.0, msg='10x', x0= 0.25, y0= 0.72
	process_one_sfr, "tests/Sb10xSb10xhs_e_32", lcolor=0, lthick= 1.0, msg='10x', x0= 0.25, y0= 0.72
	process_one_sfr, "tests/Sb10xSb10xhs_e_64", lcolor=0, lthick= 1.0, msg='10x', x0= 0.25, y0= 0.72
	process_one_sfr, "tests/Sb10xSb10xhs_e_128", lcolor=0, lthick= 1.0, msg='10x', x0= 0.25, y0= 0.72
	process_one_sfr, "tests/Sb10xSb10xhs_e_256", lcolor=0, lthick= 1.0, msg='10x', x0= 0.25, y0= 0.72
endif




; -------------------------------------

; Sbc models

;process_one_sfr, "isolated/Sbc", lcolor=100, lthick= 2.0, msg='iso/Sbc', x0= 0.70, y0= 0.88
;process_one_sfr, "Sbc", lcolor=50, lthick= 2.0, msg='Sbc', x0= 0.70, y0= 0.88
;process_one_sfr, "Sbc201", lcolor=150, lthick= 2.0, msg='Sbc201', x0= 0.70, y0= 0.84

;xmax = 2.5 ,  xmin = 0,  ymax = 1000,  ymin = 0.8

;process_one_sfr, "Sbc/Sbc201_10x_wBH", lcolor=50, lthick= 2.0, msg='Sbc201_10x_wBH', x0= 0.22, y0= 0.92
;process_one_sfr, "Sbc/Sbc201_10x_woBH", lcolor=150, lthick= 2.0, msg='Sbc201_10x_woBH', x0= 0.22, y0= 0.88
;process_one_sfr, "Sbc/Sbc201_q1", lcolor=100, lthick= 2.0, msg='Sbc201_q1', x0= 0.22, y0= 0.84
;process_one_sfr, "Sbc/Sbc201_q1_woBH", lcolor=20, lthick= 2.0, msg='Sbc201_q_woBH1', x0= 0.22, y0= 0.80
;process_one_sfr, "Sbc/Sbccut201", lcolor=70, lthick= 2.0, msg='Sbccut201', x0= 0.22, y0= 0.76
;process_one_sfr, "Sbc/Sbccut201_woBH", lcolor=200, lthick= 2.0, msg='Sbccut201_woBH', x0= 0.22, y0= 0.72

;process_one_sfr, "Sbc/Sbc_10x_wBH", lcolor=100, lthick= 2.0, msg='Sbc_10x_wBH', x0= 0.72, y0= 0.40
;process_one_sfr, "Sbc/Sbc_10x_wBHq", lcolor=200, lthick= 2.0, msg='Sbc_10x_wBHq', x0= 0.72, y0= 0.36
;process_one_sfr, "Sbc/Sbc_q1", lcolor=180, lthick= 2.0, msg='Sbc_q1', x0= 0.72, y0= 0.32
;process_one_sfr, "Sbc/Sbccut", lcolor=140, lthick= 2.0, msg='Sbccut', x0= 0.72, y0= 0.28

;process_one_sfr, "Sbc/Sbcm1201", lcolor=30, lthick= 2.0, msg='Sbcm1201', x0= 0.22, y0= 0.92
;process_one_sfr, "Sbc/Sbcm1201_woBH", lcolor=60, lthick= 2.0, msg='Sbcm1201_woBH', x0= 0.22, y0= 0.88
;process_one_sfr, "Sbc/Sbcm2201", lcolor=90, lthick= 2.0, msg='Sbcm2201', x0= 0.22, y0= 0.84
;process_one_sfr, "Sbc/Sbcm2201_woBH", lcolor=120, lthick= 2.0, msg='Sbcm2201_woBH', x0= 0.22, y0= 0.80
;process_one_sfr, "Sbc/Sbcm3201", lcolor=150, lthick= 2.0, msg='Sbcm3201', x0= 0.22, y0= 0.76
;process_one_sfr, "Sbc/Sbcm3201_woBH", lcolor=180, lthick= 2.0, msg='Sbcm3201_woBH', x0= 0.22, y0= 0.72

;process_one_sfr, "Sbc/Sbcm1_q1", lcolor=180, lthick= 2.0, msg='Sbcm1_q1', x0= 0.72, y0= 0.40
;process_one_sfr, "Sbc/Sbcm2_q1", lcolor=130, lthick= 2.0, msg='Sbcm2_q1', x0= 0.72, y0= 0.36
;process_one_sfr, "Sbc/Sbcm3_q1", lcolor=100, lthick= 2.0, msg='Sbcm3_q1', x0= 0.72, y0= 0.32
;process_one_sfr, "Sbc/Sbcm1", lcolor=20, lthick= 2.0, msg='Sbcm1', x0= 0.72, y0= 0.28
;process_one_sfr, "Sbc/Sbcm2", lcolor=50, lthick= 2.0, msg='Sbcm2', x0= 0.72, y0= 0.24
;process_one_sfr, "Sbc/Sbcm3", lcolor=80, lthick= 2.0, msg='Sbcm3', x0= 0.72, y0= 0.20



; extended gas disk versus not

;process_one_sfr, "Sbc/Sbccut201", lcolor=150, lthick= 2.0, msg=' ', x0= 0.22, y0= 0.76
;process_one_sfr, "Sbc/Sbcm2201", lcolor=150, lthick= 2.0, msg=' ', x0= 0.22, y0= 0.84

;process_one_sfr, "Sbc/Sbc201_q1", lcolor=0, lthick= 2.0, msg=' ', x0= 0.22, y0= 0.84
;process_one_sfr, "Sbc/Sbcm1201", lcolor=0, lthick= 2.0, msg=' ', x0= 0.22, y0= 0.92
;process_one_sfr, "Sbc/Sbcm3201", lcolor=0, lthick= 2.0, msg=' ', x0= 0.22, y0= 0.76


; -------------------------------------

; Alternate SFR laws   rho_SFR= rho_gas ^ N

;process_one_sfr, "isolated/d3_N0", lcolor=30, lthick= 2.0, msg='d3_N0', x0= 0.75, y0= 0.86
;process_one_sfr, "isolated/d3_N0a", lcolor=60, lthick= 2.0, msg='d3_N0a', x0= 0.75, y0= 0.82
;process_one_sfr, "isolated/d3_N0b", lcolor=90, lthick= 2.0, msg='d3_N0b', x0= 0.75, y0= 0.78
;process_one_sfr, "isolated/d3", lcolor=0, lthick= 2.0, msg='d3', x0= 0.75, y0= 0.90
;process_one_sfr, "isolated/d3_N1.5", lcolor=120, lthick= 2.0, msg='d3_N1.5', x0= 0.75, y0= 0.74
;process_one_sfr, "isolated/d3_N2", lcolor=150, lthick= 2.0, msg='d3_N2', x0= 0.75, y0= 0.70
;process_one_sfr, "isolated/d3_N2a", lcolor=180, lthick= 2.0, msg='d3_N2a', x0= 0.75, y0= 0.66
;process_one_sfr, "isolated/d3_N2b", lcolor=210, lthick= 2.0, msg='d3_N2b', x0= 0.75, y0= 0.62

;process_one_sfr, "isolated/d3_N0", lcolor=30, lthick= 2.0, msg='N=1.0 (d3_N0)', x0= 0.60, y0= 0.86
;process_one_sfr, "isolated/d3_N1.5", lcolor=120, lthick= 2.0, msg='N=1.5 (d3_N1.5)', x0= 0.60, y0= 0.81
;process_one_sfr, "isolated/d3_N2", lcolor=150, lthick= 2.0, msg='N=2.0 (d3_N2)', x0= 0.60, y0= 0.76

;process_one_sfr, "isolated/d3_N0c", lcolor=30, lthick= 2.0, msg='N=1.0 (d3_N0c)', x0= 0.60, y0= 0.86
;process_one_sfr, "isolated/d3_N1.5", lcolor=120, lthick= 2.0, msg='N=1.5 (d3_N1.5)', x0= 0.60, y0= 0.81
;process_one_sfr, "isolated/d3_N2c", lcolor=150, lthick= 2.0, msg='N=2.0 (d3_N2c)', x0= 0.60, y0= 0.76

;process_one_sfr, "ds/vc3vc3e_2", lcolor=0, lthick= 2.0, msg='vc3vc3e_2', x0= 0.60, y0= 0.86
;process_one_sfr, "ds/vc3vc3e", lcolor=50, lthick= 2.0, msg='vc3vc3e', x0= 0.60, y0= 0.82
;process_one_sfr, "altsf/vc3vc3e_N1", lcolor=120, lthick= 2.0, msg='N=1.0 (vc3vc3e_N1)', x0= 0.60, y0= 0.78
;process_one_sfr, "altsf/vc3vc3e_N2", lcolor=150, lthick= 2.0, msg='N=2.0 (vc3vc3e_N2)', x0= 0.60, y0= 0.74



; -------------------------------------

; Test Odyssey

;process_one_sfr, "ds/vc3vc3e_2", lcolor=0, lthick= 2.0, msg='vc3vc3e_2', x0= 0.60, y0= 0.90
;process_one_sfr, "ds/vc3vc3e", lcolor=200, lthick= 2.0, msg='vc3vc3e', x0= 0.60, y0= 0.86
;process_one_sfr, "odyssey/vc3vc3e", lcolor=50, lthick= 5.0, msg='odyssey/vc3vc3e', x0= 0.60, y0= 0.82
;process_one_sfr, "odyssey/vc3vc3e_1", lcolor=150, lstyle= 1, lthick= 1.0, msg='odyssey/vc3vc3e_1', x0= 0.60, y0= 0.78


; -------------------------------------

; Test Gadget

;process_one_sfr, "ds/vc3vc3e_2", lcolor=0, lthick= 2.0, msg='vc3vc3e_2', x0= 0.60, y0= 0.90
;process_one_sfr, "ds/vc3vc3e", lcolor=200, lthick= 2.0, msg='vc3vc3e', x0= 0.60, y0= 0.86
;process_one_sfr, "ds/vc3vc3e_3", lcolor=50, lthick= 2.0, msg='vc3vc3e_3', x0= 0.60, y0= 0.82
;process_one_sfr, "ds/vc3vc3e_4", lcolor=100, lthick= 2.0, msg='vc3vc3e_4', x0= 0.60, y0= 0.78
;process_one_sfr, "ds/vc3vc3e_5", lcolor=150, lthick= 2.0, msg='vc3vc3e_5', x0= 0.60, y0= 0.74



; -------------------------------------

; gaseous balls

; group A
;fruns= ["grpA_1","grpA_2","grpA_3","grpA_4","grpA_5", "grpA_no","grpA_no_1"]
;lbls= ["1e5","1e2","1e1","1e9","1e8", "grpA_no","grpA_no_1"]
;fruns= ["grpA_no_1","grpA_3","grpA_2","grpA_1","grpA_5","grpA_4"]
;lbls= ['none','1e1','1e2','1e5','1e8','1e9']
;msg= ' '



; -------------------------------------

;process_one_sfr, "ds/d3h1", lcolor=50, lthick= 3.0, msg=' ', show_specific_times= [0.9]

; -------------------------------------

; five to one mergers

;fruns= ["fivetoone/vc3V"]
;msg= ' '



; -------------------------------------

; local group

;fruns= ["localgroup/v2", "localgroup/v2_noigm"]
;fruns= ["localgroup/v3", "localgroup/v3_noigm"]
;lbls= ["!6The Local Group", "Isolated MW + M31"]
msg= ' '


; -------------------------------------

; Avishai's direct hit

;process_one_sfr, "crazy", lcolor=150, lthick= 2.0, msg='direct hit sim', x0= 0.20, y0= 0.88
msg= ' '





;--------------------------------------
;--------------------------------------

device, /close


end





;
;
;
;
;=====================================================================
;
;
;
;
;
;
;
;
;
;
;
;--------------------------------------------------------------------



;==============================================================================
;==============================================================================






pro sfr_4, junk, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		h=h


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_4, junk, filename=filename, /h"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=14.0, newysize=14.0
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
if keyword_set(cumulative) then yaxistitle="!6Mass (10!e10!n h!E-1!NM!D!9n!6!N)"
if keyword_set(gasmass) then yaxistitle="!6Gas Mass (10!e10!n h!E-1!NM!D!9n!6!N)"

xaxistitle = "!6Time (!8h!6!E-1!N Gyr)"
;xaxistitle = "!6Time (Gyr)"  & h= 0.7

;xmax = 14.0
;xmax = 12.0
;xmax = 10.0
;xmax = 7.5
;xmax = 6.0
;xmax = 5.0
;xmax = 4.25
;xmax = 4.0
xmax = 3.0
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
;ymax = 400
ymax = 300
;ymax = 250.0
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
;ymin= 0.01
ymin= 0.002
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


;-----------------------------------------
;
;   Load the panels
;
;-----------------------------------------
;
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    1     |    2     |
;     |          |          |
;     |          |          |
;     |---------------------|
;     |          |          |
;     |          |          |
;     |    3     |    4     |
;     |          |          |
;     |          |          |
;     |---------------------|
;
;

x0= 0.14
xs= 0.5*(0.98-x0)
x1= x0+xs
x2= x0+xs+xs

y0= 0.10
ys= 0.5*(0.98-y0)
y1= y0+ys
y2= y0+ys+ys



;  1
; ---
!p.position= [x0,y1,x1,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, $
	xtickformat='(a1)', ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
process_and_plot_sfr_history, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0
process_and_plot_sfr_history, "sbw/sb6", lcolor=220, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0
process_and_plot_sfr_history, "sbw/sb2", lcolor=180, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0
process_and_plot_sfr_history, "sbw/sb1", lcolor=140, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_and_plot_sfr_history, "sbw/sb3", lcolor=100, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_and_plot_sfr_history, "sbw/sb4", lcolor=60, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_and_plot_sfr_history, "sbw/sb5", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x1-0.15, y2-0.06, '!7g!6=0.005', /normal, charthick=3, size=1.33, color=0



;  2
; ---
!p.position= [x1,y1,x2,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

fload_newcolortable, 4
process_and_plot_sfr_history, "vc3vc3e_no", lcolor= 200, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_and_plot_sfr_history, "sbw/sb10", lcolor=200, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x2-0.15, y2-0.06, '!7g!6=0.5', /normal, charthick=3, size=1.33, color=0



;  3
; ---
!p.position= [x0,y0,x1,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
process_and_plot_sfr_history, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0
process_and_plot_sfr_history, "sbw/sb8", lcolor=200, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0
process_and_plot_sfr_history, "sbw/sb12", lcolor=140, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_and_plot_sfr_history, "sbw/sb13", lcolor=80, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_and_plot_sfr_history, "sbw/sb14", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x1-0.15, y1-0.06, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0




;  4
; ---
!p.position= [x1,y0,x2,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytickformat='(a1)'

fload_newcolortable, 4
process_and_plot_sfr_history, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0
process_and_plot_sfr_history, "sbw/sb11", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x2-0.15, y1-0.06, '!7g!6=5.0', /normal, charthick=3, size=1.33, color=0




;--------------------------------------
device, /close


end













;======================================================================








