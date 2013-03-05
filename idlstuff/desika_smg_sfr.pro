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
;xmax = 3.5
;xmax = 3.0
;xmax = 2.8
;xmax = 2.4
;xmax = 2.0
;xmax = 1.5
;xmax = 1.4
;xmax = 1.3
xmax = 1.0
;xmax = 0.5
xmin = 0

;ymax= 8e+5
;ymax= 8e+4
ymax= 2e+4
;ymax = 8000
;ymax = 4000
;ymax = 2500
;ymax = 1750
;ymax = 1500
;ymax = 1000
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
ymin= 1.0
;ymin = 0  & ynotlog= 1
;ymin = 0.8
;ymin = 0.4
;ymin = 0.1
;ymin= 0.07
;ymin= 0.01
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

;process_one_sfr, "vc3vc3b", lcolor=40,  lthick= 2.0, msg='b', x0= 0.85, y0= 0.78
;process_one_sfr, "vc3vc3c", lcolor=80,  lthick= 2.0, msg='c', x0= 0.85, y0= 0.75
;process_one_sfr, "vc3vc3d", lcolor=120, lthick= 2.0, msg='d', x0= 0.85, y0= 0.72
;process_one_sfr, "vc3vc3e", lcolor=160, lthick= 2.0, msg='e', x0= 0.85, y0= 0.69
;process_one_sfr, "vc3vc3f", lcolor=200, lthick= 2.0, msg='f', x0= 0.85, y0= 0.66

;process_one_sfr, "vc3vc3g", lcolor=40,  lthick= 2.0, msg='g', x0= 0.85, y0= 0.78
;process_one_sfr, "vc3vc3h", lcolor=80,  lthick= 2.0, msg='h', x0= 0.85, y0= 0.75
;process_one_sfr, "vc3vc3i", lcolor=120, lthick= 2.0, msg='i', x0= 0.85, y0= 0.72
;process_one_sfr, "vc3vc3j", lcolor=160, lthick= 2.0, msg='j', x0= 0.85, y0= 0.69
;process_one_sfr, "vc3vc3k", lcolor=200, lthick= 2.0, msg='k', x0= 0.85, y0= 0.66

;process_one_sfr, "vc3vc3l", lcolor=40,  lthick= 2.0, msg='l', x0= 0.85, y0= 0.78
;process_one_sfr, "vc3vc3m", lcolor=80,  lthick= 2.0, msg='m', x0= 0.85, y0= 0.75
;process_one_sfr, "vc3vc3n", lcolor=120, lthick= 2.0, msg='n', x0= 0.85, y0= 0.72
;process_one_sfr, "vc3vc3o", lcolor=160, lthick= 2.0, msg='o', x0= 0.85, y0= 0.69
;process_one_sfr, "vc3vc3p", lcolor=200, lthick= 2.0, msg='p', x0= 0.85, y0= 0.66

;--------------------------------------------------------------------------------
;process_one_sfr, "ds/vc3vc3e", lcolor=150, lthick= 2.0, msg='e', x0= 0.85, y0= 0.85
;process_one_sfr, "ds/vc3vc3h", lcolor=0, lthick= 2.0, msg='h', x0= 0.85, y0= 0.80
;--------------------------------------------------------------------------------

;process_one_sfr, "vc3vc3l", lcolor=50, lthick= 2.0, msg='l', x0= 0.25, y0= 0.32
;process_one_sfr, "vc3vc3l_bg", lcolor=150, lthick= 2.0, msg='l_bg', x0= 0.25, y0= 0.28

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

;process_one_sfr, "bs/b3e", lcolor=50, lthick= 2.0, msg='b3e', x0= 0.25, y0= 0.32
;process_one_sfr, "vc3vc3e_gf8", lcolor=150, lthick= 2.0, msg='vc3vc3e_gf8', x0= 0.25, y0= 0.29
;process_one_sfr, "vc3vc3e_gf8a", lcolor=100, lthick= 2.0, msg='vc3vc3e_gf8a', x0= 0.25, y0= 0.26

;msg=' '


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
;process_one_sfr, "bs/b2e", lcolor=200, lthick= 2.0, msg='e', x0= 0.25, y0= 0.32
;process_one_sfr, "bs/b2h", lcolor=150, lthick= 2.0, msg='h', x0= 0.25, y0= 0.28
;process_one_sfr, "bs/b2f", lcolor=100, lthick=2.0, msg='f', x0= 0.25, y0= 0.24
;process_one_sfr, "bs/b2k", lcolor=50, lthick=2.0, msg='k', x0= 0.25, y0= 0.20

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
;process_one_sfr, "As/A3", lcolor=150, lthick= 2.0, msg='z= 0.0', x0= 0.80, y0= 0.63
;process_one_sfr, "z3/a3", lcolor=50, lthick=2.0, msg='z= 3.0', x0= 0.80, y0= 0.67
;process_one_sfr, "z3/a4.5", lcolor=50, lthick=2.0, msg='z= 3.0', x0= 0.80, y0= 0.67
;process_one_sfr, "bs/b2e", lcolor=140, ctab= 1, lthick=2.0, msg='11.7', x0= 0.80, y0= 0.71



; z3 galaxies (brant v. mine)
; ------------------------------
;xyouts, 0.70, 0.92, 'Log(M!Dtot!N/M!D!9n!6!N)', /normal, charthick=3, size=1.33, color=0
;process_one_sfr, "brantsruns/full_model/z3/v3g2q1z3_v3g2q1z3_p", lcolor=0, lthick= 2.0, msg='v3g2q1z3', x0= 0.60, y0= 0.85
;process_one_sfr, "z3/b3h", lcolor=150, lthick= 2.0, msg='z3/b3h', x0= 0.60, y0= 0.80
;process_one_sfr, "z3/b3h2", lcolor=100, lthick= 2.0, msg='z3/b3h2', x0= 0.60, y0= 0.75
;process_one_sfr, "z3/b3h3", lcolor=200, lthick= 2.0, msg='z3/b3h3', x0= 0.60, y0= 0.70
;process_one_sfr, "z3/b3e", lcolor=50, lthick=2.0, msg='z3/b3e', x0= 0.60, y0= 0.87

;process_one_sfr, "z3/b3h2", lcolor=100, lthick= 2.0, msg='z3/b3h2 - raw', x0= 0.60, y0= 0.85, /raw
;process_one_sfr, "z3/b3h2", lcolor=150, lthick= 2.0, msg='z3/b3h2', x0= 0.60, y0= 0.75

; 80% gas
; set xmax= 202 above
;process_one_sfr, "brantsruns/full_model/z3/v3g2q1z3_v3g2q1z3_p", lcolor=0, lthick= 2.0, msg='v3g2q1z3', x0= 0.60, y0= 0.85
;process_one_sfr, "z3/b3h3", lcolor=200, lthick= 2.0, msg='z3/b3h3', x0= 0.60, y0= 0.80
;process_one_sfr, "z3/b3e2", lcolor=150, lthick= 2.0, msg='z3/b3e2', x0= 0.60, y0= 0.75

; set xmax= 600 above
;process_one_sfr, "brantsruns/full_model/z3/v4g2q1z3_v4g2q1z3_p", lcolor=0, lthick= 2.0, msg='v4g2q1z3', x0= 0.60, y0= 0.85
;process_one_sfr, "z3/b4h", lcolor=200, lthick= 2.0, msg='z3/b4h', x0= 0.60, y0= 0.80
;process_one_sfr, "z3/b4e", lcolor=150, lthick= 2.0, msg='z3/b4e', x0= 0.60, y0= 0.75

; set xmax= 2000 above
;process_one_sfr, "brantsruns/full_model/z3/v5g2q1z3_v5g2q1z3_p", lcolor=0, lthick= 2.0, msg='v5g2q1z3', x0= 0.60, y0= 0.85
;process_one_sfr, "z3/b5h", lcolor=200, lthick= 2.0, msg='z3/b5h', x0= 0.60, y0= 0.80
;process_one_sfr, "z3/b5e", lcolor=150, lthick= 2.0, msg='z3/b5e', x0= 0.60, y0= 0.75

;process_one_sfr, "brantsruns/full_model/z3/v6g2q1z3_v6g2q1z3_p", lcolor=0, lthick= 2.0, msg='v6g2q1z3', x0= 0.60, y0= 0.85
;process_one_sfr, "z3/b6h", lcolor=200, lthick= 2.0, msg='z3/b6h', x0= 0.60, y0= 0.80
;process_one_sfr, "z3/b6e", lcolor=150, lthick= 2.0, msg='z3/b6e', x0= 0.60, y0= 0.85
;process_one_sfr, "z3/b6f", lcolor=100, lthick= 2.0, msg='z3/b6f', x0= 0.60, y0= 0.80
;process_one_sfr, "z3/b6h", lcolor=50, lthick= 2.0, msg='z3/b6h', x0= 0.60, y0= 0.75
;process_one_sfr, "z3/b6k", lcolor=200, lthick= 2.0, msg='z3/b6k', x0= 0.60, y0= 0.70

;process_one_sfr, "z3/b6k", lcolor=200, lthick= 2.0, msg='z3/b6k', x0= 0.30, y0= 0.45
;process_one_sfr, "z3/b6k1", lcolor=175, lthick= 2.0, msg='z3/b6k1', x0= 0.30, y0= 0.42
;process_one_sfr, "z3/b6k2", lcolor=150, lthick= 2.0, msg='z3/b6k2', x0= 0.30, y0= 0.39
;process_one_sfr, "z3/b6k3", lcolor=125, lthick= 2.0, msg='z3/b6k3', x0= 0.30, y0= 0.36
;process_one_sfr, "z3/b6k4", lcolor=100, lthick= 2.0, msg='z3/b6k4', x0= 0.30, y0= 0.33
;process_one_sfr, "z3/b6k5", lcolor=75, lthick= 2.0, msg='z3/b6k5', x0= 0.30, y0= 0.30
;process_one_sfr, "z3/b6k6", lcolor=50, lthick= 2.0, msg='z3/b6k6', x0= 0.30, y0= 0.27
;process_one_sfr, "z3/b6k7", lcolor=25, lthick= 2.0, msg='z3/b6k7', x0= 0.30, y0= 0.24

;
; new z=3 SMG models
;
;process_one_sfr, "z3/b6e", lcolor=200, lthick= 2.0, msg='z3b6e', x0= 0.30, y0= 0.45
;process_one_sfr, "z3/b6v1e", lcolor=150, lthick= 2.0, msg='z3b6v1e', x0= 0.30, y0= 0.40
;process_one_sfr, "z3/b6v2e", lcolor=100, lthick= 2.0, msg='z3b6v2e', x0= 0.30, y0= 0.35
;process_one_sfr, "z3/b6v3e_a", lcolor=50, lthick= 2.0, msg='z3b6v3e', x0= 0.30, y0= 0.30

process_one_sfr, "z3/b5e", lcolor=0, lthick= 5.0, msg='z3/b5e', x0= 0.70, y0= 0.91
process_one_sfr, "z3/b5e_2", lcolor=200, lthick= 2.0, msg='z3/b5e_2', x0= 0.70, y0= 0.88
;process_one_sfr, "z3/b5e_3", lcolor=170, lthick= 2.0, msg='z3/b5e_3', x0= 0.70, y0= 0.83
;process_one_sfr, "z3/b5e_4", lcolor=140, lthick= 2.0, msg='z3/b5e_4', x0= 0.70, y0= 0.80
;process_one_sfr, "z3/b5e_5", lcolor=110, lthick= 2.0, msg='z3/b5e_5', x0= 0.70, y0= 0.77
;process_one_sfr, "z3/b5e_6", lcolor=80, lthick= 2.0, msg='z3/b5e_6', x0= 0.70, y0= 0.74

process_one_sfr, "z3/b5e_7", lcolor=50, lthick= 2.0, msg='z3/b5e_7', x0= 0.70, y0= 0.83
;process_one_sfr, "z3/b5e_7", lcolor=50, lthick= 2.0, msg='z3/b5e_7', x0= 0.70, y0= 0.74
;process_one_sfr, "z3/b5e_8", lcolor=150, lthick= 2.0, msg='z3/b5e_8', x0= 0.70, y0= 0.71

;process_one_sfr, "z3/b5e_9", lcolor=50, lthick= 2.0, msg='z3/b5e_9', x0= 0.70, y0= 0.74
;process_one_sfr, "z3/b5e_10", lcolor=150, lthick= 2.0, msg='z3/b5e_10', x0= 0.70, y0= 0.71

;process_one_sfr, "z3/b5e_11", lcolor=50, lthick= 2.0, msg='z3/b5e_11', x0= 0.70, y0= 0.74
;process_one_sfr, "z3/b5e_12", lcolor=150, lthick= 2.0, msg='z3/b5e_12', x0= 0.70, y0= 0.71

process_one_sfr, "z3/b5e_13", lcolor=100, lthick= 2.0, msg='z3/b5e_13', x0= 0.70, y0= 0.80
;process_one_sfr, "z3/b5e_13", lcolor=50, lthick= 2.0, msg='z3/b5e_13', x0= 0.70, y0= 0.74
;process_one_sfr, "z3/b5e_14", lcolor=150, lthick= 2.0, msg='z3/b5e_14', x0= 0.70, y0= 0.71

process_one_sfr, "z3/b5e_15", lcolor=130, lthick= 2.0, msg='z3/b5e_15', x0= 0.70, y0= 0.77
process_one_sfr, "z3/b5e_16", lcolor=160, lthick= 2.0, msg='z3/b5e_16', x0= 0.70, y0= 0.74
process_one_sfr, "z3/b5e_17", lcolor=240, lthick= 2.0, msg='z3/b5e_17', x0= 0.70, y0= 0.71

;--------------
;process_one_sfr, "z3/b5e", lcolor=0, lthick= 5.0, msg='z3/b5e', x0= 0.70, y0= 0.91, gasmass=gasmass
;process_one_sfr, "z3/b5e_2", lcolor=200, lthick= 2.0, msg='z3/b5e_2', x0= 0.70, y0= 0.88, gasmass=gasmass
;process_one_sfr, "z3/b5e_7", lcolor=50, lthick= 2.0, msg='z3/b5e_7', x0= 0.70, y0= 0.83, gasmass=gasmass
;process_one_sfr, "z3/b5e_13", lcolor=100, lthick= 2.0, msg='z3/b5e_13', x0= 0.70, y0= 0.80, gasmass=gasmass
;process_one_sfr, "z3/b5e_15", lcolor=130, lthick= 2.0, msg='z3/b5e_15', x0= 0.70, y0= 0.77, gasmass=gasmass
;process_one_sfr, "z3/b5e_16", lcolor=160, lthick= 2.0, msg='z3/b5e_16', x0= 0.70, y0= 0.74, gasmass=gasmass
;process_one_sfr, "z3/b5e_17", lcolor=240, lthick= 2.0, msg='z3/b5e_17', x0= 0.70, y0= 0.71, gasmass=gasmass
;--------------

;process_one_sfr, "z3/b5v1e", lcolor=150, lthick= 2.0, msg='z3/b5v1e', x0= 0.70, y0= 0.86
;process_one_sfr, "z3/b5v2e", lcolor=250, lthick= 2.0, msg='z3/b5v2e', x0= 0.70, y0= 0.83
;process_one_sfr, "z3/b5v3e", lcolor=220, lthick= 2.0, msg='z3/b5v3e', x0= 0.70, y0= 0.81
;process_one_sfr, "z3/b5v4e", lcolor=190, lthick= 2.0, msg='z3/b5v4e', x0= 0.70, y0= 0.78
;process_one_sfr, "z3/b5v4ea", lcolor=170, lthick= 2.0, msg='z3/b5v4ea', x0= 0.70, y0= 0.75
;process_one_sfr, "z3/b5v4ea", lcolor=120, lthick= 2.0, msg='z3/b5v4ea', x0= 0.70, y0= 0.75
;process_one_sfr, "z3/b5v9e", lcolor=150, lthick= 2.0, msg='z3/b5v9e', x0= 0.70, y0= 0.72
;process_one_sfr, "z3/b5v10e", lcolor=130, lthick= 2.0, msg='z3/b5v10e', x0= 0.70, y0= 0.69
;process_one_sfr, "z3/b5v11e", lcolor=110, lthick= 2.0, msg='z3/b5v11e', x0= 0.70, y0= 0.66
;process_one_sfr, "z3/b5v11ea", lcolor=90, lthick= 2.0, msg='z3/b5v11ea', x0= 0.70, y0= 0.63
;process_one_sfr, "z3/b5v11ea", lcolor=50, lthick= 2.0, msg='z3/b5v11ea', x0= 0.70, y0= 0.63
;process_one_sfr, "z3/b5v12e", lcolor=70, lthick= 2.0, msg='z3/b5v12e', x0= 0.70, y0= 0.60
;process_one_sfr, "z3/b5v13e", lcolor=50, lthick= 2.0, msg='z3/b5v13e', x0= 0.70, y0= 0.57
;process_one_sfr, "z3/b5v14e", lcolor=30, lthick= 2.0, msg='z3/b5v14e', x0= 0.70, y0= 0.54

;process_one_sfr, "brantsruns/full_model/z3/v6g2q1z3_v6g2q1z3_p", lcolor=50, lthick= 2.0, msg='v6g2q1z3', x0= 0.60, y0= 0.85, lstyle= 1
;process_one_sfr, "brantsruns/full_model/z3/v6g2q0z3_v6g2q0z3_p", lcolor=50, lthick= 2.0, msg='v6g2q0z3', x0= 0.60, y0= 0.81
;process_one_sfr, "brantsruns/full_model/z3/v6g1q1z3_v6g1q1z3_p", lcolor=150, lthick= 2.0, msg='v6g1q1z3', x0= 0.60, y0= 0.77, lstyle= 1
;process_one_sfr, "brantsruns/full_model/z3/v6g1q0z3_v6g1q0z3_p", lcolor=150, lthick= 2.0, msg='v6g1q0z3', x0= 0.60, y0= 0.73

;process_one_sfr, "brantsruns/full_model/z3/v5g2q1z3_v5g2q1z3_p", lcolor=50, lthick= 2.0, msg='v5g2q1z3', x0= 0.60, y0= 0.85, lstyle= 1
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

; -------------------------------------


;process_one_sfr, "minor/min_0", lcolor=50, lthick= 2.0, msg='min_0', x0= 0.70, y0= 0.92
;process_one_sfr, "minor/min_30", lcolor=80, lthick= 2.0, msg='min_30', x0= 0.70, y0= 0.88
;process_one_sfr, "minor/min_90", lcolor=110, lthick= 2.0, msg='min_90', x0= 0.70, y0= 0.84
;process_one_sfr, "minor/min_150", lcolor=140, lthick= 2.0, msg='min_150', x0= 0.70, y0= 0.80
;process_one_sfr, "minor/min_180", lcolor=170, lthick= 2.0, msg='min_180', x0= 0.70, y0= 0.76

;process_one_sfr, "minor/iso_Sb", lcolor=150, lthick= 2.0, msg='iso_Sb', x0= 0.70, y0= 0.92

;process_one_sfr, "minor/lr_min_0", lcolor=100, lthick= 2.0, msg='lr_min_0', x0= 0.70, y0= 0.88
;process_one_sfr, "minor/bd1_lr_min_0", lcolor=150, lthick= 2.0, msg='bd1_lr_min_0', x0= 0.70, y0= 0.84
;process_one_sfr, "minor/bd2_lr_min_0", lcolor=200, lthick= 2.0, msg='bd2_lr_min_0', x0= 0.70, y0= 0.80


;process_one_sfr, "minor/min_30", lcolor=150, lthick= 2.0, msg=' ', x0= 0.70, y0= 0.88, h=h



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



; -------------------------------------

; gurtina's lmc

;process_one_sfr, "gurtina/lmcq0", lcolor=50, lthick= 2.0, msg='q=0', x0= 0.70, y0= 0.88
;process_one_sfr, "gurtina/lmcq005", lcolor=200, lthick= 2.0, msg='q=0.05', x0= 0.70, y0= 0.84
;process_one_sfr, "gurtina/lmcq01", lcolor=100, lthick= 2.0, msg='q=0.1', x0= 0.70, y0= 0.80
;process_one_sfr, "gurtina/lmcq03", lcolor=150, lthick= 2.0, msg='q=0.3', x0= 0.70, y0= 0.76




; -------------------------------------

; Sbc models

;process_one_sfr, "isolated/Sbc", lcolor=100, lthick= 2.0, msg='iso/Sbc', x0= 0.70, y0= 0.88
;process_one_sfr, "Sbc", lcolor=50, lthick= 2.0, msg='Sbc', x0= 0.70, y0= 0.88
;process_one_sfr, "Sbc201", lcolor=150, lthick= 2.0, msg='Sbc201', x0= 0.70, y0= 0.84




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

; gaseous balls

; group A
;fruns= ["grpA_1","grpA_2","grpA_3","grpA_4","grpA_5", "grpA_no","grpA_no_1"]
;lbls= ["1e5","1e2","1e1","1e9","1e8", "grpA_no","grpA_no_1"]
;fruns= ["grpA_no_1","grpA_3","grpA_2","grpA_1","grpA_5","grpA_4"]
;lbls= ['none','1e1','1e2','1e5','1e8','1e9']
;msg= ' '




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







    


;===============================================================================
;
;
;      |----------------------|
;      |                      |
;      |                      |
;      |                      |
;      |      SFR             |
;      |                      |
;      |                      |
;      |                      |
;      |----------------------|
;
;
;
;
;======================================================================

pro sfr, frun, filename=filename, quiet=quiet
;pro junk

if not keyword_set(frun) then begin
;if not keyword_set(junk) then begin
   print, "  "
   print, "sfr, frun, filename=filename, quiet=quiet"
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename=frun+"/sfr.eps"

; -------------------------------------------------
   
   
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
;setup_plot_stuff, 'ps', filename=filename, colortable= 0


; ---------------------------------------------------
; ---------------------------------------------------
;
;  Print the Shit
;
; ---------------------------------------------------
; ---------------------------------------------------


yaxistitle="!6SFR (M!D!9n!6!N Yr!E-1!N)"
xaxistitle = "!6Time (Gyr/h)"


open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

xmax= max(sfrtime)
xmin= min(sfrtime)

ymax= max(sfrsfr) * 1.2
ymin= min(sfrsfr(where(sfrsfr gt 0.0))) * 0.2
;ymin= 0.002

; physical units
if keyword_set(h) then begin
        xaxistitle = "Time (Gyr)"
        h = fload_cosmology('h')

	sfrtime = sfrtime / h
        xmax = xmax / h
        xmin = xmin / h
endif


;---------------------------

!p.position= [0.18, 0.15, 0.98, 0.98]
;!p.font= 0

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        /ylog, ytickforma='exp_label', $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



;-------------------------------------------\

; plot sfr history

oplot, sfrtime, sfrsfr, psym=-3, color= 0, linestyle= 0, thick= 5.0



;-------------------------------------------

t_merg= fload_bh_mergertime(frun)
sfr_max_raw= max(sfrsfr)
t_sfr_max= sfrtime(where(sfrsfr eq sfr_max_raw))


;do_meansf= 0
do_meansf= 1
if do_meansf eq 1 then begin

	; initial period
	idx= where(sfrtime lt t_merg-0.2)
	meansfr_initial= mean(sfrsfr(idx))
	x= [min(sfrtime), t_merg-0.2]
	y= [meansfr_initial, meansfr_initial]
	oplot, x, y, psym=-3, linestyle=1, color= 150, thick=8.0

	; starburst
	idx= where((sfrtime ge t_merg-0.2) and (sfrtime le t_merg+0.2))
	meansfr_sb= mean(sfrsfr(idx))
	x= [t_merg-0.2, t_merg+0.2]
	y= [meansfr_sb, meansfr_sb]
	oplot, x, y, psym=-3, linestyle=1, color= 150, thick=8.0

	; final period
	idx= where(sfrtime gt t_merg+0.2)
	if idx(0) ne -1 then begin
		meansfr_final= mean(sfrsfr(idx))
		x= [t_merg+0.2,max(sfrtime)]
		y= [meansfr_final, meansfr_final]
		oplot, x, y, psym=-3, linestyle=1, color= 150, thick=8.0
	endif else begin
		meansfr_final= 1e-10
	endelse

endif


;do_sb_fit= 0
do_sb_fit= 1        ; first passage starburst
if do_sb_fit eq 1 then begin
        idx= where((sfrtime gt 0.1) and (sfrtime lt t_merg-0.2))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)

	sfr_max_first_peak= max(sfr)
	t_first_peak= time(where(sfr eq sfr_max_first_peak))

        fitandoverplot_gaussian, time, sfr, $
                        SFR_max= SFR_max_fp, width=width_fp, SBtime=SBtime_fp, quiet=quiet
endif


;do_sb_fit= 0
do_sb_fit= 1        ; final merger starburst
if do_sb_fit eq 1 then begin
        ; ----------------
        ; use full data
        idx= where((sfrtime gt t_merg-0.2) and (sfrtime lt t_merg+0.2))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)
        ; ----------------
        ; use smoothed data
        ;idx= where((timem gt t_merg-0.8) and (timem lt t_merg+0.8))
        ;time= timem(idx)
        ;sfr= sfrm(idx)
        ; ----------------

        fitandoverplot_gaussian, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime, quiet=quiet
endif

do_sb_fit= 0
;do_sb_fit= 1
if do_sb_fit eq 1 then begin
        idx= where((sfrtime gt t_merg-0.2) and (sfrtime lt t_merg+0.2))
        time= sfrtime(idx)
        sfr= sfrsfr(idx)
        fitandoverplot_gaussian_pluszeropt, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime
endif


;do_expdecay_fit= 0
do_expdecay_fit= 1
if do_expdecay_fit eq 1 then begin
        idx= where(sfrtime gt t_merg)
        ;idx= where((sfrtime gt SBtime) and (sfrtime lt SBtime+1.3))
        time= sfrtime(idx)
        time= time- time(0)
        sfr= sfrsfr(idx)
        fitandoverplot_expdecay, time, sfr, SBtime=t_merg, tau=tau, decay_const=decay_const
endif

do_fwhm_fit= 0
;do_fwhm_fit= 1
if do_fwhm_fit eq 1 then begin
        ; ----------------
        ; use full data
        idx= where((sfrtime gt t_merg-0.5) and (sfrtime lt t_merg+0.5))
        time= sfrtime
        sfr= sfrsfr
        ; ----------------
        ; use smoothed data
        ;time= timem
        ;sfr= sfrm
        ; ----------------
        findandoverplot_fwhm, time, sfr, $
                        ;SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, quiet=quiet
                        ;SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime, /quiet
                        SFR_max= SFR_max, fwhm=fwhm, SBtime=SBtime
endif



;--------------------------------------
;--------------------------------------



openw, 1, frun+'/sfrfits.txt', ERROR=err

printf, 1, "#   "
printf, 1, "#    Fits to the Star Formation History  "
printf, 1, "#   (all quantities are in Gadget units) "
printf, 1, "#   "
printf, 1, "#   "

printf, 1, "t_merg        ", t_merg
printf, 1, "sfr_max_raw   ", sfr_max_raw
printf, 1, "t_sfr_max     ", t_sfr_max
printf, 1, "sfr_max_fp    ", sfr_max_first_peak
printf, 1, "t_fp          ", t_first_peak
printf, 1, " "

printf, 1, "Mean SFRs (0 -> t_merg - 0.2 Myr/h -> t_merg + 0.2 Myr/h -> end)"
printf, 1, "meansfr_init  ", meansfr_initial
printf, 1, "meansfr_sb    ", meansfr_sb
printf, 1, "meansfr_final ", meansfr_final
printf, 1, " "

printf, 1, "Gaussian Fit to First Passage Starburst"
printf, 1, "SFR_max_fp    ", SFR_max_fp
printf, 1, "SBtime_fp     ", SBtime_fp
printf, 1, "width_fp      ", width_fp
printf, 1, " "

printf, 1, "Gaussian Fit to Starburst at t_merg"
printf, 1, "SFR_max       ", SFR_max
printf, 1, "SBtime        ", SBtime
printf, 1, "width         ", width
printf, 1, " "

printf, 1, "Exp. Fit to Post-SB SF"
printf, 1, "tau           ", tau
printf, 1, "decay_const   ", decay_const
printf, 1, " "

close, 1



;--------------------------------------
;--------------------------------------

device, /close


end















;===============================================================================


pro process_one_sfr, frun, lcolor=lcolor, lstyle=lstyle, lthick=lthick, $
				ctab=ctab, msg=msg, y0=y0, x0=x0, lpsym=lpsym, $
				cumulative=cumulative, $
				gasmass=gasmass, $
				gasfraction=gasfraction, $
				normalize=normalize, h=h, $
				multiby=multiby, $
				raw=raw


	if not keyword_set(lthick) then lthick= 1.0
	if not keyword_set(lstyle) then lstyle= 0
	if not keyword_set(lcolor) then lcolor= 0
	if not keyword_set(lpsym) then lpsym= -3


	if keyword_set(ctab) then begin
                        loadct, ctab
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0
	endif


	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	if keyword_set(raw) then begin
		open_raw_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac
	endif else begin
		open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac
	endelse


	; physical units
	if keyword_set(h) then sfrtime = sfrtime / h
	;if i eq 1 then sfrtime = sfrtime/(0.7)    ; did this for comparing h and noh

	if keyword_set(normalize) then sfrsfr=sfrsfr/sfrsfr(0)

	if keyword_set(multiby) then sfrsfr=sfrsfr * multiby


	    lthick= 4.0 * lthick

	if keyword_set(cumulative) then begin
		oplot, sfrtime, sfrmfs, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
		goto, moveon
	endif


	if keyword_set(gasmass) then begin
		oplot, sfrtime, sfrgasmass, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
		goto, moveon
	endif


	if keyword_set(gasfraction) then begin
		gasf= gasfraction * gasfrac
		oplot, sfrtime, gasf, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick
		goto, moveon
	endif


	; default option
	oplot, sfrtime, sfrsfr, psym=lpsym, color= lcolor, linestyle= lstyle, thick= lthick


moveon:


if not keyword_set(y0) then y0= 0.92
if not keyword_set(x0) then x0= 0.75

if keyword_set(msg) then begin
    xyouts, x0, y0, msg, /normal, charthick=3.0, size=1.33, color= lcolor
endif else begin
    xyouts, x0, y0, frun, /normal, charthick=3.0, size=1.33, color= lcolor
endelse



do_sb_fit= 0
;do_sb_fit= 1
if do_sb_fit eq 1 then begin
	t_merg= fload_bh_mergertime('/odyssey/tcox/'+frun)
	idx= where((sfrtime gt t_merg-0.3) and (sfrtime lt t_merg+0.3))
	time= sfrtime(idx)
	sfr= sfrsfr(idx)
	fitandoverplot_gaussian, time, sfr, $
			SFR_max= SFR_max, width=width, SBtime=SBtime
endif

do_expdecay_fit= 0
;do_expdecay_fit= 1
if do_expdecay_fit eq 1 then begin
	;idx= where(sfrtime gt SBtime)
	idx= where((sfrtime gt SBtime) and (sfrtime lt SBtime+1.3))
	time= sfrtime(idx)
	time= time- time(0)
	sfr= sfrsfr(idx)
	fitandoverplot_expdecay, time, sfr, SBtime=SBtime
endif


do_specific_point= 0
;do_specific_point= 1
if do_specific_point eq 1 then begin
	specific_point= 1.205
	idx=where(sfrtime ge specific_point)
	oplot, [sfrtime[idx[0]]], [sfrsfr[idx[0]]], psym= 4, color= 0, thick=5.0, symsize=1.6

	specific_point= 1.145
	idx=where(sfrtime ge specific_point)
	oplot, [sfrtime[idx[0]]], [sfrsfr[idx[0]]], psym= 4, color= 0, thick=5.0, symsize=1.6

	specific_point= 1.115
	idx=where(sfrtime ge specific_point)
	oplot, [sfrtime[idx[0]]], [sfrsfr[idx[0]]], psym= 4, color= 0, thick=5.0, symsize=1.6
endif

end







;======================================================================






pro open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	    ; get sfr data (from file)
	    ; -------------------------
	    loopnum= 0
	    repeat begin
	    ;sfrfile= '/data/tcox/sfr/'+fload_getid(fruns[i])+'.sfr'        ; twopiter
	    ;sfrfile= '/home/tcox/data/sfr/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid4/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid2/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/data6/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard

                ;if strmid(frun,1,4) eq 'raid' then begin
                if strmid(frun,1,7) eq 'odyssey' then begin
                        spawn, '/bin/ls '+frun+'/*.sfr ',result
                endif else begin
			case loopnum of
			  0: datadir='/n/home/tcox/data'
			  1: datadir='/n/circelfs/hernquist_lab'
			  2: datadir='/n/home/tcox'
			  3: datadir='/n/home'
			  4: datadir='/data6'
			  5: datadir='/data7'
			  else: break
			endcase

			;sfrfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'

			nfrun= fload_getid(frun)
			spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/*.sfr ',result
		endelse

		sfrfile=strcompress(result[0],/remove_all)



		get_lun, unit
		openr, unit, sfrfile, ERROR=err
		close, unit
		free_lun, unit

		if (err NE 0) then begin
		   ERR= 0
		   foundit= 0
		endif else begin
		   foundit= 1
		endelse

		loopnum= loopnum+1

	    endrep until (foundit eq 1)

	    if foundit eq 0 then begin
		msg= 'Problem: sfr file no found  file=[data]/tcox/'+fload_getid(frun)+'/'+fload_getid(frun)+'.sfr'
	        print, "  "
	        print, msg
	        print, "  "
		return
	    endif

	    ; read to open
	    print, "opening: ",sfrfile
	    sfrdata= read_ascii(sfrfile)
	    sfrtime= sfrdata.field1[0,*]
	    sfrsfr= sfrdata.field1[1,*]
	    sfrmfs= sfrdata.field1[3,*]    ; amount of new stars
	    gasfrac= sfrdata.field1[4,*]
	    n_cols= n_elements(sfrdata.field1[*,0])
	    n_rows= n_elements(sfrdata.field1[0,*])
	    finalnewstarmass= sfrmfs[n_rows-1]
	    ; gas mass
	    if n_cols ge 6 then sfrgasmass= sfrdata.field1[6,*] else sfrgasmass=[20.0,20.0-finalnewstarmass]
	    sfrmsg= ''


            t1per= 0
            idx=where(sfrtime ge 1.0)
            if idx(0) ne -1 and n_cols gt 5 then begin
                gm= sfrgasmass[0]-sfrgasmass(idx(0))
                t1per= 100.0*gm/sfrgasmass[0]
            endif

            t4per= 0
            idx=where(sfrtime ge 4.0)
            if idx(0) ne -1 and n_cols gt 5 then begin
                gm= sfrgasmass[0]-sfrgasmass(idx(0))
                t4per= 100.0*gm/sfrgasmass[0]
            endif

            maxsfr= max(sfrsfr)
            idx= where(sfrsfr eq maxsfr)
            maxsfrtime= sfrtime(idx)

            sfraftertmerg= 0.0
            smaftertmerg= 0.0
            taftertmerg= 0.0
            tmerger= 0.0
            if n_elements(tmerg) gt 1 then begin
                idx= where(sfrtime ge tmerg[i]+0.2)
                if idx(0) eq -1 then begin
                        sfraftertmerg= sfrsfr[n_rows-1]
                        smaftertmerg= sfrmfs[n_rows-1]
                        tmerger= sfrtime[n_rows-1]
                        taftertmerg= tmerger
                endif else begin
                        sfraftertmerg= sfrsfr[idx(0)]
                        smaftertmerg= sfrmfs[idx(0)]
                        tmerger= tmerg[i]
                        taftertmerg= sfrtime[idx(0)]
                endelse
            endif

            print, "-------------------------------------"
            print, "maximum sfr= ", maxsfr, " occurs at ", maxsfrtime
            print, "average sfr= ", total(sfrsfr)/n_rows
            print, "  "
            print, "      tmerg= ",tmerger,'  and 200 Myr later= ',taftertmerg
            print, "        sfr= ", sfraftertmerg,'  200 Myr after t_merger'
            print, "         sm= ", smaftertmerg,'  200 Myr after t_merger'
            n_mfs= n_elements(sfrgasmass)
            print, "-------------------------------------"
            print, "original gas mass   =", sfrgasmass[0]
            print, "remnant gas mass    =", sfrgasmass[n_mfs-1]
            print, "gas consumed        =", sfrgasmass[0]-sfrgasmass[n_mfs-1]
            print,"                      ", 100.0*(sfrgasmass[0]-sfrgasmass[n_mfs-1])/sfrgasmass[0],'  %'
            print,"                      ", t1per,' % after 1 Gyr'
            print,"                      ", t4per,' % after 4 Gyr'
            print, "-------------------------------------"


end





;======================================================================






pro open_raw_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac

	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	    ; get sfr data (from file)
	    ; -------------------------
	    loopnum= 0
	    repeat begin
	    ;sfrfile= '/data/tcox/sfr/'+fload_getid(fruns[i])+'.sfr'        ; twopiter
	    ;sfrfile= '/home/tcox/data/sfr/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid4/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/raid2/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard
	    ;sfrfile= '/data6/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'    ; harvard

                if strmid(frun,1,4) eq 'raid' then begin
                        spawn, '/bin/ls '+frun+'/sfr.txt ',result
                endif else begin
			case loopnum of
			  0: datadir='/raid4'
			  1: datadir='/home'
			  2: datadir='/raid12'
			  3: datadir='/raid2'
			  4: datadir='/data6'
			  5: datadir='/data7'
			  else: break
			endcase

			;sfrfile= datadir+'/tcox/'+fload_getid(fruns[i])+'/'+fload_getid(fruns[i])+'.sfr'

			nfrun= fload_getid(frun)
			spawn, '/bin/ls '+datadir+'/tcox/'+nfrun+'/sfr.txt ',result
		endelse

		sfrfile=strcompress(result[0],/remove_all)



		get_lun, unit
		openr, unit, sfrfile, ERROR=err
		close, unit
		free_lun, unit

		if (err NE 0) then begin
		   ERR= 0
		   foundit= 0
		endif else begin
		   foundit= 1
		endelse

		loopnum= loopnum+1

	    endrep until (foundit eq 1)

	    if foundit eq 0 then begin
		msg= 'Problem: sfr file no found  file=[data]/tcox/'+fload_getid(frun)+'/sfr.txt'
	        print, "  "
	        print, msg
	        print, "  "
		return
	    endif

	    ; read to open
	    print, "opening: ",sfrfile
	    sfrdata= read_ascii(sfrfile)
	    sfrtime= sfrdata.field1[0,*]
	    sfrmfs= sfrdata.field1[1,*]
	    sfrsfr_gadu= sfrdata.field1[2,*]    ; amount of new stars
	    sfrsfr= sfrdata.field1[3,*]
	    sfrmfs_spawned= sfrdata.field1[4,*]
	    n_cols= n_elements(sfrdata.field1[*,0])
	    n_rows= n_elements(sfrdata.field1[0,*])

            maxsfr= max(sfrsfr)
            idx= where(sfrsfr eq maxsfr)
            maxsfrtime= sfrtime(idx)

            print, "-------------------------------------"
            print, "maximum sfr= ", maxsfr, " occurs at ", maxsfrtime
            print, "average sfr= ", total(sfrsfr)/n_rows
            print, "-------------------------------------"

	    N_sample= 5
            print, "-------------------------------------"
            print, "-------------------------------------"
	    print, "   RESAMPLING N_sample= ", N_sample
	    sfrsfr= transpose(sfrsfr)
	    sfrtime= transpose(sfrtime)
	    n= n_elements(sfrsfr)
	    n_sub= n/N_sample
	    n_extra= n-N_sample*n_sub
	    for i=0,(N_sample-n_extra-1) do begin
		sfrsfr= [sfrsfr, sfrsfr[n-1]]
		sfrtime= [sfrtime, sfrtime[n-1]]
	    endfor
	    sfrsfr= rebin(sfrsfr, n_sub+1)
	    sfrtime= rebin(sfrtime, n_sub+1)
            print, "-------------------------------------"
            print, "-------------------------------------"

end





;--------------------------------------------------------------------



pro read_sfrfits_file, frun, $
		t_merg, sfr_max, t_sfr_max, sfr_max_fp, t_sfr_max_fp, $
		meansfr_initial, meansfr_sb, meansfr_final, $
		gSFR_max_fp, gSBtime_fp, gwidth_fp, $
		gSFR_max, gSBtime, gwidth, $
		tau, decay_const


;sfrfitsfile= '/raid4/tcox/'+frun+'/sfrfits.txt'
sfrfitsfile= frun+'/sfrfits.txt'

openr, 1, sfrfitsfile
junk=''
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_merg= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & sfr_max= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_sfr_max= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & sfr_max_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & t_sfr_max_fp= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & meansfr_initital= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & meansfr_sb= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & meansfr_final= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSFR_max_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSBtime_fp= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gwidth_fp= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSFR_max= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gSBtime= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & gwidth= float(tempjunk(1))
readf, 1, junk

readf, 1, junk
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & tau= float(tempjunk(1))
readf, 1, junk  & tempjunk=strsplit(junk,/extract,count=count) & decay_const= float(tempjunk(1))
readf, 1, junk

close, 1


end










;==============================================================================
;==============================================================================


pro sfr_test, junk, $
		filename=filename, $
		cumulative=cumulative, $
		gasmass=gasmass, $
		ynotlog=ynotlog, $
		h=h


if not keyword_set(junk) then begin
   print, "  "
   print, "sfr_test, junk, filename=filename, /h, "
   print, "           /cumulative, /gasmass, /ynotlog"
   print, "  "
   print, "  "
   print, "  WARNING: cumulative and gasmass do not currently work!  "
   print, "  "
   print, "  default filename: sfrtest.eps"
   print, "  "
   print, "  "
   print, "  "
   return
endif

if not keyword_set(filename) then filename='sfrtest.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4
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
;xmax = 3.5
xmax = 3.0
;xmax = 2.8
;xmax = 2.4
;xmax = 2.0
;xmax = 1.5
;xmax = 1.3
;xmax = 1.0
;xmax = 0.5
xmin = 0

;ymax= 8e+5
;ymax= 8e+4
;ymax = 8000
;ymax = 2500
;ymax = 1750
;ymax = 1500
;ymax = 1000
;ymax = 600
;ymax = 400
;ymax = 300
;ymax = 250.0
;ymax = 180
;ymax = 120
;ymax = 100
;ymax = 90
;ymax = 80
;ymax = 60
;ymax = 50.0
;ymax = 40.0
ymax = 20.0
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
ymin = 0  & ynotlog= 1
;ymin = 0.1
;ymin= 0.07
;ymin= 0.01
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



frun_orig= "ssft/iso2a"



;-------------------------------------------
;   Load the sfr.txt file
;-------------------------------------------
frun= frun_orig
process_one_sfr, frun, lcolor=50, lthick= 2.0, msg='from sfr.txt file', x0= 0.60, y0= 0.83, /raw




;-------------------------------------------
;   Load the *.sfr file
;-------------------------------------------
frun= frun_orig
process_one_sfr, frun, lcolor=150, lthick= 2.0, msg='from parse file', x0= 0.60, y0= 0.88






;--------------------------------------
;--------------------------------------

device, /close


end

    









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
process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0
process_one_sfr, "sbw/sb6", lcolor=220, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb2", lcolor=180, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb1", lcolor=140, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb3", lcolor=100, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb4", lcolor=60, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb5", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x1-0.15, y2-0.06, '!7g!6=0.005', /normal, charthick=3, size=1.33, color=0



;  2
; ---
!p.position= [x1,y1,x2,y2]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtickformat='(a1)', ytickformat='(a1)'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_no", lcolor= 200, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb10", lcolor=200, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x2-0.15, y2-0.06, '!7g!6=0.5', /normal, charthick=3, size=1.33, color=0



;  3
; ---
!p.position= [x0,y0,x1,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytitle=yaxistitle, ytickformat='exp_label'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0
process_one_sfr, "sbw/sb8", lcolor=200, ctab= 1, lthick= 2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb12", lcolor=140, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb13", lcolor=80, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0
process_one_sfr, "sbw/sb14", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x1-0.15, y1-0.06, '!7g!6=2.0', /normal, charthick=3, size=1.33, color=0




;  4
; ---
!p.position= [x1,y0,x2,y1]

plot, [1.0],[1.0], psym=-3, xrange=[xmin,xmax], yrange=[ymin,ymax], color= 0, $
        xstyle=1, ystyle=1, xcharsize=1.50, ycharsize=1.50, /nodata, /ylog, $
        xthick=4.0, ythick=4.0, charthick=3.0, /noerase, $
        xtitle=xaxistitle, ytickformat='(a1)'

fload_newcolortable, 4
process_one_sfr, "vc3vc3e_no", lcolor=200, lthick=2.0, msg=' ', x0= 0.0, y0=0.0
process_one_sfr, "sbw/sb11", lcolor=20, ctab= 1, lthick=2.0, msg=' ', x0= 0.0, y0= 0.0

xyouts, x2-0.15, y1-0.06, '!7g!6=5.0', /normal, charthick=3, size=1.33, color=0




;--------------------------------------
device, /close


end













;======================================================================








pro sfr_at_snaptimes, frun, h=h


if not keyword_set(frun) then begin
   print, "  "
   print, "sfr_at_snaptimes, frun, /h"
   print, "  "
   print, "  "
   return
endif


;-------------------------------------------
;-------------------------------------------


; determine the number of
; snapshots in frun directory
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snap* | wc ",result
nsnaps=long(result[0])




; ---------------------------------------------------


time= fltarr(nsnaps)

sfr_inst= fltarr(nsnaps)
sfr_avg= fltarr(nsnaps)
sfr_snap= fltarr(nsnaps)


ymin = 0 


; physical units
if keyword_set(h) then begin
	h = fload_cosmology('h')
endif


; ---------------------------------------------------


	;-------------------------------------------
	;   Get SFR rate from txt - for each file
	;-------------------------------------------
	open_sfr_file, frun, sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass, gasfrac



	    ; physical units
	    if keyword_set(h) then sfrtime = sfrtime / h
	    ;if i eq 1 then sfrtime = sfrtime/(0.7)    ; did this for comparing h and noh


; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ;ok=fload_snapshot_bh(frun,i,/nopot_in_snap)
        ok=fload_snapshot_bh(frun,i,/skip_center)


	sfr_snap[i]= total(fload_gas_sfr(1))

        ; what time is it?
        time[i]= fload_time(1)

	idx_gtr_snaptime= where(sfrtime ge (time[i]-0.0001))
	idx_curr_snaptime= idx_gtr_snaptime[0]


	; instantaneous sfr
	sfr_inst[i]= -1
	if idx_curr_snaptime(0) ne -1 then sfr_inst[i]= sfrsfr(idx_curr_snaptime)


	; set time window
	dt= 0.01
	idx_window= where((sfrtime ge (time[i]-dt)) and (sfrtime le (time[i]+dt)))
	sfr_avg[i]= -1
	if idx_window(0) ne -1 then sfr_avg[i]= mean(sfrsfr(idx_window))

	print, "T= ", time[i], sfr_inst[i], sfr_avg[i]

endfor



;--------------------------------------
;--------------------------------------


; SFR(snaptimes) FILE
; --------------------
openw, 1, frun+'/sfr_snaptimes.txt', ERROR=err

printf, 1, "#   sfr_snaptimes.txt"
printf, 1, "# "
printf, 1, "#           inst.       avg     from snap   "
printf, 1, "# time       sfr        sfr        sfr      "
printf, 1, "# (Gyr)    (Mo /Yr)   (Mo /Yr)   (Mo /Yr)   "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"    ",3(F8.3,"   "))', $
                time[i], sfr_inst[i], sfr_avg[i], sfr_snap[i]
endfor
close, 1




; ----------------------------------------

print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            "
print, "---------------------------------"
print, "---------------------------------"





end








;====================================================================



; Fit Gaussian to Starburst
;-----------------------------
pro fitandoverplot_gaussian, time, sfr, $
			SFR_max= SFR_max, width=width, SBtime=SBtime, $
			quiet=quiet

        x_tofit= time
        y_tofit= sfr
        weight_tofit= 1.0 + 0.0*sfr
        ;weight_tofit= 1.0 + 0.1*sfr


        ;  parameters
        ;---------------- 
	; p[0] = Normalization, i.e., the SFR_max
        ; p[1] = Width of the dist.; sigma
        ; p[2] = mean of the dist., i.e., time of SFR_max

        ; initial guess
        ;guess= [100.0,0.2,1.5]
	normz= max(sfr)
        idd= where(sfr eq normz)
        meanz= time(idd[0])
        wid= max(time)/3.0
        guess= [normz, wid, meanz]
print, guess

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
        ; Normalization is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; width is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001
        ; mean is greater than 0
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.001


        ; markwardt mpfit procedure
        sb_result = MPFITFUN('func_gaussian', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        print, "--------------------------"
        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        chi2= strcompress(string(redchi2),/remove_all)
        chi2= strmid(chi2,0,4)

        SFR_max= sb_result[0]
        print, "SFR_max= ", SFR_max
        sfrmaxlbl= strcompress(string(SFR_max),/remove_all)
        sfrmaxlbl= strmid(sfrmaxlbl,0,5)

        width= sb_result[1]
        print, "Starburst Width= ", width," Gyr"
        wlbl= strcompress(string(width),/remove_all)
        wlbl= strmid(wlbl,0,4)

        SBtime= sb_result[2]
        print, "Time at SFR_max= ", SBtime," Gyr"
        sbtm= strcompress(string(SBtime),/remove_all)
        sbtm= strmid(sbtm,0,4)


        ; overplot fit
        ; ----------
        ;x=min(time) - 0.2 + (0.2+max(time)-min(time))*findgen(100)/100.0
        ;x= SBtime + 1.5*(findgen(100)/100.0 - 0.5)
        x= SBtime + 4.0*(findgen(100)/100.0 - 0.5)
        y= func_gaussian(x,sb_result)

	idx= where(y gt 0.3*min(sfr))
	y= y(idx)
	x= x(idx)

	if not keyword_set(quiet) then begin
        	oplot, x, y, psym=-3, linestyle= 2, color= 50, thick= 6.0
	endif


        if not keyword_set(quiet) then begin
                ;xyouts, 0.22, 0.75, 'Gaussian ('+chi2+')', size=1.2, color=0, /normal
                ;xyouts, 0.23, 0.71, 'SFRmax= '+sfrmaxlbl, size=1.0, color=0, /normal
                ;xyouts, 0.23, 0.67, 'sigma=  '+wlbl, size=1.0, color=0, /normal
                ;xyouts, 0.23, 0.63, 'sbtime= '+sbtm, size=1.0, color=0, /normal
                ;xyouts, 0.65, 0.84, 'Gaussian', size=1.6, color=0, /normal
                ;xyouts, 0.65, 0.80, '!7r!6= '+wlbl+' Gyr', size=1.5, color=0, /normal
        endif

        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end




;====================================================================


; Fit Gaussian to Starburst
;-----------------------------
pro fitandoverplot_gaussian_pluszeropt, time, sfr, $
                        SFR_max= SFR_max, width=width, SBtime=SBtime

        x_tofit= time
        y_tofit= sfr
        weight_tofit= 1.0 + 0.0*sfr


        ;  parameters
        ;---------------- 
        ; p[0] = Normalization, i.e., the SFR_max
        ; p[1] = Width of the dist.; sigma
        ; p[2] = mean of the dist., i.e., time of SFR_max

        ; initial guess
        guess= [100.0,0.2,1.5,1.0]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
        ; Normalization is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; width is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.01
        ; mean is greater than 0
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.001
        ; zero point is greater than 0
        pi[3].limited(0) = 1
        pi[3].limits(0) = 0.1



        ; markwardt mpfit procedure
        sb_result = MPFITFUN('func_gaussian_pluszeropt', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        ;chi2= strcompress(string(redchi2),/remove_all)
        ;chi2= strmid(chi2,0,4)

        SFR_max= sb_result[0]
        print, "SFR_max= ", SFR_max
        ;rho0lbl= strcompress(string(rho0lbl),/remove_all)
        ;rho0lbl= strmid(rho0lbl,0,5)

        width= sb_result[1]
        print, "Starburst Width= ", width," Gyr"
        ;rslbl= strcompress(string(rslbl),/remove_all)
        ;rslbl= strmid(rslbl,0,5)

        SBtime= sb_result[2]
        print, "Time at SFR_max= ", SBtime," Gyr"

        zeropt= sb_result[3]
        print, "Zero Pt= ", zeropt," M_solar/Yr"


        ; overplot fit
        ; ----------
        x=min(time) - 0.2 + (0.2+max(time)-min(time))*findgen(100)/100.0
        y= func_gaussian_pluszeropt(x,sb_result)
        oplot, x, y, psym=-3, linestyle= 0, color= 200, thick=2.0

        ;xyouts, 0.65, 0.75, 'NFW', /normal, size= 1.0, color=150
        ;xyouts, 0.65, 0.71, 'r!Ds!N='+rslbl, /normal, size= 1.0, color=150
        ;xyouts, 0.77, 0.71, ', c='+clbl, /normal, size= 1.0, color=150


        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end










;====================================================================




; Fit Exponential to Fading Starburst
;--------------------------------------
pro fitandoverplot_expdecay_0, time, sfr, $
			SFR_SB= SFR_SB, tau=tau, SBtime=SBtime

        x_tofit= time
        y_tofit= sfr
        weight_tofit= 1.0 + 0.0*sfr


        ;  parameters
        ;---------------- 
	; p[0] = Normalization, i.e., the SFR_max
        ; p[1] = Width of the dist.; sigma
        ; p[2] = mean of the dist., i.e., time of SFR_max

        ; initial guess
        guess= [100.0,1.0]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)
        ; SFR_SB is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; tau is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001


        ; markwardt mpfit procedure
        tau_result = MPFITFUN('func_exp', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)

        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        ;chi2= strcompress(string(redchi2),/remove_all)
        ;chi2= strmid(chi2,0,4)

        SFR_SB= tau_result[0]
	print, "SFR_SB= ", SFR_SB
        ;rho0lbl= strcompress(string(rho0lbl),/remove_all)
        ;rho0lbl= strmid(rho0lbl,0,5)

        tau= tau_result[1]
	print, "SB tau= ", tau," Gyr"
        ;rslbl= strcompress(string(rslbl),/remove_all)
        ;rslbl= strmid(rslbl,0,5)


        ; overplot fit
        ; ----------
        x=time
        y= func_exp(x,tau_result)
        x=x+SBtime
        oplot, x, y, psym=-3, linestyle= 0, color= 100

        ;xyouts, 0.65, 0.75, 'NFW', /normal, size= 1.0, color=150
        ;xyouts, 0.65, 0.71, 'r!Ds!N='+rslbl, /normal, size= 1.0, color=150
        ;xyouts, 0.77, 0.71, ', c='+clbl, /normal, size= 1.0, color=150


        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end





;====================================================================




; Fit Exponential to Fading Starburst
;--------------------------------------
pro fitandoverplot_expdecay, time, sfr, $
                        SFR_SB= SFR_SB, tau=tau, SBtime=SBtime, decay_const=decay_const, $
                        quiet=quiet

        x_tofit= time
        y_tofit= sfr
        ;weight_tofit= 1.0 + 0.0*sfr
        weight_tofit= 0.1*sfr + 1.0


        ;  parameters
        ;---------------- 
        ; p[0] = Normalization, i.e., the SFR_max ( - p[2])
        ; p[1] = Exp. decay, i.e., tau
        ; p[2] = constant, so it decays to this, rather than 0

        ; initial guess
        guess= [100.0,1.0,0.1]

        ; constraints
        pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
        ; SFR_SB is greater than 0
        pi[0].limited(0) = 1
        pi[0].limits(0) = 0.001
        ; tau is greater than 0
        pi[1].limited(0) = 1
        pi[1].limits(0) = 0.001
        ; const is greater than 0
        pi[2].limited(0) = 1
        pi[2].limits(0) = 0.00001


        ; markwardt mpfit procedure
        tau_result = MPFITFUN('func_exp_sf', x_tofit, y_tofit, weight_tofit, guess, $
                                PARINFO=pi, BESTNORM=bestnorm, DOF=dof)


        redchi2= bestnorm/dof
        print, "Reduced Chi^2 = ", redchi2
        ;chi2= strcompress(string(redchi2),/remove_all)
        ;chi2= strmid(chi2,0,4)

        SFR_SB= tau_result[0]
        print, "SFR_SB= ", SFR_SB, " M_solar/Yr"
        ;rho0lbl= strcompress(string(rho0lbl),/remove_all)
        ;rho0lbl= strmid(rho0lbl,0,5)

        tau= tau_result[1]
        print, "SB tau= ", tau," Gyr"
        taulbl= strcompress(string(tau),/remove_all)
        taulbl= strmid(taulbl,0,4)

        decay_const= tau_result[2]
        print, "decay const= ", decay_const," M_solar/Yr"
        ;rslbl= strcompress(string(rslbl),/remove_all)
        ;rslbl= strmid(rslbl,0,5)


        ; overplot fit
        ; ----------
        x=time
        y= func_exp_sf(x,tau_result)
        x=x+SBtime
        if not keyword_set(quiet) then begin
                oplot, x, y, psym=-3, linestyle= 0, color= 200, thick= 6.0
        endif

        ;xyouts, 0.65, 0.75, 'NFW', /normal, size= 1.0, color=150
        ;xyouts, 0.65, 0.71, 'r!Ds!N='+rslbl, /normal, size= 1.0, color=150
        ;xyouts, 0.77, 0.71, ', c='+clbl, /normal, size= 1.0, color=150


        if not keyword_set(quiet) then begin
                ;xyouts, 0.65, 0.71, 'Tau Model', color=0, size=1.5, /normal
                ;xyouts, 0.65, 0.66, '!7s!6!DSF!N= '+taulbl+' Gyr', color=0, size=1.5, /normal
        endif

        ; reduced chi squared
        ;xyouts, 0.45, 0.26, '!7v!3!E2!N='+chi2_nfw, /normal, size= 1.0, color=150

end





;====================================================================




