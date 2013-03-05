;------------------------------------------------------------------------
;
;    Timing Info
;
;
;
;
;------------------------------------------------------------------------


pro ti, junk


if not keyword_set(junk) then begin
   print, "  "
   print, " ti, junk"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='ti.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; -----------------------------

xmax= 75
xmin= 0
xaxistitle='!6Processors'

ymax = 600.0
ymin = 0.0
yaxistitle= '!6Time Steps'



; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
	;/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; -----------------------------------------------


; SAURON (PGadget2)
proc= [4, 8, 16, 32, 64]
tsteps= [128, 246, 342, 374, 407]
oplot, proc, tsteps, psym=-2, color=50, thick=6.0, symsize=1.3
oplot, [23.0,25.0], [170.0,170.0], psym=-2, color=50, thick=3.0
xyouts, 28.0, 160.0, 'Sauron', size=1.2, color=50, /data, charthick=4.0
;xyouts, 28.0, 160.0, 'PGadget2 (mpich1)', size=1.2, color=50, /data, charthick=4.0


; Microway test cluster
;   1node= 4 2.2 GHz Opteron, 3 Gbyte DDR RAM
;   Mellanox SDR Infiniband
;
proc= [4, 8, 16, 32, 64]
tsteps= [158, 262, 358, 512, 561]
oplot, proc, tsteps, psym=-5, color=100, thick=6.0, symsize=1.3
oplot, [23.0,25.0], [130.0,130.0], psym=-5, color=100, thick=3.0
xyouts, 28.0, 120.0, 'Microway Test Cluster', size=1.2, color=100, /data, charthick=4.0



; -----------------------------------------------

doit= 0
if doit eq 1 then begin
	; SAURON (PGadget3 - mpich1)
	proc= [8, 16, 32, 64]
	tsteps= [278, 474, 268, 263]
	oplot, proc, tsteps, psym=-2, color=20, thick=6.0, symsize=1.3
	oplot, [23.0,25.0], [145.0,145.0], psym=-2, color=20, thick=3.0
	xyouts, 28.0, 135.0, 'PG3 (mpich1, MD=8)', size=1.2, color=20, /data, charthick=4.0

	; SAURON (PGadget3 - mpich1)
	proc= [16, 32, 64]
	tsteps= [318, 312, 314]
	oplot, proc, tsteps, psym=-2, color=110, thick=6.0, symsize=1.3
	oplot, [23.0,25.0], [120.0,120.0], psym=-2, color=110, thick=3.0
	xyouts, 28.0, 110.0, 'PG3 (mpich1, MD=4)', size=1.2, color=110, /data, charthick=4.0
endif


; -----------------------------------------------

doit= 0
if doit eq 1 then begin
	; SAURON (PGadget3 - mpich2)
	proc= [4, 8, 16, 32, 64]
	tsteps= [98, 178, 194, 258, 193]
	oplot, proc, tsteps, psym=-2, color=150, thick=6.0, symsize=1.3
	oplot, [23.0,25.0], [95.0,95.0], psym=-2, color=150, thick=3.0
	xyouts, 28.0, 85.0, 'PG3 (mpich2, MD=4)', size=1.2, color=150, /data, charthick=4.0

	; SAURON (PGadget3 - mpich2)
	proc= [4, 8, 16, 32, 64]
	tsteps= [118, 198, 258, 322, 310]
	oplot, proc, tsteps, psym=-2, color=180, thick=6.0, symsize=1.3
	oplot, [23.0,25.0], [70.0,70.0], psym=-2, color=180, thick=3.0
	xyouts, 28.0, 60.0, 'PG3 (mpich2, MD=8)', size=1.2, color=180, /data, charthick=4.0

	; SAURON (PGadget3 - mpich2)
	proc= [4, 8, 16, 32, 64]
	tsteps= [122, 210, 351, 378,346]
	oplot, proc, tsteps, psym=-2, color=220, thick=6.0, symsize=1.3
	oplot, [23.0,25.0], [35.0,35.0], psym=-2, color=220, thick=3.0
	xyouts, 28.0, 35.0, 'PG3 (mpich2, MD=16)', size=1.2, color=220, /data, charthick=4.0


	xyouts, 40.0, 520.0, 'Sauron', size=1.7, color=0, /data, charthick=4.0
endif

; -----------------------------------------------




; done
; -----
device, /close


end
   




;--------------------------------------------------------------------------




pro ti2, junk


if not keyword_set(junk) then begin
   print, "  "
   print, " ti2, junk"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='ti.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; -----------------------------

xmax= 275
;xmax= 128
;xmax= 90
xmin= 0
xaxistitle='!6Number of Cores'

ymax = 1.1
;ymin = 0.1
ymin = 0.05
yaxistitle= '!6Normalized speed (time)'



; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
	;/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; -----------------------------------------------
; **  DONE **
;
; test # 1   SbSbhs_e_xx   (no letter)
;
; Sb x Sbhs (1x)
; MULTIPLEDOMAINS= 8
; ptile= 4
;
;proc = [8, 16, 24, 32, 64, 128]
std_proc = [8, 16, 24, 32, 48, 64]      ; > 64 never really went anywhere
std_cputime_02= [  9277.26,   5064.88,   3398.12,   2641.59,   1891.64,   4491.89]
std_cputime_10= [ 69063.64,  38625.93,  25536.35,  20760.91,  15572.93,  33972.72]
std_cputime_20= [300905.72, 206482.58, 155754.32, 159670.44, 135239.13, 152719.80]
std_cputime_25= [435838.68, 360715.88, 270741.31, 260830.85, 206180.07, 207720.38]

;
;  things scale really well all the way to 32 processors at times 0.2 and 1.0,
;  and have a very fast scaling.  At times 2.0 and 2.5, the scaling is still 
;  decent to 24 processors, but with a smaller speed-up factor.
;
;

oplot, std_proc, std_cputime_02/std_cputime_02[0], psym=-2, color=50, thick=1.0, symsize=1.3, linestyle= 1
oplot, std_proc, std_cputime_10/std_cputime_10[0], psym=-2, color=50, thick=1.0, symsize=1.3, linestyle= 1
oplot, std_proc, std_cputime_20/std_cputime_20[0], psym=-2, color=50, thick=1.0, symsize=1.3, linestyle= 1
oplot, std_proc, std_cputime_25/std_cputime_25[0], psym=-2, color=50, thick=6.0, symsize=1.3
oplot, [138.0,152.0], [0.74,0.74], psym=-2, color=50, thick=6.0, symsize= 1.3
xyouts, 162.0, 0.72, '1x, ptile=4, MD=8', size=1.0, color=50, /data, charthick=2.0



; -----------------------------------------------
; ** DONE **
;
; test # 2   SbSbhs_e_xxh   (i.e., the "h" runs)
;
; Sb x Sbhs (1x)
; MULTIPLEDOMAINS= 32
; ptile= 1
;
;proc = [8, 16, 24, 32, 64, 128]
;h_proc = [8, 16, 24, 32, 64]
h_proc = [8, 16, 24, 32, 48, 64]      ; > 64 never really went anywhere
;h_cputime_02= [  7464.32,   4053.45,   2862.22,   5552.00,   3646.47]
;h_cputime_10= [ 54309.08,  29890.75,  21390.72,  39667.79,  25706.13]
;h_cputime_20= [240507.07, 136280.64,  93968.58, 153828.08, 106208.63]
;h_cputime_25= [356581.44, 199628.14, 135809.30, 212786.79, 147450.36]
; used the newer version of SbSbhs_e_32h2 rather then _32h
h_cputime_02= [  7464.32,   4053.45,   2862.22,   2180.00,   2246.75,   3646.47]
h_cputime_10= [ 54309.08,  29890.75,  21390.72,  16521.15,  15872.68,  25706.13]
h_cputime_20= [240507.07, 136280.64,  93968.58,  73386.01,  65397.35, 106208.63]
h_cputime_25= [356581.44, 199628.14, 135809.30, 106928.35,  92869.02, 147450.36]

;
;  it's hard to tell what's going on until things finish
;
;

;oplot, h_proc, h_cputime_02/h_cputime_02[0], psym=-4, color=50, thick=6.0, symsize=1.3
;oplot, h_proc, h_cputime_10/h_cputime_10[0], psym=-4, color=0, thick=6.0, symsize=1.3
;oplot, h_proc, h_cputime_20/h_cputime_20[0], psym=-4, color=100, thick=6.0, symsize=1.3
;oplot, h_proc, h_cputime_25/h_cputime_25[0], psym=-4, color=150, thick=6.0, symsize=1.3
oplot, h_proc, h_cputime_02/std_cputime_02[0], psym=-4, color=150, thick=1.0, symsize=1.0, linestyle= 1
oplot, h_proc, h_cputime_10/std_cputime_10[0], psym=-4, color=150, thick=1.0, symsize=1.0, linestyle= 1
oplot, h_proc, h_cputime_20/std_cputime_20[0], psym=-4, color=150, thick=1.0, symsize=1.0, linestyle= 1
;oplot, h_proc, h_cputime_25/std_cputime_25[0], psym=-4, color=150, thick=1.0, symsize=1.0, linestyle= 1
oplot, h_proc, h_cputime_25/std_cputime_25[0], psym=-4, color=150, thick=6.0, symsize=1.3
oplot, [138.0,153.0], [0.66,0.66], psym=-4, color=150, thick=3.0
;xyouts, 80.0, 0.67, "'h'", size=1.2, color=150, /data, charthick=4.0
xyouts, 162.0, 0.64, '1x, ptile=1, MD=32', size=1.0, color=150, /data, charthick=2.0




; -----------------------------------------------
; ** DONE **
;
; test # 3   SbSbhs_e_xxe   (i.e., the "e" runs)
;
; Sb x Sbhs (1x)
; MULTIPLEDOMAINS= 32
; ptile= 4
;
;proc = [8, 16, 24, 32, 64, 128]
e_proc = [8, 16, 24, 32, 48, 64]      ; > 64 never really went anywhere
e_cputime_02= [  8164.96,   4294.16,   3714.53,   2313.07,   2059.15,   5751.27]
e_cputime_10= [ 60084.50,  31770.45,  27904.18,  17951.63,  15873.53,  55677.91]
e_cputime_20= [247980.04, 144415.11, 113622.76,  80104.37,  67711.99, 175289.97]
e_cputime_25= [359354.51, 210617.96, 162497.10, 115747.66,  96510.53, 235107.14]

;
;  it's hard to tell what's going on until things finish
;
;

;oplot, e_proc, e_cputime_02/e_cputime_02[0], psym=-7, color=50, thick=6.0, symsize=1.3
;oplot, e_proc, e_cputime_10/e_cputime_10[0], psym=-7, color=0, thick=6.0, symsize=1.3
;oplot, e_proc, e_cputime_20/e_cputime_20[0], psym=-7, color=100, thick=6.0, symsize=1.3
;oplot, e_proc, e_cputime_25/e_cputime_25[0], psym=-7, color=200, thick=6.0, symsize=1.3
oplot, e_proc, e_cputime_02/std_cputime_02[0], psym=-7, color=200, thick=1.0, symsize=1.0, linestyle= 1
oplot, e_proc, e_cputime_10/std_cputime_10[0], psym=-7, color=200, thick=1.0, symsize=1.0, linestyle= 1
oplot, e_proc, e_cputime_20/std_cputime_20[0], psym=-7, color=200, thick=1.0, symsize=1.0, linestyle= 1
;oplot, e_proc, e_cputime_25/std_cputime_25[0], psym=-7, color=200, thick=1.0, symsize=1.0, linestyle= 1
oplot, e_proc, e_cputime_25/std_cputime_25[0], psym=-7, color=200, thick=6.0, symsize=1.3
oplot, [138.0,153.0], [0.58,0.58], psym=-7, color=200, thick=3.0
;xyouts, 80.0, 0.67, "'h'", size=1.2, color=150, /data, charthick=4.0
xyouts, 162.0, 0.56, '1x, ptile=4, MD=32', size=1.0, color=200, /data, charthick=2.0





; -----------------------------------------------

;
; test # 4   Sb10xSb10xhs_e_xx   (std 10x runs)
;
; Sb x Sbhs (10x)
; MULTIPLEDOMAINS= 16
; ptile= 4
;
x10_proc =      [     8,           16,         32,         64,        128,       256]
x10_cputime_01= [ 152918.53,   81138.69,   44641.83,   24798.18,   16941.55,   14163.95]
x10_cputime_02= [ 311497.69,  166081.18,   91339.13,   50973.63,   34910.39,   28704.65]
x10_cputime_04= [ 715934.79,  387660.26,  214472.24,  120580.09,   82715.60,   88534.72]
x10_cputime_05= [ 935796.85,  506735.39,  281360.50,  159739.00,  109996.22,  126598.18]
x10_cputime_08= [1850739.11, 1005061.62,  565402.37,  320372.22,  272869.75,  279728.83]
x10_cputime_10= [2588878.56, 1422376.24,  805324.80,  450098.48,  436857.15,  389244.45]
x10_cputime_12= [3800000.00, 2079809.51, 1194729.29,  794054.36,  681055.97,  539273.93]
x10_cputime_13= [6900000.00, 3625145.76, 2155907.58, 1412557.23, 1075303.32,  805542.88]
;x10_cputime_14= [          ,           ,           , 2652506.78, 1802440.00, 1343610.19]
;x10_cputime_15= [          ,           ,           , 3861780.55, 2500183.85, 1889303.26]
;x10_cputime_16= [          ,           ,           ,           , 3181586.73, 2468708.69]
;x10_cputime_17= [          ,           ,           ,           , 3839236.66, 2976649.09]
;x10_cputime_18= [          ,           ,           ,           ,           , 3474604.69]
;x10_cputime_19= [          ,           ,           ,           ,           , 3954831.17]
;x10_cputime_20= [          ,           ,           ,           ,           , 
;x10_cputime_25= [          ,           ,           ,           ,           , 

;
;  it's hard to tell what's going on until things finish
;
;

oplot, x10_proc, x10_cputime_01/x10_cputime_01[0], psym=-5, color=100, thick=1.0, symsize=1.0, linestyle= 1
oplot, x10_proc, x10_cputime_02/x10_cputime_02[0], psym=-5, color=100, thick=1.0, symsize=1.0, linestyle= 1
oplot, x10_proc, x10_cputime_04/x10_cputime_04[0], psym=-5, color=100, thick=1.0, symsize=1.0, linestyle= 1
;oplot, x10_proc, x10_cputime_05/x10_cputime_05[0], psym=-5, color=100, thick=1.0, symsize=1.0, linestyle= 1
oplot, x10_proc, x10_cputime_05/x10_cputime_05[0], psym=-5, color=100, thick=1.0, symsize=1.0, linestyle= 1
oplot, x10_proc, x10_cputime_08/x10_cputime_08[0], psym=-5, color=100, thick=1.0, symsize=1.0, linestyle= 1
oplot, x10_proc, x10_cputime_10/x10_cputime_10[0], psym=-5, color=100, thick=6.0, symsize=1.3
oplot, x10_proc, x10_cputime_12/x10_cputime_12[0], psym=-5, color=120, thick=1.0, symsize=1.8, linestyle= 2
oplot, x10_proc, x10_cputime_13/x10_cputime_13[0], psym=-5, color=120, thick=1.0, symsize=1.8, linestyle= 2
;oplot, x10_proc, x10_cputime_20/x10_cputime_20[0], psym=-5, color=150, thick=6.0, symsize=1.3
;oplot, x10_proc, x10_cputime_25/x10_cputime_25[0], psym=-5, color=100, thick=6.0, symsize=1.3
oplot, [138.0,153.0], [0.51,0.51], psym=-5, color=100, thick=3.0
;xyouts, 80.0, 0.47, "'h'", size=1.2, color=150, /data, charthick=4.0
xyouts, 162.0, 0.50, '10x, ptile=4, MD=16', size=1.0, color=100, /data, charthick=2.0
;oplot, [138.0,153.0], [0.48,0.48], psym=-5, color=100, thick=3.0
;;xyouts, 80.0, 0.47, "'h'", size=1.2, color=150, /data, charthick=4.0
;xyouts, 162.0, 0.47, '10x, ptile=4, MD=16', size=1.0, color=100, /data, charthick=2.0





; -----------------------------------------------

;
; test # 5   Sb5x1Sb5x1_e_xx   (5x runs, but have higher baryon resolution)
;
; Sb x Sb (50x)
; MULTIPLEDOMAINS= 
; ptile= 4
;
;x5_proc =      [      8,       16,         32,       48,           64,          96]
;x5_cputime_01= [   750000.,  370000.,  185626.43,   40693.16,   29823.21,   69032.99]
;x5_cputime_02= [  1655000.,  827000.,  413490.61,   99755.66,   71763.83,  155446.12]
;x5_cputime_04= [  4100000., 2050000., 1024670.37,  306799.75,  204051.45,  354537.48]
;x5_cputime_05= [ 15400000., 7690000., 3843000.00, 1087880.33,  766518.02]
x5_proc =      [       32,       48,           64,          96]
x5_cputime_01= [    185626.43,   40693.16,   29823.21,   69032.99]
x5_cputime_02= [    413490.61,   99755.66,   71763.83,  155446.12]
x5_cputime_04= [   1024670.37,  306799.75,  204051.45,  354537.48]
x5_cputime_05= [   3843000.00, 1087880.33,  766518.02]
;x5_cputime_08= [          ,         ,           ,           ,  2456522.19]
;x5_cputime_10= [
;x5_cputime_12= [
;x5_cputime_13= [
;x5_cputime_14= [
;x5_cputime_15= [
;x5_cputime_20= [
;x5_cputime_25= [

;
;  it's hard to tell what's going on until things finish
;
;

;oplot, x5_proc, (1./4.) * x5_cputime_01/x5_cputime_01[0], psym=-5, color=0, thick=1.0, symsize=1.0, linestyle= 1
;oplot, x5_proc, (1./4.) * x5_cputime_02/x5_cputime_02[0], psym=-5, color=0, thick=1.0, symsize=1.0, linestyle= 1
;oplot, x5_proc, (1./4.) * x5_cputime_04/x5_cputime_04[0], psym=-5, color=0, thick=1.0, symsize=1.0, linestyle= 1
;oplot, x5_proc, x5_cputime_05/x5_cputime_05[0], psym=-5, color=100, thick=1.0, symsize=1.0, linestyle= 1
;oplot, x5_proc, x5_cputime_05/x5_cputime_05[0], psym=-5, color=100, thick=1.0, symsize=1.0, linestyle= 1
;oplot, x5_proc, x5_cputime_08/x5_cputime_08[0], psym=-5, color=100, thick=1.0, symsize=1.0, linestyle= 1
;oplot, x5_proc, x5_cputime_10/x5_cputime_10[0], psym=-5, color=100, thick=6.0, symsize=1.3
;oplot, x5_proc, x5_cputime_12/x5_cputime_12[0], psym=-5, color=120, thick=1.0, symsize=1.8, linestyle= 2
;oplot, x5_proc, x5_cputime_13/x5_cputime_13[0], psym=-5, color=120, thick=1.0, symsize=1.8, linestyle= 2
;oplot, x5_proc, x5_cputime_20/x5_cputime_20[0], psym=-5, color=150, thick=6.0, symsize=1.3
;oplot, x5_proc, x5_cputime_25/x5_cputime_25[0], psym=-5, color=100, thick=6.0, symsize=1.3
;oplot, [138.0,153.0], [0.51,0.51], psym=-5, color=0, thick=3.0
;xyouts, 162.0, 0.50, '5x, ptile=4, MD=24', size=1.0, color=0, /data, charthick=2.0





; -----------------------------------------------


x= 8.0 + 3.0 * findgen(80)
y= x[0]/x
oplot, x, y, psym=-3, linestyle= 1, color= 0, thick= 3.0
oplot, [135.0,157.0], [0.84,0.84], psym=-3, linestyle= 1, color=0, thick=3.0
xyouts, 162.0, 0.81, 'perfect scaling', size=1.2, color=0, /data, charthick=1.0


; done
; -----
device, /close


end
   




;--------------------------------------------------------------------------
;--------------------------------------------------------------------------




pro ti3, junk


if not keyword_set(junk) then begin
   print, "  "
   print, " ti3, junk"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='ti3.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


; -----------------------------
; -----------------------------

xmax= 27.5
xmin= 0
;xaxistitle='!6Particle Number (x 10!E6!N)'
xmax= 3.5
xaxistitle='!6Gas Particle Number (x 10!E6!N)'

ymax = 285.0
ymin = 0.2
yaxistitle= '!6Normalized speed (time)'



; ------------------
; Plot this up
; ------------------
!p.position= [0.18, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
	;/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=4.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata


; -----------------------------------------------
; **  DONE **
;
; test # 1   SbSbhs_e_xx   (no letter)
;
; Sb x Sbhs (1x)
; MULTIPLEDOMAINS= 8
; ptile= 4
;
std_N= 2.5  ; in millions
std_Ngas= 0.09
std_proc = [8, 16, 24, 32, 48, 64]      ; > 64 never really went anywhere
std_cputime_02= [  9277.26,   5064.88,   3398.12,   2641.59,   1891.64,   4491.89]
std_cputime_05= [ 29205.93,  15992.18,  10767.04,   8569.88,   6234.52,  14263.62]
std_cputime_10= [ 69063.64,  38625.93,  25536.35,  20760.91,  15572.93,  33972.72]
std_cputime_20= [300905.72, 206482.58, 155754.32, 159670.44, 135239.13, 152719.80]
std_cputime_25= [435838.68, 360715.88, 270741.31, 260830.85, 206180.07, 207720.38]




; -----------------------------------------------
; ** DONE **
;
; test # 2   SbSbhs_e_xxh   (i.e., the "h" runs)
;
; Sb x Sbhs (1x)
; MULTIPLEDOMAINS= 32
; ptile= 1
;
h_N= 2.5
h_Ngas= 0.09
h_proc = [8, 16, 24, 32, 48, 64]      ; > 64 never really went anywhere
h_cputime_02= [  7464.32,   4053.45,   2862.22,   2180.00,   2246.75,   3646.47]
h_cputime_10= [ 54309.08,  29890.75,  21390.72,  16521.15,  15872.68,  25706.13]
h_cputime_20= [240507.07, 136280.64,  93968.58,  73386.01,  65397.35, 106208.63]
h_cputime_25= [356581.44, 199628.14, 135809.30, 106928.35,  92869.02, 147450.36]




; -----------------------------------------------
; ** DONE **
;
; test # 3   SbSbhs_e_xxe   (i.e., the "e" runs)
;
; Sb x Sbhs (1x)
; MULTIPLEDOMAINS= 32
; ptile= 4
;
e_N= 2.5
e_Ngas= 0.09
e_proc = [8, 16, 24, 32, 48, 64]      ; > 64 never really went anywhere
e_cputime_02= [  8164.96,   4294.16,   3714.53,   2313.07,   2059.15,   5751.27]
e_cputime_10= [ 60084.50,  31770.45,  27904.18,  17951.63,  15873.53,  55677.91]
e_cputime_20= [247980.04, 144415.11, 113622.76,  80104.37,  67711.99, 175289.97]
e_cputime_25= [359354.51, 210617.96, 162497.10, 115747.66,  96510.53, 235107.14]





; -----------------------------------------------

;
; test # 4   Sb10xSb10xhs_e_xx   (std 10x runs)
;
; Sb x Sbhs (10x)
; MULTIPLEDOMAINS= 16
; ptile= 4
;
x10_N= 25.0
x10_Ngas= 0.9
x10_proc =      [     8,           16,         32,         64,        128,       256]
x10_cputime_01= [ 152918.53,   81138.69,   44641.83,   24798.18,   16941.55,   14163.95]
x10_cputime_02= [ 311497.69,  166081.18,   91339.13,   50973.63,   34910.39,   28704.65]
x10_cputime_04= [ 715934.79,  387660.26,  214472.24,  120580.09,   82715.60,   88534.72]
x10_cputime_05= [ 935796.85,  506735.39,  281360.50,  159739.00,  109996.22,  126598.18]
x10_cputime_08= [1850739.11, 1005061.62,  565402.37,  320372.22,  272869.75,  279728.83]
x10_cputime_10= [2588878.56, 1422376.24,  805324.80,  450098.48,  436857.15,  389244.45]
x10_cputime_12= [3800000.00, 2079809.51, 1194729.29,  794054.36,  681055.97,  539273.93]
x10_cputime_13= [6900000.00, 3750000.00, 2155907.58, 1412557.23, 1075303.32,  805542.88]
;x10_cputime_14= [          ,           ,           , 2652506.78, 1802440.00, 1343610.19]
;x10_cputime_15= [          ,           ,           ,           , 2500183.85, 1889303.26]
;x10_cputime_20= [        ,          ,
;x10_cputime_25= [        ,          ,





; -----------------------------------------------

;
; test # 5   Sb5x1Sb5x1_e_xx   (5x runs, but have higher baryon resolution)
;
; Sb x Sb (5x)
; MULTIPLEDOMAINS= 
; ptile= 4
;
x5_N= 15.5
x5_Ngas= 2.4
x5_proc =      [       32,       48,           64,          96]
x5_cputime_01= [    185626.43,   40693.16,   29823.21,   69032.99]
x5_cputime_02= [    413490.61,   99755.66,   71763.83,  155446.12]
x5_cputime_04= [   1024670.37,  306799.75,  204051.45,  354537.48]
x5_cputime_05= [   3843000.00, 1087880.33,  766518.02]
;x5_cputime_08= [          ,         ,           ,           ,  2456522.19]
;x5_cputime_10= [
;x5_cputime_12= [
;x5_cputime_13= [
;x5_cputime_14= [
;x5_cputime_15= [
;x5_cputime_20= [
;x5_cputime_25= [




; -----------------------------------------------

;
; test # 6   Sb3x1Sb3x1_e_xx   (3x runs, but have higher baryon resolution)
;
; Sb x Sb (3x)
; MULTIPLEDOMAINS= 
; ptile= 4
;
x3_N= 7.75
x3_Ngas= 1.2
x3_proc =      [       16,       32]
x3_cputime_01= [    56816.24,   15746.54]
x3_cputime_02= [   114889.31,   32608.83]
x3_cputime_04= [   310310.76,  100403.15]
x3_cputime_05= [   604256.19,  269435.52]
x3_cputime_08= [  1338375.47,  754717.81]
;x3_cputime_10= [            , 1050014.21]
;x3_cputime_12= [
;x3_cputime_13= [
;x3_cputime_14= [
;x3_cputime_15= [
;x3_cputime_20= [
;x3_cputime_25= [




; -----------------------------------------------

;
; test # 8   Sb2x1Sb2x1_e_xx   (2x runs, but have higher baryon resolution)
;
; Sb x Sb (2x)
; MULTIPLEDOMAINS= 16
; ptile= 4
;
x2_N= 3.875
x2_Ngas= 0.6
x2_proc =      [       16,       32]
x2_cputime_01= [    34018.13,    21365.89]
x2_cputime_02= [    64107.56,    44478.66]
x2_cputime_04= [   198057.18,   128082.39]
x2_cputime_05= [   344500.31,   256982.59]
x2_cputime_08= [   566349.53,   576504.16]
x2_cputime_10= [   711302.79,   749820.32]
x2_cputime_12= [   853430.87,   957158.14]
x2_cputime_13= [   941672.66,  1079820.00]
;x2_cputime_14= [  1263520.61, ]
;x2_cputime_15= [
;x2_cputime_20= [
;x2_cputime_25= [





; -----------------------------------------------



Ns= [std_N, x2_N, x3_N, x5_N, x10_N]
Ngass= [std_Ngas, x2_Ngas, x3_Ngas, x5_Ngas, x10_Ngas]

;
; t= 0.2
; 16 proc, except 5x
;
Ns= [std_N, x2_N, x3_N, x5_N, x10_N]
Ns= Ngass
speed= [std_cputime_02[1], x2_cputime_02[0], x3_cputime_02[0], x5_cputime_02[0], x10_cputime_02[1]]
relspeed= speed/speed[0]
oplot, Ns, relspeed, psym= 4, linestyle= 0, color= 10, thick= 3.0, symsize= 2.0

;
; t= 0.2
; 32 proc
;
Ns= [std_N, x2_N, x3_N, x5_N, x10_N]
Ns= Ngass
speed= [std_cputime_02[3], x2_cputime_02[1], x3_cputime_02[1], x5_cputime_02[0], x10_cputime_02[2]]
relspeed= speed/speed[0]
oplot, Ns, relspeed, psym= 5, linestyle= 0, color= 50, thick= 3.0, symsize= 2.0

;
; t= 0.2   (use e and h instead of std)
; 32 proc
;
Ns= [std_N, x2_N, x3_N, x5_N, x10_N]
Ns= Ngass
speed= [e_cputime_02[3], x2_cputime_02[1], x3_cputime_02[1], x5_cputime_02[0], x10_cputime_02[2]]
relspeed= speed/speed[0]
oplot, Ns, relspeed, psym= 5, linestyle= 0, color= 80, thick= 3.0, symsize= 2.0
speed= [h_cputime_02[3], x2_cputime_02[1], x3_cputime_02[1], x5_cputime_02[0], x10_cputime_02[2]]
relspeed= speed/speed[0]
oplot, Ns, relspeed, psym= 6, linestyle= 0, color= 120, thick= 3.0, symsize= 2.0

;
; t= 0.5   (use e and h instead of std)
; 32 proc
;
Ns= [std_N, x2_N, x3_N, x5_N, x10_N]
Ns= Ngass
speed= [std_cputime_05[3], x2_cputime_05[1], x3_cputime_05[1], x5_cputime_05[0], x10_cputime_05[2]]
relspeed= speed/speed[0]
oplot, Ns, relspeed, psym= 5, linestyle= 0, color= 140, thick= 3.0, symsize= 2.0

;
; t= 1.0  (for those that have it)
; 32 proc
;
Ns= [std_N, x2_N, x10_N]
Ns= [std_Ngas, x2_Ngas, x10_Ngas]
speed= [std_cputime_10[3], x2_cputime_02[1], x10_cputime_10[2]]
relspeed= speed/speed[0]
oplot, Ns, relspeed, psym= 1, linestyle= 0, color= 180, thick= 3.0, symsize= 2.0
speed= [e_cputime_10[3], x2_cputime_02[1], x10_cputime_10[2]]
relspeed= speed/speed[0]
oplot, Ns, relspeed, psym= 2, linestyle= 0, color= 220, thick= 3.0, symsize= 2.0
speed= [h_cputime_10[3], x2_cputime_02[1], x10_cputime_10[2]]
relspeed= speed/speed[0]
oplot, Ns, relspeed, psym= 7, linestyle= 0, color= 200, thick= 3.0, symsize= 2.0





; -----------------------------------------------


; For N
; --------
x= 2.0 + 0.333 * findgen(80)   ; OK for N
y= x * alog(x)
x= x - 1.8
oplot, x, y, psym=-3, linestyle= 1, color= 0, thick= 3.0
oplot, [11.0,14.0], [0.90,0.90], psym=-3, linestyle= 1, color=0, thick=3.0
;xyouts, 15.0, 0.8, 'N LogN ', size=1.5, color=0, /data, charthick=2.0
xyouts, 2.2, 1.9, 'N LogN ', size=1.5, color=0, /data, charthick=2.0


;  For Ngas
; -----------
x= 0.1 + 0.2 * findgen(20)    ; OK for Ngas
y= x -1.0
oplot, x, y, psym=-3, linestyle= 1, color= 50, thick= 3.0
oplot, [11.0,14.0], [0.90,0.90], psym=-3, linestyle= 2, color=50, thick=3.0
xyouts, 2.2, 0.8, 'N', size=1.5, color=50, /data, charthick=2.0



; done
; -----
device, /close


end
   




;--------------------------------------------------------------------------



