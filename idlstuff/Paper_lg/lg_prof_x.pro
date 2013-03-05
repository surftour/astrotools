;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Projected Gas Profiles
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------






; --------------------------
;  Gas Profile - projected
; --------------------------
pro plot_gas_profile, frun, sendto, smoothlen=smoothlen, filename=filename

if not keyword_set(sendto) then sendto='ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_gas_profile"
   print, "  "
   return
endif


if fload_npart(0) eq 0 then begin
	print, " "
	print, " plot_gas_profile didn't find any gas particles"
	print, " "
	return
endif

setup_plot_stuff, sendto, filename=filename



; -----------------------------
; set up constants
; -----------------------------

bins = 50

;xislog= 1
xislog= 0      ; manually set /xlog

if xislog eq 1 then begin
    xmax = alog10(300)
    xmin = 0.0
    xaxistitle= 'Log R (kpc)'
endif else begin
    xmax = 300
    xmin = 1.0
    xaxistitle= 'R (h!E-1!Nkpc)'
endelse

ymax = 1e+4
ymin = 1e-4

;yaxistitle= 'Log !4R!3(R) !N(M!D!9n!3!Npc!E-2!N)'
yaxistitle= '!4R!3(R) !N(M!D!9n!3!Npc!E-2!N)'


smoothlen2= alog10(smoothlen)

;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
	/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




process_and_plot_onesimgasdensity, xmax, xmin, bins, frun, snapnum, /includeerr


; ---------------
; print extras
; -------------
x=[smoothlen2,smoothlen2]
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0

x0= 0.65
y0= 0.86
;xyouts, x0, y0, fload_fid(1), /normal, charthick= 2.0, size= 1.7, color= 0

y0= 0.80
xyouts, x0, y0, 'Gas Profile', /normal, charthick= 2.0, size= 1.7, color= 0



; Done
;-------

if (sendto EQ 'ps') then device, /close



end














; --------------------------
;  Gas Profile - projected
; --------------------------
pro plot_gas_profile_comparison, junk, smoothlen=smoothlen, filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_gas_profile"
   print, "  "
   return
endif

if not keyword_set(filename) then filename="gas_g_prof.eps"
if not keyword_set(smoothlen) then smoothlen=0.1


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 1



; -----------------------------
; set up constants
; -----------------------------

bins = 50

;xislog= 1
xislog= 0      ; manually set /xlog

if xislog eq 1 then begin
    xmax = alog10(300)
    xmin = 0.0
    xaxistitle= 'Log R (kpc)'
endif else begin
    xmax = 1000
    ;xmax = 300
    xmin = 1.0
    xaxistitle= 'R (h!E-1!Nkpc)'
endelse

ymax = 1e+4
ymin = 1e-4

;yaxistitle= 'Log !4R!3(R) !N(M!D!9n!3!N pc!E-2!N)'
yaxistitle= '!4R!3(R) !N(M!D!9n!3!N pc!E-2!N)'


smoothlen2= alog10(smoothlen)

;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
	/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




;process_and_plot_onesimgasdensity, xmax, xmin, bins, frun, snapnum, /includeerr
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc3vd3e", 25, linecolor= 150, /includeerr, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc3vd3e", 30, linecolor= 150, xislog=xislog


; orientations with/without BH
; -------------------------------
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc3vc3f", 25, linecolor= 150, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc3vc3h", 25, linecolor= 150, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc3vc3j", 25, linecolor= 150, xislog=xislog

;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc3vc3f_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc3vc3h_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc3vc3j_no", 25, linecolor= 50, xislog=xislog
	

; different size halos
; ----------------------
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc4vc4a", 25, linecolor= 150, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc5vc5a", 25, linecolor= 150, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc6vc6a", 25, linecolor= 150, xislog=xislog

;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc4vc4_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc5vc5_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/vc6vc6_no", 25, linecolor= 50, xislog=xislog



; different gas fractions
; ----------------------
process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/gfs/vc3vc3u_k", 30, linecolor= 200, xislog=xislog
process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/gfs/vc3vc3v_k", 30, linecolor= 170, xislog=xislog
process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/gfs/vc3vc3w_k", 30, linecolor= 140, xislog=xislog
process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/gfs/vc3vc3x2_k", 30, linecolor= 110, xislog=xislog
process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/gfs/vc3vc3y_k", 30, linecolor= 80, xislog=xislog
process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid4/tcox/gfs/vc3vc3z_k", 30, linecolor= 50, xislog=xislog



; group halo
; -----------
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 0, linecolor= 200, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 1, linecolor= 180, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 2, linecolor= 160, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 3, linecolor= 140, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 4, linecolor= 120, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 5, linecolor= 100, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 10, linecolor= 80, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 20, linecolor= 60, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 30, linecolor= 40, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 40, linecolor= 20, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 50, linecolor= 0, xislog=xislog

;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA_1", 0, linecolor= 200, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA_1", 1, linecolor= 180, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA_1", 2, linecolor= 160, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA_1", 3, linecolor= 140, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA_1", 4, linecolor= 120, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA_1", 5, linecolor= 100, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA_1", 10, linecolor= 80, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA_1", 20, linecolor= 60, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA_1", 30, linecolor= 40, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA_1", 40, linecolor= 20, xislog=xislog
;process_and_plot_onesimgasdensity, xmax, xmin, bins, "/raid2/tcox/grpA", 50, linecolor= 0, xislog=xislog


; ---------------
; print extras
; -------------
x=[smoothlen2,smoothlen2]
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0

x0= 0.65
y0= 0.82
;xyouts, x0, y0, fload_timelbl(1,2), /normal, charthick= 2.0, size= 1.3, color= 0
;xyouts, x0, y0-0.05, 'black hole', /normal, charthick= 2.0, size= 1.3, color= 150
;xyouts, x0, y0-0.10, 'no black hole', /normal, charthick= 2.0, size= 1.3, color= 50



; Done
;-------
device, /close



end







; process one profile
;----------------------
pro process_and_plot_onesimgasdensity, xmax, xmin, bins, frun, snapnum, $
					linecolor=linecolor, $
					xislog=xislog, $
					includeerr=includeerr					


ok=fload_snapshot_bh(frun, snapnum)


; ------------------------
;  determine gas profile
; ------------------------

x=fload_gas_xyz('x')
y=fload_gas_xyz('y')
z=fload_gas_xyz('z')
m=fload_gas_mass(1)

; use this to weight the measurement
xray=fload_gas_xray_luminosity(1,/diffuse_hotgas)
;xraylum= xray/1.0d+22
xraylum= xray*0.0 + 1.0  ; - should be the same as no weighting


r_log_gas = fltarr(bins)
mass_sd_avg = fltarr(bins)
mass_sd_1sig = fltarr(bins)


tempxmin=xmin
tempxmax=xmax

if xislog eq 0 then begin
  tempxmin=alog10(xmin)
  tempxmax=alog10(xmax)
endif

process_averages_formanyprojections, x, y, z, m, $
				tempxmin, tempxmax, bins, /x_is_log, $
				mass_sd_avg, mass_sd_1sig, r_log_gas, $
				y_weighting= xraylum

if xislog eq 0 then begin
  r_log_gas= 10^(r_log_gas)
endif



mass_sd_avg= mass_sd_avg*1.0e+4
mass_sd_1sig= mass_sd_1sig*1.0e+4

idx= where(mass_sd_avg gt 0)
if idx(0) ne -1 then begin
    mass_sd_avg= mass_sd_avg(idx)
    mass_sd_1sig= mass_sd_1sig(idx)
    r_log_gas= r_log_gas(idx)
endif

idx= where(mass_sd_avg-mass_sd_1sig lt 0)
if idx(0) ne -1 then begin
    mass_sd_1sig(idx)= mass_sd_avg(idx)-0.0000001
endif

oplot, r_log_gas, mass_sd_avg, psym=-3, linestyle=0, color= linecolor, thick=4.0

if keyword_set(includeerr) then begin
   oplot, r_log_gas, mass_sd_avg+mass_sd_1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
   oplot, r_log_gas, mass_sd_avg-mass_sd_1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
endif



end






;===============================================================================
;===============================================================================

; --------------------------
;  The whole thing
; --------------------------
pro plot_whole_thing, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "plot_whole_think, junk"
   print, "  "
   return
endif

filename='gasZT.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=18, newysize=22




; ----------------------------- 
; set up constants      
; -----------------------------


;frun1="/raid4/tcox/vc3vc3e_2"
;frun2="/raid4/tcox/vc3vc3e_no"
;snapnum1= 107
;snapnum2= snapnum1
frun1="/raid4/tcox/vc3vc3h_2"
frun2="/raid4/tcox/vc3vc3h_no"
snapnum1= 107
snapnum2= 30
;frun1="/raid4/tcox/vc3bvc3b"
;frun2="/raid4/tcox/vc3bvc3b_no"
;snapnum= 30
          
bins = 20
          
;xmax = 3
;xmin = -1
xmax = 350
xmin = 1.0

;xaxistitle= 'R (!13h!3!E-1!Nkpc)'
xaxistitle= 'R (kpc)'

;---------------------------
;  Print it
;---------------------------
;
; Produce a multi-figure
; such as below.
;
; --------
; |      |
; |  #1  |
; |      |
; --------
; |      |
; |  #2  |
; |      |
; --------
;
;

x0= 0.18 & x1= 0.99

y0= 0.10 & y1= 0.545 & y2= 0.99



; plot #1
; -------
; Temperature
;
yaxistitle= '!6Temperature (keV)'
ymax = 0.8
ymin = 0.02
generate_axes, x0, y1, x1, y2, xmax, xmin, ymax, ymin, yaxistitle

process_and_plot_onetempprof, xmax, xmin, bins, frun1, snapnum1, linecolor= 50

process_and_plot_onetempprof, xmax, xmin, bins, frun2, snapnum2, linecolor= 150


; want these on the top graph, not the bottom
symsize= 1.5
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=6.0
oplot, [20.0, 28.4], [0.52,0.52], thick=5.0, psym=-8, color= 50, linestyle=2
;xyouts, 0.67, 0.93, 'black hole (RS)', /normal, charthick=3.0, size=1.33, color=50
xyouts, 0.67, 0.93, '!6BH (RS)', /normal, charthick=5.0, size=1.8, color=50


symsize= 1.5
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, [18.0, 24.4], [0.405,0.405], thick=3.0, psym=-2, color= 150
oplot, [20.0, 28.4], [0.40,0.40], thick=5.0, psym=-8, color= 150
;xyouts, 0.67, 0.90, 'no black hole', /normal, charthick=3.0, size=1.33, color=150
xyouts, 0.67, 0.90, '!6std (RS)', /normal, charthick=5.0, size=1.8, color=150


; plot #2
; -------
; gas mass
yaxistitle= '!6Metallicity (Z!D!9n!3!N)'
ymax = 2.0
ymin = 0.02
generate_axes, x0, y0, x1, y1, xmax, xmin, ymax, ymin, yaxistitle, xaxistitle, /yexplabel

process_and_plot_onegasmet, bins, xmax, xmin, frun1, snapnum1, linecolor= 50

process_and_plot_onegasmet, bins, xmax, xmin, frun2, snapnum2, linecolor= 150



; -----------------
;  Plot Extras
; -----------------

;xyouts, x0+0.03, y2-0.04, 'black hole', /normal, charthick=3.5, size= 1.4, color= 150
;xyouts, x0+0.03, y2-0.07, 'no black hole', /normal, charthick=3.5, size= 1.4, color= 50

; want these on the top graph, not the bottom
;symsize= 1.0
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
;oplot, [18.0, 24.4], [0.55,0.55], thick=3.0, psym=-8, color= 50, linestyle=2
;xyouts, 0.67, 0.93, 'black hole (RS)', /normal, charthick=3.0, size=1.33, color=50

;oplot, [18.0, 24.4], [0.40,0.40], thick=3.0, psym=-2, color= 150
;xyouts, 0.67, 0.90, 'no black hole', /normal, charthick=3.0, size=1.33, color=150


;--------------------------------------
;--------------------------------------

device, /close



end





; ------------------------------------------------------------------




; ----------------------
; generate axes and such
; ----------------------
pro generate_axes, x0, y0, x1, y1, $
                xmax, xmin, ymax, ymin, $
                yaxistitle, xaxistitle, $
                linear=linear, $
                yexplabel=yexplabel

!p.position= [x0, y0, x1, y1]

if keyword_set(linear) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, $
        xcharsize=1.50, ycharsize=1.50, xthick=5.0, ythick=5.0, charthick=6.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
        return
endif

if keyword_set(yexplabel) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /xlog, /ylog, $
        xcharsize=2.10, ycharsize=2.10, xthick=5.0, ythick=5.0, charthick=6.0, $
        ;xtickformat='(a1)', ytickformat='exp_label', ytitle=yaxistitle, /nodata, /noerase
        xtitle=xaxistitle, ytickformat='exp_label', ytitle=yaxistitle, /nodata, /noerase
        return
endif

if keyword_set(xaxistitle) then begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, /xlog, $
        xcharsize=2.10, ycharsize=2.10, xthick=5.0, ythick=5.0, charthick=6.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase
endif else begin
        plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, /xlog, $
        xcharsize=2.10, ycharsize=2.10, xthick=5.0, ythick=5.0, charthick=6.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
endelse


end










;===============================================================================
;===============================================================================



; --------------------------
;  Gas Temperature Profile
; --------------------------
pro gastemperature_profile, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "gastemperature_profile, junk"
   print, "  "
   return
endif

filename="gas_t_prof.eps"

setup_plot_stuff, 'ps', filename=filename



; -----------------------------
; set up constants
; -----------------------------

bins = 50

;xislog= 1
xislog= 0      ; manually set /xlog

if xislog eq 1 then begin
    xmax = alog10(300)
    xmin = 0.0
    xaxistitle= 'Log R (kpc)'
endif else begin
    xmax = 300
    xmin = 1.0
    xaxistitle= 'R (h!E-1!Nkpc)'
endelse

;ymax = 1e+9
;ymin = 1e+4
ymax = 8.0
ymin = 0.02

;yaxistitle= 'Log Temperature (K)'
yaxistitle= 'Temperature (keV)'

smoothlen2= alog10(smoothlen)


;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
	/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
	;ytickformat='exp_label', $
        /nodata


;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3f", 25, linecolor= 150, xislog=xislog
;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3h", 25, linecolor= 150, xislog=xislog
;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3j", 25, linecolor= 150, xislog=xislog

;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3f_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3h_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3j_no", 25, linecolor= 50, xislog=xislog

;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc4vc4a", 25, linecolor= 150, xislog=xislog
;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc5vc5a", 25, linecolor= 150, xislog=xislog
;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc6vc6a", 25, linecolor= 150, xislog=xislog

;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc4vc4_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc5vc5_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc6vc6_no", 25, linecolor= 50, xislog=xislog

process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3e", 22, linecolor= 150, xislog=xislog
process_and_plot_onetempprof, xmax, xmin, bins, "/raid2/tcox/vc3HGe", 22, linecolor= 50, xislog=xislog
;process_and_plot_onetempprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3h", 22, linecolor= 150, xislog=xislog
;process_and_plot_onetempprof, xmax, xmin, bins, "/raid2/tcox/vc3HGh", 22, linecolor= 50, xislog=xislog



; print extras
; -------------
x=[smoothlen2,smoothlen2]
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0

x0= 0.65
y0= 0.82
;xyouts, x0, y0, fload_timelbl(1,2), /normal, charthick= 2.0, size= 1.3, color= 0
;xyouts, x0, y0-0.05, 'black hole', /normal, charthick= 2.0, size= 1.3, color= 150
;xyouts, x0, y0-0.10, 'no black hole', /normal, charthick= 2.0, size= 1.3, color= 50
xyouts, x0, y0-0.05, 'std', /normal, charthick= 2.0, size= 1.3, color= 150
xyouts, x0, y0-0.10, 'with halo gas', /normal, charthick= 2.0, size= 1.3, color= 50


; Done
;-------
device, /close

end







; procedure to process single snapshots
; -------------------------------------------

pro process_and_plot_onetempprof, xmax, xmin, bins, frun, snapnum, $
				includeerr=incluederr, $
				linecolor=linecolor, $
				xislog=xislog

if not keyword_set(xislog) then xislog= 0

ok=fload_snapshot_bh(frun,snapnum)

x=fload_gas_xyz('x')
y=fload_gas_xyz('y')
z=fload_gas_xyz('z')
x= x/0.7
y= y/0.7
z= z/0.7
temp=fload_gas_temperature(1,/keV)
;logtemp= alog10(temp)
logtemp= temp

; use this to weight the metallicity measurement
xray=fload_gas_xray_luminosity(1,/diffuse_hotgas)
xraylum= xray/1.0d+22
;xraylum= xray*0.0 + 1.0  ; - should be the same as no weighting



r_log = fltarr(bins)
temp_avg = fltarr(bins)
temp_1sig = fltarr(bins)

if xislog eq 0 then begin
    xmin_temp= alog10(xmin)
    xmax_temp= alog10(xmax)
endif


process_averages_formanyprojections, x, y, z, logtemp, $
				xmin_temp, xmax_temp, bins, /x_is_log, /average, $
				temp_avg, temp_1sig, r_log, $
				y_weighting= xraylum

if xislog eq 0 then begin
  r_log= 10^(r_log)
endif


idx= where(temp_avg gt 0)
if idx(0) ne -1 then begin
	temp_avg= temp_avg(idx)
	temp_1sig= temp_1sig(idx)
	r_log= r_log(idx)
endif

; added some stuff to match the other plots
; -------------------------------------------
;oplot, r_log, temp_avg, psym=-3, linestyle=0, color= linecolor, thick=4.0
if linecolor eq 150 then begin
	symsize= 1.5
	usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	oplot, r_log, temp_avg, psym=-8, linestyle=0, color= linecolor, thick=5.0
	;oplot, r_log, temp_avg, psym=-2, linestyle=0, color= linecolor, thick=3.0
endif

if linecolor eq 50 then begin
	symsize= 1.5
	usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=3.0
	;oplot, r_log, xray_avg, psym=-3, linestyle=0, color= linecolor, thick=4.0
	oplot, r_log, temp_avg, psym=-8, linestyle=2, color= linecolor, thick=5.0
endif


if keyword_set(includeerr) then begin
   oplot, r_log, temp_avg+temp_1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
   oplot, r_log, temp_avg-temp_1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
endif



end





;=================================================================================





;======================================================================================
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;======================================================================================


; --------------------------
; xxxxxxxxxxxxxxxxxxxxxxxxxx
; --------------------------
;
;  Gas Metallicity Profile
;
; --------------------------
; xxxxxxxxxxxxxxxxxxxxxxxxxx
; --------------------------

pro plot_gasmetallicity_profile_comparison, junk, filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_gasmetallicity_profile"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='gas_z_mets.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 50

;xislog= 1
xislog= 0      ; manually set /xlog

if xislog eq 1 then begin
    xmax = alog10(300)
    xmin = 0.0
    xaxistitle= 'Log R (kpc)'
endif else begin
    xmax = 300
    xmin = 1.0
    xaxistitle= 'R (h!E-1!Nkpc)'
endelse

;ymax = 0.5
;ymin = -2.0
ymax = 3.0
ymin = 0.005

;yaxistitle= 'Log Metallicity (Z!D!9n!3!N)'
yaxistitle= 'Metallicity (Z!D!9n!3!N)'



;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
	/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=2.0, $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
	ytickformat='exp_label', $
        /nodata




; ------------------------
;  determine gas profile
; ------------------------

;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc3vc3f", 25, linecolor= 150, xislog=xislog
;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc3vc3h", 25, linecolor= 150, xislog=xislog
;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc3vc3j", 25, linecolor= 150, xislog=xislog

;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc3vc3f_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc3vc3h_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc3vc3j_no", 25, linecolor= 50, xislog=xislog
	
;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc4vc4a", 25, linecolor= 150, xislog=xislog
;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc5vc5a", 25, linecolor= 150, xislog=xislog
;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc6vc6a", 25, linecolor= 150, xislog=xislog

;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc4vc4_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc5vc5_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc6vc6_no", 25, linecolor= 50, xislog=xislog

;process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc3vc3e", 22, linecolor= 150, xislog=xislog
;process_and_plot_onegasmet, bins, xmax, xmin, "/raid2/tcox/vc3HGe", 22, linecolor= 50, xislog=xislog
process_and_plot_onegasmet, bins, xmax, xmin, "/raid4/tcox/vc3vc3h", 22, linecolor= 150, xislog=xislog
process_and_plot_onegasmet, bins, xmax, xmin, "/raid2/tcox/vc3HGh", 22, linecolor= 50, xislog=xislog



; print extras
; -------------

x0= 0.25
y0= 0.32
;xyouts, x0, y0, fload_timelbl(1,2), /normal, charthick= 2.0, size= 1.3, color= 0
;xyouts, x0, y0-0.05, 'black hole', /normal, charthick= 2.0, size= 1.3, color= 150
;xyouts, x0, y0-0.10, 'no black hole', /normal, charthick= 2.0, size= 1.3, color= 50
xyouts, x0, y0-0.05, 'std', /normal, charthick= 2.0, size= 1.3, color= 150
xyouts, x0, y0-0.10, 'with halo gas', /normal, charthick= 2.0, size= 1.3, color= 50


; Done
;-------

device, /close



end






; process one profile at a time
; -------------------------------
pro process_and_plot_onegasmet, bins, xmax, xmin, frun, snapnum, $
				linecolor=linecolor, $
				linepsym=linepsym, $
				includeerr=lincludeerr, $
				xislog=xislog

if not keyword_set(xislog) then xislog= 0

ok=fload_snapshot_bh(frun,snapnum)

x=fload_gas_xyz('x')
y=fload_gas_xyz('y')
z=fload_gas_xyz('z')
x= x/0.7
y= y/0.7
z= z/0.7
zmets=fload_gas_metallicity(1)
zmets=zmets/0.02
;logzmets= alog10(zmets)
logzmets= zmets

; use this to weight the metallicity measurement
xray=fload_gas_xray_luminosity(1,/diffuse_hotgas)
xraylum= xray/1.0d+22
;xraylum= xray*0.0 + 1.0  ; - should be the same as no weighting


r_log = fltarr(bins)
z_avg = fltarr(bins)
z_1sig = fltarr(bins)

if xislog eq 0 then begin
  xmax_temp=alog10(xmax)
  xmin_temp=alog10(xmin)
endif

process_averages_formanyprojections, x, y, z, logzmets, $
				xmin_temp, xmax_temp, bins, /x_is_log, /average, $
				z_avg, z_1sig, r_log, $
				y_weighting= xraylum

if xislog eq 0 then begin
  r_log= 10^(r_log)
endif

idx= where(z_avg gt 0)
if idx(0) ne -1 then begin
	z_avg= z_avg(idx)
	z_1sig= z_1sig(idx)
	r_log= r_log(idx)
endif

if not keyword_set(linepsym) then linepsym= 3

; added some stuff to match the other plots
; -------------------------------------------
;oplot, r_log, z_avg, psym=-linepsym, linestyle=0, color= linecolor, thick=4.0
if linecolor eq 150 then begin
	symsize= 1.5
	usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	oplot, r_log, z_avg, psym=-8, linestyle=0, color= linecolor, thick=5.0
	;oplot, r_log, z_avg, psym=-2, linestyle=0, color= linecolor, thick=3.0
endif

if linecolor eq 50 then begin
	symsize= 1.5
	usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=3.0
	;oplot, r_log, xray_avg, psym=-3, linestyle=0, color= linecolor, thick=4.0
	oplot, r_log, z_avg, psym=-8, linestyle=2, color= linecolor, thick=5.0
endif


if keyword_set(includeerr) then begin
   oplot, r_log, z_avg+z_1sig, psym=-linepsym, linestyle=1, color= linecolor, thick=2.0
   oplot, r_log, z_avg-z_1sig, psym=-linepsym, linestyle=1, color= linecolor, thick=2.0
endif

end







;================================================================================


; ------------------------------------------
; xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; ------------------------------------------
;
;
;  X-ray Luminosity Profile
;
;
; ------------------------------------------
; xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; ------------------------------------------
pro lgx, junk, filename=filename, smoothlen=smoothlen


if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: lgx"
   print, "  "
   return
endif

frun='/raid4/tcox/localgroup/v5'
snapnum=105

;if not keyword_set(filename) then filename=frun+'/profx.eps'
if not keyword_set(filename) then filename=frun+'/profx2.eps'
if not keyword_set(smoothlen) then smoothlen=0.1

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------

bins = 20

;xislog= 1
xislog= 0      ; manually set /xlog

if xislog eq 1 then begin
    xmax = alog10(300)
    xmin = 0.0
    xaxistitle= 'Log R (kpc)'
endif else begin
    xmax = 500
    xmin = 0.5
    ;xaxistitle= 'R (!13h!3!E-1!Nkpc)'
    xaxistitle= '!6R (kpc)'
endelse

ymax = 39.8
ymin = 34.0

;yaxistitle= 'Log mu!DX!N (ergs s!E-1!N kpc!E-2!N)'
yaxistitle= '!6Log Surface Brightness (ergs s!E-1!N kpc!E-2!N)'



;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.98, 0.98]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
	/xlog, $
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





process_and_plot_onexrayprof, bins, xmax, xmin, frun, snapnum, linecolor= 150, xislog=xislog


; print extras
; -------------

x0= 0.65
y0= 0.82
;xyouts, x0, y0, fload_timelbl(1,2), /normal, charthick= 2.0, size= 1.3, color= 0
;xyouts, x0, y0-0.05, 'black hole', /normal, charthick= 2.0, size= 1.3, color= 150
;xyouts, x0, y0-0.10, 'no black hole', /normal, charthick= 2.0, size= 1.3, color= 50
;xyouts, x0, y0-0.05, 'std', /normal, charthick= 2.0, size= 1.3, color= 150
;xyouts, x0, y0-0.10, 'with halo gas', /normal, charthick= 2.0, size= 1.3, color= 50


;symsize= 1.2
;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
;oplot, [22.0, 32.4], [38.21,38.21], thick=3.0, psym=-8, color= 50, linestyle=2
;xyouts, 0.65, 0.87, 'black hole (RS)', /normal, charthick=3.0, size=1.33, color=50
;xyouts, 0.69, 0.87, 'BH (RS)', /normal, charthick=3.0, size=1.33, color=50

;symsize= 1.0
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, [22.0, 32.4], [37.75,37.75], thick=3.0, psym=-8, color= 150
;xyouts, 0.65, 0.82, 'no black hole', /normal, charthick=3.0, size=1.33, color=150
;xyouts, 0.69, 0.82, 'std (RS)', /normal, charthick=3.0, size=1.33, color=150


xyouts, 0.8, 0.9, '(e)', /normal, size=1.8, color= 0, charthick=3.0


; Done
;-------
device, /close

end







;
;  determine gas profile   - for one galaxy
; -------------------------------------------
pro process_and_plot_onexrayprof, bins, xmax, xmin, frun, snapnum, $
				linecolor=linecolor, $
				linepsym=linepsym, $
				includeerr=includeerr, $
				xislog=xislog, $
				fitbeta=fitbeta

ok=fload_snapshot_bh(frun,snapnum)


x=fload_gas_xyz('x')
y=fload_gas_xyz('y')
z=fload_gas_xyz('z')

x= x/0.7
y= y/0.7
z= z/0.7

;xray=fload_gas_xray_luminosity(1)
xray=fload_gas_xray_luminosity(1,/diffuse_hotgas)
print, "Brem total=", total(xray)

; -----------------------------------------------
; Alternate methods: Raymond-Smith
;
; need to .run time_hotgas for this
;
;do_rs_xray= 1
do_rs_xray= 0
if do_rs_xray eq 1 then begin
   load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum
   xray= soft_xray_lum

   print, "RS total=",total(xray)

   ; if we use RS, we need to make the following cut
   rho= fload_gas_rho(1)
   temp= float(fload_gas_temperature(1))
   idx=where((rho lt 0.000854924) and (temp ge 1.0e+5))
   if idx(0) ne -1 then begin
	x= x(idx)
	y= y(idx)
	z= z(idx)
	print, n_elements(idx), " out of ", n_elements(rho), " used "
	print, "n_xray= ", n_elements(xray)
   endif
endif
; -----------------------------------------------

; for averaging (don't know why)
;idx= where(xray le 0.0)
;xray(idx)= 1.0
;xray= alog10(xray)
;xray(idx)= 0.0

; for surface density - smaller numbers seem to work better
xray= xray/1.0d+22


r_log = fltarr(bins)
xray_avg = fltarr(bins)
xray_1sig = fltarr(bins)

if xislog eq 0 then begin
  xmax_temp=alog10(xmax)
  xmin_temp=alog10(xmin)
endif

process_averages_formanyprojections, x, y, z, xray, $
				;xmin, xmax, bins, /x_is_log, /average, $
				xmin_temp, xmax_temp, bins, /x_is_log, $
				xray_avg, xray_1sig, r_log

if xislog eq 0 then begin
  r_log= 10^(r_log)
endif


idx= where(xray_avg gt 0)
if idx(0) ne -1 then begin
	xray_avg= xray_avg(idx)
	xray_1sig= xray_1sig(idx)
	r_log= r_log(idx)
endif

; back to ergs/s/kpc^2
xray_avg= xray_avg*1.0d+22
xray_1sig= xray_1sig*1.0d+22

; fit a beta profile
if keyword_set(fitbeta) then begin
	print, "Fit Beta Profile"
	;weight=1.0/xray_1sig
	weight=xray_1sig
	;weight(*)= 1.0
	if linecolor eq 150 then r_cutoff= 15
	if linecolor eq 50 then r_cutoff= 50
	if r_cutoff gt 0 then begin
	    idx=where(r_log le r_cutoff)
	    r_log_tofit= r_log(idx)
	    xray_avg_tofit= xray_avg(idx)
	    weight_tofit= weight(idx)
	    print, "Cutting this off at r_cutoff= ", r_cutoff
	endif else begin
	    r_log_tofit= r_log
	    xray_avg_tofit= xray_avg
	    weight_tofit= weight
	endelse

	fitandoverplot_betaprofile, r_log_tofit, xray_avg_tofit, weight_tofit, /ylogaxis
endif

; y is log axis
xray_avg_p1sig= alog10(xray_avg+xray_1sig)
xray_avg_m1sig= alog10(xray_avg-xray_1sig)
xray_avg= alog10(xray_avg)
print, min(xray_avg), max(xray_avg)
print, '0=',xray_avg[0],r_log[0]
print, '1=',xray_avg[1],r_log[1]


if linecolor eq 150 then begin
	symsize= 1.0
	usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
	;oplot, r_log, xray_avg, psym=-3, linestyle=0, color= linecolor, thick=4.0
	;oplot, r_log, xray_avg, psym=-2, linestyle=0, color= linecolor, thick=3.0
	oplot, r_log, xray_avg, psym=-8, linestyle=0, color= linecolor, thick=3.0
endif

if linecolor eq 50 then begin
	symsize= 1.2
	usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=3.0
	;oplot, r_log, xray_avg, psym=-3, linestyle=0, color= linecolor, thick=4.0
	oplot, r_log, xray_avg, psym=-8, linestyle=2, color= linecolor, thick=3.0
endif

if (linecolor ne 150) and (linecolor ne 50) then begin
	oplot, r_log, xray_avg, psym=-3, linestyle=1, color= linecolor, thick=3.0
endif

if keyword_set(includeerr) then begin
   oplot, r_log, xray_avg_p1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
   oplot, r_log, xray_avg_m1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
endif



end








;=============================================================================











;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Surface Density
;     ---------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------


; --------------------------
; this is my master program which takes in
; x, y and z coordinates, plus some quantity
; and randomly projects this in many 
; different directions, averaging the results.
; it returns the average and 1 sigma 
; dispersion.
; --------------------------

; result_avg, result_1sig, xaxisvals - all need to be predefined as fltarr(bins)
;
; note that the default is to use the surface density routine,
; and if we set the /average flag, it will average instead
;
;

pro process_averages_formanyprojections, x, y, z, quantitytoaverage, $
					quantitytoaverage_y=quantitytoaverage_y, $
					quantitytoaverage_z=quantitytoaverage_z, $
					center=center, $
					xmin, xmax, bins, $
					result_avg, result_1sig, xaxisvals, $
					x_is_log=x_is_log, x_is_devac=x_is_devac, $
					average=average, $
					y_weighting= y_weighting


;n_projections= 100    ; maybe do 1000 at some point
;n_projections= 40
n_projections= 10
seed= 154L

; 2D array which stores projected quantity for
; each projection
temp_average= fltarr(n_projections+1,bins)

if keyword_set(center) then c=center else c=[0,0,0]

original_quantitytoaverage= quantitytoaverage
if keyword_set(quantitytoaverage_y) and keyword_set(quantitytoaverage_z) then begin
	qx= original_quantitytoaverage
	qy= quantitytoaverage_y
	qz= quantitytoaverage_z
endif


; ---------------------------
;  begin projection loop
; ---------------------------
for i=0,n_projections do begin

        rdphi= randomu(seed)
        rdtheta= randomu(seed)

        ; in radians
	theta= rdtheta*!PI
	phi= rdphi*2*!PI


        ; rotate
        rot_x= (x-c[0])*(cos(theta)*cos(phi)) + (y-c[1])*(cos(theta)*sin(phi)) + (z-c[2])*sin(theta)
        rot_y= -(x-c[0])*sin(theta) + (y-c[1])*cos(theta)
        ;rot_z= x*(sin(theta)*cos(phi)) - y*(sin(theta)*sin(phi)) + z*cos(theta)
        ;rot_vx= vx*(cos(theta)*cos(phi)) + vy*(cos(theta)*sin(phi)) + vz*sin(theta)
        ;rot_vy= -vx*sin(theta) + vy*cos(theta)
        ;rot_vz= vx*(sin(theta)*cos(phi)) - vy*(sin(theta)*sin(phi)) + vz*cos(theta)

	rot_r= sqrt(rot_x*rot_x + rot_y*rot_y)

	if keyword_set(quantitytoaverage_y) and keyword_set(quantitytoaverage_z) then begin
	    quantitytoaverage= qx*(sin(theta)*cos(phi)) - qy*(sin(theta)*sin(phi)) + qz*cos(theta)
	endif

	if keyword_set(x_is_log) then rot_r=alog10(rot_r)
	if keyword_set(x_is_devac) then rot_r=rot_r^(0.25)

	; calculate the profile for this projection
	rad= fltarr(bins)
	thisquant= fltarr(bins)

	if not keyword_set(average) then begin
		process_sd_profile, rot_r, quantitytoaverage, bins, xmax, xmin, rad, thisquant, $
				x_is_log=x_is_log, x_is_devac=x_is_devac, $
				y_weighting=y_weighting
	endif else begin
		process_avg_profile, rot_r, quantitytoaverage, bins, xmax, xmin, rad, thisquant, $
								y_weighting=y_weighting
	endelse

	temp_average[i,*]= thisquant

	xaxisvals= rad
	
        if (i mod 100) eq 0 then print, "i= ",i

endfor

for i=0,bins-1 do begin
	result_moment= moment(temp_average[*,i])
	result_avg[i]= result_moment[0]
	result_1sig[i]= sqrt(result_moment[1])
endfor


end


; --------------------------
;  compute a radial profile
;
;  Note that this procedure is
;  completely general, all it
;  needs is x_var, y_var, the
;  number of bins along with
;  the max and min of the x range
;  and it returns the cumulative
;  x and y
; --------------------------
pro process_sd_profile, x_var, y_var, bins, xmax, xmin,  r_nlog, mass_sd_nlog, $
					x_is_log=x_is_log, x_is_devac=x_is_devac, $
					y_weighting=y_weighting

        ;r_l= alog10(radius)    ; actually surface density is linear
	r_l= x_var

        binsize = float((xmax-xmin))/bins

        for i=1,bins do begin
           lg_r = i*binsize + xmin
           sm_r = (i-1)*binsize + xmin
           r_nlog(i-1) = 0.5*(lg_r + sm_r)

           idx= where((r_l GE sm_r) AND (r_l LT lg_r))
	   weight= 1.0
           if idx(0) lt 0 then begin 
		m= [0.0]
	   endif else begin
		m= y_var(idx)
		if keyword_set(y_weighting) then begin
		   m= y_var(idx)*y_weighting(idx)
		   weight= total(y_weighting(idx))/n_elements(y_weighting(idx))
		endif
	   endelse
           ;thisarea= !PI*((10^(lg_r))^2 - (10^(sm_r))^2)       ; (kpc/h)^2
	   thisarea= !PI*((lg_r)^2 - (sm_r)^2)       ; (kpc/h)^2
	   if keyword_set(x_is_log) then thisarea= !PI*((10^(lg_r))^2 - (10^(sm_r))^2)
	   if keyword_set(x_is_devac) then thisarea= !PI*((lg_r)^8 - (sm_r)^8)
           ;mass_sd_nlog(i-1)= total(m)*1e10/thisarea                     ; h Msolar kpc-2
	   mass_sd_nlog(i-1)= total(m)/(thisarea*weight)          ; Gadget units - h 10^10 Msolar kpc-2
        endfor

end



; -------------------------------
;  compute average profile
;
;  This is a general procedure to
; average the y variables along the
; x range.
; -------------------------------
pro process_avg_profile, x_var, y_var, bins, xmax, xmin,  r_nlog, mass_avg_nlog, $
					y_weighting=y_weighting

	r_l= x_var

        binsize = float((xmax-xmin))/bins

        for i=1,bins do begin
           lg_r = i*binsize + xmin
           sm_r = (i-1)*binsize + xmin
           r_nlog(i-1) = 0.5*(lg_r + sm_r)

           idx= where((r_l GE sm_r) AND (r_l LT lg_r))
           if idx(0) lt 0 then begin
		m= [0.0]
		weight= 1.0
	   endif else begin
		m= y_var(idx)
		weight= n_elements(idx)
		if keyword_set(y_weighting) then begin
		   m= y_var(idx)*y_weighting(idx)
		   weight= total(y_weighting(idx))
		endif
	   endelse
	   mass_avg_nlog(i-1)= total(m)/weight
        endfor

end






;================================================================================







; ==================================
;   Surface Density (averaged)
; ==================================
pro plot_surface_density, frun, sendto, filename=filename, $
				convert_to_cm2=convert_to_cm2, $
				smoothlen=smoothlen

if not keyword_set(sendto) then sendto='ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_surface_density"
   print, "  "
   return
endif


setup_plot_stuff, sendto, filename=filename

; -----------------------------
; set up constants
; -----------------------------

bins = 100

; x is radius in kpc/h
;xmax = 20.0
;xmax = 40.0
xmax = 100.0
xmin = 0.05

; y is surface density in 10^10 Msolar / kpc^2 = 10^4 Msolar/pc^2
;ymax = 1e+4
;ymin = 1e-1
ymax = 6.0
ymin = -1.0


;xtit="R (h!E-1!Nkpc)"
xtit="R (kpc)"
;ytit="Log !4R!3(r) !N(h!E1!N M!D!9n!3!Npc!E-2!N)"
ytit="Log !4R!3(R) !N(M!D!9n!3!Npc!E-2!N)"


; convert to nH/cm2 rather than msolar/pc2
; -----------------------------------------

;SurfDen_Unit_Conversion= 20.290480       ; convert to nH/cm2 (msolar/pc2 = 1.952d20 nH/cm2)
SurfDen_Unit_Conversion= 20.096910       ; convert to nH/cm2 (msolar/pc2 = 1.25d20 nH/cm2)

if keyword_set(convert_to_cm2) then begin
        ytit="Log !4R!3(R) (cm!E-2!N)"

	ymax = 24.0
	ymin = 19.0
endif


; ------------------
; Plot this up
; ------------------
!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, $
        ;/ylog, $
	/xlog, $
        color= 0, $
        xrange=[xmin,xmax], $
        yrange=[ymin,ymax], $
        xstyle=1, ystyle=1, $
	xcharsize=1.50, ycharsize=1.50, $
	xthick=4.0, ythick=4.0, $
	charthick=3.0, $
	xtitle=xtit, $
	ytitle=ytit, $
        /nodata




        ; total (baryonic) profile
        ; --------------------------
        ;a=fload_baryon_xyz('rxy')
        ;c=fload_baryon_mass(1)

        ;r_xy = fltarr(bins)
        ;mass_sd = fltarr(bins)

        ;process_sd_profile, a, c, bins, xmax, xmin, r_xy, mass_sd

        ;print, "Surface Density: (max/min)   ",max(mass_sd), min(mass_sd)
        ;print, "Radius:          (max/min)   ",max(r_xy), min(r_xy)


	; total (all star) profile
	; --------------------------
	x=fload_allstars_xyz('x')
	y=fload_allstars_xyz('y')
	z=fload_allstars_xyz('z')
        m=fload_allstars_mass(1)


	process_and_plot_surfacedensity, x, y, z, m, xmin, xmax, bins, linecolor=0, /includeerr


	; old stellar disk profile
	; ------------------------
	if fload_npart(2) GT 0 then begin
		x=fload_disk_xyz('x')
		y=fload_disk_xyz('y')
		z=fload_disk_xyz('z')
		m=fload_disk_mass(1)

		process_and_plot_surfacedensity, x, y, z, m, xmin, xmax, bins, linecolor=50
	endif


	; old stellar bulge profile
	; -------------------------
        if fload_npart(3) GT 0 then begin
                x=fload_bulge_xyz('x')
		y=fload_bulge_xyz('y')
		z=fload_bulge_xyz('z')
                m=fload_bulge_mass(1)

                process_and_plot_surfacedensity, x, y, z, m, xmin, xmax, bins, linecolor=100
        endif

        ; new stars profile
        ; -------------------
        if fload_npart(4) GT 0 then begin
                x=fload_newstars_xyz('x')
		y=fload_newstars_xyz('y')
		z=fload_newstars_xyz('z')
                m=fload_newstars_mass(1)

		process_and_plot_surfacedensity, x, y, z, m, xmin, xmax, bins, linecolor=150
        endif


; --------------------
; plot random extras
; --------------------

if keyword_set(smoothlen) then begin
	x=[smoothlen,smoothlen]
	y=[ymin,ymax]
	oplot, x, y, linestyle=1, color= 0
endif

; draw HI cutoff
if keyword_set(convert_to_cm2) then begin
	x= findgen(10)*xmax
	x= [xmin,x]
	y= 19.6 + 0.0*x
	oplot, x, y, psym=-3, linestyle= 1, thick=2.0, color= 0
endif

;xyouts, 0.25, 0.87, fload_fid(1), size=1.7, color= 0, /normal
;xyouts, 0.45, 0.87, 'Sc', size=2.9, color= 0, /normal, charthick= 3.0
;xyouts, 0.45, 0.87, 'Sbc', size=2.9, color= 0, /normal, charthick= 3.0

xyouts, 0.75, 0.90, fload_fid(1), /normal, size= 1.5, color= 0, charthick=3.0
xyouts, 0.75, 0.85, fload_timelbl(1,2), /normal, size= 1.5, color= 0, charthick=3.0

; done
; -----
if (sendto EQ 'ps') then device, /close


end







;----------------------------------------
;  Does most of the hard work
;----------------------------------------
pro process_and_plot_surfacedensity, x, y, z, m, $
			xmin, xmax, bins, $
			linecolor=linecolor, includeerr=includeerr, $
			do_linefit=do_linefit, convert_to_cm2=convert_to_cm2

	r_log = fltarr(bins)
	mass_sd_avg = fltarr(bins)
	mass_sd_1sig = fltarr(bins)

	tempxmin=alog10(xmin)
	tempxmax=alog10(xmax)
	process_averages_formanyprojections, x, y, z, m, $
					tempxmin, tempxmax, bins, /x_is_log, $
					mass_sd_avg, mass_sd_1sig, r_log

	idx= where(mass_sd_avg gt 0)
	if idx(0) ne -1 then begin
		mass_sd_avg= mass_sd_avg(idx)
		mass_sd_1sig= mass_sd_1sig(idx)
		r_log= r_log(idx)
	endif

	mass_sd_p1sig= alog10(mass_sd_avg+mass_sd_1sig) + 4.0
	mass_sd_m1sig= alog10(mass_sd_avg-mass_sd_1sig) + 4.0    ; this can give an error when sig > avg.

	;mass_sd_p1sig= alog10(mass_sd_avg)+alog10(mass_sd_1sig) + 4.0
        ;mass_sd_m1sig= alog10(mass_sd_avg)-alog10(mass_sd_1sig) + 4.0
	mass_sd= alog10(mass_sd_avg) + 4.0     ; change to m_solar pc-2
	
	; we're choosing to do /xlog rather than log the coordinate
	; ourselves
	r_log= 10^(r_log)


	; ------------------
	;  Fit a line
	; ------------------

	; note we are computing a line fit to these values
	; after they have been converted to log base 10, but
	; we will want natural log for Rd's and hence we
	; get this extra factor of log10(e)=0.43429 in here.
	;
	if keyword_set(do_linefit) then begin
		r_to_fit= r_log
		m_to_fit= mass_sd

		rd_disk= linfit(r_to_fit, m_to_fit, chisq=chisq, /double)
		print, "Fitting line to surface density"
		print, "  Sigma_0= ", 10^rd_disk[0]
		print, "       Rd= ",-(0.43429)/rd_disk[1]

		x=(findgen(10) * 1000.0) - 10.0
		y= rd_disk[0] + rd_disk[1]*x
		oplot, x, y, psym=-3, linestyle= 2, thick=2, color= 0
		rdlbl= strcompress(string(-(0.43429)/rd_disk[1]),/remove_all)
		rdlbl= strmid(rdlbl,0,4)
		rdlbl= 'R!Dd!N='+rdlbl+' kpc'
		xyouts, x0, y0, rdlbl, /normal, charthick= 1.7, size= 1.2, color= 0

	endif


	; convert to nH/cm2 rather than msolar/pc2
	; -----------------------------------------
	;SurfDen_Unit_Conversion= 20.290480       ; convert to nH/cm2 (msolar/pc2 = 1.952d20 nH/cm2)
	SurfDen_Unit_Conversion= 20.096910       ; convert to nH/cm2 (msolar/pc2 = 1.25d20 nH/cm2)

	if keyword_set(convert_to_cm2) then begin
		mass_sd= mass_sd + SurfDen_Unit_Conversion
	endif


	; plot the total
	; ---------------
	oplot, r_log, mass_sd, psym=-3, linestyle=0, color= linecolor, thick=4.0

	if keyword_set(includeerr) then begin
		oplot, r_log, mass_sd_p1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
		oplot, r_log, mass_sd_m1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
	endif


end














