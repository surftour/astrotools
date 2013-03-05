;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   Compiled Gas Profiles
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------






; --------------------------
;  Gas Profile - projected
; --------------------------
pro plot_2d_gas_profile, frun, sendto, smoothlen=smoothlen, filename=filename

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
pro plot_2d_gas_profile_comparison, junk, smoothlen=smoothlen, filename=filename

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
    ;xmax = 1000
    ;xmax = 300
    xmax= 100.0
    xmin = 1.0
    ;xaxistitle= '!6R (h!E-1!Nkpc)'
    xaxistitle= '!6R (kpc)'
endelse

ymax = 1e+3
ymin = 1e-3

;yaxistitle= 'Log !4R!3(R) !N(M!D!9n!3!N pc!E-2!N)'
yaxistitle= '!4R!6(R) !N(M!D!9n!6!N pc!E-2!N)'


smoothlen2= alog10(smoothlen)

;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        /ylog, $
	;/xlog, $
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

process_and_plot_2d_onesimgasdensity, xmax, xmin, bins, "data1/Sbc/Sbc_10x_wBH", 0, linecolor= 0, xislog=xislog
process_and_plot_2d_onesimgasdensity, xmax, xmin, bins, "data1/Sbc/Sbccut", 0, linecolor= 50, xislog=xislog
process_and_plot_2d_onesimgasdensity, xmax, xmin, bins, "data1/Sbc/Sbcm1", 0, linecolor= 150, xislog=xislog
process_and_plot_2d_onesimgasdensity, xmax, xmin, bins, "data1/Sbc/Sbcm2", 0, linecolor= 175, xislog=xislog
process_and_plot_2d_onesimgasdensity, xmax, xmin, bins, "data1/Sbc/Sbcm3", 0, linecolor= 200, xislog=xislog


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
pro process_and_plot_2d_onesimgasdensity, xmax, xmin, bins, frun, snapnum, $
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

;process_averages_formanyprojections, x, y, z, m, $   ; no longer called this
process_prof_frommanyprojections, x, y, z, m, $
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




; --------------------------
;  Gas Temperature Profile
; --------------------------
pro plot_2d_gastemperature_profile, frun, sendto, $
				smoothlen=smoothlen, $
				filename=filename

if not keyword_set(sendto) then sendto='ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_gastemperature_profile"
   print, "  "
   return
endif


if fload_npart(0) eq 0 then begin
	print, " "
	print, " plot_gastemperature_profile didn't find any gas particles"
	print, " "
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

;ymax = 1e+9
;ymin = 1e+4
ymax = 2.0
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

process_and_plot_onetempprof, xmax, xmax, bins, frun, snapnum, /includeerr, linecolor=150


; print extras
; -------------
x=[smoothlen2,smoothlen2]
y=[ymin,ymax]
oplot, x, y, linestyle=1, color= 0

xyouts, 0.65, 0.86, fload_fid(1), /normal, charthick= 2.0, size= 1.3, color= 0
xyouts, 0.65, 0.82, fload_timelbl(1,2), /normal, charthick= 2.0, size= 1.3, color= 0


; Done
;-------
if (sendto EQ 'ps') then device, /close



end






; --------------------------
;  Gas Temperature Profile
; --------------------------
pro plot_gastemperature_profile_comparison, junk, $
				smoothlen=smoothlen, $
				filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_gastemperature_profile"
   print, "  "
   return
endif

if not keyword_set(filename) then filename="gas_t_prof.eps"
if not keyword_set(smoothlen) then smoothlen=0.1


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

ok=fload_snapshot_bh(frun,snapnum)

x=fload_gas_xyz('x')
y=fload_gas_xyz('y')
z=fload_gas_xyz('z')
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

oplot, r_log, temp_avg, psym=-3, linestyle=0, color= linecolor, thick=4.0

if keyword_set(includeerr) then begin
   oplot, r_log, temp_avg+temp_1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
   oplot, r_log, temp_avg-temp_1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
endif



end





;=================================================================================


; --------------------------
; xxxxxxxxxxxxxxxxxxxxxxxxxx
; --------------------------
;
;  Gas Entropy Profile
;
; --------------------------
; xxxxxxxxxxxxxxxxxxxxxxxxxx
; --------------------------
pro plot_gasentropy_profile, frun, sendto, filename=filename

if not keyword_set(sendto) then sendto='ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_gasentropy_profile"
   print, "  "
   return
endif


if fload_npart(0) eq 0 then begin
	print, " "
	print, " plot_gasentropy_profile didn't find any gas particles"
	print, " "
endif

initialize_plotinfo, 1
setup_plot_stuff, sendto, filename=filename, colortable=4



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
ymax = 5.0
ymin = 0.0

yaxistitle= 'Log Entropy (keV cm!E2!N)'



;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
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

process_and_plot_oneentprof, xmax, xmin, bins, frun, snapnum, linecolor= 150, xislog=xislog


; print extras
; -------------

xyouts, 0.65, 0.86, fload_fid(1), /normal, charthick= 2.0, size= 1.3, color= 0
xyouts, 0.65, 0.82, fload_timelbl(1,2), /normal, charthick= 2.0, size= 1.3, color= 0



; Done
;-------

if (sendto EQ 'ps') then device, /close



end





; --------------------------
;  Gas Entropy Profile
; --------------------------
pro plot_gasentropy_profile_comparison, junk, $
				smoothlen=smoothlen, $
				filename=filename

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_gasentropy_profile"
   print, "  "
   return
endif

if not keyword_set(filename) then filename="gas_e_prof.eps"
if not keyword_set(smoothlen) then smoothlen=0.1

if fload_npart(0) eq 0 then begin
	print, " "
	print, " plot_gasentropy_profile didn't find any gas particles"
	print, " "
endif

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

bins = 50

ymax = 1e+5
ymin = 1e+1
;ymax = 5.0
;ymin = 0.0

;yaxistitle= 'Log Entropy (keV cm!E2!N)'
yaxistitle= 'Entropy (keV cm!E2!N)'



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

;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3f", 25, linecolor= 150, xislog=xislog
;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3h", 25, linecolor= 150, xislog=xislog
;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3j", 25, linecolor= 150, xislog=xislog

;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3f_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3h_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3j_no", 25, linecolor= 50, xislog=xislog

;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc4vc4a", 25, linecolor= 150, xislog=xislog
;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc5vc5a", 25, linecolor= 150, xislog=xislog
;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc6vc6a", 25, linecolor= 150, xislog=xislog

;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc4vc4_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc5vc5_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc6vc6_no", 25, linecolor= 50, xislog=xislog

process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3e", 22, linecolor= 150, xislog=xislog
process_and_plot_oneentprof, xmax, xmin, bins, "/raid2/tcox/vc3HGe", 22, linecolor= 50, xislog=xislog
;process_and_plot_oneentprof, xmax, xmin, bins, "/raid4/tcox/vc3vc3h", 22, linecolor= 150, xislog=xislog
;process_and_plot_oneentprof, xmax, xmin, bins, "/raid2/tcox/vc3HGh", 22, linecolor= 50, xislog=xislog


; print extras
; -------------

x0= 0.65
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










; procedure to process single snapshots
; -------------------------------------------

pro process_and_plot_oneentprof, xmax, xmin, bins, frun, snapnum, $
				includeerr=incluederr, $
				linecolor=linecolor, $
				xislog=xislog

ok=fload_snapshot_bh(frun,snapnum)


x=fload_gas_xyz('x')
y=fload_gas_xyz('y')
z=fload_gas_xyz('z')
ent=fload_gas_entropy(1)
;logent= alog10(ent)
logent= ent

; use this to weight the metallicity measurement
xray=fload_gas_xray_luminosity(1,/diffuse_hotgas)
xraylum= xray/1.0d+22
;xraylum= xray*0.0 + 1.0  ; - should be the same as no weighting


r_log = fltarr(bins)
ent_avg = fltarr(bins)
ent_1sig = fltarr(bins)


if xislog eq 0 then begin
  xmax_temp=alog10(xmax)
  xmin_temp=alog10(xmin)
endif


process_averages_formanyprojections, x, y, z, logent, $
				xmin_temp, xmax_temp, bins, /x_is_log, /average, $
				ent_avg, ent_1sig, r_log, $
				y_weighting= xraylum


idx= where(ent_avg gt 0)
if idx(0) ne -1 then begin
	ent_avg= ent_avg(idx)
	ent_1sig= ent_1sig(idx)
	r_log= r_log(idx)
endif

if xislog eq 0 then begin
  r_log= 10^(r_log)
endif

print, "ent max/min= ", max(ent_avg), min(ent_avg)

oplot, r_log, ent_avg, psym=-3, linestyle=0, color=linecolor, thick=4.0

if keyword_set(includeerr) then begin
   oplot, r_log, ent_avg+ent_1sig, psym=-3, linestyle=1, color=linecolor, thick=2.0
   oplot, r_log, ent_avg-ent_1sig, psym=-3, linestyle=1, color=linecolor, thick=2.0
endif


end






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

pro plot_gasmetallicity_profile, frun, snapnum, filename=filename

if not keyword_set(frun) or not keyword_set(snapnum) then begin
   print, "  "
   print, "PROBLEM: plot_gasmetallicity_profile"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='gas_z_mets.eps'

if fload_npart(0) eq 0 then begin
	print, " "
	print, " plot_gasmetallicity_profile didn't find any gas particles"
	print, " "
endif

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

;ymax = 0.5
;ymin = -2.0
ymax = 3.0
;ymin = 0.001
ymin = 0.1

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

process_and_plot_onegasmet, bins, xmax, xmin, frun, snapnum, linecolor= 150, /includeerr, xislog=xislog
;process_and_plot_onegasmet, bins, xmax, xmin, frun, snapnum, linecolor= 150, xislog=xislog

	

; print extras
; -------------

xyouts, 0.65, 0.86, fload_fid(1), /normal, charthick= 2.0, size= 1.3, color= 0
xyouts, 0.65, 0.82, fload_timelbl(1,2), /normal, charthick= 2.0, size= 1.3, color= 0



; Done
;-------

device, /close



end











; --------------------------
;  Gas Metallicity Profile
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

ok=fload_snapshot_bh(frun,snapnum)

x=fload_gas_xyz('x')
y=fload_gas_xyz('y')
z=fload_gas_xyz('z')
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

oplot, r_log, z_avg, psym=-linepsym, linestyle=0, color= linecolor, thick=4.0

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

pro plot_xraylum_profile, frun, sendto, filename=filename

if not keyword_set(sendto) then sendto='ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "PROBLEM: plot_xraylum_profile"
   print, "  "
   return
endif


if fload_npart(0) eq 0 then begin
	print, " "
	print, " plot_xraylum_profile didn't find any gas particles"
	print, " "
endif

setup_plot_stuff, sendto, filename=filename



; -----------------------------
; set up constants
; -----------------------------

bins = 100

xmax = 2.8
xmin = -0.2
ymax = 40.0
ymin = 30.0

xaxistitle= 'Log R (kpc)'
yaxistitle= 'Log L!DX!N (ergs s!E-1!N)'



;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
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



process_and_plot_onexrayprof, bins, xmax, xmin, frun, snapnum, linecolor=150, /includeerr


; print extras
; -------------

xyouts, 0.65, 0.86, fload_fid(1), /normal, charthick= 2.0, size= 1.3, color= 0
xyouts, 0.65, 0.82, fload_timelbl(1,2), /normal, charthick= 2.0, size= 1.3, color= 0



; Done
;-------

if (sendto EQ 'ps') then device, /close

end





; ----------------------------
;  X-ray Luminosity Profile
; ----------------------------
pro plot_xraylum_profile_comparison, junk, filename=filename, $
					smoothlen=smoothlen

if not keyword_set(junk) then begin
   print, "  "
   print, "PROBLEM: plot_xraylum_profile"
   print, "  "
   return
endif

if not keyword_set(filename) then filename='gasxprof.eps'
if not keyword_set(smoothlen) then smoothlen=0.1

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

ymax = 40.0
ymin = 30.0

yaxistitle= 'Log mu!DX!N (ergs s!E-1!N kpc!E-2!N)'



;---------------------------
;  Print it
;---------------------------

!p.position= [0.2, 0.15, 0.95, 0.95]

plot, [0], [0], psym=-3, linestyle=0, $
        ;/ylog, $
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




;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc3vd3e", 25, linecolor= 150, /includeerr, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc3vd3e", 30, linecolor= 150, xislog=xislog

;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc3vc3f", 25, linecolor= 150, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc3vc3h", 25, linecolor= 150, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc3vc3j", 25, linecolor= 150, xislog=xislog

;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc3vc3f_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc3vc3h_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc3vc3j_no", 25, linecolor= 50, xislog=xislog

;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc4vc4a", 25, linecolor= 150, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc5vc5a", 25, linecolor= 150, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc6vc6a", 25, linecolor= 150, xislog=xislog

;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc4vc4_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc5vc5_no", 25, linecolor= 50, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc6vc6_no", 25, linecolor= 50, xislog=xislog

process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc3vc3e", 22, linecolor= 150, xislog=xislog
process_and_plot_onexrayprof, bins, xmax, xmin, "/raid2/tcox/vc3HGe", 22, linecolor= 50, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid4/tcox/vc3vc3h", 22, linecolor= 150, xislog=xislog
;process_and_plot_onexrayprof, bins, xmax, xmin, "/raid2/tcox/vc3HGh", 22, linecolor= 50, xislog=xislog

; print extras
; -------------

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







;
;  determine gas profile   - for one galaxy
; -------------------------------------------
pro process_and_plot_onexrayprof, bins, xmax, xmin, frun, snapnum, $
				linecolor=linecolor, $
				linepsym=linepsym, $
				includeerr=lincludeerr, $
				xislog=xislog

ok=fload_snapshot_bh(frun,snapnum)


x=fload_gas_xyz('x')
y=fload_gas_xyz('y')
z=fload_gas_xyz('z')
;xray=fload_gas_xray_luminosity(1)
xray=fload_gas_xray_luminosity(1,/diffuse_hotgas)

; for averaging (don't know why)
;idx= where(xray le 0.0)
;xray(idx)= 1.0
;xray= alog10(xray)
;xray(idx)= 0.0

; for surface density
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

xray_avg= xray_avg*1.0d+22
xray_avg= alog10(xray_avg)
print, min(xray_avg), max(xray_avg)
print, '0=',xray_avg[0],r_log[0]
print, '1=',xray_avg[1],r_log[1]

oplot, r_log, xray_avg, psym=-3, linestyle=0, color= linecolor, thick=4.0

if keyword_set(includeerr) then begin
   oplot, r_log, xray_avg+xray_1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
   oplot, r_log, xray_avg-xray_1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
endif



end








;=============================================================================









