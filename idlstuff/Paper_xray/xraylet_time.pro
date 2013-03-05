;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;   3D Gas Profiles
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------





;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     This is our 3D figure for the *revised* 
;     X-ray letter.
;     -------------------------------------------
;  
;  
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------




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

filename='proptime.eps'
                
initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=16, newysize=20



                
; ----------------------------- 
; set up constants      
; -----------------------------


;frun1="/raid4/tcox/vc3vc3e_2"
;frun2="/raid4/tcox/vc3vc3e_no"
frun1="/raid4/tcox/vc4vc4a"
frun2="/raid4/tcox/vc4vc4_no"
snapnum= 30
                        
bins = 30
                
xmax = 3.0
xmin = 0.0

xaxistitle= 'Time (!13h!e-1!n!3Gyr)'




;---------------------------
;  Print it
;---------------------------
;
; Produce a multi-figure
; such as below.
;
; --------   --------
; |      |   |      |
; |  #1  |   |  #2  |
; |      |   |      |
; --------   --------
; |      |   |      |
; |  #3  |   |  #4  |
; |      |   |      |
; --------   --------
; |      |   |      |
; |  #5  |   |  #6  |
; |      |   |      |
; --------   --------
;
;

x0= 0.12 & x1= 0.495 & x2= 0.615 & x3= 0.99

y0= 0.09 & y1= 0.39 & y2= 0.69 & y3= 0.99



; plot #1
; -------
; L_B
;
yaxistitle= 'L!DB!N (L!DB,!9n!3!N)'
ymax = 1.0e+12
ymin = 7.0e+10
generate_axes, x0, y2, x1, y3, xmax, xmin, ymax, ymin, yaxistitle, /yexplabel

processandplot_onelb, frun1, snapnum, xmin, xmax, bins, linecolor= 150

processandplot_onelb, frun2, snapnum, xmin, xmax, bins, linecolor= 50


; plot #2
; -------
; L_X / L_B 
;
yaxistitle= 'L!DX!N/L!DB!N (erg s!e-1!n L!DB,!9n!3!N!E-1!N)'
ymax = 7.0d+30
ymin = 6.0d+28
generate_axes, x2, y2, x3, y3, xmax, xmin, ymax, ymin, yaxistitle, /yexplabel

processandplot_lxlb, frun1, snapnum, xmin, xmax, bins, linecolor= 150

processandplot_lxlb, frun2, snapnum, xmin, xmax, bins, linecolor= 50




; plot #3
; -------
; gas (and L_X weighted) mass
;
yaxistitle= 'Mass (10!E10!N M!D!9n!3!N)'
ymax = 4.5
ymin = 0.0
generate_axes, x0, y1, x1, y2, xmax, xmin, ymax, ymin, yaxistitle, /linear

processandplot_hotgas, frun1, snapnum, xmin, xmax, bins, linecolor= 150, /mass

processandplot_hotgas, frun2, snapnum, xmin, xmax, bins, linecolor= 50, /mass



; plot #4
; -------
; sigma (km/sec)
;
yaxistitle= '!7r!3 (km sec!E-1!N) '
ymax = 325.0
ymin = 125.0
generate_axes, x2, y1, x3, y2, xmax, xmin, ymax, ymin, yaxistitle, /linear

processandplot_sigma, frun1, snapnum, xmin, xmax, bins, linecolor= 150

processandplot_sigma, frun2, snapnum, xmin, xmax, bins, linecolor= 50




; plot #5
; -------
; temperature (keV)
;
yaxistitle= 'Temperature (keV)'
ymax = 1.5
ymin = 0.08
generate_axes, x0, y0, x1, y1, xmax, xmin, ymax, ymin, yaxistitle, xaxistitle

processandplot_hotgas, frun1, snapnum, xmin, xmax, bins, linecolor= 150, /tempkeV

processandplot_hotgas, frun2, snapnum, xmin, xmax, bins, linecolor= 50, /tempkeV




; plot #6
; -------
; metallicity (solar)
;
yaxistitle= 'Metallicity (Z!D!9n!3!N) '
ymax = 7.0
ymin = 0.002
generate_axes, x2, y0, x3, y1, xmax, xmin, ymax, ymin, yaxistitle, xaxistitle

processandplot_hotgas, frun1, snapnum, xmin, xmax, bins, linecolor= 150, /mets

processandplot_hotgas, frun2, snapnum, xmin, xmax, bins, linecolor= 50, /mets




; -----------------
;  Plot Extras
; -----------------


xyouts, x0+0.03, y3-0.04, 'black hole', /normal, charthick=3.5, size= 1.4, color= 150
xyouts, x0+0.03, y3-0.07, 'no black hole', /normal, charthick=3.5, size= 1.4, color= 50


;--------------------------------------
;--------------------------------------

device, /close



end









; ------------------------------------------------------------------






; --------------------------
;  The gas phases
; --------------------------
pro plot_gasphases, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "plot_gasphases, junk"
   print, "  "
   return                       
endif                           

filename='gasproptime.eps'
                
initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=13, newysize=12
setup_plot_stuff, 'ps', filename=filename, colortable= 4



                
; ----------------------------- 
; set up constants      
; -----------------------------


;frun="/raid4/tcox/vc3vc3e_2"
;frun="/raid4/tcox/vc3vc3e_no"
;frun="/raid4/tcox/vc3vc3h_2"
frun="/raid4/tcox/vc3vc3h_no"
;frun="/raid4/tcox/vc4vc4a"
;frun="/raid4/tcox/vc4vc4_no"
snapnum= 30
                        
bins = 30
                
xmax = 4.25
xmin = 0.0

;xaxistitle= 'Time (!13h!e-1!n!3Gyr)'
xaxistitle= '!6Time (Gyr)'




;---------------------------
;  Print it
;---------------------------
;

x0= 0.18 & x1= 0.98
y0= 0.15 & y1= 0.98




; plot #3
; -------
; gas (and L_X weighted) mass
;
;yaxistitle= '!6Mass (10!E10!N M!D!9n!3!N)'
yaxistitle= '!6Mass (M!D!9n!3!N)'
;ymax = 3.25
;ymin = 0.0
;generate_axes, x0, y0, x1, y1, xmax, xmin, ymax, ymin, yaxistitle, xaxistitle, /linear
ymax = 5.25e+10
ymin = 0.001e+10
generate_axes, x0, y0, x1, y1, xmax, xmin, ymax, ymin, yaxistitle, xaxistitle, /yexplabel

processandplot_hotgas, frun, snapnum, xmin, xmax, bins, linecolor= 0, /totmass
;processandplot_hotgas, frun, snapnum, xmin, xmax, bins, linecolor= 50, /cmass
;processandplot_hotgas, frun, snapnum, xmin, xmax, bins, linecolor= 50, /hmass
processandplot_hotgas, frun, snapnum, xmin, xmax, bins, linecolor= 50, /xmass
processandplot_hotgas, frun, snapnum, xmin, xmax, bins, linecolor= 150, /wmass




; -----------------
;  Plot Extras
; -----------------


xyouts, 0.60, 0.90, '!6black hole', /normal, charthick=5.5, size= 1.6, color= 0
;xyouts, 0.60, 0.90, '!6no black hole', /normal, charthick=5.5, size= 1.6, color= 0

xyouts, 0.63, 0.82, 'Total', /normal, charthick=3.5, size= 1.4, color= 0
;xyouts, 0.63, 0.80, 'Cold', /normal, charthick=3.5, size= 1.4, color= 50
;xyouts, 0.63, 0.75, 'Hot', /normal, charthick=3.5, size= 1.4, color= 100
xyouts, 0.63, 0.70, 'X-ray', /normal, charthick=3.5, size= 1.4, color= 50
xyouts, 0.63, 0.51, 'Wind', /normal, charthick=3.5, size= 1.4, color= 150


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
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
	return
endif

if keyword_set(yexplabel) then begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        ;xtickformat='(a1)', ytickformat='exp_label', ytitle=yaxistitle, /nodata, /noerase
        xtitle=xaxistitle, ytickformat='exp_label', ytitle=yaxistitle, /nodata, /noerase
	return
endif

if keyword_set(xaxistitle) then begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, $
        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase
endif else begin
	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
	xrange=[xmin,xmax], yrange=[ymin,ymax], /ylog, $
	xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase
endelse


end









;==================================================================================


; --------------------------
;
;  hotgas.txt 
;
; --------------------------
pro processandplot_hotgas, frun, snapnum, xmin, xmax, bins, $
					linecolor=linecolor, $
					totmass=totmass, mets=mets, $
					cmass=cmass, hmass=hmass, $
					wmass=wmass, xmass=xmass, $
					tempkeV=tempkeV

read_wind_file, frun, time, mass_gtx, mass_egt, wind_ke, wind_rke, wind_z
read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

time= time/0.7
gas_tot= gas_tot*1.0e+10/0.7
gas_cold= gas_cold*1.0e+10/0.7
gas_hot= gas_hot*1.0e+10/0.7
mass_xraygas= mass_xraygas*1.0e+10/0.7
mass_egt= mass_egt*1.0e+10/0.7

; mass
;
if keyword_set(totmass) then begin
	oplot, time, gas_tot, psym=-3, linestyle=0, color= linecolor, thick=4.0
endif
if keyword_set(hmass) then begin
	oplot, time, gas_hot, psym=-3, linestyle=1, color= linecolor, thick=4.0
endif
if keyword_set(cmass) then begin
	oplot, time, gas_cold, psym=-3, linestyle=1, color= linecolor, thick=4.0
endif
if keyword_set(wmass) then begin
	oplot, time, mass_egt, psym=-2, linestyle=0, color= linecolor, thick=4.0
endif
if keyword_set(xmass) then begin
	oplot, time, mass_xraygas, psym=-1, linestyle=0, color= linecolor, thick=4.0
endif


; temp
;
if keyword_set(tempkeV) then begin
	oplot, time, temp_keV_X, psym=-3, linestyle=0, color= linecolor, thick=4.0
endif

; mets
;
if keyword_set(mets) then begin
	mets=time
	mets(*)= 0.1
	oplot, time, mets, psym=-3, linestyle=0, color= linecolor, thick=4.0
endif


end





; --------------------------
;
;  sigma.txt
; 
; --------------------------
pro processandplot_sigma, frun, snapnum, xmin, xmax, bins, $
					linecolor=linecolor

read_sigma_file, frun, time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
                                Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, $
                                gsigavg, gsigerr

oplot, time, Asigavg, psym=-3, linestyle=0, color= linecolor, thick=4.0

end




; --------------------------
;
;  Read colors.txt
; 
; --------------------------
pro processandplot_onelb, frun, snapnum, xmin, xmax, bins, $
                                        linecolor=linecolor

read_colors_file, frun, time, mags

;bolo= mags[0,*]
Bmag= mags[2,*]
;Kmag= mags[8,*]
Lum_B= 10^(-0.4*(Bmag-5.51))
;print, Lum_B

oplot, time, Lum_B, psym=-3, linestyle=0, color= linecolor, thick=4.0

end




pro processandplot_lxlb, frun, snapnum, xmin, xmax, bins, $
                                        linecolor=linecolor

read_colors_file, frun, time, mags
read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h


;bolo= mags[0,*]
Bmag= mags[2,*]
;Kmag= mags[8,*]
Lum_B= -0.4*(Bmag-5.51)

y= 10^(xray_rs_s-Lum_B)

;oplot, time, xray_rs_s-Lum_B, psym=-3, linestyle=0, color= linecolor, thick=4.0
oplot, time, y, psym=-3, linestyle=0, color= linecolor, thick=4.0

end













