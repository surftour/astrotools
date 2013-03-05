;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;       Histograms of Wind Properties
;     ------------------------------------------------
;
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------


;
;   Velocities
; --------------------------------
pro vw_hist, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "vw_hist, junk"
   print, "  "
   return
endif

filename='vwhist.eps'


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius, linear scale
xaxistitle= "!6Velocity (km s!E-1!N)"
xmax = 1000.0     ; also set these in the process_vwhist procedure
xmin = 0.0

; number (histogram)
;yaxistitle= "!6V!Dr!N (km s!E-1!N)"
yaxistitle= ' '
ymax = 1.2
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.05
y0= 0.15
x1= 0.95
y1= 0.98

!p.position= [x0, y0, x1, y1]

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
	ytickformat='(a1)', $
	;ytitle=yaxistitle, $
	/nodata   ;, /noerase


; -----------------------------------------------

;snapnum= 8
;snapnum= 49
;snapnum= 78
;snapnum= 81

;process_vwhist, "/raid4/tcox/vc3vc3e_no", snapnum, oplotit= 50
;xyouts, 0.70, 0.85, "e_no", size=1.5, color=50, /normal, charthick=3.0

;process_vwhist, "/raid4/tcox/vc3vc3e_2", snapnum, oplotit= 0
;xyouts, 0.70, 0.80, "e_2", size=1.5, color=0, /normal, charthick=3.0

;process_vwhist, "/raid4/tcox/sbw/sb10", snapnum, oplotit= 100
;xyouts, 0.70, 0.75, "sb10", size=1.5, color=100, /normal, charthick=3.0

;process_vwhist, "/raid4/tcox/sbw/sb10BH", snapnum, oplotit= 150
;xyouts, 0.70, 0.70, "sb10BH", size=1.5, color=150, /normal, charthick=3.0

;xyouts, 0.70, 0.90, fload_timelbl(1,2,/noteq), /normal, color= 0, charthick=3.0, size=1.5


; -----------------------------------------------


process_vwhist, "/raid4/tcox/sb10_mass/d1e", 5, oplotit= 10
xyouts, 0.70, 0.85, "d1e", size=1.5, color=10, /normal, charthick=3.0

process_vwhist, "/raid4/tcox/sb10_mass/d2e", 5, oplotit= 50
xyouts, 0.70, 0.80, "d2e", size=1.5, color=50, /normal, charthick=3.0

process_vwhist, "/raid4/tcox/sb10_mass/d4e", 5, oplotit= 100
xyouts, 0.70, 0.75, "d4e", size=1.5, color=100, /normal, charthick=3.0

process_vwhist, "/raid4/tcox/sb10_mass/d5e", 5, oplotit= 200
xyouts, 0.70, 0.70, "d5e", size=1.5, color=200, /normal, charthick=3.0

; -----------------------------------------------


; print extras
; -------------

; done
; -----
device, /close


end







; ------------------
; compute histogram
; ------------------
pro process_vwhist, frun, snapnum, oplotit= oplotit


xmax= 1000.0
xmin= 0.0

ok=fload_snapshot_bh(frun,snapnum)

rad= fload_gas_xyz('r')
vtot= fload_gas_v('mag')
vr= fload_gas_v('r')


windvelocities= vtot
print, "wind vel max/min= ", max(windvelocities), min(windvelocities)

; -----------------------------------------------

temp= process_histogram(windvelocities, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
print, min(temp), max(temp)

; -----------------------------------------------


end





;========================================================================================


;
;   Metallicities
; --------------------------------
pro z_hist, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "z_hist, junk"
   print, "  "
   return
endif

filename='zhist.eps'


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius, linear scale
xaxistitle= "!6Log Metallicity (Z!D!9n!6!N)"
xmax = 1.0     ; also set these in the process_vwhist procedure
xmin = -4.0

; number (histogram)
;yaxistitle= "!6V!Dr!N (km s!E-1!N)"
yaxistitle= ' '
ymax = 1.2
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.05
y0= 0.15
x1= 0.95
y1= 0.98

!p.position= [x0, y0, x1, y1]

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
	ytickformat='(a1)', $
	;ytitle=yaxistitle, $
	/nodata   ;, /noerase


; -----------------------------------------------

;snapnum= 8
;snapnum= 49
;snapnum= 78
;snapnum= 81
snapnum= 107

;process_zhist, "/raid4/tcox/vc3vc3e_no", snapnum, oplotit= 50
;xyouts, 0.70, 0.85, "e_no", size=1.5, color=50, /normal, charthick=3.0

process_zhist, "/raid4/tcox/vc3vc3e_2", snapnum, oplotit= 150, xmax=xmax, xmin=xmin
xyouts, 0.10, 0.80, "BH-driven wind", size=1.5, color=150, /normal, charthick=3.0

process_zhist, "/raid4/tcox/sbw/sb10", snapnum, oplotit= 50, xmax=xmax, xmin=xmin
xyouts, 0.10, 0.75, "SB-driven wind", size=1.5, color=50, /normal, charthick=3.0

;process_zhist, "/raid4/tcox/sbw/sb10BH", snapnum, oplotit= 150
;xyouts, 0.70, 0.70, "sb10BH", size=1.5, color=150, /normal, charthick=3.0

;xyouts, 0.70, 0.90, fload_timelbl(1,2,/noteq), /normal, color= 0, charthick=3.0, size=1.5


; -----------------------------------------------



; print extras
; -------------

; done
; -----
device, /close


end







; ------------------
; compute histogram
; ------------------
pro process_zhist, frun, snapnum, oplotit= oplotit, xmax=xmax, xmin=xmin


ok=fload_snapshot_bh(frun,snapnum)

        ; gas properties
        ; -----------------
        ;rad= fload_gas_xyz('r')
        gas_mass= fload_gas_mass(1)

        ;e= fload_gas_energy()
        ke= fload_gas_energy(1,/kinetic)
        ;ke_radial= fload_gas_energy(1,/radial_kinetic)
        pe= fload_gas_energy(1,/potential)
        the= fload_gas_energy(1,/thermal)
        e= ke + pe + the
        ;temp= fload_gas_temperature(1)
        ;keV= fload_gas_temperature(1,/keV)
        entropy= fload_gas_entropy(1)

        ; not really used at current, and slightly
        ; outdated
        ;xr= fload_gas_xray_luminosity(1)

        metallicity= fload_gas_metallicity(1)


        ; define wind
        ; --------------

	;r_wind= 50.0
        ; greater than r_wind
        ;gtr_idx= where(rad gt r_wind)
        ;if gtr_idx(0) ne -1 then mass_gtx[i]= total(gas_mass(gtr_idx)) else mass_gtx[i]= 0

        ; positive energy
        egt0_idx= where(e gt 0.0)
        ;if egt0_idx(0) ne -1 then mass_egt0[i]= total(gas_mass(egt0_idx)) else mass_egt0[i]= 0

        ; wind index
        ;wind_idx= gtr_idx
        wind_idx= egt0_idx



windz= metallicity(wind_idx) / 0.02
print, "wind z max/min= ", max(windz), min(windz)
print, "wind z mean = ", mean(windz)
print, "wind z median = ", median(windz)

idx=where(windz le 0)
if idx(0) ne -1 then windz(idx)= 1.0d-9
windz= alog10(windz)

; -----------------------------------------------

;if oplotit eq 150 then mannorm = 250.0
;if oplotit eq 50 then mannorm = 800.0

;temp= process_histogram(windz, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=1000.0)
;temp= process_histogram(windz, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=mannorm)
temp= process_histogram(windz, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
print, min(temp), max(temp)

; -----------------------------------------------


end





;========================================================================================



;
;   Entropy
; --------------------------------
pro ent_hist, junk

if not keyword_set(junk) then begin
   print, "  "
   print, "ent_hist, junk"
   print, "  "
   return
endif

filename='enthist.eps'


initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=1
setup_plot_stuff, 'ps', filename=filename, colortable=4



; -----------------------------
; set up constants
; -----------------------------


; radius, linear scale
xaxistitle= "!6Log Entropy (keV cm!E2!N)"
xmax = 4.0     ; also set these in the process_vwhist procedure
xmin = -1.0

; number (histogram)
;yaxistitle= "!6V!Dr!N (km s!E-1!N)"
yaxistitle= ' '
ymax = 1.2
ymin = 0.0



; ------------------
; Plot this up
; ------------------
x0= 0.05
y0= 0.15
x1= 0.95
y1= 0.98

!p.position= [x0, y0, x1, y1]

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
	ytickformat='(a1)', $
	;ytitle=yaxistitle, $
	/nodata   ;, /noerase


; -----------------------------------------------

;snapnum= 8
;snapnum= 49
;snapnum= 78
;snapnum= 81
snapnum= 107

;process_enthist, "/raid4/tcox/vc3vc3e_no", snapnum, oplotit= 50
;xyouts, 0.70, 0.85, "e_no", size=1.5, color=50, /normal, charthick=3.0

process_enthist, "/raid4/tcox/vc3vc3e_2", snapnum, oplotit= 150, xmax=xmax, xmin=xmin
xyouts, 0.10, 0.80, "BH-driven wind", size=1.5, color=150, /normal, charthick=3.0

process_enthist, "/raid4/tcox/sbw/sb10", snapnum, oplotit= 50, xmax=xmax, xmin=xmin
xyouts, 0.10, 0.75, "SB-driven wind", size=1.5, color=50, /normal, charthick=3.0

;process_enthist, "/raid4/tcox/sbw/sb10BH", snapnum, oplotit= 150, xmax=xmax, xmin=xmin
;xyouts, 0.10, 0.70, "sb10BH", size=1.5, color=150, /normal, charthick=3.0

;xyouts, 0.70, 0.90, fload_timelbl(1,2,/noteq), /normal, color= 0, charthick=3.0, size=1.5


; -----------------------------------------------



; print extras
; -------------

; done
; -----
device, /close


end







; ------------------
; compute histogram
; ------------------
pro process_enthist, frun, snapnum, oplotit= oplotit, xmax=xmax, xmin=xmin


ok=fload_snapshot_bh(frun,snapnum)

        ; gas properties
        ; -----------------
        ;rad= fload_gas_xyz('r')
        gas_mass= fload_gas_mass(1)

        ;e= fload_gas_energy()
        ke= fload_gas_energy(1,/kinetic)
        ;ke_radial= fload_gas_energy(1,/radial_kinetic)
        pe= fload_gas_energy(1,/potential)
        the= fload_gas_energy(1,/thermal)
        e= ke + pe + the
        ;temp= fload_gas_temperature(1)
        ;keV= fload_gas_temperature(1,/keV)
        entropy= fload_gas_entropy(1)

        ; not really used at current, and slightly
        ; outdated
        ;xr= fload_gas_xray_luminosity(1)

        metallicity= fload_gas_metallicity(1)


        ; define wind
        ; --------------

	;r_wind= 50.0
        ; greater than r_wind
        ;gtr_idx= where(rad gt r_wind)
        ;if gtr_idx(0) ne -1 then mass_gtx[i]= total(gas_mass(gtr_idx)) else mass_gtx[i]= 0

        ; positive energy
        egt0_idx= where(e gt 0.0)
        ;if egt0_idx(0) ne -1 then mass_egt0[i]= total(gas_mass(egt0_idx)) else mass_egt0[i]= 0

        ; wind index
        ;wind_idx= gtr_idx
        wind_idx= egt0_idx



windent= entropy(wind_idx)
print, "wind ent max/min= ", max(windent), min(windent)
print, "wind ent mean = ", mean(windent)
print, "wind ent median = ", median(windent)

idx=where(windent le 0)
if idx(0) ne -1 then windent(idx)= 1.0d-9
windent= alog10(windent)

; -----------------------------------------------

;if oplotit eq 150 then mannorm = 250.0
;if oplotit eq 50 then mannorm = 800.0

;temp= process_histogram(windent, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=1000.0)
;temp= process_histogram(windent, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=mannorm)
temp= process_histogram(windent, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
print, min(temp), max(temp)

; -----------------------------------------------


end





;========================================================================================


