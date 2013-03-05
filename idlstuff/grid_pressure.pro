

;-----------------------------------------------------------------------------------

pro p_hist, junk



if not keyword_set(junk) then begin
        print, " "
        print, " p_hist, junk"
        print, " "
        print, " "
        return 
endif




; -------------------
filename='phist.eps'



initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


; -----------------------------
; set up constants
; -----------------------------


; alignement angle
xaxistitle= "!6Log Pressure (K cm!E-3!N)"
xmax =  12.0
xmin = -7.0

; number (histogram)
yaxistitle= ' '
;ymax = 11.0
ymax = 1.05
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






; =====================
; =====================


snapnum= 52
frun= "/raid4/tcox/ds/vc3vc3e_2"

ok= fload_snapshot_bh(frun,snapnum)

; 
; pressure from the snapshot
oplotit= 150
gamma_minus_1= 5./3. - 1
rho = fload_gas_rho(1) * 674.59    ; cm^-3
;u = fload_gas_u(1)
temperature = fload_gas_temperature(1)    ; K
;pressure = gamma_minus_1 * rho * u
pressure = alog10(gamma_minus_1 * rho * temperature)
print, "N= ", n_elements(pressure), min(pressure), max(pressure)
;temp= process_histogram(pressure, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=13000.0)
temp= process_histogram(pressure, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, normalization=normalization)
print, min(temp), max(temp)
xyouts, 0.12, 0.90, "snapshot", size=1.2, color=oplotit, /normal, charthick=3.0


r= fload_gas_xyz('r')
idx=where(r lt 9.0)
temp= process_histogram(pressure(idx), xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=normalization)
print, min(temp), max(temp)
xyouts, 0.13, 0.86, "(dashed is < 9 kpc)", size=1.2, color=oplotit, /normal, charthick=3.0


gridfile= frun+"/desikagrids/griddata_52_50x50x50.txt"
;
; pressure from the grid
oplotit= 50
read_desika_grid, gridfile, Grid, Velocity, Density, Temp, StellarMass
rho = 10^Density              ; cm^-3
temperature = 10^Temp            ; K
pressure = alog10(gamma_minus_1 * rho * temperature)
print, "N= ", n_elements(pressure), min(pressure), max(pressure)
weighting = 10^(Density)
;temp= process_histogram(pressure, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=13000.0)
;temp= process_histogram(pressure, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
temp= process_histogram(pressure, weighting=weighting, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
print, min(temp), max(temp)
xyouts, 0.12, 0.81, '50!E3!N grid', size=1.2, color=oplotit, /normal, charthick=3.0




gridfile= frun+"/desikagrids/griddata_52_100x100x100.txt"
;
; pressure from the grid
oplotit= 100
read_desika_grid, gridfile, Grid, Velocity, Density, Temp, StellarMass
rho = 10^Density              ; cm^-3
temperature = 10^Temp            ; K
pressure = alog10(gamma_minus_1 * rho * temperature)
print, "N= ", n_elements(pressure), min(pressure), max(pressure)
weighting = 10^(Density)
;temp= process_histogram(pressure, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=13000.0)
;temp= process_histogram(pressure, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
temp= process_histogram(pressure, weighting=weighting, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
print, min(temp), max(temp)
xyouts, 0.12, 0.76, '100!E3!N grid', size=1.2, color=oplotit, /normal, charthick=3.0




; -----------------------------------------------------------------------------



; -------------
;  Done
; -------------

device, /close



end









;-----------------------------------------------------------------------------------

pro rho_hist, junk



if not keyword_set(junk) then begin
        print, " "
        print, " rho_hist, junk"
        print, " "
        print, " "
        return 
endif




; -------------------
filename='rhist.eps'



initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


; -----------------------------
; set up constants
; -----------------------------


; alignement angle
xaxistitle= "!6Log Density (cm!E-3!N)"
xmax =  5.0
xmin = -9.0

; number (histogram)
yaxistitle= ' '
;ymax = 11.0
ymax = 1.05
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






; =====================
; =====================


snapnum= 52
frun= "/raid4/tcox/ds/vc3vc3e_2"

ok= fload_snapshot_bh(frun,snapnum)

; 
; rho from the snapshot
rho = alog10(fload_gas_rho(1) * 674.59)    ; cm^-3
print, "N= ", n_elements(rho), min(rho), max(rho)

oplotit= 151
temp= process_histogram(rho, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, normalization=normalization)
print, min(temp), max(temp)
xyouts, 0.12, 0.90, "snapshot", size=1.2, color=oplotit, /normal, charthick=3.0

;
; this was to make sure the "weighting" was working OK
;oplotit= 200
;weighting= rho*0.0 + 1.0
;temp= process_histogram(rho, weighting=weighting, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, normalization=normalization)
;print, min(temp), max(temp)
;xyouts, 0.12, 0.90, "snapshot - testing", size=1.2, color=oplotit, /normal, charthick=3.0
;
;


oplotit= 150
r= fload_gas_xyz('r')
idx=where(r lt 9.0)
;temp= process_histogram(rho(idx), xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
temp= process_histogram(rho(idx), xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=normalization)
print, min(temp), max(temp)
xyouts, 0.13, 0.86, "(dashed is < 9 kpc)", size=1.2, color=oplotit, /normal, charthick=3.0




;gridfile= frun+"/desikagrids/griddata_52_50x50x50.txt"
gridfile= frun+"/desikagrids/griddata_52.txt"
;
; rho from the grid
oplotit= 50
read_desika_grid, gridfile, Grid, Velocity, Density, Temp, StellarMass
rho = Density              ; cm^-3
weighting = 10^(Density)
print, "N= ", n_elements(rho), min(rho), max(rho)
;temp= process_histogram(rho, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=13000.0)
;temp= process_histogram(rho, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
temp= process_histogram(rho, weighting=weighting, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
print, min(temp), max(temp)
xyouts, 0.12, 0.81, '50!E3!N grid', size=1.2, color=oplotit, /normal, charthick=3.0



;
;gridfile= frun+"/desikagrids/griddata_52_100x100x100.txt"
;;
;; rho from the grid
;oplotit= 100
;read_desika_grid, gridfile, Grid, Velocity, Density, Temp, StellarMass
;weighting= 10^(Density)
;rho = Density              ; cm^-3
;print, "N= ", n_elements(rho), min(rho), max(rho)
;;temp= process_histogram(rho, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=13000.0)
;;temp= process_histogram(rho, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
;temp= process_histogram(rho, weighting=weighting, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
;print, min(temp), max(temp)
;xyouts, 0.12, 0.76, '100!E3!N grid', size=1.2, color=oplotit, /normal, charthick=3.0
;



gridfile= frun+"/griddata_new_52.txt"
;
; rho from the grid
oplotit= 100
read_new_desika_grid, gridfile, Grid, Velocity, Density_Cold, Temp_Cold, Density_Hot, Temp_Hot, SFR, StellarMass
rho = Density_Hot              ; cm^-3
weighting = 10^(Density_Hot)
print, "N= ", n_elements(rho), min(rho), max(rho)
;temp= process_histogram(rho, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=13000.0)
;temp= process_histogram(rho, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
temp= process_histogram(rho, weighting=weighting, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
print, min(temp), max(temp)
xyouts, 0.12, 0.76, 'New 50!E3!N grid', size=1.2, color=oplotit, /normal, charthick=3.0





; -----------------------------------------------------------------------------



; -------------
;  Done
; -------------

device, /close



end









;-----------------------------------------------------------------------------------

pro temp_hist, junk



if not keyword_set(junk) then begin
        print, " "
        print, " temp_hist, junk"
        print, " "
        print, " "
        return 
endif




; -------------------
filename='thist.eps'



initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


; -----------------------------
; set up constants
; -----------------------------


; alignement angle
xaxistitle= "!6Log Temperature (K)"
xmax =  8.0
xmin = 1.0

; number (histogram)
yaxistitle= ' '
;ymax = 11.0
ymax = 1.05
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






; =====================
; =====================


snapnum= 52
frun= "/raid4/tcox/ds/vc3vc3e_2"

ok= fload_snapshot_bh(frun,snapnum)

; 
; temperature from the snapshot
temperature = alog10(fload_gas_temperature(1))    ; K
temperature_cold= alog10(fload_gas_temperature_multi(1,/cold))  ; K
print, "N= ", n_elements(temperature), min(temperature), max(temperature)

;oplotit= 151
;temp= process_histogram(temperature, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
;temp= process_histogram(temperature, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, normalization=normalization)
;temp= process_histogram(temperature_cold, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, normalization=normalization)
;print, min(temp), max(temp)
;xyouts, 0.12, 0.90, "snapshot (cold)", size=1.2, color=oplotit, /normal, charthick=3.0

oplotit= 150
r= fload_gas_xyz('r')
idx=where(r lt 9.0)
;temp= process_histogram(temperature(idx), xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=normalization)
temp= process_histogram(temperature_cold(idx), xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=normalization)
print, min(temp), max(temp)
xyouts, 0.12, 0.86, "cold gas within 9 kpc (dashed)", size=1.2, color=oplotit, /normal, charthick=3.0
;xyouts, 0.13, 0.86, "(dashed is < 9 kpc)", size=1.2, color=oplotit, /normal, charthick=3.0


oplotit= 200
r= fload_gas_xyz('r')
idx=where(r lt 9.0)
temp= process_histogram(temperature(idx), xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=normalization)
print, min(temp), max(temp)
xyouts, 0.12, 0.90, "snapshot within 9 kpc, no multi-breakdown", size=1.2, color=oplotit, /normal, charthick=3.0



;-----------------------------------------------------

;gridfile= frun+"/desikagrids/griddata_52_50x50x50.txt"
gridfile= frun+"/desikagrids/griddata_52.txt"
;
; temperature from the grid
oplotit= 50
read_desika_grid, gridfile, Grid, Velocity, Density, Temp, StellarMass
;
; v0
temperature = Temp            ; K
print, "N= ", n_elements(temperature), min(temperature), max(temperature)
;
; v1
weighting = 10^(Density)
;weighting = 10^(Density) / total(10^(Density))
;temperature = weighting * 10^Temp      ; do a density/mass weighting
;print, "N= ", n_elements(temperature), min(temperature), max(temperature)
;temperature= alog10(temperature)
;
;temp= process_histogram(temperature, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=13000.0)
temp= process_histogram(temperature, weighting=weighting, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
print, min(temp), max(temp)
xyouts, 0.12, 0.80, '50!E3!N grid', size=1.2, color=oplotit, /normal, charthick=3.0




do_100= 0
;do_100= 1
if do_100 eq 1 then begin
	gridfile= frun+"/desikagrids/griddata_52_100x100x100.txt"
	;
	; temperature from the grid
	oplotit= 100
	read_desika_grid, gridfile, Grid, Velocity, Density, Temp, StellarMass
	temperature = Temp            ; K
	print, "N= ", n_elements(temperature), min(temperature), max(temperature)
	weighting = 10^(Density)
	;temp= process_histogram(temperature, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit, mannorm=13000.0)
	;temp= process_histogram(temperature, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
	temp= process_histogram(temperature, weighting=weighting, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
	print, min(temp), max(temp)
	xyouts, 0.12, 0.75, '100!E3!N grid', size=1.2, color=oplotit, /normal, charthick=3.0
endif



;do_new= 0
do_new= 1
if do_new eq 1 then begin
	gridfile= frun+"/griddata_new_52.txt"
	;
	; temperature from the grid
	oplotit= 100
	read_new_desika_grid, gridfile, Grid, Velocity, Density_Cold, Temp_Cold, Density_Hot, Temp_Hot, SFR, StellarMass
	temperature = Temp_Hot              ; K
	print, "N= ", n_elements(temperature), min(temperature), max(temperature)
	weighting = 10^(Density_Hot)
	temp= process_histogram(temperature, weighting=weighting, xmax=xmax, xmin=xmin, levels=100, oplotit=oplotit)
	print, min(temp), max(temp)
	xyouts, 0.12, 0.75, 'New 50!E3!N grid', size=1.2, color=oplotit, /normal, charthick=3.0
endif





; -----------------------------------------------------------------------------



; -------------
;  Done
; -------------

device, /close



end







