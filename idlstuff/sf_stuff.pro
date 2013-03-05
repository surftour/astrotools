; ------------------------------------------------------------------------------------------
;
; returns the gas surface density in gadget units (10^10 m_solar / kpc^2)
;
; r= radial array

function initialize_gas_surface_density, r, M_gas, R_d

Sigma_0= M_gas / (2.0 * !PI * R_d * R_d)
print, "Disk Sigma_0= ", Sigma_0
return, Sigma_0 * exp(-r/R_d)

end


; ------------------------------------------------------------------------------------------
;
; returns the sfr surface density in gadget units
; so 10^10 Msolar / Gyr / kpc^2 = 10 Msolar / Yr / kpc^2
;

function sfr_surface_density, r, Gas_sd, N=N

; std sfr law
SFR_Normalization= 2.5d-4                      ; M_solar Yr-1  kpc-2
;print, "SFR Normalization= ", SFR_Normalization
;print, "SFR Gas Law N= ", N
;gas_sd= gas_surface_density(r) * 1.0d+4        ; convert to m_solar pc^-2
gas_sd_i= Gas_sd * 1.0d+4                        ; convert to m_solar pc^-2
sfr_sd= SFR_Normalization * (gas_sd_i^N)
sfr_sd= sfr_sd / 10.0                           ; convert to gadget units


; density cut-off
sd_cut_off= 10.0    ; 10 M_solar pc^-2
idx= where(gas_sd_i LE sd_cut_off)
if idx(0) ne -1 then sfr_sd(idx)= 1.0d-12


return, sfr_sd

end


; ------------------------------------------------------------------------------------------
;
; analytical SFR for an exponential disk and std sfr law
;

function analytic_SFR, M_gas, R_d, N=N

if not keyword_set(N) then N= 3./2.
SFR_Normalization= 2.5d-4                      ; M_solar Yr-1  kpc-2

Sigma_0= M_gas / (2.0 * !PI * R_d * R_d)
Sigma_0= Sigma_0 * 1.0d+4                      ; convert to m_solar pc^-2

return, 2.*!PI*SFR_Normalization*(Sigma_0^N)*(R_d/N)^2.

end


; ------------------------------------------------------------------------------------------
;
; find radius equal to critical surface density
;

function manual_func_exp, x, p
	;print, x, p
        return, p - exp(-x)
end

function sdthresh_radius, M_gas, R_d

Sigma_Thresh= 10.0                    ;  m_solar pc^-2

Sigma_0= M_gas / (2.0 * !PI * R_d * R_d)
Sigma_0= Sigma_0 * 1.0d+4                      ; convert to m_solar pc^-2

if Sigma_Thresh GE Sigma_0 then return, -1

pp= Sigma_Thresh / Sigma_0
print, "Sigma_Thresh / Sigma_0= ", pp

;init_guess= [0, 1.0, 10.0]
;x= fx_root(init_guess, 'func_exponential_stripped', /double)
;x= fx_root_tj(init_guess, 'manual_func_exp', PARINFO=pp, /double)

rmin= 0.0
rmax= 10.0
x= zbrent_tj(rmin, rmax, FUNC_NAME='manual_func_exp', PARINFO=pp)

return, x*R_d

end




; ==========================================================================================
;



pro disk_sfr_time, junk

if not keyword_set(junk) then begin
	print, " "
	print, " disk_sfr_time, junk"
	print, " "
	print, " "
	return
endif

;------------------------------
;   (gadget) system of units
;       mass - 10^10 Msolar
;       length - kpc
;       time - Gyr
;

M_gas= 3.12387
R_d= 2.75

N= 3./2.


;------------------------------

print, "Expected Initial SFR= ", analytic_SFR(M_gas, R_d)

;------------------------------
;  Time Integration Parameters
;   (units are Gyr)
;
T_max= 2.0
dt= 0.1
N_time= long(T_max/dt)
;
time= fltarr(N_time)
sfr_time= fltarr(N_time)
mass_time= fltarr(N_time)
;
for ii= 0L, N_time-1 do time[ii]= 0.5 * (ii*dt + (ii+1)*dt)
;
;------------------------------
;  Radial Integration Parameters
;   (units are kpc)
;
R_max= 50.0
dr= 0.01
N_r= long(R_max/dr)
;
R_i= fltarr(N_r)
for ii= 0L, N_r-1 do R_i[ii]=  0.5 * (ii*dr + (ii+1)*dr)
;
area= 2. * !PI * R_i * dr
;
;M_gas_r= fltarr(N_r)
Gas_sd= initialize_gas_surface_density(R_i, M_gas, R_d)
M_gas_r= Gas_sd * 2. * !PI * R_i * dr
print, "Total Gas Mass= ", total(M_gas_r)
;

;------------------------------
;------------------------------
;

for i= 0L, N_time-1 do begin

   sfr_sd= sfr_surface_density(R_i, Gas_sd, N=N)

   SFR_r= sfr_sd * area        ; SFR in gadget units (10^10 msolar / Gyr)
   SF_r= sfr_sd * dt * area
   M_gas_r= M_gas_r - SF_r
   Gas_sd= Gas_sd - (SF_r/area)

   ; keep trackRof global quantities
   sfr_time[i]= 10.0 * total(SFR_r)     ; converts gadget -> msolar/yr
   mass_time[i]= total(M_gas_r)

   print, "i= ", i, time[i], sfr_time[i], total(SF_r), mass_time[i]

endfor


;------------------------------
;------------------------------
;


plot_sfr_v_time, time, sfr_time
plot_mass_v_time, [0,time], [M_gas,mass_time]


end




;=========================================================================


pro plot_sfr_v_time, time, sfr

filename= 'disk_sfr.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6SFR (M!D!9n!6!N yr!E-1!N)"
xaxistitle = "!6Time (Gyr)"

xmax= 2.0
xmin= 0.0

ymax= 120.0
;ymax= 12.0
;ymin= 0.0
ymin= 0.01

;---------------------------

!p.position= [0.16, 0.14, 0.98, 0.98]

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


oplot, time, sfr, psym= -2, color= 50, linestyle= 0, thick= 4.0

;--------------------------------------

;   Get SFR rate from txt 
open_sfr_file, "vc3c", sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
;open_sfr_file, "vc3vc3e_2", sfrtime, sfrsfr, sfrmfs, finalnewstarmass, sfrgasmass
sfrtime = sfrtime / 0.7
;sfrsfr=sfrsfr/2.0
oplot, sfrtime, sfrsfr, psym=-3, color= 150, linestyle= 0, thick= 3.0

;--------------------------------------

device, /close

end




;=========================================================================



pro plot_mass_v_time, time, mass

filename= 'disk_mass.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4


;--------------------------------------
;  Print the Shit
;--------------------------------------

yaxistitle="!6Gas Mass (10!e10!n M!D!9n!6!N)"
xaxistitle = "!6Time (Gyr)"

xmax= 2.0
xmin= 0.0

ymax= 3.5
;ymin= 0.0
ymin= 0.01

;---------------------------

!p.position= [0.16, 0.14, 0.98, 0.98]

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


oplot, time, mass, psym= -2, color= 50, linestyle= 0, thick= 4.0

;--------------------------------------

device, /close


end




;=========================================================================





