pro Xwhatup, junk


if not keyword_set(junk) then begin
        print, " "
        print, " Xwhatup, junk"
        print, " "
        print, " (everything is fixed from within)"
        print, " "
        print, " "
        return
endif

; this assumes they are ordered
;  0 through


;frun="/raid4/tcox/vc3vc3h"
frun="/raid4/tcox/vc3vc3h_2"
;frun="/raid4/tcox/vc3vc3h_no"

snap=107

; open snapshot
ok=fload_snapshot_bh(frun,snap)


; determine gas masses
; ---------------------
gmass= fload_gas_mass(1)
print, "Total Gas Mass= ", total(gmass)
print, "Total Hot Gas Mass= ", total(fload_gas_mass(1,/hot))
print, "Total Cold Gas Mass= ", total(fload_gas_mass(1,/cold))
print, "Total SF Gas Mass= ", total(fload_gas_mass(1,/sf))


; entropy  (keV cm2)
; ------------------
print, "Average Entropy= ", fload_gas_entropy(1,/averageit)


; diffuse gas
; X-ray luminosity (ergs s-1)
;  (bremstauhlung)
; ----------------------------
xr_hg= fload_gas_xray_luminosity(1,/diffuse_hotgas)
xrtot= total(xr_hg)
print, "Xray Luminosity (Brem) = ", alog10(xrtot)

idx= where(xr_hg gt 0.0)
print, "Mass of X-ray Gas= ", total(gmass(idx))

;xr= fload_gas_xray_luminosity(1,/sf_hotgas)
;xrtot= total(xr)
;if xrtot gt 0.0 then xray_sf[i]= alog10(xrtot) else xray_sf[i]= 0.0


; raymond-smith x-ray luminosity
;  (uses z)
; -------------------------------
;load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum, $
;	xr_weighted_temp=xr_weighted_temp, $
;	xr_weighted_z=xr_weighted_z
;xrayz_soft_tot= total(soft_xray_lum)
;if xrayz_soft_tot gt 0.0 then xrayz_soft[i]= alog10(xrayz_soft_tot) else xrayz_soft[i]= 0.0
;xrayz_hard_tot= total(hard_xray_lum)
;if xrayz_hard_tot gt 0.0 then xrayz_hard[i]= alog10(xrayz_hard_tot) else xrayz_hard[i]= 0.0

;load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum, /zero_metallicity
;soft_xray_lum_z0= 0
;hard_xray_lum_z0= 0
;xrayz0_soft_tot= total(soft_xray_lum_z0)
;if xrayz0_soft_tot gt 0.0 then xrayz0_soft[i]= alog10(xrayz0_soft_tot) else xrayz0_soft[i]= 0.0
;xrayz0_hard_tot= total(hard_xray_lum_z0)
;if xrayz0_hard_tot gt 0.0 then xrayz0_hard[i]= alog10(xrayz0_hard_tot) else xrayz0_hard[i]= 0.0


; gas temperature
;  * mass weighted (keV and K)
;  * x-ray emission weighted (keV)
; -----------------------------------
;keV[i]= fload_gas_temperature_kT(1,/averageit,/keV)
;keV[i]= fload_gas_temperature(1,/averageit,/keV)
;Temp[i]= fload_gas_temperature(1,/averageit)
keVi= fload_gas_temperature(1,/keV)
Tempi= fload_gas_temperature(1)

; mass weighted
keV_moment= total(keVi*gmass)
print, "Mass weighted average Temperature (keV)= ", keV_moment/total(gmass)

; weighted by the diffuse HG emission
keV_moment= total(keVi*xr_hg)
print, "Lx weighted average Temperature (keV)=", keV_moment/total(xr_hg)


; mass weighted - in K
Temp_moment=total(Tempi*gmass)
print, "Mass weighted Temperature (K)=", Temp_moment/total(gmass)



; Densities
; ----------
rho=fload_gas_rho(1)
print, "Average LX rho= ", mean(rho(idx))
rho_moment= total(rho*xr_hg)
print, "X-ray weighted Average rho= ", rho_moment/total(xr_hg)


; X-ray weighted metallicity
; --------------------------
;gasz= fload_gas_metallicity(1)
;z_moment= total(gasz*xr_hg)
;if total(xr_hg) gt 0 then z_Xray[i]= z_moment/total(xr_hg) else z_Xray[i]= 0.0




; ----------------------------------------
        
print, "---------------------------------"
print, "---------------------------------"
print, "         DONE - enjoy            " 
print, "---------------------------------"
print, "---------------------------------"
        




end





; if we don't want to do all snaps, then
; we use the following to just compile
; one snapshot info (such as the final one)

pro onetime_hotgas, frun, snapnum

if not keyword_set(frun) then begin
        print, " "
        print, " onetime_hotgas, frun, snapnum"
        print, " "
        return
endif

; open snapshot
ok=fload_snapshot_bh(frun,snapnum)

; what time is it?
ttime= fload_time(1)
print, "time= ", ttime


; determine gas masses
; ---------------------
gmass= fload_gas_mass(1)
print, "total gas mass= ", total(gmass)
print, "  hot gas mass= ", total(fload_gas_mass(1,/hot))
print, " cold gas mass= ", total(fload_gas_mass(1,/cold))
print, "   sf gas mass= ", total(fload_gas_mass(1,/sf))


; entropy  (keV cm2)
; ------------------
entropy= fload_gas_entropy(1,/averageit)
entropy= alog10(entropy)
print, "average entropy= ", entropy


; diffuse gas
; X-ray luminosity (ergs s-1)
;  (bremstauhlung)
; ----------------------------
xr_hg= fload_gas_xray_luminosity(1,/diffuse_hotgas)
xray= total(xr_hg)
xray= alog10(xray)
print, "Xray lum, diffuse hot gas= ", total(xr_hg)
xray_sf= alog10(total(fload_gas_xray_luminosity(1,/sf_hotgas)))
print, "Xray lum, sf gas= ", total(fload_gas_xray_luminosity(1,/sf_hotgas))


; raymond-smith x-ray luminosity
;  (uses z)
; -------------------------------
load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum
xrayz_soft_tot= total(soft_xray_lum)
xrayz_soft= alog10(xrayz_soft_tot)
print, "Xray lum, RS soft= ", xrayz_soft_tot
xrayz_hard_tot= total(hard_xray_lum)
xrayz_hard= alog10(xrayz_hard_tot)
print, "Xray lum, RS hard= ", xrayz_hard_tot

;load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum, /zero_metallicity
print, "WARNING: not calculating Raymond & Smith zero metallicity X-ray emission"
hard_xray_lum=10.0
soft_xray_lum=10.0
xrayz0_soft_tot= total(soft_xray_lum)
xrayz0_soft= alog10(xrayz0_soft_tot)
print, "Xray lum, RS(z0) soft= ", xrayz0_soft_tot
xrayz0_hard_tot= total(hard_xray_lum)
xrayz0_hard= alog10(xrayz0_hard_tot)
print, "Xray lum, RS(z0) hard= ", xrayz0_hard_tot


; gas temperature
;  * mass weighted (keV and K)
;  * x-ray emission weighted (keV)
; -----------------------------------
;keV[i]= fload_gas_temperature_kT(1,/averageit,/keV)
;keV[i]= fload_gas_temperature(1,/averageit,/keV)
;Temp[i]= fload_gas_temperature(1,/averageit)
keVi= fload_gas_temperature(1,/keV)
Tempi= fload_gas_temperature(1)
keV_moment= total(keVi*gmass)
print, "Temp (mass wt.)= ", keV_moment/total(gmass)
keV_moment= total(keVi*xr_hg)
if total(xr_hg) gt 0 then keV_Xray= keV_moment/total(xr_hg) else keV_Xray= 0.0
print, "Temp (xray wt.)= ", keV_Xray
Temp_moment=total(Tempi*gmass)
print, "Temp (K)= ",  Temp_moment/total(gmass)


; X-ray weighted metallicity
; --------------------------
gasz= fload_gas_metallicity(1)
z_moment= total(gasz*xr_hg)
if total(xr_hg) gt 0 then z_Xray= z_moment/total(xr_hg) else z_Xray= 0.0


print, " "
print, " --------------------------------- "
print, " "

; Xrays FILE
; -----------
openw, 1, frun+'/xrays.txt', ERROR=err
        
printf, 1, "#   xrays.txt"
printf, 1, "# "
printf, 1, "#        xr wt.            dif_hg     sf_hg      rs_s      rs_h    rs_z0_s   rs_z0_h  "
printf, 1, "# time   <temp>  log(ent) log(xray) log(xray) log(xray) log(xray) log(xray) log(xray) "               
printf, 1, "# (Gyr) (keV/mp)(keV cm2)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s) "
printf, 1, FORMAT= '(F6.3," ",F8.5,"  ", 7(F8.3,"  "))', $
                ttime, keV_Xray, entropy, xray, xray_sf, $
                        xrayz_soft, xrayz_hard, xrayz0_soft, xrayz0_hard
close, 1


end




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




; ----------------------------
;  Read xrays.txt file
; ----------------------------
pro read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, $
				xray, xray_sf, $
				xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

hgasfile= frun+'/xrays.txt'

spawn, "wc "+hgasfile,result
lines=long(result)
lines=lines(0)-5
if lines GT 0 then hgas_data= fltarr(11,lines)

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

; with new xrays.txt file (i.e. it has X-ray emission weighted Z)
z_X= hgas_data[2,*]
mass_xraygas= hgas_data[3,*]
entropy= hgas_data[4,*]
xray= hgas_data[5,*]
xray_sf= hgas_data[6,*]
xray_rs_s= hgas_data[7,*]
xray_rs_h= hgas_data[8,*]
xray_rs0_s= hgas_data[9,*]
xray_rs0_h= hgas_data[10,*]



end





; -------------------------------------------------------------------------------------------






;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro plot_xraylum_vs_time_comparison, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_xraylum_vs_time, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='xraylum.eps'

initialize_plotinfo, 1
;setup_plot_stuff, 'ps', filename=filename, colortable=4
setup_plot_stuff, 'ps', filename=filename, colortable=0



yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
xaxistitle = "Time (Gyr)"
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
frun="/raid4/tcox/vc3bvc3b"
;frun="/raid4/tcox/vc3avc3a_3"
;frun="/raid4/tcox/vc3ba_2_3"
read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

time= time/0.7
symsize= 1.0
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
oplot, time, xray, thick=3.0, psym=-8, color= 50
;oplot, time, xray, thick=3.0, psym=-2, color= 50
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
oplot, time, xray_rs_s, thick=3.0, psym=-8, color= 50, linestyle= 2
;oplot, time, xray_rs_s, thick=3.0, psym=-3, color= 50, linestyle= 1

;xyouts, 0.65, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=50
usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
oplot, [2.2, 2.4], [41.25,41.25], thick=3.0, psym=-8, color= 50
;oplot, [2.2, 2.4], [40.65,40.65], thick=3.0, psym=-8, color= 50
xyouts, 0.65, 0.82, 'black hole', /normal, charthick=3.0, size=1.33, color=50
;xyouts, 0.82, 0.87, '1e3', /normal, charthick=3.0, size=1.33, color=50

usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
oplot, [2.2, 2.4], [41.55,41.55], thick=3.0, psym=-8, color= 50, linestyle=2
xyouts, 0.65, 0.87, 'black hole (RS)', /normal, charthick=3.0, size=1.33, color=50


; Galaxy 2
;----------
;frun="/raid4/tcox/As/A2"
;frun="/raid4/tcox/vc2vc2"
;frun="pool/vc3vc3_wBH"
frun="/raid4/tcox/vc3bvc3b_no"
;frun="/raid4/tcox/vc3avc3a_2"
;frun="/raid4/tcox/vc3ba_2_2"
read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

time= time/0.7
symsize= 1.0
;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
;oplot, time, xray, thick=3.0, psym=-8, color= 150
oplot, time, xray, thick=3.0, psym=-2, color= 150
usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
;oplot, time, xray_rs_s, thick=3.0, psym=-8, color= 150, linestyle= 1
;oplot, time, xray_rs_s, thick=3.0, psym=-3, color= 150, linestyle= 1

;xyouts, 0.65, 0.85, fload_getid(frun), /normal, charthick=1, size=1.33, color=150
;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
oplot, [2.2, 2.4], [40.95,40.95], thick=3.0, psym=-2, color= 150
;oplot, [2.2, 2.4], [40.40,40.40], thick=3.0, psym=-2, color= 150
xyouts, 0.65, 0.77, 'no black hole', /normal, charthick=3.0, size=1.33, color=150
;xyouts, 0.82, 0.82, '1e4', /normal, charthick=3.0, size=1.33, color=150


; Galaxy 3
;----------
;frun="/raid4/tcox/As/A3"
;frun="/raid4/tcox/vc3vc3"
;frun="/raid4/tcox/vc3avc3a_1"
;frun="/raid4/tcox/vc3ba_2_1"
;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
;;read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
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
;;read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
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
;;read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
;oplot, time/0.7, xray, thick=3.0, psym=-6, color= 0
;oplot, [2.2, 2.4], [41.08,41.08], thick=3.0, psym=-6, color= 150
;xyouts, 0.82, 0.67, '1e7', /normal, charthick=3.0, size=1.33, color=0


; Galaxy 6
;----------
;frun="/raid4/tcox/vc6vc6a"
;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
;;read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
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






;--------------------------------------
;  Plot multiple ones
;----------------------------------------
pro plot_xraydivB_vs_time, junk


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


yaxistitle = "Log L!DX!N / L!DB!N"
xaxistitle = "Time (Gyr)"

;xmax = max(time)
xmax = 3.0
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
frun="pool/vc3bvc3b"
read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

read_colors_file, frun, time,  mags
Mb= mags[2,*]
Lb= -0.4*(Mb-5.51)

; xray and Lb are both in log here
xrayB= xray - Lb

; make it log
;idx= where(xrayB le 0)
;if idx(0) eq -1 then begin
;	xrayB= alog10(xrayB)
;endif else begin
;	xrayB(idx)= 1e10
;	xrayB= alog10(xrayB)
;	xrayB(idx)= 0.0
;endelse


oplot, time, xrayB, thick=3.0, psym=-2, color= 50

xyouts, 0.7, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=50


; Galaxy 2
;----------
;frun="pool/vc3vc3_wBH"
frun="pool/vc3bvc3b_wBH"
read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

read_colors_file, frun, time,  mags
Mb= mags[2,*]
Lb= -0.4*(Mb-5.51)

xrayB= xray-Lb


oplot, time, xrayB, thick=3.0, psym=-6, color= 150

xyouts, 0.7, 0.85, fload_getid(frun), /normal, charthick=1, size=1.33, color=150



device, /close




end




;--------------------------------------
;  Plot X-ray Luminosity
;----------------------------------------
pro plot_xraylum_vs_time, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_xraylum_vs_time, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='xraylum.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
xaxistitle = "Time (Gyr)"
;xmax = max(time)
xmax = 4.25
;xmax = 3.0
xmin = 0
ymax = 45.0
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



; Galaxy 1
;----------
;frun="/raid4/tcox/vc3bvc3b"
;frun="/raid4/tcox/vc3vc3"
frun="/raid4/tcox/localgroup/v2"

read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

time=time/0.7

oplot, time, xray, thick=3.0, psym=-2, color= 50
;oplot, time, xray_sf, thick=3.0, psym=-2, color= 150, linestyle= 2
oplot, time, xray_rs_s, thick=3.0, psym=-2, color= 0, linestyle= 2
oplot, time, xray_rs0_s, thick=3.0, psym=-2, color= 0, linestyle= 1

;xyouts, 0.65, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=50


include_agn= 0
if include_agn eq 1 then begin
   ; now include the AGN component
   ; -----------------------------
   get_blackhole_data, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd

   bhtime= bhtime/0.7

   ;
   ; ------------------
   ; cgs units, cm, sec, g                        
   L_solar= 3.9d+33                               ; in ergs per sec
   cc= 2.9979d+10                                 ; speed of light in cm/sec
   convert_sunyr_gsec= 6.30428d+25                ; convert msun/yr -> g/sec

   boloL_ergss= 0.1*bh_mdot_sunyr*cc*cc*convert_sunyr_gsec
   boloL_sun= boloL_ergss / L_solar
   print, "L_sun    max/min ", max(boloL_sun), min(boloL_sun)
   print, "L_erg s  max/min ", max(boloL_ergss), min(boloL_ergss)
   boloL_ergss_log= alog10(boloL_ergss)
   ;oplot, bhtime, boloL_ergss_log, thick=3.0, psym=-3, color= 0

   ;bolo_corr_f= BH_hardXlum_inlog(ar_sunyr)       ; L/L_hard
   bolo_corr_f= BH_hardXlum_inlog(boloL_sun)       ; L/L_hard
   ;print, "f_hard= ",bolo_corr_f
   L_hard= boloL_ergss / bolo_corr_f
   print, "Lx_hard  max/min ", max(L_hard), min(L_hard)
   L_hard= alog10(L_hard)
   oplot, bhtime, L_hard, thick=3.0, psym=-3, color=100, linestyle=2

   ;bolo_corr_f= BH_softXlum_inlog(ar_sunyr)       ; L/L_soft
   bolo_corr_f= BH_softXlum_inlog(boloL_sun)       ; L/L_soft
   ;print, "f_soft= ",bolo_corr_f
   L_soft= boloL_ergss / bolo_corr_f
   print, "Lx_soft  max/min ", max(L_soft), min(L_soft)
   L_soft= alog10(L_soft)
   oplot, bhtime, L_soft, thick=3.0, psym=-3, color=150, linestyle=2


   ;xyouts, 0.65, 0.85, 'AGN, bolometric', /normal, charthick=1, size=1.33, color=50
   oplot, [2.18, 2.42], [44.20,44.20], thick=3.0, psym=-3, color= 100
   ;xyouts, 0.65, 0.87, 'AGN, hard', /normal, charthick=3, size=1.33, color=100
   xyouts, 0.65, 0.87, 'AGN, 2-8 keV', /normal, charthick=3, size=1.33, color=100
   oplot, [2.18, 2.42], [43.65,43.65], thick=3.0, psym=-3, color= 150
   ;xyouts, 0.65, 0.82, 'AGN, soft', /normal, charthick=3, size=1.33, color=150
   xyouts, 0.65, 0.82, 'AGN, 0.5-2 keV', /normal, charthick=3, size=1.33, color=150
endif


;xyouts, 0.65, 0.73, 'SF Gas (dashed)', /normal, charthick=2, size=1.33, color=150
oplot, [2.2, 2.4], [43.08,43.08], thick=3.0, psym=-2, color= 50
xyouts, 0.65, 0.77, 'diffuse gas', /normal, charthick=3, size=1.33, color=50


; draw arrow for merger time
;timemerge=1.05/0.7
;arrow, timemerge, 37.2, timemerge, 38.2, COLOR=0, THICK=3.0, hthick=3.0, /data



device, /close


end










;----------------------------------------------------------------------------------------


;===================================
;
;  X-ray Luminosity (Hot Gas + BH)
;
;===================================

pro xraylum, frun, filename=filename, $
		h=h, bhmsg=bhmsg

if not keyword_set(bhmsg) then bhmsg= ''
if not keyword_set(sendto) then sendto= 'ps'
if not keyword_set(frun) then begin
   print, "  "
   print, "xraylum, frun, /h, bhmsg=bhmsg, filename=filename"
   print, "  "
   print, "  need: .run bh"
   print, "  "
   return
endif


;-------------------------------------
;   Get BH Info from fid.bh file
;-------------------------------------


if keyword_set(h) then begin
	xaxistitle = "Time (Gyr)"
	h = fload_cosmology('h')
	bhtime = bhtime / h
	bh_mass = bh_mass / h
	bh_totalmass = bh_totalmass / h
endif

;---------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4

xaxistitle = "Time (Gyr)"
yaxistitle = "Log Luminosity (ergs s!E-1!N)"

;xmax = max(time)
xmax = 3.0
xmin = 0
ymax = 42.0
ymin = 37.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        xstyle=1, $
	ystyle=1, $
	;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
	;ytickformat='exp_label', $
	xtitle=xaxistitle, $
	ytitle=yaxistitle, $
	/nodata



; Hot Gas: Galaxy
;-------------------
;frun="pool/vc3vc3"
;frun="pool/vc3vc3_wBH"
frun="/raid4/tcox/vc3bvc3b"
frun="/raid4/tcox/vc3bvc3b_no"
;frun="pool/vc3bvc3b_wBH"
read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

oplot, time, xray, thick=3.0, psym=-2, color= 50
;oplot, time, xray_sf, thick=3.0, psym=-2, color= 150, linestyle= 2
oplot, time, xray_rs_s, thick=3.0, psym=-5, color= 0, linestyle= 2
oplot, time, xray_rs0_s, thick=3.0, psym=-3, color= 150, linestyle= 1

xyouts, 0.65, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=50

snaptime= time


; Plot BH X-ray Luminosity
; -------------------------
plot_bh_too= 0
if plot_bh_too eq 1 then begin
   get_blackhole_data, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd
   idx= where(bh_num eq 1)
   if idx(0) ne -1 then time_merge= bhtime(idx(0))
   ar_sunyr= bh_mdot_sunyr

   ; cgs units, cm, sec, g
   L_solar= 3.9d+33                               ; in ergs per sec
   cc= 2.9979d+10                                 ; speed of light in cm/sec
   ;convert_sunyr_gsec= 6.446d+27                  ; convert msun/yr -> g/sec
   convert_sunyr_gsec= 6.30428d+25                  ; convert msun/yr -> g/sec
   boloL_ergss= 0.1*ar_sunyr*cc*cc*convert_sunyr_gsec
   boloL_sun= boloL_ergss / L_solar
   print, "L_sun    max/min ", max(boloL_sun), min(boloL_sun)
   print, "L_erg s  max/min ", max(boloL_ergss), min(boloL_ergss)
   boloL_ergss_log= alog10(boloL_ergss)
   boloL_ergss_log_avg= resample_array(snaptime, bhtime, boloL_ergss_log)
   oplot, snaptime, boloL_ergss_log_avg, thick=3.0, psym=-3, color= 0

   ;bolo_corr_f= BH_hardXlum_inlog(ar_sunyr)       ; L/L_hard
   bolo_corr_f= BH_hardXlum_inlog(boloL_sun)       ; L/L_hard
   ;print, "f_hard= ",bolo_corr_f
   L_hard= boloL_ergss / bolo_corr_f
   print, "Lx_hard  max/min ", max(L_hard), min(L_hard)
   L_hard= alog10(L_hard)
   L_hard_avg= resample_array(snaptime, bhtime, L_hard)
   oplot, snaptime, L_hard_avg, thick=3.0, psym=-3, color=50, linestyle=2

   ;bolo_corr_f= BH_softXlum_inlog(ar_sunyr)       ; L/L_soft
   bolo_corr_f= BH_softXlum_inlog(boloL_sun)       ; L/L_soft
   ;print, "f_soft= ",bolo_corr_f
   L_soft= boloL_ergss / bolo_corr_f
   print, "Lx_soft  max/min ", max(L_soft), min(L_soft)
   L_soft= alog10(L_soft)
   L_soft_avg= resample_array(snaptime, bhtime, L_soft)
   oplot, snaptime, L_soft_avg, thick=3.0, psym=-3, color=150, linestyle=2
endif



; plot raymond-smith X-ray lum
; -----------------------------
plot_rs_too= 1
if plot_rs_too eq 1 then begin
	;ok=fload_snapshot_bh(frun,30)
	;load_raymondsmith_lums, hard, soft
	;oplot, time, xray, thick=3.0, psym=-2, color= 50
endif



; plot extras
; ------------
xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0
if bhmsg NE '' then xyouts, 0.25, 0.80, bhmsg, /normal, charthick=3.0, size=1.33, color=0

timemerge= 1.1
if timemerge gt 0 then begin
        arrow, timemerge, 1.0e-4, timemerge, 10^(-3.5), COLOR=0, THICK=3.0, hthick=3.0, /data
endif


device, /close


end





;  Do boxcar average so many points
;  matches the small times.
;----------------------------------
function resample_array, sntime, manytime, manypts

print, "N_smalltime= ", n_elements(sntime)
print, "N_manytime=  ", n_elements(manytime)

avgpts= fltarr(n_elements(sntime))

for i=0, n_elements(sntime)-1 do begin

	idx= where(manytime ge sntime[i])
	avgpts[i]= manypts[idx(0)]

endfor

return, avgpts

end




; -----------------------------------------------------------------------------






;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;       L_x   -   T_x  
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro plot_Lx_Tx, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_Lx_Tx, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='LxTx.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "Log L!DX!N (ergs s!E-1!N)"
xaxistitle = "Temperature (keV)"
xmax = 4.0
;xmin = 0.10
xmin = 0.04
ymax = 43.5
;ymin = 38.6
ymin = 38.5


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
	/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; vc1's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc1vc1", 3, /tx, yevolfac=2.0/3.0
;readandplot_lx_and_else, "/raid4/tcox/vc1vc1", 7, /tx

; vc2's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc2vc2", 3, /tx, yevolfac=0.5
;readandplot_lx_and_else, "/raid4/tcox/vc2vc2", 8, /tx

; vc3b with various masses
;----------------------------
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_3", 1, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_2", 1, /tx
readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_1", 3, /tx, yevolfac=2.0/3.0
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_1", 1, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_4", 1, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_5", 1, /tx

; vc3 with varous bh seed masses
;--------------------------------
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_3", 2, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_2", 2, /tx
readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_1", 3, /tx, yevolfac=7.0/6.0
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_1", 2, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_4", 2, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_6", 2, /tx

; vc4's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc4vc4a", 3, /tx, yevolfac=4.0/3.0
;readandplot_lx_and_else, "/raid4/tcox/vc4vc4a", 3, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc4avc4a", 3, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc4vc4_no", 3, /tx

; vc5's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc5vc5a", 3, /tx, yevolfac=2.0
;readandplot_lx_and_else, "/raid4/tcox/vc5vc5a", 4, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc5avc5a", 4, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc5vc5_no", 4, /tx

; vc6's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc6vc6a", 3, /tx, yevolfac=13.0/6.0
;readandplot_lx_and_else, "/raid4/tcox/vc6vc6a", 5, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc6avc6a", 5, /tx
;readandplot_lx_and_else, "/raid4/tcox/vc6vc6_no", 5, /tx

; a's
; -----
;readandplot_lx_and_else, "/raid4/tcox/As/A1", 6, /tx
;readandplot_lx_and_else, "/raid4/tcox/As/A2", 6, /tx
;readandplot_lx_and_else, "/raid4/tcox/As/A3", 6, /tx
;readandplot_lx_and_else, "/raid4/tcox/As/A4", 6, /tx
;readandplot_lx_and_else, "/raid4/tcox/As/A5", 6, /tx
;



; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0


; Put actual data on there
; -------------------------
read_osullivan_03, loglx, loglb, tempx
;symsize= 0.2
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, tempx, loglx, psym=8, color= 0
oplot, tempx, loglx, psym=7, color= 0   ;, symsize=0.5


; a little key
x0=0.06  &  x1= 0.28
y0=42.6  &  y1= 43.0
;xyouts, 0.7, 0.28, 'OFP (2001)', color= 0, charthick=2.0, /normal
;oplot, [x0+0.13], [y0+0.75], psym=8, color= 0
xyouts, 0.31, 0.83, 'OPC (2003)', color= 0, charthick=2.0, /normal
oplot, [x0+0.01], [y0+0.20], psym=7, color= 0
oplot, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], psym=-3, color=0



; Lx_T fit
; ----------
x= [0.01,100.0]
y= 42.45 + 4.8*(alog10(x))

oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0

; draw arrow for merger time
;timemerge=1.05/0.7
;arrow, timemerge, 36.8, timemerge, 37.5, COLOR=0, THICK=3.0, hthick=3.0, /data


; helpful info
; --------------

; bh seed mass
;xyouts, 0.13, 42.6, 'BH seed mass', color= 0, charthick=1.2, /data
;arrow, 0.13, 42.5, 0.13, 42.0, COLOR=0, THICK=3.0, hthick=3.0, /data


;xyouts, 0.22, 0.22, 'Without discrete source contribution', /normal, color= 0, size=1.4


device, /close




end






;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;       L_x   -   L_b  
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro plot_Lx_Lb, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_Lx_Lb, junk"
	print, " "
	print, " "
	return
endif

filename='LxLb.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "Log L!DX!N (ergs s!E-1!N)"
xaxistitle = "Log L!DB!N (L!D!9n!3!N)"
xmax = 12.5
xmin = 8.5
ymax = 43.5
ymin = 37.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	;/ylog, $
	;/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; vc1's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc1vc1", 3, /lb, yevolfac=0.66
;readandplot_lx_and_else, "/raid4/tcox/vc1vc1", 7, /lb

; vc2's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc2vc2", 3, /lb, yevolfac=0.5
;readandplot_lx_and_else, "/raid4/tcox/vc2vc2", 8, /lb

; vc3b with various masses
;----------------------------
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_3", 1, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_2", 1, /lb
readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_1", 3, /lb, yevolfac=0.66
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_1", 1, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_4", 1, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_5", 1, /lb

; vc3 with varous bh seed masses
;--------------------------------
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_3", 2, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_2", 2, /lb
readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_1", 3, /lb, yevolfac=7.0/6.0
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_1", 2, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_4", 2, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_6", 2, /lb



; vc4's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc4vc4a", 3, /lb, yevolfac=4.0/3.0
;readandplot_lx_and_else, "/raid4/tcox/vc4vc4a", 3, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc4avc4a", 3, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc4vc4_no", 3, /lb

; vc5's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc5vc5a", 3, /lb, yevolfac=2.0
;readandplot_lx_and_else, "/raid4/tcox/vc5vc5a", 4, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc5avc5a", 4, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc5vc5_no", 4, /lb

; vc6's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc6vc6a", 3, /lb, yevolfac=13.0/6.0
;readandplot_lx_and_else, "/raid4/tcox/vc6vc6a", 5, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc6avc6a", 5, /lb
;readandplot_lx_and_else, "/raid4/tcox/vc6vc6_no", 5, /lb

; a's
; -----
;readandplot_lx_and_else, "/raid4/tcox/As/A1", 6, /lb
;readandplot_lx_and_else, "/raid4/tcox/As/A2", 6, /lb
;readandplot_lx_and_else, "/raid4/tcox/As/A3", 6, /lb
;readandplot_lx_and_else, "/raid4/tcox/As/A4", 6, /lb
;readandplot_lx_and_else, "/raid4/tcox/As/A5", 6, /lb


; evolutionary tracks
; --------------------
;readandplot_lx_and_else_time, "/raid4/tcox/vc3rem_vc3rem", /lb
;readandplot_lx_and_else_time, "/raid4/tcox/vc3rem_vc2", /lb
;readandplot_lx_and_else_time, "/raid4/tcox/vc3rem_vc1", /lb
;readandplot_lx_and_else_time, "/raid4/tcox/vc3rem_vc1a", /lb
;readandplot_lx_and_else_time, "/raid4/tcox/vc3rem_vc1b", /lb
;readandplot_lx_and_else_time, "/raid4/tcox/vc3rem_vc1c", /lb


; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0


; Put actual data on there
; -------------------------
read_osullivan_01, loglb, loglx, ttype
symsize= 0.2
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
oplot, loglb, loglx, psym=8, color= 0
;oplot, loglb, loglx, psym=7, color= 0, symsize= 0.5

read_osullivan_03, loglx, loglb, tempx
oplot, loglb, loglx, psym=7, color= 0


; a little key
x0=10.9  &  x1= 12.15
y0=37.4  &  y1= 38.4
xyouts, 0.7, 0.28, 'OFP (2001)', color= 0, charthick=2.0, /normal
oplot, [x0+0.13], [y0+0.75], psym=8, color= 0
xyouts, 0.7, 0.22, 'OPC (2003)', color= 0, charthick=2.0, /normal
oplot, [x0+0.13], [y0+0.25], psym=7, color= 0
oplot, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], psym=-3, color=0





; Lx_Lb fit
; ----------
; fit for all E's (OFP '01)
;x= [1.0,10.0,100.0]
;y= 17.98 + 2.17*x
;oplot, x, y, psym=-3, linestyle= 3, thick=2.0, color= 0

; bright E galaxy slope (OPC '03)
x= [1.0,10.0,100.0]
y= 11.9 + 2.7*x
oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0



; discrete source Lx  (Ciotti et al. 1991)
x= [1.0,10.0,100.0]
y= 29.45 + x
oplot, x, y, psym=-3, linestyle= 1, thick=2.0, color= 0


; helpful info
; --------------

; bh seed mass
;xyouts, 8.8, 42.6, 'BH seed mass', color= 0, charthick=1.2, /data
;arrow, 8.8, 42.5, 8.8, 42.0, COLOR=0, THICK=3.0, hthick=3.0, /data


;xyouts, 0.22, 0.22, 'Without discrete source contribution', /normal, color= 0, size=1.4

device, /close


end




;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------
;
;       Sigma   -   T_x  
;
;----------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;----------------------------------------

pro plot_sigma_Tx, junk


if not keyword_set(junk) then begin
	print, " "
	print, " plot_sigma_Tx, junk"
	print, " "
	print, " "
	return
endif

;filename=frun+'/sigma.eps'
filename='SigmaTx.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4



;yaxistitle = "Log X-ray Luminosity (ergs s!E-1!N)"
yaxistitle = "!7r!3 (km s!E-1!N)"
xaxistitle = "Temperature (keV)"
xmax = 4.0
xmin = 0.04
ymax = 500.0
;ymin = 38.6
ymin = 50.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        xstyle=1, $
        ystyle=1, $
	/ylog, $
	/xlog, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata



; vc1's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc1vc1", 3, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc1vc1", 7, /sig

; vc2's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc2vc2", 3, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc2vc2", 8, /sig

; vc3b with various masses
;----------------------------
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_3", 1, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_2", 1, /sig
readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_1", 3, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_1", 1, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_4", 1, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc3ba_2_5", 1, /sig

; vc3 with varous bh seed masses
;--------------------------------
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_3", 2, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_2", 2, /sig
readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_1", 3, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_1", 2, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_4", 2, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc3avc3a_6", 2, /sig

; vc4's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc4vc4a", 3, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc4vc4a", 3, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc4avc4a", 3, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc4vc4_no", 3, /sig

; vc5's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc5vc5a", 3, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc5vc5a", 4, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc5avc5a", 4, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc5vc5_no", 4, /sig

; vc6's
;--------
readandplot_lx_and_else, "/raid4/tcox/vc6vc6a", 3, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc6vc6a", 5, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc6avc6a", 5, /sig
;readandplot_lx_and_else, "/raid4/tcox/vc6vc6_no", 5, /sig

; a's
; -----
;readandplot_lx_and_else, "/raid4/tcox/As/A1", 6, /sig
;readandplot_lx_and_else, "/raid4/tcox/As/A2", 6, /sig
;readandplot_lx_and_else, "/raid4/tcox/As/A3", 6, /sig
;readandplot_lx_and_else, "/raid4/tcox/As/A4", 6, /sig
;readandplot_lx_and_else, "/raid4/tcox/As/A5", 6, /sig
;



; extras
; ---------
;xyouts, 0.2, 0.92, 'bulge', /normal, charthick=3.0, size=1.33, color= 0


; Put actual data on there
; -------------------------
read_osullivan_03, loglx, loglb, tempx
read_osullivan_03b, sigma
;symsize= 0.2
;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
;oplot, tempx, loglx, psym=8, color= 0
oplot, tempx, sigma, psym=7, color= 0   ;, symsize=0.5


; a little key
x0=0.06  &  x1= 0.28
y0=335.0  &  y1= 390.0
;xyouts, 0.7, 0.28, 'OFP (2001)', color= 0, charthick=2.0, /normal
;oplot, [x0+0.13], [y0+0.75], psym=8, color= 0
xyouts, 0.31, 0.83, 'OPC (2003)', color= 0, charthick=2.0, /normal
oplot, [x0+0.01], [y0+25.0], psym=7, color= 0
oplot, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], psym=-3, color=0



; Sigma_T Beta_spec=1 (I think?)
; -------------------
x= [0.001,10.0]
y= 10^(0.480426*alog10(x) + alog10(300.0))

oplot, x, y, psym=-3, linestyle= 2, thick=2.0, color= 0



; helpful info
; --------------

; bh seed mass
;xyouts, 0.13, 42.6, 'BH seed mass', color= 0, charthick=1.2, /data
;arrow, 0.13, 42.5, 0.13, 42.0, COLOR=0, THICK=3.0, hthick=3.0, /data


;xyouts, 0.22, 0.22, 'Without discrete source contribution', /normal, color= 0, size=1.4


device, /close




end










; genereric plotting program
; ----------------------------
pro readandplot_lx_and_else, frun, pointselection, msg, y0, $
				tx=tx, $
				lb=lb, $
				sig=sig, $
				yevolfac=yevolfac

    thispsym= -3
    thiscolor= 150   ;red 

    ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
    read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

    ; get last index
    lstidx= n_elements(time)-1  ; take last

    ; default y-axis
    ;yval= xray[lstidx]
    yval= xray_rs_s[lstidx]

    ; T_x
    ; ----
    if keyword_set(tx) or keyword_set(sig) then begin
	xval= temp_keV_X[lstidx]
    endif

    ; L_b
    ; -----
    read_colors_file, frun, time,  mags
    lstidx= n_elements(time)-1
    Mb= mags[2,*]
    Lumb= -0.4*(Mb-5.51)

    if keyword_set(lb) then begin
	xval= Lumb[lstidx]
    endif

    ; Sigma
    ; ------
    if keyword_set(sig) then begin
	read_sigma_file, frun, time, Asigxy, Asigxz, Asigyz, Asigavg, Asigerr, $
                                Bsigxy, Bsigxz, Bsigyz, Bsigavg, Bsigerr, $
                                gsigavg, gsigerr
	slstidx= n_elements(time)-1
	yval= Asigavg[slstidx]
    endif

    ; add point source component to L_x
    ;L_x_discrete = 29.5 + Lumb[lstidx]
    ;yval= yval - 38.0
    ;L_x_discrete= L_x_discrete - 38.0
    ;yval= (10^(yval)) + (10^(L_x_discrete))
    ;yval= alog10(yval) + 38.0


;  open (or filled) square
; --------------------------
if pointselection eq 1 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 150
endif

;  x's
; -----
if pointselection eq 2 then begin
        symsel= 7
        symcolor= 100
endif

;  open circle
; --------------
if pointselection eq 3 then begin
        symsize= 1.5
        usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
        ;usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0
        symsel= 8
        symcolor= 150
endif

;  open triangle
; ----------------
if pointselection eq 4 then begin
        symsize= 1.0
        ;usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0, /fill
        usersym,symsize*[-1,0,1,-1], symsize*[-1,1,-1,-1], thick=4.0
        symsel= 8
        symcolor= 0
endif

;  asterisk
; ----------
if pointselection eq 5 then begin
        symsel= 2
        symcolor= 200
endif

;  filled square
; --------------------------
if pointselection eq 6 then begin
        symsize= 1.0
        usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0,/fill
        ;usersym,symsize*[-1,-1,1,1,-1],symsize*[-1,1,1,-1,-1],thick=4.0
        symsel= 8
        symcolor= 20
endif

;  diamond
; ----------
if pointselection eq 7 then begin
        symsel= 4 
        symcolor= 120
endif

;  triangle (open)
; -----------------
if pointselection eq 8 then begin
        symsel= 5 
        symcolor= 180
endif

;  plus sign
; ----------
if pointselection eq 9 then begin
        symsel= 1 
        symcolor= 170
endif



    oplot, [xval], [yval], thick=3.0, psym=symsel, color= symcolor

    ;oplot, [2.2, 2.4], [41.48,41.48], thick=3.0, psym=-2, color= 50
    if not keyword_set(msg) then msg= ' '
    if not keyword_set(y0) then y0= 0.0
    xyouts, 0.73, y0, msg, /normal, charthick=3.0, size=1.33, color=thiscolor

if not keyword_set(sig) then begin
   ; draw arrow for pt. source contribution
   L_x_discrete = 29.5 + Lumb[lstidx]
   yval_wpts= yval - 38.0
   L_x_discrete= L_x_discrete - 38.0
   yval_wpts= (10^(yval_wpts)) + (10^(L_x_discrete))
   yval_wpts= alog10(yval_wpts) + 38.0
   arrow, [xval], [yval], [xval], [yval_wpts], COLOR=100, THICK=3.0, hthick=3.0, /data

   yval=yval_wpts

   ; draw arrow for 5 Gyr evolution
   xevolfac=0.0
   if not keyword_set(yevolfac) then yevolfac= 0.0
   yevolfac=alog10(1+yevolfac)
   if keyword_set(lb) then xevolfac=-1.0*alog10(2.0)
   arrow, [xval], [yval], [xval+xevolfac], [yval+yevolfac], COLOR=symcolor, THICK=3.0, hthick=3.0, /data

endif


end




; genereric plotting program
;  - only this one does time
;    evolution
; ----------------------------
pro readandplot_lx_and_else_time, frun, msg, y0, $
				tx=tx, $
				lb=lb

    thispsym= -3
    thiscolor= 150   ;red 

    ;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
    read_xrays_file, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h


    ;yval= xray
    yval= xray_rs_s

    ; T_x
    ; ----
    if keyword_set(tx) then begin
	xval= temp_keV_X
    endif

    ; L_b
    ; -----
    read_colors_file, frun, time,  mags
    lstidx= n_elements(time)-1
    Mb= mags[2,*]
    Lumb= -0.4*(Mb-5.51)

    if keyword_set(lb) then begin
	xval= Lumb
    endif

    ; add point source component to L_x
    ; ------------------------------------
    ;L_x_discrete = 29.5 + Lumb
    ;yval= yval - 38.0
    ;L_x_discrete= L_x_discrete - 38.0
    ;yval= (10^(yval)) + (10^(L_x_discrete))
    ;yval= alog10(yval) + 38.0


    oplot, xval, yval, thick=3.0, psym=thispsym, color= thiscolor

    ;oplot, [2.2, 2.4], [41.48,41.48], thick=3.0, psym=-2, color= 50
    if not keyword_set(msg) then msg= ' '
    if not keyword_set(y0) then y0= 0.0
    xyouts, 0.73, y0, msg, /normal, charthick=3.0, size=1.33, color=thiscolor

end



; -----------------------------------------------------------------------------------




;=========================================================
;
;  Use Raymond & Smith ('77, I think) to determine
;  the X-ray luminosity.
;
;
;=========================================================

pro load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum, $
					zero_metallicity=zero_metallicity, $
					xr_weighted_temp=xr_weighted_temp, $
					xr_weighted_z=xr_weighted_z


ngas= long(fload_npart(0))
;Ttime= float(fload_time(1))
Redshift= double(0.0)

GasMetallicity = (1/0.02)*fload_gas_metallicity(1)
GasNe = fload_gas_nume(1)
GasMass = 1.0e+10*fload_gas_mass(1)  ; in solar masses
GasTemp = float(fload_gas_temperature(1))   ; in K (for some reason this comes back in double)
GasTemp_keV= fload_gas_temperature(1,/keV)   ; not used, except to be passed back
GasHsml = fload_gas_hsml(1)


rho= fload_gas_rho(1)

; units
;hub= 0.7
hub= 1.0
ProtonMass=         1.6726e-24
UnitLength_in_cm=   3.085678d+21
UnitMass_in_g=      1.989d+43
UnitDensity_in_cgs= UnitMass_in_g/(UnitLength_in_cm^3)
print, UnitDensity_in_cgs

; hydrogen number density (nH cm-3)
GasHIRho= rho*0.76*(hub*hub)*UnitDensity_in_cgs/ProtonMass
GasHIRho= float(GasHIRho)

; electron number density (ne cm-3)
GasNeRho= float(GasNe*GasHIRho)

idx=where((rho lt 0.000854924) and (GasTemp ge 1.0e+5))
if idx(0) ne -1 then begin
    print, n_elements(idx), " out of ",ngas," particles will have diffuse X-ray emission."
    ngas= n_elements(idx)
    GasMetallicity= GasMetallicity(idx)
    GasNeRho= GasNeRho(idx)
    GasHIRho= GasHIRho(idx)
    GasMass= GasMass(idx)
    GasTemp= GasTemp(idx)
    GasTemp_keV= GasTemp_keV(idx)
    GasHsml= GasHsml(idx)
endif else begin
    ngas= 0
endelse


; for testing purposes
;    ngas= 10L
;    GasMetallicity= GasMetallicity(5600:5609)
;    GasNeRho= GasNeRho(5600:5609)
;    GasHIRho= GasHIRho(5600:5609)
;    GasMass= GasMass(5600:5609)
;    GasTemp= GasTemp(5600:5609)
;    GasHsml= GasHsml(5600:5609)
;

if keyword_set(zero_metallicity) then begin
	GasMetallicity(*)= 0.0
endif


;
; new incarnation of code doesn't
; like the zero metallicities
;
idx=where(GasMetallicity le 0.0)
if idx(0) ne -1 then begin
    GasMetallicity(idx)= 1.0e-5
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;       Find luminosities
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;NGas= long(1)
;GasMetallicity= [mean(GasMetallicity)]
;GasNeRho= [mean(GasNeRho)]
;GasHIRho= [mean(GasHIRho)]
;;GasMass= [1.0e+10]
;GasMass= total(GasMass)
;GasTemp= [mean(GasTemp)]
;GasHsml= [mean(GasHsml)]

help, NGas, GasMetallicity, GasNeRho, GasHIRho, GasMass, GasTemp, GasHsml
print, NGas
;print, "Z= ", GasMetallicity
;print, "Rho_Ne= ", GasNeRho
;print, "Rho_HI= ",GasHIRho
;print, "Mass= ",GasMass
;print, "Temp= ",GasTemp
;print, "Hsml= ",GasHsml


if(ngas gt 0) then begin

        soft_xray_lum = fltarr(ngas)
        hard_xray_lum = fltarr(ngas)

        ;S = CALL_EXTERNAL('/home/brant/code/idl/raymond_smith/raymond_smith', $
        S = CALL_EXTERNAL('/home/tcox/Tools/C-Routines_for_IDL/RaymondSmith/raymond_smith', $
                'raymond_smith', $
                Ngas, $
                Redshift, $
                GasMetallicity, $
                GasNeRho, $
		GasHIRho, $
                GasMass, $
                GasTemp, $
                GasHsml, $
                soft_xray_lum, $
                hard_xray_lum)

        ;LUMINOSITIES are in h^-1 units

endif else begin
        print,'No gas, no x-ray luminosity.'
	hard_xray_lum= [0]
	soft_xray_lum= [0]
endelse


xr_weighted_temp= 0.0
xr_weighted_z= 0.0
if total(soft_xray_lum) gt 0.0 then begin
	xr_weighted_temp= total(GasTemp_keV*soft_xray_lum)/total(soft_xray_lum)
	xr_weighted_z= total(GasMetallicity*soft_xray_lum)/total(soft_xray_lum)
endif

; brant returns this in solar luminosities
; multiply by 3.989e33 to get ergs/sec


hard_xray_lum = hard_xray_lum*3.989d33
soft_xray_lum = soft_xray_lum*3.989d33
print, "Total hard= ",total(hard_xray_lum)," erg sec-1"
print, "Total soft= ",total(soft_xray_lum)," erg sec-1"


end












;========================================
;
;   Read the O'Sullivan data
;
;========================================

; This is the O'Sullivan, Forbes & Ponman
; 2001, MNRAS, 328, 461 data:
; A catalogue and analysis of L_X of early-type galaxies 
;
pro read_osullivan_01, loglb, loglx, ttype


; this first file has spaces in the name which 
; screws things up, use  the second one.
;osullivanfile= '/home/tcox/osullivan_01_all.txt'

osullivanfile= '/home/tcox/OSullivan/osullivan_01_all2.txt'


;
;  Format of this file is:
;#
;#
;#
;#  Name		D	Log LB	 	Log LX	Source	T
;# 		(Mpc)	(LB)	 	(erg s1)	 	 
;ESO10114	30.12	9.93*	<	41.02	B	3.0
;ESO1074	38.89	10.22	<	40.94	B	4.0
; etc...
;

spawn, "wc "+osullivanfile,result
lines=long(result)
datalines=lines(0)-5
loglb= fltarr(datalines)
loglx= fltarr(datalines)
ttype= fltarr(datalines)

openr, 1, osullivanfile
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
	readf, 1, junk
	;print, junk
	tempjunk= strsplit(junk,/extract,count=count)
	name= tempjunk(0)
	;distance= float(tempjunk(1))
	loglb(i)= float(tempjunk(2))   ; doesn't need a trap for asterisk
	if count eq 7 then begin
		idontknow= tempjunk(3)
		loglx(i)= float(tempjunk(4))
		source= tempjunk(5)
		ttype(i)= float(tempjunk(6))
	endif else begin
		loglx(i)= float(tempjunk(3))
		source= tempjunk(4)
		ttype(i)= float(tempjunk(5))
	endelse

endfor

close, 1



end




;========================================
;========================================
; This is the O'Sullivan, Ponman & Collins
; 2003, MNRAS, 340, 1375 data:
; X-ray scaling properties of early-type galaxies
;
pro read_osullivan_03, loglx, loglb, tempx


; first we'll read the '01 data, to load the L_B
; values
osullivanfile01= '/home/tcox/OSullivan/osullivan_01_all2.txt'

spawn, "wc "+osullivanfile01,result
lines=long(result)
datalines01=lines(0)-5
name01= strarr(datalines01)
loglb01= fltarr(datalines01)
loglx01= fltarr(datalines01)
ttype01= fltarr(datalines01)

openr, 1, osullivanfile01
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
	readf, 1, junk
	tempjunk= strsplit(junk,/extract,count=count)
	name01(i)= tempjunk(0)
	loglb01(i)= float(tempjunk(2))   ; doesn't need a trap for asterisk
	if count eq 7 then begin
		loglx01(i)= float(tempjunk(4))
		ttype01(i)= float(tempjunk(6))
	endif else begin
		loglx01(i)= float(tempjunk(3))
		ttype01(i)= float(tempjunk(5))
	endelse

endfor

close, 1





; next we'll actually read the 03 data, and fix
; the values to send back
osullivanfile03= '/home/tcox/OSullivan/osullivan_03_lxfixed.txt'

spawn, "wc "+osullivanfile03,result
lines=long(result)
datalines03=lines(0)-5
loglb= fltarr(datalines03)
loglx= fltarr(datalines03)
tempx= fltarr(datalines03)

openr, 1, osullivanfile03
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
        readf, 1, junk
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
	idx= where(name eq name01)
	if idx(0) ne -1 then begin
		loglb(i)= loglb01(idx)
		;nh(i)= float(tempjunk(1))                 ; 10^21 cm^2
		tempx(i)= float(tempjunk(2))               ; keV
		;metallicity(i)= float(tempjunk(3))        ; in solar
		loglx(i)= float(tempjunk(4))
	endif else begin
		print, "PROBLEM: can't find ",name
	endelse

endfor

close, 1




end




;========================================
;========================================
; This is the O'Sullivan, Ponman & Collins
; 2003, MNRAS, 340, 1375 data (informational):
; X-ray scaling properties of early-type galaxies
;
pro read_osullivan_03b, sigma


; we read the 03 data, note that this is mainly
; informational data, so sigma, R_e, etc.  L_x, T_x,
; are all enclosed above.
osullivanfile03= '/home/tcox/OSullivan/osullivan_03_infofixed.txt'

spawn, "wc "+osullivanfile03,result
lines=long(result)
datalines03=lines(0)-5
sigma= fltarr(datalines03)

openr, 1, osullivanfile03
junk=''

; read the header
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk
readf, 1, junk

; read the data
for i=0,lines(0)-6 do begin
        readf, 1, junk
        tempjunk= strsplit(junk,/extract,count=count)
        name= tempjunk(0)
	sigma(i)= float(tempjunk(1))                 ; km sec^-1
	;distance(i)= float(tempjunk(3))             ; Mpc
	;re(i)= float(tempjunk(4))                   ; arcmin
	;ttype(i)= float(tempjunk(5))
	;environment(i)= float(tempjunk(6))

endfor

close, 1




end




;========================================
;========================================





pro sansome, junk

frun="/raid4/tcox/vc3vc3h_2"

;read_hotgas_file, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_xrays_file, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, $
				xray, xray_sf, $
				xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

read_colors_file, frun, time,  mags

Mb= mags[2,*]
Lb= -0.4*(Mb-5.51)

Mk= mags[8,*]
Lk= -0.4*(Mb-3.33)

time= time/0.7
nsnaps= n_elements(time)

; ----------------------------------------

openw, 1, 'sansome.txt', ERROR=err

printf, 1, "#            rs_s"
printf, 1, "# time   log(xray)  log(L_b)  log(L_k)"
printf, 1, "# (Gyr)   (ergs/s)   (solar)   (solar)"
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3,"   ", 3(F8.3,"  "))', $
                time[i], xray_rs_s[i], Lb[i], Lk[i]
endfor
close, 1


end




