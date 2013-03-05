pro time_hotgas, frun


if not keyword_set(frun) then begin
        print, " "
        print, " time_hotgas, frun"
        print, " "
        ;print, " needs: .run determine_lums"
        print, " "
        print, " "
        print, " "
        return
endif

; this assumes they are ordered
;  0 through


; determine the number of
; snapshots in frun directory
;spawn, "ls "+frun+"/snapshot* | wc ",result
;spawn, "/usr/bin/ls "+frun+"/snap* | wc ",result
spawn, "/bin/ls "+frun+"/snapshot* | wc ",result
nsnaps=long(result[0])
print, "found: ", nsnaps, " snapshot files"



time= fltarr(nsnaps)

entropy=      fltarr(nsnaps)
keV=          fltarr(nsnaps)
keV_Xray=     fltarr(nsnaps)
Temp=         fltarr(nsnaps)
xray_gasmass= fltarr(nsnaps)
xray=         fltarr(nsnaps)
xray_sf=      fltarr(nsnaps)
xrayz_hard=   fltarr(nsnaps)
xrayz_soft=   fltarr(nsnaps)
xrayz0_hard=  fltarr(nsnaps)
xrayz0_soft=  fltarr(nsnaps)
z_Xray=       fltarr(nsnaps)
mass_hotgas=  fltarr(nsnaps)
mass_coldgas= fltarr(nsnaps)
mass_sfgas=   fltarr(nsnaps)
mass_totalgas=fltarr(nsnaps)


;t_merg= fload_mergertime(frun)

; ----------------------------------------
; This part loops through the snapshots
; and compiles sigma and BH information
; ----------------------------------------

for i=0,nsnaps-1 do begin

        print, "--------------------------------------"


        ; open snapshot
        ;ok=fload_snapshot(frun,snapnum)
        ok=fload_snapshot_bh(frun,i)
        ;ok=fload_snapshot_bh(frun,i,/nopot_in_snap)    ; --> makes the rest of this impossible


        ; what time is it?
        time[i]= fload_time(1)
        TTime= float(time[i])


	; determine gas masses
	; ---------------------
	gmass= fload_gas_mass(1)
	mass_totalgas[i]= total(gmass)

	mass_hotgas[i]= total(fload_gas_mass(1,/hot))
	mass_coldgas[i]= total(fload_gas_mass(1,/cold))
	mass_sfgas[i]= total(fload_gas_mass(1,/sf))


	; entropy  (keV cm2)
	; ------------------
	entropy[i]= fload_gas_entropy(1,/averageit)
	if entropy[i] gt 0.0 then entropy[i]= alog10(entropy[i]) else entropy[i]= 0.0


	; diffuse gas
	; X-ray luminosity (ergs s-1)
	;  (bremstauhlung)
	; ----------------------------
	xr_hg= fload_gas_xray_luminosity(1,/diffuse_hotgas)
	xrtot= total(xr_hg)
	if xrtot gt 0.0 then xray[i]= alog10(xrtot) else xray[i]= 0.0

	idx= where(xr_hg gt 0.0)
	if idx(0) ne -1 then xray_gasmass[i]= total(gmass(idx)) else xray_gasmass[i]= 0.0

	xr= fload_gas_xray_luminosity(1,/sf_hotgas)
	xrtot= total(xr)
	if xrtot gt 0.0 then xray_sf[i]= alog10(xrtot) else xray_sf[i]= 0.0


	; raymond-smith x-ray luminosity
	;  (uses z)
	; -------------------------------
	;load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum, $
	;				xr_weighted_temp=xr_weighted_temp, $
	;				xr_weighted_z=xr_weighted_z
	soft_xray_lum= 0
	hard_xray_lum= 0
	xr_weighted_temp=0
	xr_weighted_z=0
	xrayz_soft_tot= total(soft_xray_lum)
	if xrayz_soft_tot gt 0.0 then xrayz_soft[i]= alog10(xrayz_soft_tot) else xrayz_soft[i]= 0.0
	xrayz_hard_tot= total(hard_xray_lum)
	if xrayz_hard_tot gt 0.0 then xrayz_hard[i]= alog10(xrayz_hard_tot) else xrayz_hard[i]= 0.0

	;load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum, /zero_metallicity
	soft_xray_lum_z0= 0
	hard_xray_lum_z0= 0
	xrayz0_soft_tot= total(soft_xray_lum_z0)
	if xrayz0_soft_tot gt 0.0 then xrayz0_soft[i]= alog10(xrayz0_soft_tot) else xrayz0_soft[i]= 0.0
	xrayz0_hard_tot= total(hard_xray_lum_z0)
	if xrayz0_hard_tot gt 0.0 then xrayz0_hard[i]= alog10(xrayz0_hard_tot) else xrayz0_hard[i]= 0.0


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
	keV[i]= keV_moment/total(gmass)

	; weighted by the diffuse HG emission
	keV_moment= total(keVi*xr_hg)
	if total(xr_hg) gt 0 then keV_Xray[i]= keV_moment/total(xr_hg) else keV_Xray[i]= 0.0

	temp_kev= keV_Xray[i]
	; weighted by the RS X-ray calculated emission  - doesn't work because size(soft_xray_lum) ne size(anything else)
	;keV_moment= total(keVi*soft_xray_lum)
	;if total(soft_xray_lum) gt 0 then keV_Xray[i]= keV_moment/total(soft_xray_lum) else keV_Xray[i]= 0.0
	keV_Xray[i]= xr_weighted_temp
print, "Temp comparison: ", temp_kev, " (diffuseHG weighted)   ", keV_Xray[i], " (RS soft weighted)"
	keV_Xray[i]= temp_kev
print, "Selecting diffuseHG weighting for temperature!"

	; mass weighted - in K
	Temp_moment=total(Tempi*gmass)
	Temp[i]= Temp_moment/total(gmass)


	; X-ray weighted metallicity
	; --------------------------
	gasz= fload_gas_metallicity(1)
	z_moment= total(gasz*xr_hg)
	if total(xr_hg) gt 0 then z_Xray[i]= z_moment/total(xr_hg) else z_Xray[i]= 0.0

	temp_z= z_Xray[i]
print, "Z comparison: ", z_Xray[i]," (diffuseHG weighted)   ", xr_weighted_z, " (RS soft weighted)"
	; RS calculation
	z_Xray[i]= xr_weighted_z
	z_Xray[i]= temp_z
print, "Selecting diffuseHG weighting for temperature!"

endfor

; ----------------------------------------
; write to files

; HOTGAS FILE
; -----------
openw, 1, frun+'/hotgas.txt', ERROR=err
        
printf, 1, "#   hotgas.txt"     
printf, 1, "# "
printf, 1, "#        xr wt.   mass wt     mass wt           "
printf, 1, "# time   <temp>    <temp>      <temp>  log(ent)  gas mass    hot       cold       sf "               
printf, 1, "# (Gyr) (keV/mp) (keV/mp)        (K)  (keV cm2)  (10^10 msolar ---->)   "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3," ",2(F8.5,"  "),F11.2," ", 5(F8.5,"  "))', $
                time[i], keV_Xray[i], keV[i], Temp[i], entropy[i], $
                        mass_totalgas[i], mass_hotgas[i], mass_coldgas[i], mass_sfgas[i]
endfor
close, 1


; Xrays FILE
; -----------
openw, 1, frun+'/xrays.txt', ERROR=err
        
printf, 1, "#   xrays.txt"
printf, 1, "# "
printf, 1, "#        xr wt.   xr wt.    xr gas              dif_hg     sf_hg      rs_s      rs_h   rs_z0_s   rs_z0_h  "
printf, 1, "# time   <temp>     Z        mass    log(ent) log(xray) log(xray) log(xray) log(xray) log(xray) log(xray) "
printf, 1, "# (Gyr) (keV/mp)(Z_solar) (10^10 m) (keV cm2)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s) "
for i=0,nsnaps-1 do begin
        printf, 1, FORMAT= '(F6.3," ", 4(F8.5,"  "), 6(F8.3,"  "))', $
                time[i], keV_Xray[i], z_Xray[i], xray_gasmass[i], entropy[i], xray[i], xray_sf[i], $
                        xrayz_soft[i], xrayz_hard[i], xrayz0_soft[i], xrayz0_hard[i]
endfor
close, 1


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
;ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap)

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
;xr_hg= fload_gas_xray_luminosity(1,/diffuse_hotgas)
xr_hg= fload_gas_xray_luminosity(1,/diffuse_hotgas,minkeVTemp=0.5)
xray= total(xr_hg)
xray= alog10(xray)
print, "Xray lum, diffuse hot gas= ", total(xr_hg)

xr= fload_gas_xray_luminosity(1,/sf_hotgas)
xrtot= total(xr)
if xrtot gt 0.0 then xray_sf= alog10(xrtot) else xray_sf= 0.0
print, "Xray lum, sf gas= ", xray_sf

idx= where(xr_hg gt 0.0)
if idx(0) ne -1 then xray_gasmass= total(gmass(idx)) else xray_gasmass= 0.0


; raymond-smith x-ray luminosity
;  (uses z)
; -------------------------------
;load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum
print, "WARNING: not calculating Raymond & Smith X-ray emission"   &   hard_xray_lum=10.0    & soft_xray_lum=10.0
xrayz_soft_tot= total(soft_xray_lum)
xrayz_soft= alog10(xrayz_soft_tot)
print, "Xray lum, RS soft= ", xrayz_soft_tot
xrayz_hard_tot= total(hard_xray_lum)
xrayz_hard= alog10(xrayz_hard_tot)
print, "Xray lum, RS hard= ", xrayz_hard_tot

;load_raymondsmith_xraylum, hard_xray_lum, soft_xray_lum, /zero_metallicity
print, "WARNING: not calculating Raymond & Smith zero metallicity X-ray emission"   & hard_xray_lum=10.0   & soft_xray_lum=10.0
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
;stop


print, " "
print, " --------------------------------- "
print, " "


; Xrays FILE
; -----------
openw, 1, frun+'/xrays.txt', ERROR=err

printf, 1, "#   xrays.txt"
printf, 1, "# "
printf, 1, "#        xr wt.   xr wt.    xr gas              dif_hg     sf_hg      rs_s      rs_h   rs_z0_s   rs_z0_h  "
printf, 1, "# time   <temp>     Z        mass    log(ent) log(xray) log(xray) log(xray) log(xray) log(xray) log(xray) "
printf, 1, "# (Gyr) (keV/mp)(Z_solar) (10^10 m) (keV cm2)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s) "
printf, 1, FORMAT= '(F6.3," ", 4(F8.5,"  "), 6(F8.3,"  "))', $
                ttime, keV_Xray, z_Xray, xray_gasmass, entropy, xray, xray_sf, $
                        xrayz_soft, xrayz_hard, xrayz0_soft, xrayz0_hard
close, 1




end




; -------------------------------------------------------------------------------------------------




; -------------------------------------------------------------------------------------------
; -------------------------------------------------------------------------------------------




;------------------------------------------------------
;  Scripts from Phil Hopkins related to BH luminosity
;------------------------------------------------------




;--------------------------------------------------------------
; Calculate BH luminosity/magnitude in different bands, given
;   an accretion rate (assumed to be in M_solar/yr)
;
;  Adapted this from Phil's script to calculate the
; soft and hard X-ray luminosity from the bolometric luminosity
;
; I think it returns it in L_solar (why not ergs sec-1)
;
;--------------------------------------------------------------
function BH_hardXlum_inlog, bolometric_lum
  L0 = alog10(bolometric_lum)                   ; Log(L) for magnitudes, conversions, etc.
print, min(L0), max(L0)
  L1 = L0 - 12.0                                ; (script L of Marconi et al. 2004)

  ; Marconi et al. (2004) relations for bolometric & x-ray luminosity
  ;logL_hard_xray = L0 - (1.54 + 0.24*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)        ; 2-10 keV
  logL_hard_xray = (1.54 + 0.24*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)        ; 2-10 keV

  L_hardX = 10^(logL_hard_xray)                 ; bolometric correction (L/L_hard)
print, min(L_hardX), max(L_hardX)

  return, L_hardX
end




function BH_softXlum_inlog, bolometric_lum
  L0 = alog10(bolometric_lum)                   ; Log(L) for magnitudes, conversions, etc.
print, min(L0), max(L0)
  L1 = L0 - 12.0                                ; (script L of Marconi et al. 2004)
  
  ; Marconi et al. (2004) relations for bolometric & x-ray luminosity
  ;logL_soft_xray = L0 - (1.65 + 0.22*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)        ; 0.5-2 keV
  logL_soft_xray = (1.65 + 0.22*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)        ; 0.5-2 keV
  
  L_softX = 10^(logL_soft_xray)                 ; bolometric correction (L/L_soft)
print, min(L_softX), max(L_softX)
  
  return, L_softX
end



function BH_Bband_correction, bolometric_lum
  L0 = alog10(bolometric_lum)                   ; Log(L) for magnitudes, conversions, etc.
print, min(L0), max(L0)
  L1 = L0 - 12.0                                ; (script L of Marconi et al. 2004)
  
  ; Marconi et al. (2004) relations for bolometric & x-ray luminosity
  logL_B= (0.80 - 0.067*L1 + 0.017*L1*L1 - 0.0023*L1*L1*L1)       ; B-band
  
  L_B = 10^(logL_B)                 ; bolometric correction (L/L_soft)
  
  return, L_B
end





;--------------------------------------------------------------
; Calculate BH luminosity/magnitude in different bands, given
;   an accretion rate (assumed to be in M_solar/yr)
;--------------------------------------------------------------
pro band_luminosity, bolometric_lum, L_hardX, L_B
  L0 = alog10(bolometric_lum)   ; Log(L) for magnitudes, conversions, etc.
  L1 = L0 - 12.0                                ; (script L of Marconi et al. 2004)
  
  ; Marconi et al. (2004) relations for bolometric & x-ray luminosity
  logL_soft_xray = L0 - (1.65 + 0.22*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)        ; 0.5-2 keV
  logL_hard_xray = L0 - (1.54 + 0.24*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)        ; 2-10 keV
  log_nuB_LB     = L0 - (0.80 - 0.067*L1 + 0.017*L1*L1 - 0.0023*L1*L1*L1)       ; B-band
  nuB = 6.818d14                ; assumes B-band at 0.44 microns
  L_B = (0.2) * (10^log_nuB_LB)         ; assumes spectral resolution about 5 (rough estimate only of B)


  ; Marconi alpha_OX (Zamorani 82, emperical Vignali, Brandt, & Scheider 03)
                R = 2.605       ; [defined as -log(nu(2500 Ang) / nu(2 keV)) ]
  logL_UV = (1.85*R + logL_soft_xray)                                           ; 2500 Angstroms
  
  L_hardX = 10^(logL_hard_xray)
  
  return
end





;--------------------------------------------------------------------------
; Function to calculate the luminosity (L_solar per Hz) 
;   at a given frequency f [GHz], from a given bolometric luminosity
;   (Complete broken-power-law model intrinsic QSO spectrum up to 10keV)
;--------------------------------------------------------------------------
function dL_df, L_bol, f
  L0 = alog10(L_bol)                    ; Log(L) for magnitudes, conversions, etc.
  L1 = L0 - 12.0                                ; (script L of Marconi et al. 2004)
  
  ; Marconi et al. (2004) relations for bolometric & x-ray luminosity
  logL_soft_xray = L0 - (1.65 + 0.22*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)        ; 0.5-2 keV
  logL_hard_xray = L0 - (1.54 + 0.24*L1 + 0.012*L1*L1 - 0.0015*L1*L1*L1)        ;  2-10 keV
  
  ; Now assume a photon index (photons/keV ~ E^(-Gamma)) Gamma=1.9
  ;             (Marconi refs. George et al. 1998, Perola et al. 2002, Nandra & Pounds 1994)
  Gamma = 1.9
    integration_factor = (1.0/(2.0-Gamma)) * ((10.0/2.0)^(2.0-Gamma) - 1.0)
    L_hard_xray = 10.0^logL_hard_xray   ; 2-10 keV band luminosity in L_solar
    f_2keV = 2.0 * (2.41775d17)                 ; frequency of 2keV in Hz
  dL_df_2keV = L_hard_xray / (integration_factor * f_2keV)


  ; Now can calculate other xray luminosities up to 10keV (at larger energies, the 
  ;   Compton reflection hump comes in, needs to be modeled)
  f_05keV = 0.25 * f_2keV               ; frequency [Hz] at 0.5 keV
  f_10keV = 5.0 * f_2keV                ; frequency [Hz] at 10 keV 
  if ((f_05keV LE f) AND (f_10keV GE f)) then DL_DF = dL_df_2keV * ((f/f_2keV)^(1 - Gamma))
  if (f GT f_10keV) then DL_DF = 0.0    ; Sorry, not incorporated yet! (need Compton modeling)
  dL_df_05keV = dL_df_2keV * ((f_05keV/f_2keV)^(1 - Gamma))
  
  ; Now assume a photon index (photons/keV ~ E^(-Gamma)) Gamma=1.9
  ;             (Marconi refs. George et al. 1998, Perola et al. 2002, Nandra & Pounds 1994)
  Gamma = 1.9
    integration_factor = (1.0/(2.0-Gamma)) * ((10.0/2.0)^(2.0-Gamma) - 1.0)
    L_hard_xray = 10.0^logL_hard_xray   ; 2-10 keV band luminosity in L_solar
    f_2keV = 2.0 * (2.41775d17)                 ; frequency of 2keV in Hz
  dL_df_2keV = L_hard_xray / (integration_factor * f_2keV)


  ; Now can calculate other xray luminosities up to 10keV (at larger energies, the 
  ;   Compton reflection hump comes in, needs to be modeled)
  f_05keV = 0.25 * f_2keV               ; frequency [Hz] at 0.5 keV
  f_10keV = 5.0 * f_2keV                ; frequency [Hz] at 10 keV 
  if ((f_05keV LE f) AND (f_10keV GE f)) then DL_DF = dL_df_2keV * ((f/f_2keV)^(1 - Gamma))
  if (f GT f_10keV) then DL_DF = 0.0    ; Sorry, not incorporated yet! (need Compton modeling)
  dL_df_05keV = dL_df_2keV * ((f_05keV/f_2keV)^(1 - Gamma))


  ; Using Marconi alpha_OX (Zamorani 82, emperical Vignali, Brandt, & Scheider 03)
  ;   which relates flux at 2500 Angstroms to that at 2keV  by:
  ;    a_OX = - LOG(L_2500 / L_2keV) / LOG(f_2500 / f_2keV)   and
  ;    a_OX = - 0.11*LOG(L_2500) + 1.85
  ;  (from Vignali, Brandt, & Schneider 2003, 
  ;     which needs everything in units of ergs/s/Hz
  R = -2.605    ; [defined as log(nu(2500 Ang) / nu(2 keV)) ]
  logL_solar_ergs = 33.591              ; Log of solar luminosity in ergs/s
  log_dLdf2kev_ergs = alog10(dL_df_2kev) + logL_solar_ergs 
  logdL_df_2500_ergs = (1.85*R + log_dLdf2kev_ergs)/(1.0 + 0.11*R)
  logdL_df_2500 = logdL_df_2500_ergs - logL_solar_ergs
  dL_df_2500 = 10.0^logdL_df_2500


  ; Now use the SDSS, etc. average optical power law L ~ f^(-0.44)
  ;   from 1micrometer to 1300 Angstroms  (Vanden Berk et al. 2001)
  f_micrometer = 2.998d14               ; frequency in [Hz] at 1 micrometer
  f_2500ang    = 1.199d15               ; frequency in [Hz] at 2500 Angstroms
  f_1300ang    = 2.306d15               ; frequency in [Hz] at 1300 Angstroms
  alpha = -0.44                                 ; spectral index in this range
  if ((f_micrometer LE f) AND (f_1300ang GE f)) then DL_DF = dL_df_2500 * ((f/f_2500ang)^(alpha))
  dL_df_micrometer = dL_df_2500 * ((f_micrometer/f_2500ang)^(alpha))
  dL_df_1300 = dL_df_2500 * ((f_1300ang/f_2500ang)^(alpha))



  ; For wavelengths larger than a micrometer, assume the spectrum cuts off with 
  ;   alpha = 2  (as in Rayleigh-Jeans tail of blackbody)
  alpha = 2.0
  if (f LT f_micrometer) then DL_DF = dL_df_micrometer * ((f/f_micrometer)^(alpha))


  ; Assume flat over small break 1300-1200 Angstroms
  f_1300ang    = 2.306d15               ; frequency in [Hz] at 1300 Angstroms
  f_1200ang    = 2.498d15               ; frequency in [Hz] at 1200 Angstroms
  if ((f_1300ang LE f) AND (f_1200ang GE f)) then DL_DF = dL_df_1300
  dL_df_1200 = dL_df_1300


  ; From 1200 - 500 Angstroms we have alpha = -1.76 
  ;     (Telfer et al. 2002, Vanden Berk et al. 2001)
  f_500ang    = 5.996d15                ; frequency in [Hz] at 500 Angstroms
  alpha = -1.76
  if ((f_1200ang LE f) AND (f_500ang GE f)) then DL_DF = dL_df_1200 * ((f/f_1200ang)^(alpha))
  dL_df_500 = dL_df_1200 * ((f_500ang/f_1200ang)^(alpha))


  ; Now, finally, we use the values of dL_df at 500 Angstroms and 0.5 keV
  ;   to create a power law connecting the two frequencies and project over it
  alpha = alog10(dL_df_05keV/dL_df_500) / alog10(f_05keV/f_500ang) 
  if ((f_500ang LE f) AND (f_05keV GE f)) then DL_DF = dL_df_500 * ((f/f_500ang)^(alpha))
  
  ;print, alog10(f_micrometer*dL_df_micrometer)
  ;print, alog10(f_2500ang*dL_df_2500)
  ;print, alog10(f_1300ang*dL_df_1300)
  ;print, alog10(f_1200ang*dL_df_1200)
  ;print, alog10(f_500ang*dL_df_500)
  ;print, alog10(f_05keV*dL_df_05keV)
  ;print, alog10(f_2keV*dL_df_2keV)
  ;checked out, matches well with Marconi et al. (2004) model spectrum
  
  return, DL_DF
end
    





; -------------------------------------------------------------------------------------------
; -------------------------------------------------------------------------------------------





;  Compile the X-ray emission, by component,
;   for one run and write it to a file.
;-------------------------------------------
pro compile_all_xrayinfo, frun, nobh=nobh


if not keyword_set(frun) then begin
	print, " "
	print, " compile_all_xrayinfo, frun, nobh=nobh"
	print, " "
	print, " "
	return
endif


; ----------------------

read_file_xrays, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, $
				xray, xray_sf, $
				xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

ntime= n_elements(time)

xray_sfism= xray_sf
xray_diff_brem= xray
xray_diff_soft= xray_rs_s
xray_diff_hard= xray_rs_h


if not keyword_set(nobh) then begin

   ; now include the AGN component
   ; -----------------------------
   open_blackhole_txt, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, bh_totalmass, bh_mdot_edd

   ; ------------------
   ; cgs units, cm, sec, g                        
   L_solar= 3.9d+33                               ; in ergs per sec
   cc= 2.9979d+10                                 ; speed of light in cm/sec
   convert_sunyr_gsec= 6.30428d+25                ; convert msun/yr -> g/sec

   boloL_ergss= 0.1*bh_mdot_sunyr*cc*cc*convert_sunyr_gsec
   boloL_sun= boloL_ergss / L_solar
   print, "L_sun    max/min ", max(boloL_sun), min(boloL_sun)
   print, "L_erg s  max/min ", max(boloL_ergss), min(boloL_ergss)

   ;bolo_corr_f= BH_hardXlum_inlog(ar_sunyr)       ; L/L_hard
   bolo_corr_f= BH_hardXlum_inlog(boloL_sun)       ; L/L_hard
   ;print, "f_hard= ",bolo_corr_f
   L_hard= boloL_ergss / bolo_corr_f
   print, "Lx_hard  max/min ", max(L_hard), min(L_hard)
   L_hard= alog10(L_hard)

   ;bolo_corr_f= BH_softXlum_inlog(ar_sunyr)       ; L/L_soft
   bolo_corr_f= BH_softXlum_inlog(boloL_sun)       ; L/L_soft
   ;print, "f_soft= ",bolo_corr_f
   L_soft= boloL_ergss / bolo_corr_f
   print, "Lx_soft  max/min ", max(L_soft), min(L_soft)
   L_soft= alog10(L_soft)


   ;-------------
   xray_agn_soft=    fltarr(ntime)
   xray_agn_hard=    fltarr(ntime)

   for i= 0, ntime-1 do begin

	timeidx= where(bhtime ge time[i])

	xray_agn_soft[i]= L_soft(timeidx(0))
	xray_agn_hard[i]= L_hard(timeidx(0))
   endfor

endif else begin

   xray_agn_soft=    fltarr(ntime)
   xray_agn_hard=    fltarr(ntime)

   xray_agn_soft(*)= 0.0
   xray_agn_hard(*)= 0.0

endelse



; ---------------------


; Xrays FILE
; -----------
print, " "
print, " writing file:  xrays_bycomponent.txt"
print, " "

openw, 1, frun+'/xrays_bycomponent.txt', ERROR=err
        
printf, 1, "#   xrays_bycomponent.txt"
printf, 1, "# "
printf, 1, "#             sfism   diff_brem diff_soft diff_hard  agn_soft  agn_hard "
printf, 1, "#   time    log(xray) log(xray) log(xray) log(xray) log(xray) log(xray) "               
printf, 1, "# (Gyr/h)    (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s)  (ergs/s) "
for i=0,ntime-1 do begin
	printf, 1, FORMAT= '(F8.3,"    ", 6(F8.3,"  "))', $
                time[i], xray_sfism[i], xray_diff_brem[i], xray_diff_soft[i], xray_diff_hard[i], xray_agn_soft[i], xray_agn_hard[i]

endfor

close, 1

end




; -------------------------------------------------------------------------------------------------






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
read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

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
read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

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
;read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
;;read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
;oplot, time/0.7, xray, thick=3.0, psym=-6, color= 220
;oplot, [2.2, 2.4], [41.08,41.08], thick=3.0, psym=-6, color= 150
;xyouts, 0.82, 0.77, '1e5', /normal, charthick=3.0, size=1.33, color=220


; Galaxy 4
;----------
;frun="/raid4/tcox/As/A4"
;frun="/raid4/tcox/vc4vc4a"
;frun="/raid4/tcox/vc3avc3a_4"
;frun="/raid4/tcox/vc3ba_2_4"
;read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
;;read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
;oplot, time/0.7, xray, thick=3.0, psym=-6, color= 100
;oplot, [2.2, 2.4], [41.08,41.08], thick=3.0, psym=-6, color= 150
;xyouts, 0.82, 0.72, '1e6', /normal, charthick=3.0, size=1.33, color=100


; Galaxy 5
;----------
;frun="/raid4/tcox/As/A5"
;frun="/raid4/tcox/vc5vc5a"
;frun="/raid4/tcox/vc3avc3a_6"
;frun="/raid4/tcox/vc3ba_2_5"
;read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
;;read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
;oplot, time/0.7, xray, thick=3.0, psym=-6, color= 0
;oplot, [2.2, 2.4], [41.08,41.08], thick=3.0, psym=-6, color= 150
;xyouts, 0.82, 0.67, '1e7', /normal, charthick=3.0, size=1.33, color=0


; Galaxy 6
;----------
;frun="/raid4/tcox/vc6vc6a"
;read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
;;read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h
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
read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

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
read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

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

read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

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
read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

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
read_data_osullivan_03, loglx, loglb, tempx
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
read_data_osullivan_01, loglb, loglx, ttype
symsize= 0.2
usersym,symsize*cos(findgen(49)/49*2*!pi),symsize*sin(findgen(49)/49*2*!pi), thick=4.0, /fill
oplot, loglb, loglx, psym=8, color= 0
;oplot, loglb, loglx, psym=7, color= 0, symsize= 0.5

read_data_osullivan_03, loglx, loglb, tempx
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
read_data_osullivan_03, loglx, loglb, tempx
read_data_osullivan_03b, sigma
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

    ;read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
    read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h

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

    ;read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
    read_file_xrays, frun, time, temp_keV_X, z_X, entropy, xray, xray_sf, xray_rs_s, xray_rs_h, xray_rs0_s, xray_rs0_h


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
;  Use Raymond & Smith ('77) to determine
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

        S = CALL_EXTERNAL('/n/home/tcox/Tools/C-Routines_for_IDL/RaymondSmith/raymond_smith', $
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
;========================================





pro sansome, junk

frun="/raid4/tcox/vc3vc3h_2"

;read_file_hotgas, frun, time, temp_keV_X, temp_keV, temp_K, entropy, gas_tot, gas_hot, gas_cold, gas_sf
read_file_xrays, frun, time, temp_keV_X, z_X, mass_xraygas, entropy, $
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




