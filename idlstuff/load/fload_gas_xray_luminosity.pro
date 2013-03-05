function fload_gas_xray_luminosity, dummy, $
			diffuse_hotgas=diffuse_hotgas, $
			sf_hotgas=sf_hotgas, $
			minkeVTemp=minkeVTemp

    ; returns x-ray luminosity in ergs sec-1
    ;
    ; this formula comes from Eke, Navarro & Frenk '98

    COMMON GalaxyHeader

    if dummy lt 0 or npart(0) le 0 then return, 0

    PROTONMASS = 1.6726d-24
    Brem_Normalization= 1.2d-24


    ; load gas data
    mu= fload_gas_mu(1)
    ;keV= fload_gas_temperature_keV(1,/keV)
    keV= fload_gas_temperature(1,/keV)
    m= fload_gas_mass(1)
    m= m*1.989d+43                    ; convert to g
    rho= fload_gas_rho(1)
    nume= fload_gas_nume(1)


    ; ---------------------------
    ;  determine hot gas density
    ; ---------------------------
    Ngas= long(npart(0))
    simtime = float(time)
    hubbleparam= 0.7
    
    simrho= rho
    simne= nume

    ; need these to match the run
    ; ----------------------------
    tstar = 4.5
    tSN   = float(3.0e+8)
    fEVP  = 3000.0
    fSN   = 0.1
    comoving = 0.0

    print, "t_star=    ",tstar
    print, "tSN=       ",tSN
    print, "fEVP=      ",fEVP
    print, "fSN=       ",fSN
    ;print, "coming=    ",comoving

  
    cold_mass_fraction   = fltarr(Ngas)     ;mass fraction
    cold_volume_fraction = fltarr(Ngas)     ;volume fraction
    hot_temperature      = fltarr(Ngas)     ;hot phase temp




     ; do this first, or else the recalculated hot-phase
     ; density if very likely going to satify this and
     ; be set to zero
     ; --------------------------------------------------
     if keyword_set(sf_hotgas) then begin
        ; if we're doing sf_hotgas only, then
        ; set stuff below rho_crit to zero
        idx= where(rho lt 0.000854924)
        if idx(0) ne -1 then begin
           m(idx)= 0.0
        endif
     endif


    ; if we're just doing diffuse_hotgas then don't
    ; worry about the sf multiphase model, just set
    ; all this gas to zero
    ; ----------------------------------------------
    if keyword_set(diffuse_hotgas) then begin
        ; if we're doing diffuse_hotgas, then everything
        ; above rho_crit doesn't contribute to flux
        idx= where(rho ge 0.000854924)
        if idx(0) ne -1 then begin
           m(idx)= 0.0
        endif
    endif else begin
         ; Brant's external code to calculate the cooling and
         ; return the calculated mass fraction
         ; ----------------------------------------------------
	 spawn, 'echo $TJHOME', result
	 homedir= strcompress(result,/remove_all)
	 libfile= homedir+'/Tools/C-Routines_for_IDL/ComputeMultiphase/mphase'
	 S = CALL_EXTERNAL(libfile, $
                'multiphase_info', $
                Ngas, $
                tstar, $
                tSN, $
                fEVP, $
                fSN, $
                hubbleparam, $
                simtime, $
                comoving, $
                simrho, $
                simne, $
                cold_mass_fraction, $
                cold_volume_fraction, $
                hot_temperature )

	idx= where(hot_temperature gt 0.0)
	if idx(0) ne -1 then begin
	   rho(idx)= (1-cold_mass_fraction(idx))*rho(idx)/(1-cold_volume_fraction(idx))
	   keV(idx)= hot_temperature(idx)*8.617d-5*0.001
	   m(idx)= (1-cold_mass_fraction(idx))*m(idx)
	   mu(idx)= 0.6
	endif
    endelse


    ; ditch low density, cold gas
    ; always do this, but it will be
    ; redundant if sf_hotgas is set
    ; --------------------------------
    ; only gas over a certain temperature contributes
    ;    (not clear what the best cutoff is?)
    ;idx= where(keV lt 0.01)    ; 10^5 K
    ;idx= where(keV lt 0.016)    ; 10^5.2 K (primordial plasma is fully ionized above this)
    idx= where((keV lt 0.016) and (rho lt 0.000854924))    ; 10^5.2 K (primordial plasma is fully ionized above this)
    ;idx= where(keV lt 0.1)    ; 10^6 K
    if idx(0) ge 0 then m(idx)= 0.0



    ;
    ; add temperature floor
    ;
    if keyword_set(minkeVTemp) then begin
	idx= where(keV lt minkeVTemp)
	if idx(0) ne -1 then m(idx)= 0.0
    endif


    ; actually, we'll self-consistently calculate
    ; the hot density for star-forming gas
    ; --------------------------------------------
    ;; make sure star forming gas isn't contributing
    ;idx= where(rho gt 0.00085)
    ;if idx(0) ge 0 then begin
    ;	m(idx)= 0.0
    ;endif
    ;idx= where(hot_temperature gt 0.0)
    ;if idx(0) ne -1 then begin
    ;	rho(idx)= (1-cold_mass_fraction(idx))*rho(idx)/(1-cold_volume_fraction(idx))
    ;	keV(idx)= hot_temperature(idx)*8.617d-5*0.001
    ;endif


    MeanWeight= mu * PROTONMASS       ; mu better be ~0.6 for all the contributing gas
    density= rho * 6.76991d-22        ; convert to g/cm3

    ; bolometric luminosity per particle
    xray_lum= Brem_Normalization * (m/MeanWeight) * (density/MeanWeight) * sqrt(keV)

    return, xray_lum
	
end


