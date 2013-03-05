function fload_gas_temperature_multi, dummy, $
			hot=hot, $
			cold=cold

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
    TempK= fload_gas_temperature(1)
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




    ;
    ; worry about the sf multiphase model
    ; ----------------------------------------------


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



	; parses out multiphase temps
	; -------------------------------

	; return hot temperature
	; 
	; this is a bit dubious, and should not be used.
	; it currently does not take into consideration the 
	; actual particle temperature!
	;
	if keyword_set(hot) then begin
		hot_temps= TempK
		idx= where(hot_temperature gt 0.0)   ; selects multiphase parts
		hot_temps(idx)= hot_temperature(idx)
		return, hot_temps
	endif



	; default is to return cold gas temperatures
	cold_temperature= TempK
	idx= where(cold_mass_fraction gt 0.0)        ; selects multiphase parts
	if idx(0) ne -1 then cold_temperature(idx)= 1000.0    ; cold phase is always at 1000K

	return, cold_temperature


end


