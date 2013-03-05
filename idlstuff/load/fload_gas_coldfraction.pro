function fload_gas_coldfraction, dummy, idtofollow=idtofollow

    COMMON GalaxyHeader
    COMMON GasData
    COMMON OtherData
    COMMON CoolingData


    if(npart(0) le 0) then return, [0]

    Ngas= npart(0)
    simtime = float(time)
    hubbleparam= 0.7;
    
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
                TSN, $
                fEVP, $
                fSN, $
                hubbleparam, $
                simtime, $
                comoving, $
                simrho, $
                simne, $
                cold_mass_fraction, $
                cold_volume_fraction, $
		hot_temperature)



    ; grab a specific id numbers xyz
    if keyword_set(idtofollow) then begin
        if dummy NE 1 then begin

                gas_ids= id(0:npart(0)-1)
                idx= where(gas_ids EQ idtofollow)
                if idx(0) LE 0 or n_elements(idx) NE 1 then begin
                        print, "Hey, what's up with idx?  printing"
                        print, idx
                        return, 0.0
                endif

                return, cold_mass_fraction(idx)
        endif
    endif


    if dummy EQ 1 then return, cold_mass_fraction


end


