function fload_gas_temperature_kT, dummy, $
		averageit=averageit, keV=keV

    COMMON GalaxyHeader

; returns each gas particles temperature in kT, actually the
; units are eV.   there is an option to produce the average
; not the mass weighted average, just the average.  perhaps
; we'll adjust this at some point


    if dummy lt 0 or npart(0) le 0 then return, 0

    g_gamma= 5.0/3.0
    g_minus_1= g_gamma-1.0
    PROTONMASS = 1.6726e-24
    BoltzMann_ergs= 1.3806d-16
    BoltzMann_eV= 8.617d-5
    UnitMass_in_g    =        1.989d43 ;  1.0e10 solar masses
    UnitEnergy_in_cgs=        1.989d53

    ; note gadget units of energy/mass = 1e10 ergs/g


    ; load gas data
    mu= fload_gas_mu(1)
    u= fload_gas_u(1)

    MeanWeight= mu*PROTONMASS

    kT= MeanWeight/BoltzMann_ergs * g_minus_1 * u * 1e+10 * BoltzMann_eV

    if keyword_set(keV) then kT= kT*0.001

    ; if we just want an average
    if keyword_set(averageit) then begin
	kT_moment= moment(kT)
	return, kT_moment[0]
    endif


    ; return on a perparticle basis
    ; -----------------------------
    return, kT
	


end


