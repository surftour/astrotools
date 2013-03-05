function fload_1gal_gas_temperature, dummy, startid, numpart, averageit=averageit, keV=keV
				


    ; returns gas particles temperature in Kelvin
    ;
    ;

    COMMON GalaxyHeader
    COMMON OtherData


    if dummy lt 0 or npart(0) le 0 then return, 0

    g_gamma= 5.0/3.0
    g_minus_1= g_gamma-1.0
    PROTONMASS = 1.6726e-24
    BoltzMann_ergs= 1.3806d-16
    BoltzMann_eV= 8.617d-5
    UnitMass_in_g    =        1.989d43 ;  1.0e10 solar masses
    UnitEnergy_in_cgs=        1.989d53

    ; note gadget units of energy/mass = 1e10 ergs/g,
    ; this comes into the formula below


    ; load gas data
    mu= fload_gas_mu(1)
    u= fload_gas_u(1)

    MeanWeight= mu*PROTONMASS

    Temp= MeanWeight/BoltzMann_ergs * g_minus_1 * u * 1e+10


    ; do we want units of keV  (0.001 factor converts from eV to keV)
    if keyword_set(keV) then begin
	Temp= Temp * BoltzMann_eV * 0.001
    endif

    ; if we just want an average
    if keyword_set(averageit) then begin
	Temp_moment= moment(Temp)
	return, Temp_moment[0]
    endif



    ;gid = id(0:npart(0)-1)
    gid = fload_gas_id(1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))



    ; return on a perparticle basis
    ; -----------------------------
    return, Temp(idx)
	


end


