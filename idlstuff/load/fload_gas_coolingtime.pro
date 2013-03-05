function fload_gas_coolingtime, Metallicity, averageit=averageit

    ; returns entropy in keV cm2

    ;COMMON GalaxyHeader
    ;COMMON GasData
    ;COMMON FeedbackData

    g_gamma= 5.0/3.0
    g_minus_1= g_gamma-1.0
    m_p= 1.6726d-24         ; proton mass (in cgs, i.e. g)
    k_b= 1.3806e-16         ; boltzmann constant (in cgs)
    t_cgs= 3.08568e+16          ; convert time to seconds (cgs units)
    UnitEnergy_in_cgs=        1.989d+53
    UnitMass_in_g    =        1.989d+43
    UnitDensity_in_cgs =      6.7699112d-22


    ; load gas quantities
    ; --------------------
    u= fload_gas_u(1)
    rho= fload_gas_rho(1)


    ; determine cooling rate from sd files
    ; -------------------------------------
    ;ok=load_cooling_table('zero')
    ;ok=fload_cooling_table('-10')           ; load cooling file
    ok=fload_cooling_table(Metallicity)           ; load cooling file
    ;temp= u*fload_gas_mu(1)/0.012381322    ; temp in Kelvin
    temp= fload_gas_temperature(1)         ; temp in Kelvin
    coolrate= fload_lambdan_array(temp)    ; cooling rate (ergs cm3 s-1)
    zeroidx= where(coolrate le 0)
    if zeroidx(0) ne -1 then coolrate(zeroidx)= 1.0


    ; convert to physical units
    ; --------------------------
    rho_cgs= rho*UnitDensity_in_cgs         ; g/cm3
    u_cgs= u*UnitEnergy_in_cgs/UnitMass_in_g  ; in ergs/g per unit volume


    cooltime= u_cgs*m_p*m_p/(rho_cgs*coolrate)     ; cooling time (s)
    cooltime= cooltime/3.08568d+16         ; cooling time (Gadget Units - Gyr)
    if zeroidx(0) ne -1 then cooltime(zeroidx)= 0.0


    ; average it
    if keyword_set(averageit) then begin
	ct_moment= moment(cooltime)
	return, ct_moment[0]
    endif
    
    return, cooltime

end


