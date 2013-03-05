function fload_gas_entropy, dummy, averageit=averageit

    ; returns entropy in keV cm2

    COMMON GalaxyHeader
    COMMON GasData
    COMMON FeedbackData

    g_gamma= 5.0/3.0
    g_minus_1= g_gamma-1.0
    PROTONMASS = 1.6726e-24

    mu= fload_gas_mu(1)
    MeanWeight= mu * PROTONMASS

    ;keV= fload_gas_temperature_kT(1,/keV)
    keV= fload_gas_temperature(1,/keV)

    rho=fload_gas_rho(1)
    density= rho * 6.76991e-22    ; convert to g/cm3
    density= density / MeanWeight  ; convert to cm-3 (i.e. take out g)


    ; actually calculate the entropy
    ;s= keV / (rho^(g_minus_1))
    s= keV / (density^(g_minus_1))


    ; average it
    if keyword_set(averageit) then begin
	s_moment= moment(s)
	return, s[0]
    endif
    
    return, s

end


