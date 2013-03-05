function fload_gas_rho, dummy, idtofollow=idtofollow, $
		cold=cold, hot=hot, TurbulentFactor=TurbulentFactor, $
		rho_crit=rho_crit

    COMMON GalaxyHeader
    COMMON GasData
    COMMON OtherData


    if not keyword_set(TurbulentFactor) then TurbulentFactor= 1.0

    ; calculate multiphase density?
    ;--------------------------------
    if keyword_set(cold) or keyword_set(hot) then begin

	;  Contants
	; ---------------------------
	HYDROGEN_MASSFRAC= 0.76
	GAMMA_MINUS1= (5./3.) - 1.0
	BOLTZMANN=   1.3806d-16
	PROTONMASS=  1.6726d-24
	UnitMass_in_g= 1.989d+43
	UnitDensity_in_cgs = 6.76991d-22 
	UnitEnergy_in_cgs = 1.989d+53

	;  Cold Phase Temperature
        ; ---------------------------
	TempClouds    = 1000.0           ; in K
	meanweight = 4. / (1. + 3. * HYDROGEN_MASSFRAC)
	EgySpecCold = 1. / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * TempClouds
	u_cold = EgySpecCold * UnitMass_in_g / UnitEnergy_in_cgs   ;  now in gadget units 


	; Critical Overdensity for SF
	; ----------------------------
	PhysDensThresh= 0.000854924
	if keyword_set(rho_crit) then PhysDensThresh= rho_crit


	;  
	; ---------------------------
	Ngas= long(npart(0))
	simtime = float(time)
	hubbleparam= 0.7

	simrho= fload_gas_rho(1)
	simne= fload_gas_nume(1)
	simu= fload_gas_u(1)

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


	; this is the x from SH03  (eq. 17, rho_cold/rho)
	x= cold_mass_fraction
	;print, "x (mass fraction)= ", x[10318:10320]

	; fix cold "temperature"
	u_cold= u*0.0 + u_cold*TurbulentFactor
	idx= where(u_cold gt 0.75*u)
	if idx(0) ne -1 then u_cold(idx)= 0.75*u(idx)

	; we assume that the total energy (u) is a 
	; mass weighted average, so
	;  rho*u =  rho_hot*u_hot + rho_cold*u_cold
	;or
	; u = (rho_hot/rho)u_hot  +  (rho_cold/rho)u_cold
	; u = (1-x)u_hot  +  (x)u_cold
	;u_hot = (u - x*u_cold) / (1.0-x)
	u_hot = (u - x*u_cold) / (1.0-x)

	; ********
	; not sure if this is correct, but
	; add this to trap for all cold particles
	idx= where(u_hot lt 1.1*u_cold) 
        if idx(0) ne -1 then u_hot(idx) = 1.1*u_cold
	; ********

	density_ratio = u_cold / u_hot
	;print, "rho_hot / rho_cold = ", density_ratio[10318:10320]

	; this is volume_cold / volume_hot
	volume_filling_factor= (x/(1.0-x)) * density_ratio
	;print, "cold volume_filling_factor= ", volume_filling_factor[10318:10320]

	; cold_volume_fraction is
	;   volume_cold / volume_total
	;print, "cold_volume_fraction= ", cold_volume_fraction[10318:10320]


	; what's original rho?
	;print, "rho= ", rho[10318:10320]

	; ***  tj's method ***
	;rho_cold= rho * u / (2.0 * u_cold)
	rho_cold= rho * u /  u_cold
	rho_hot = rho_cold * density_ratio
	;print, "---- TJ's methods ----"
	;print, "rho_cold= ", rho_cold[10318:10320]
	;print, "rho_hot= ", rho_hot[10318:10320]
	;print, " "

	; ***  phil's method  ***
	xc0= (1.0 + volume_filling_factor) * x / (1.e-10 + volume_filling_factor)
	xh0= xc0 * density_ratio
	rho_cold= rho *  xc0
	rho_hot= rho * xh0
	;print, "---- Phil's methods ----"
	;print, "rho_cold= ", rho_cold[10318:10320]
	;print, "rho_hot= ", rho_hot[10318:10320]
	;print, " "

	;print, "Pressures:"
	;print, "tot:   ",  rho[10318:10320]*u[10318:10320]
	;print, "cold:  ",  rho_cold[10318:10320]*u_cold
	;print, "hot:   ",  rho_hot[10318:10320]*u_hot[10318:10320]

	idx= where(rho lt PhysDensThresh)
	if idx(0) ne -1 then rho_hot(idx)= rho(idx)
	if idx(0) ne -1 then rho_cold(idx)= 0.0

	if keyword_set(cold) then return, rho_cold
	if keyword_set(hot) then return, rho_hot

    endif






    ;--------------------------------
    if dummy EQ 1 then return, rho

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

                return, rho(idx)
        endif
    endif


end


