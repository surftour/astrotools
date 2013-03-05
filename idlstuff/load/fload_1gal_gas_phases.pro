function fload_1gal_gas_phases, startid, numpart, sfrhothresh=sfrhothresh

    COMMON GalaxyHeader
    COMMON GasData
    COMMON SfrData
    COMMON OtherData
    COMMON FeedbackData

    if keyword_set(sfrhothresh) then begin
	;sfr_rhocutoff= sfrhothresh
	sfr_rhocutoff= sfrhothresh/10.0    ; convert Msolar/pc3 to gadget units
    endif else begin
	sfr_rhocutoff= 0.00170994      ; gadget units
    endelse

    cold_cutoff= 150.0      ; gadget units = 1.2 e4 T(K)

    if numpart eq 0 then return, [0,0,0,0]

    ; get just this galaxies information
    gid = id(0:npart(0)-1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))
    g1gasrho = rho(idx)
    g1gasu = u(idx)
    g1gastpu = tpu(idx)
    g1gasmass = mgas(idx)
    g1gasmfs = mfs(idx)

    if numpart GT 0 then begin

	sf_idx= where(g1gasrho GE sfr_rhocutoff, count_sf)
	if sf_idx(0) ne -1 then begin
	    sf_mgas= g1gasmass(sf_idx)-g1gasmfs(sf_idx)
	endif else begin
	    sf_mgas= 0.0
        endelse

	rest_idx= where(g1gasrho LT sfr_rhocutoff, count_r)
	rest_u= g1gasu(rest_idx)+g1gastpu(rest_idx)

	if((count_r+count_sf) NE n_elements(g1gasmass)) then begin
		print, "  "
		print, "PROBLEM: not all gas particles are above or below rho cutoff"
		print, "  "
		return, 0
	endif

	cold_idx= where((rest_u) LT 150.0, count_cold)
	if cold_idx(0) GE 0 then cold_mgas= g1gasmass(cold_idx)-g1gasmfs(cold_idx) else cold_mgas= [0.0]

	hot_idx= where((rest_u) GE 150, count_hot)
	if hot_idx(0) GE 0 then begin
		hot_mgas= g1gasmass(hot_idx)-g1gasmfs(hot_idx)
		turb_percentage= total(g1gastpu(hot_idx))/total(g1gasu(hot_idx)+g1gastpu(hot_idx))
	endif else begin
		hot_mgas= [0.0]
		turb_percentage= 0.0
	endelse

        if((count_cold+count_hot) NE count_r) then begin
                print, "  "
                print, "PROBLEM: not all gas particles are above or below rho cold cutoff"
                print, "  "
                return, 0
        endif

	;print, "Mtot=",total(g1gasmass-g1gasmfs)
	;print, "Mcold=",total(cold_mgas)
	;print, "Msf=",total(sf_mgas)
	;print, "Mhot=",total(hot_mgas), "(turbulent percentage:",100.0*turb_percentage,")"

	masses= [total(g1gasmass-g1gasmfs), total(cold_mgas), total(sf_mgas), turb_percentage*total(hot_mgas)]

	return, masses
    endif

end




