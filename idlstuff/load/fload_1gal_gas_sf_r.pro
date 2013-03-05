function fload_1gal_gas_sf_r, startid, numpart, center=center

    COMMON GalaxyHeader
    COMMON GasData
    COMMON FeedbackData
    COMMON OtherData
    COMMON Center


    if numpart eq 0 then return, [0]

    if keyword_set(center) then begin
        c= center
    endif else begin
        c= com
    endelse


    sfr_rhocutoff= 0.00170994      ; gadget units
    cold_cutoff= 150.0      ; gadget units = 1.2 e4 T(K)

    ; get just this galaxies information
    gid = id(0:npart(0)-1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))
    g1gasrho = rho(idx)
    g1gasu = u(idx)
    g1gastpu = tpu(idx)
    g1gasmass = mgas(idx)
    g1gasx = xgas(idx)
    g1gasy = ygas(idx)
    g1gasz = zgas(idx)


    if numpart GT 0 then begin

        sf_idx= where(g1gasrho GE sfr_rhocutoff, count_sf)

        rest_idx= where(g1gasrho LT sfr_rhocutoff, count_r)

        if((count_r+count_sf) NE n_elements(g1gasmass)) then begin
                print, "  "
                print, "PROBLEM: not all gas particles are above or below rho cutoff"
                print, "  "
                return, 0
        endif



	if sf_idx(0) GE 0 then begin
		x= g1gasx(sf_idx)
		y= g1gasy(sf_idx)
		z= g1gasz(sf_idx)
		return, sqrt((x-c(0))*(x-c(0)) + (y-c(1))*(y-c(1)) + (z-c(2))*(z-c(2)))
	endif else begin
		return, [-1]
	endelse
    endif


end


