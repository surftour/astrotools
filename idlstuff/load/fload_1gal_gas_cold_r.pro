function fload_1gal_gas_cold_r, startid, numpart, center=center

    COMMON GalaxyHeader
    COMMON GasData
    COMMON Center
    COMMON OtherData
    COMMON FeedbackData


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
        if flag_feedbacktp then begin
                rest_u= g1gasu(rest_idx)+g1gastpu(rest_idx)
        endif else begin
                rest_u= g1gasu(rest_idx)
        endelse


        if((count_r+count_sf) NE n_elements(g1gasmass)) then begin
                print, "  "
                print, "PROBLEM: not all gas particles are above or below rho cutoff"
                print, "  "
                return, 0
        endif

        cold_idx= where((rest_u) LT 150.0)

	if cold_idx(0) ne -1 then begin
	   x= g1gasx(cold_idx)
	   y= g1gasy(cold_idx)
	   z= g1gasz(cold_idx)
	   return, sqrt((x-c(0))*(x-c(0)) + (y-c(1))*(y-c(1)) + (z-c(2))*(z-c(2)))
	endif else begin
	   return, [0]
	endelse
    endif


end


