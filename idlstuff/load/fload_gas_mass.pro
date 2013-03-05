function fload_gas_mass, dummy, $
		idtofollow=idtofollow, $
		cold=cold, hot=hot, sf=sf, $
		HI=HI


    COMMON GalaxyHeader
    COMMON GasData
    COMMON FeedbackData
    COMMON OtherData


    sfr_rhocutoff= 0.00170994      ; gadget units (??)
    ;sfr_rhocutoff= 0.000854924     ; gadget units (standard value)
    cold_cutoff= 150.0      ; gadget units = 12,000 K


    if dummy EQ 1 and npart(0) gt 0 then begin
	m=mgas
        ;gas_ids= id(0:npart(0)-1)
	gas_ids= fload_gas_id(1)


	; star forming particles
	if keyword_set(sf) then begin
		sf_idx= where(rho GE sfr_rhocutoff, count_sf)
		if sf_idx(0) lt 0 then return, [0]
		return, mgas(sf_idx)
	endif


	; hot particles
	if keyword_set(hot) then begin
                lowrho_idx= where(rho LT sfr_rhocutoff, count_r)

                lowrho_u= u(lowrho_idx)
		lowrho_m= mgas(lowrho_idx)

                hot_idx= where((lowrho_u) GT 150.0)
		if hot_idx(0) lt 0 then return, [0]
                return, lowrho_m(hot_idx)
        endif


	; grab cold particles
	if keyword_set(cold) then begin
		lowrho_idx= where(rho LT sfr_rhocutoff, count_r)

		lowrho_u= u(lowrho_idx)
		lowrho_m= mgas(lowrho_idx)

		cold_idx= where((lowrho_u) LE 150.0)
		if cold_idx(0) lt 0 then return, [0]
		return, lowrho_m(cold_idx)
	endif


	; grab "HI" particles
	if keyword_set(HI) then begin
		HI_m= mgas

		;
		; SF gas is assumed to be 10% HI
		;
		sf_idx= where(rho GE sfr_rhocutoff, count_sf)
		if sf_idx(0) ne -1 then HI_m(sf_idx)= 0.1 * HI_m(sf_idx)


		; 
		; hot, low density gas is given zero mass
		;
		hot_idx= where((rho LT sfr_rhocutoff) and (u GT cold_cutoff), count_hot)
		if hot_idx(0) ne -1 then HI_m(hot_idx)= 0.0


		; 
		; cold, low density gas is all HI
		;
		; do nothing, since this should be the rest.
		return, HI_m
	endif


	; look for specific id
        if keyword_set(idtofollow) then begin
          if dummy EQ 99 then begin
                idx= where(gas_ids EQ idtofollow)
                if idx(0) EQ -1 or n_elements(idx) NE 1 then begin
                        print, "Hey, what's up with idx?  printing"
                        print, idx
                        return, [0,0,0]
                endif

                return, m(idx)
          endif
        endif

	; no other hits
	return, m

    endif else begin
	return, [0]
    endelse


end


