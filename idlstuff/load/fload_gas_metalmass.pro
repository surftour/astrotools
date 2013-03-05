function fload_gas_metallicity, dummy, idtofollow=idtofollow

    COMMON GalaxyHeader
    COMMON SfrData
    COMMON OtherData

    if flag_sfr then begin
	if dummy EQ 1 then begin
	   if flag_metals gt 0 then begin
		return, gmetals 
	   endif else begin
		return, 0
	   endelse
	endif

    	if keyword_set(idtofollow) then begin
            if dummy EQ 99 then begin

                gas_ids= id(0:npart(0)-1)
                idx= where(gas_ids EQ idtofollow)
                if idx(0) EQ -1 or n_elements(idx) NE 1 then begin
                        print, "Hey, what's up with idx?  printing"
                        print, idx
                        return, [0,0,0]
                endif

                return, gmetals(idx)
            endif
    	endif
    endif else begin
	return, 0
    endelse

end


