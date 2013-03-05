function fload_gas_mu, dummy, idtofollow=idtofollow


    COMMON GalaxyHeader
    COMMON CoolingData
    COMMON OtherData

    XH= 0.76
    yhelium= (1-XH)/(4*XH)     ; 0.0789474

    if flag_cooling eq 0 then return, 1.0

    if dummy EQ 1 then begin
	mu= (1+4*yhelium)/(1+yhelium+nume)
	return, mu
    endif

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

		mu= (1+4*yhelium)/(1+yhelium+nume(idx))
                return, mu
        endif
    endif
end


