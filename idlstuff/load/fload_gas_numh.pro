;
;  The returns the nuetral 
; hydrogen relative to the 
; total hydrogen density
;
;
function fload_gas_numh, dummy, idtofollow=idtofollow

    COMMON GalaxyHeader
    COMMON CoolingData
    COMMON OtherData

    if dummy EQ 1 then return, numh

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

                return, numh(idx)
        endif
    endif


end


