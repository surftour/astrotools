function fload_gas_tpu, dummy, idtofollow=idtofollow

    COMMON GalaxyHeader
    COMMON FeedbackData
    COMMON OtherData


    if dummy EQ 1 and npart(0) gt 0 then begin
	if flag_feedbacktp GT 0 then return, tpu else return, 0.0
    endif else begin
        return, [0]
    endelse



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

                if flag_feedbacktp GT 0 then return, tpu(idx) else return, 0.0
        endif
    endif



end


