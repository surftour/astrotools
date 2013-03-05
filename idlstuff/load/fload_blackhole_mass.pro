function fload_blackhole_mass, dummy, idtofollow=idtofollow

    COMMON GalaxyHeader
    COMMON BlackHoleData
    COMMON OtherData


    if npart(5) le 0 then begin
	return, -1
    endif

    if dummy EQ 1 then return, mbh


    if keyword_set(idtofollow) then begin
        if dummy ne 99 then begin

		startids=npart(0)+npart(1)+npart(2)+npart(3)+npart(4)
                bh_ids= id(startids:startids+npart(5)-1)
                idx= where(bh_ids EQ idtofollow)
		;print, bh_ids
		;print, idtofollow
		;print, idx
                if idx(0) EQ -1 or n_elements(idx) NE 1 then begin
                        print, "Hey, what's up with idx?  printing"
                        print, idx
                        return, [0,0,0]
                endif

                return, mbh(idx)
        endif else begin
		return, -2
	endelse
    endif


end


