function fload_newstars_z, dummy, idtofollow=idtofollow


    COMMON GalaxyHeader
    COMMON OtherData
    COMMON SfrData


    if dummy EQ 1 then begin
	if npart(4) EQ 0 then return, 0
	return, smetals
    endif


    if keyword_set(idtofollow) then begin
        if dummy EQ 99 then begin

        	sid= npart(0)+npart(1)+npart(2)+npart(3)
        	eid= npart(0)+npart(1)+npart(2)+npart(3)+npart(4)-1
        	stars_ids= id(sid:eid)
                idx= where(stars_ids EQ idtofollow)
                if idx(0) EQ -1 or n_elements(idx) NE 1 then begin
                        print, "Hey, what's up with idx?  printing"
                        print, idx
                        return, [0,0,0]
                endif

                return, smetals(idx)
        endif
    endif

end


