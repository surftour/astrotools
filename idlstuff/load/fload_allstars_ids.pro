function fload_allstars_ids, dummy

    COMMON GalaxyHeader
    COMMON OtherData

    if dummy NE -1 then begin

        if npart(2)+npart(3)+npart(4) le 0 then return, 0

        sid= npart(0)+npart(1)
        eid= npart(0)+npart(1)+npart(2)+npart(3)-1
        if eid gt sid then old_star_ids= id(sid:eid) else old_star_ids= [-1]

        sid= npart(0)+npart(1)+npart(2)+npart(3)
        eid= npart(0)+npart(1)+npart(2)+npart(3)+npart(4)-1
        if eid gt sid then new_star_ids= id(sid:eid) else new_star_ids= [-1]

	; recover same ID as original gas particle
        ;print, "HEY !!!!!!  WARNING:  turned off new star ID correction"
	if new_star_ids(0) ne -1 then begin
        	whacky_ids= where((new_star_ids lt 0) or (new_star_ids gt 1.0e7))
        	num= long(1)
        	if whacky_ids[0] ne -1 then new_star_ids(whacky_ids)= new_star_ids(whacky_ids) - ishft(num,31)

		; for the moment we'll screw with this to
		; make sure we have unique id's for all stars
		;
		;if whacky_ids[0] ne -1 then new_star_ids(whacky_ids)= new_star_ids(whacky_ids) + 3000000L
		;
		; took this off because many of the other codes are selecting
		; particle numbers by the original IDs only, and thus this makes
		; these particle not fall into this range.

		if old_star_ids[0] ne -1 then begin
			star_ids= [old_star_ids,new_star_ids]
		endif else begin
			star_ids= new_star_ids
		endelse
	endif else begin

		; we're screwed if old_star_ids[0] eq -1 here, hopefully
		; the particle number check above catched that

		star_ids= [old_star_ids]
	endelse

        return, star_ids

    endif

end



