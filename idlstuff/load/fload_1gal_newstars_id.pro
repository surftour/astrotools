function fload_1gal_newstars_id, dummy, startid, numpart


    COMMON GalaxyHeader
    COMMON OtherData


    if dummy NE -1 then begin
        if npart(4) EQ 0 then return, 0

        sid= npart(0)+npart(1)+npart(2)+npart(3)
        eid= npart(0)+npart(1)+npart(2)+npart(3)+npart(4)-1
        stars_ids= id(sid:eid)

	; recover same ID as original gas particle
	;print, "HEY !!!!!!  WARNING:  turned off new star ID correction"
	whacky_ids= where((stars_ids lt 0) or (stars_ids gt 2.0e6))
	num= long(1)
	if whacky_ids(0) ne -1 then stars_ids(whacky_ids)= stars_ids(whacky_ids) - ishft(num,31)
	;
	; ishft(num,31) comes from the bit flip operation volker
	; does in sfr_eff.c, essentially one of the left-most (?)
	; bits is flipped to one.  here we subtract that off
	;
	;

	idx = where((stars_ids ge long(startid)) and (stars_ids lt long(startid)+long(numpart)))
        return, stars_ids(idx)

    endif


end


