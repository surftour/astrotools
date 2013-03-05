function fload_bulge_id, dummy

    COMMON GalaxyHeader
    COMMON OtherData

    if dummy NE -1 then begin

        if npart(3) le 0 then return, 0

        sid= npart(0)+npart(1)+npart(2)
        eid= npart(0)+npart(1)+npart(2)+npart(3)-1

	return, id(sid:eid)

    endif

end



