function fload_halo_id, dummy

    COMMON GalaxyHeader
    COMMON OtherData

    if dummy NE -1 then begin

        if npart(1) le 0 then return, 0

        sid= npart(0)
        eid= npart(0)+npart(1)-1

	return, id(sid:eid)

    endif

end



