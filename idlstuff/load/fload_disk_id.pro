function fload_disk_id, dummy

    COMMON GalaxyHeader
    COMMON OtherData

    if dummy NE -1 then begin

        if npart(2) le 0 then return, 0

        sid= npart(0)+npart(1)
        eid= npart(0)+npart(1)+npart(2)-1

	return, id(sid:eid)

    endif

end



