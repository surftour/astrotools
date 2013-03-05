function fload_blackhole_id, dummy, n_bh=n_bh

    COMMON GalaxyHeader
    COMMON BlackHoleData
    COMMON OtherData


    if dummy ge 0 and npart(5) gt 0 then begin
	startids=npart(0)+npart(1)+npart(2)+npart(3)+npart(4)
	bh_ids= id(startids:startids+npart(5)-1)
	n_bh= n_elements(bh_ids)
	return, long(bh_ids(sort(bh_ids)))
    endif else begin
	return, [-1]
    endelse


end


