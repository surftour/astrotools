function fload_1gal_disk_xyz, dummy, startid, numpart, center=center

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON DiskData
    COMMON Center


    if numpart eq 0 then return, [0]

    if keyword_set(center) then begin
	c= center
    endif else begin
	c= com
    endelse

    ;gid = id(0:npart(0)-1)
    gid = fload_disk_id(1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))


    if dummy EQ 'x' then begin
	g1xdisk = xdisk(idx)
	return, g1xdisk-c(0)
    endif

    if dummy EQ 'y' then begin
	g1ydisk = ydisk(idx)
	return, g1ydisk-c(1)
    endif

    if dummy EQ 'z' then begin
	g1zdisk = zdisk(idx)
	return, g1zdisk-c(2)
    endif

    if dummy EQ 'phi' then begin
	g1xdisk = xdisk(idx)-c(0)
        g1ydisk = ydisk(idx)-c(1)
        phi = atan(g1ydisk,g1xdisk)
        return, phi
    endif

    if dummy EQ 'theta' then begin
	g1xdisk = xdisk(idx)-c(0)
        g1ydisk = ydisk(idx)-c(1)
        g1zdisk = zdisk(idx)-c(2)
        theta= atan(sqrt(g1xdisk*g1xdisk + g1ydisk*g1ydisk),g1zdisk)
        return, theta
    endif

    if dummy EQ 'rxy' then begin
        g1xdisk = xdisk(idx)-c(0)
        g1ydisk = ydisk(idx)-c(1)
        rxy= sqrt(g1xdisk*g1xdisk + g1ydisk*g1ydisk)
        return, rxy
    endif

    if dummy EQ 'r' then begin
	g1xdisk = xdisk(idx)-c(0)
	g1ydisk = ydisk(idx)-c(1)
	g1zdisk = zdisk(idx)-c(2)
        r= sqrt(g1xdisk*g1xdisk + g1ydisk*g1ydisk + g1zdisk*g1zdisk)
        return, r
    endif

end


