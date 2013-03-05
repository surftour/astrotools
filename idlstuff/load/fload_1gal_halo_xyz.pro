function fload_1gal_halo_xyz, dummy, startid, numpart, center=center

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON HaloData
    COMMON Center


    if numpart eq 0 then return, [0]

    if keyword_set(center) then begin
	c= center
    endif else begin
	c= com
    endelse

    ;gid = id(0:npart(0)-1)
    gid = fload_halo_id(1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))


    if dummy EQ 'x' then begin
	g1xhalo = xhalo(idx)
	return, g1xhalo-c(0)
    endif

    if dummy EQ 'y' then begin
	g1yhalo = yhalo(idx)
	return, g1yhalo-c(1)
    endif

    if dummy EQ 'z' then begin
	g1zhalo = zhalo(idx)
	return, g1zhalo-c(2)
    endif

    if dummy EQ 'phi' then begin
	g1xhalo = xhalo(idx)-c(0)
        g1yhalo = yhalo(idx)-c(1)
        phi = atan(g1yhalo,g1xhalo)
        return, phi
    endif

    if dummy EQ 'theta' then begin
	g1xhalo = xhalo(idx)-c(0)
        g1yhalo = yhalo(idx)-c(1)
        g1zhalo = zhalo(idx)-c(2)
        theta= atan(sqrt(g1xhalo*g1xhalo + g1yhalo*g1yhalo),g1zhalo)
        return, theta
    endif

    if dummy EQ 'rxy' then begin
        g1xhalo = xhalo(idx)-c(0)
        g1yhalo = yhalo(idx)-c(1)
        rxy= sqrt(g1xhalo*g1xhalo + g1yhalo*g1yhalo)
        return, rxy
    endif

    if dummy EQ 'r' then begin
	g1xhalo = xhalo(idx)-c(0)
	g1yhalo = yhalo(idx)-c(1)
	g1zhalo = zhalo(idx)-c(2)
        r= sqrt(g1xhalo*g1xhalo + g1yhalo*g1yhalo + g1zhalo*g1zhalo)
        return, r
    endif

end


