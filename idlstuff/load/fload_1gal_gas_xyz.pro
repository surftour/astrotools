function fload_1gal_gas_xyz, dummy, startid, numpart, center=center

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON GasData
    COMMON Center


    if numpart eq 0 then return, [0]

    if keyword_set(center) then begin
	c= center
    endif else begin
	c= com
    endelse

    ;gid = id(0:npart(0)-1)
    gid = fload_gas_id(1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))


    if dummy EQ 'x' then begin
	g1xgas = xgas(idx)
	return, g1xgas-c(0)
    endif

    if dummy EQ 'y' then begin
	g1ygas = ygas(idx)
	return, g1ygas-c(1)
    endif

    if dummy EQ 'z' then begin
	g1zgas = zgas(idx)
	return, g1zgas-c(2)
    endif

    if dummy EQ 'phi' then begin
	g1xgas = xgas(idx)-c(0)
        g1ygas = ygas(idx)-c(1)
        phi = atan(g1ygas,g1xgas)
        return, phi
    endif

    if dummy EQ 'theta' then begin
	g1xgas = xgas(idx)-c(0)
        g1ygas = ygas(idx)-c(1)
        g1zgas = zgas(idx)-c(2)
        theta= atan(sqrt(g1xgas*g1xgas + g1ygas*g1ygas),g1zgas)
        return, theta
    endif

    if dummy EQ 'rxy' then begin
        g1xgas = xgas(idx)-c(0)
        g1ygas = ygas(idx)-c(1)
        rxy= sqrt(g1xgas*g1xgas + g1ygas*g1ygas)
        return, rxy
    endif

    if dummy EQ 'r' then begin
	g1xgas = xgas(idx)-c(0)
	g1ygas = ygas(idx)-c(1)
	g1zgas = zgas(idx)-c(2)
        r= sqrt(g1xgas*g1xgas + g1ygas*g1ygas + g1zgas*g1zgas)
        return, r
    endif

end


