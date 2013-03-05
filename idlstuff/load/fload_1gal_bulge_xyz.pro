function fload_1gal_bulge_xyz, dummy, startid, numpart, center=center

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON BulgeData
    COMMON Center


    if numpart eq 0 then return, [0]
    if npart(3) le 0 then return, [0]

    if keyword_set(center) then begin
	c= center
    endif else begin
	c= com
    endelse

    ;gid = id(0:npart(0)-1)
    gid = fload_bulge_id(1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))


    if dummy EQ 'x' then begin
	g1xbulge = xbulge(idx)
	return, g1xbulge-c(0)
    endif

    if dummy EQ 'y' then begin
	g1ybulge = ybulge(idx)
	return, g1ybulge-c(1)
    endif

    if dummy EQ 'z' then begin
	g1zbulge = zbulge(idx)
	return, g1zbulge-c(2)
    endif

    if dummy EQ 'phi' then begin
	g1xbulge = xbulge(idx)-c(0)
        g1ybulge = ybulge(idx)-c(1)
        phi = atan(g1ybulge,g1xbulge)
        return, phi
    endif

    if dummy EQ 'theta' then begin
	g1xbulge = xbulge(idx)-c(0)
        g1ybulge = ybulge(idx)-c(1)
        g1zbulge = zbulge(idx)-c(2)
        theta= atan(sqrt(g1xbulge*g1xbulge + g1ybulge*g1ybulge),g1zbulge)
        return, theta
    endif

    if dummy EQ 'rxy' then begin
        g1xbulge = xbulge(idx)-c(0)
        g1ybulge = ybulge(idx)-c(1)
        rxy= sqrt(g1xbulge*g1xbulge + g1ybulge*g1ybulge)
        return, rxy
    endif

    if dummy EQ 'r' then begin
	g1xbulge = xbulge(idx)-c(0)
	g1ybulge = ybulge(idx)-c(1)
	g1zbulge = zbulge(idx)-c(2)
        r= sqrt(g1xbulge*g1xbulge + g1ybulge*g1ybulge + g1zbulge*g1zbulge)
        return, r
    endif

end


