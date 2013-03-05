function fload_1gal_gas_v, dummy, startid, numpart, $
			center=center, com_velocity=com_velocity

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

    if keyword_set(com_velocity) then begin
        cv= com_velocity
    endif else begin
        cv=[0.0 , 0.0 , 0.0]
    endelse


    ;gid = id(0:npart(0)-1)
    gid = fload_gas_id(1)
    idx = where((gid ge long(startid)) and (gid lt (long(startid)+long(numpart))))

    if (dummy EQ 'mag') or (dummy EQ 'tot') then begin
	g1vgas2 = (vxgas(idx)-cv(0))*(vxgas(idx)-cv(0)) + (vygas(idx)-cv(1))*(vygas(idx)-cv(1)) + (vzgas(idx)-cv(2))*(vzgas(idx)-cv(2))
	g1vgas = sqrt(g1vgas2)
	return, g1vgas
    endif

    if dummy EQ 'x' then begin
	g1vxgas = vxgas(idx)-cv(0)
	return, g1vxgas
    endif

    if dummy EQ 'y' then begin
	g1vygas = vygas(idx)-cv(1)
	return, g1vygas
    endif

    if dummy EQ 'z' then begin
	g1vzgas = vzgas(idx)-cv(2)
	return, g1vzgas
    endif

    if dummy EQ 'phi' then begin
        g1xgas = xgas(idx)-c(0)
        g1ygas = ygas(idx)-c(1)
	g1vxgas = vxgas(idx)-cv(0)
	g1vygas = vygas(idx)-cv(1)
        phi = atan(g1ygas,g1xgas)
        v_phi= -g1vxgas*sin(phi) + g1vygas*cos(phi)
        return, v_phi
    endif

    if dummy EQ 'theta' then begin
        g1xgas = xgas(idx)-c(0)
        g1ygas = ygas(idx)-c(1)
        g1zgas = zgas(idx)-c(2)
	g1vxgas = vxgas(idx)-cv(0)
	g1vygas = vygas(idx)-cv(1)
	g1vzgas = vzgas(idx)-cv(2)
        phi= atan(g1ygas,g1xgas)
	theta= atan(sqrt(g1xgas*g1xgas + g1ygas*g1ygas),g1zgas)
        v_theta = cos(phi)*cos(theta)*g1vxgas + cos(theta)*sin(phi)*g1vygas + cos(theta)*g1vzgas
        return, v_theta
    endif

    if dummy EQ 'r' then begin
        g1xgas = xgas(idx)-c(0)
        g1ygas = ygas(idx)-c(1)
        g1zgas = zgas(idx)-c(2)
        g1vxgas = vxgas(idx)-cv(0)
        g1vygas = vygas(idx)-cv(1)
        g1vzgas = vzgas(idx)-cv(2)
        phi= atan(g1ygas,g1xgas)
        theta= atan(sqrt(g1xgas*g1xgas + g1ygas*g1ygas),g1zgas)
        v_r = cos(phi)*sin(theta)*g1vxgas + sin(theta)*sin(phi)*g1vygas + cos(theta)*g1vzgas
        return, v_r
    endif

end


