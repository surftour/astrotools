function fload_1gal_allstars_v, dummy, startid, numpart, center=center, vcom=vcom

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON ParentID


    if keyword_set(center) then begin
        c= center
    endif else begin
        c=fload_1gal_center(startid,numpart)
    endelse


    if keyword_set(vcom) then begin
        cv= vcom
    endif else begin
        cv=[0.0 , 0.0 , 0.0]
    endelse


    ; first check for stars
    if (npart(2)+npart(3)+npart(4) le 0) then begin
        print, "PROBLEM"
        return, [0]
    endif


    ; first grab id's for disk and bulge particles
    gid = fload_allstars_ids(1)


    ; then grab appropriate ones for galaxy 1
    idx = where((gid ge startid) and (gid lt startid+numpart))


    if npart(2) gt 0 then begin
		allxs = [xdisk]
		allys = [ydisk]
		allzs = [zdisk]
		allvxs = [vxdisk]
		allvys = [vydisk]
		allvzs = [vzdisk]
    endif

    if npart(3) GT 0 then begin
	if n_elements(allxs) gt 0 then begin
	  allxs = [allxs, xbulge]
	  allys = [allys, ybulge]
	  allzs = [allzs, zbulge]
	  allvxs = [allvxs, vxbulge]
	  allvys = [allvys, vybulge]
	  allvzs = [allvzs, vzbulge]
	endif else begin
          allxs = [xbulge]
          allys = [ybulge]
          allzs = [zbulge]
          allvxs = [vxbulge]
          allvys = [vybulge]
          allvzs = [vzbulge]
	endelse
    endif

    if npart(4) GT 0 then begin
          allxs = [allxs, xstars]
          allys = [allys, ystars]
	  allzs = [allzs, zstars]
	  allvxs = [allvxs, vxstars]
	  allvys = [allvys, vystars]
	  allvzs = [allvzs, vzstars]
    endif




    if dummy EQ 'x' then begin
	g1vx= allvxs(idx)-cv(0)
	return, g1vx
    endif

    if dummy EQ 'y' then begin
	g1vy= allvys(idx)-cv(1)
	return, g1vy
    endif

    if dummy EQ 'z' then begin
	g1vz= allvzs(idx)-cv(2)
	return, g1vz
    endif

    if dummy EQ 'phi' then begin
        g1x= allxs(idx)-c(0)
        g1y= allys(idx)-c(1)
	g1vx= allvxs(idx)-cv(0)
	g1vy= allvys(idx)-cv(1)
        phi = atan(g1y,g1x)
        v_phi= -g1vx*sin(phi) + g1vy*cos(phi)
        return, v_phi
    endif

    if dummy EQ 'theta' then begin
        g1x= allxs(idx)-c(0)
        g1y= allys(idx)-c(1)
        g1z= allzs(idx)-c(2)
	g1vx= allvxs(idx)-cv(0)
	g1vy= allvys(idx)-cv(1)
	g1vz= allvzs(idx)-cv(2)
        phi= atan(g1y,g1x)
	theta= atan(sqrt(g1x*g1x+ g1y*g1y),g1z)
        v_theta = cos(phi)*cos(theta)*g1vx+ cos(theta)*sin(phi)*g1vy+ cos(theta)*g1vz
        return, v_theta
    endif

    if dummy EQ 'r' then begin
        g1x= allxs(idx)-c(0)
        g1y= allys(idx)-c(1)
        g1z= allzs(idx)-c(2)
        g1vx= allvxs(idx)-cv(0)
        g1vy= allvys(idx)-cv(1)
        g1vz= allvzs(idx)-cv(2)
        phi= atan(g1y,g1x)
        theta= atan(sqrt(g1x*g1x+ g1y*g1y),g1z)
        v_r = cos(phi)*sin(theta)*g1vx+ sin(theta)*sin(phi)*g1vy+ cos(theta)*g1vz
        return, v_r
    endif

end


