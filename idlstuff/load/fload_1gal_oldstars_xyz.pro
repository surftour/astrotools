function fload_1gal_oldstars_xyz, dummy, startid, numpart, center=center, $
					xy=xy, xz=xz, yz=yz

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON DiskData
    COMMON BulgeData

    if keyword_set(center) then begin
	c= center
    endif else begin
	c=fload_1gal_center(startid,numpart)
    endelse

    ; first grab id's for all baryons
    ;if (npart(2)+npart(3) gt 0) then begin
;	gid= [id(npart(0)+npart(1):N-1-npart(4))]
    ;endif else begin
;	print, "PROBLEM"
;	return, [0]
    ;endelse
    gid= fload_oldstars_id(1)


    ; then grab appropriate ones for galaxy 1
    idx = where((gid ge startid) and (gid lt startid+numpart))


    if npart(2) gt 0 then begin
	allxs = [xdisk]
	allys = [ydisk]
	allzs = [zdisk]
	allmasses = [mdisk]
    endif

    if npart(3) GT 0 then begin
	if n_elements(allmasses) gt 0 then begin
	  allxs = [allxs, xbulge]
	  allys = [allys, ybulge]
	  allzs = [allzs, zbulge]
	  allmasses = [allmasses, mbulge]
	endif else begin
          allxs = [xbulge]
          allys = [ybulge]
          allzs = [zbulge]
          allmasses = [mbulge]
	endelse
    endif



    if dummy EQ 'x' then begin
        g1x= allxs(idx)
        return, g1x-c(0)
    endif

    if dummy EQ 'y' then begin
        g1y= allys(idx)
        return, g1y-c(1)
    endif

    if dummy EQ 'z' then begin
        g1z= allzs(idx)
        return, g1z-c(2)
    endif

    if dummy EQ 'phi' then begin
        g1x= allxs(idx)-c(0)
        g1y= allys(idx)-c(1)
        phi = atan(g1y,g1x)
        return, phi
    endif

    if dummy EQ 'theta' then begin
        g1x= allxs(idx)-c(0)
        g1y= allys(idx)-c(1)
        g1z= allzs(idx)-c(2)
        theta= atan(sqrt(g1x*g1x+ g1y*g1y),g1z)
        return, theta
    endif

    if dummy EQ 'rxy' then begin
        g1x= allxs(idx)-c(0)
        g1y= allys(idx)-c(1)
        rxy= sqrt(g1x*g1x+ g1y*g1y)
        return, rxy
    endif

    if dummy EQ 'r' then begin
        g1x= allxs(idx)-c(0)
        g1y= allys(idx)-c(1)
        g1z= allzs(idx)-c(2)
        r= sqrt(g1x*g1x+ g1y*g1y+ g1z*g1z)
        return, r
    endif


    if dummy EQ 'reff' then begin
	; prepare projected r and mass
	if keyword_set(xy) then begin
		g1x= allxs(idx)-c(0)
		g1y= allys(idx)-c(1)
		rxy= sqrt(g1x*g1x + g1y*g1y)
	endif

	if keyword_set(xz) then begin
		g1x= allxs(idx)-c(0)
		g1z= allzs(idx)-c(2)
		rxy= sqrt(g1x*g1x + g1z*g1z)
	endif

	if keyword_set(yz) then begin
		g1y= allys(idx)-c(1)
		g1z= allzs(idx)-c(2)
		rxy= sqrt(g1y*g1y + g1z*g1z)
	endif

	ms = allmasses(idx)
	mtot = total(ms)
	hm = 0.5*mtot
	nidx = n_elements(ms)

	; ok, preparations done
	sorta= sort(rxy)
	r= rxy(sorta)
	m= ms(sorta)

	; find, effective radius
	n_guess = long(nidx/2.0)
	niterations= 0

	repeat begin
	    m_guess = total(m[0:n_guess])
	    dm = m_guess-hm
	    n_guess=long(n_guess - nidx*dm/mtot)
	    niterations= niterations+1
	endrep until ((abs(dm/mtot) lt 0.0001) or (niterations gt 100))

	r_eff = r[n_guess]

	return, r_eff
	
    endif

end


