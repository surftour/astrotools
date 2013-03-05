function fload_1gal_newstars_xyz, dummy, startid, numpart, center=center, $
					xy=xy, xz=xz, yz=yz, $
					useparentinfo=useparentinfo

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON NewStarData
    COMMON Center
    COMMON ParentID


    if npart(4) le 0 then return, [0]

    if keyword_set(center) then begin
	c= center
    endif else begin
	c= com
    endelse

    nsid= fload_newstars_id(1)

    if keyword_set(useparentinfo) then begin
	; now instead of grabbing new star id number's
	;  we'll grab their parent id's
	if flag_stargens eq 1 and flag_parentid eq 1 then begin
		nsid= [parentid]
	endif else begin
		nsstartid= npart(0)+npart(1)+npart(2)+npart(3)
		N= total(npart)
		nsid= id(nsstartid:N-1)
	endelse
    endif


    idx = where((nsid ge long(startid)) and (nsid lt long(startid)+long(numpart)))


    allxs= xstars
    allys= ystars
    allzs= zstars
    allmasses= mstars

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


