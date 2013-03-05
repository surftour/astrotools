function fload_blackhole_xyz, dummy, idtofollow=idtofollow, center=center


    COMMON GalaxyHeader
    COMMON BlackHoleData
    COMMON OtherData
    COMMON Center

    if keyword_set(center) then begin
        c= center
    endif else begin
        c=com
    endelse



    ; grab a specific id numbers xyz
    if keyword_set(idtofollow) then begin
	if dummy EQ 'xyz' then begin

		startids=npart(0)+npart(1)+npart(2)+npart(3)+npart(4)
		bh_ids= id(startids:startids+npart(5)-1)
		idx= where(bh_ids EQ idtofollow)
		if idx(0) LT 0 or n_elements(idx) NE 1 then begin
			print, "Hey, what's up with idx?  printing"
			print, idx
			return, [-1,-1,-1]
		endif

		x= xbh(idx)-c(0)
		y= ybh(idx)-c(1)
		z= zbh(idx)-c(2)

		;print, "bh_ids= ", bh_ids
		;print, "idtofollow= ", idtofollow
		;print, "x, y, z= ", x, y, z

		return, [x,y,z]
	endif
    endif



    if dummy EQ 'x' then return, xbh-c(0)
    if dummy EQ 'y' then return, ybh-c(1)
    if dummy EQ 'z' then return, zbh-c(2)

    if dummy EQ 'phi' then begin
        phi = atan(ybh-c(1),xbh-c(0))
        return, phi
    endif

    if dummy EQ 'theta' then begin
        theta= atan(sqrt((xbh-c(0))*(xbh-c(0)) + (ybh-c(1))*(ybh-c(1))),(zbh-c(2)))
        return, theta
    endif

    if dummy EQ 'r' then begin
        r= sqrt((xbh-c(0))*(xbh-c(0)) + (ybh-c(1))*(ybh-c(1)) + (zbh-c(2))*(zbh-c(2)))
        return, r
    endif

    if dummy EQ 'rxy' then begin
        r= sqrt((xbh-c(0))*(xbh-c(0)) + (ybh-c(1))*(ybh-c(1)))
        return, r
    endif


    ; get the effective radii
    if dummy EQ 'reff' then begin

	; xy projection is the default
	r= sqrt((xgas-c(0))*(xgas-c(0)) + (ygas-c(1))*(ygas-c(1)))
	sortidx= sort(r)
	r= r(sortidx)
	m= mgas(sortidx)

	reff= r(0:n_elements(r)/2-1)
	meff= m(0:n_elements(r)/2-1)
	print, "Reff= ",max(reff), "      Mgas=",total(mgas),"      Meff=",total(meff)
	return, max(reff)

    endif
end




