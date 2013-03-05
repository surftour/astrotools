function fload_gas_xyz, dummy, idtofollow=idtofollow, center=center


    COMMON GalaxyHeader
    COMMON GasData
    COMMON Center
    COMMON OtherData
    COMMON Center

    if keyword_set(center) then begin
        c= center
    endif else begin
        c= com
    endelse
    

    if dummy EQ 'x' then return, xgas-c(0)
    if dummy EQ 'y' then return, ygas-c(1)
    if dummy EQ 'z' then return, zgas-c(2)

    if dummy EQ 'phi' then begin
        phi = atan(ygas-c(1),xgas-c(0))
        return, phi
    endif

    if dummy EQ 'theta' then begin
        theta= atan(sqrt((xgas-c(0))*(xgas-c(0)) + (ygas-c(1))*(ygas-c(1))),(zgas-c(2)))
        return, theta
    endif

    if dummy EQ 'r' then begin
        r= sqrt((xgas-c(0))*(xgas-c(0)) + (ygas-c(1))*(ygas-c(1)) + (zgas-c(2))*(zgas-c(2)))
        return, r
    endif

    if dummy EQ 'rxy' then begin
        r= sqrt((xgas-c(0))*(xgas-c(0)) + (ygas-c(1))*(ygas-c(1)))
        return, r
    endif

    ; grab a specific id numbers xyz
    if keyword_set(idtofollow) then begin
	if dummy EQ 'xyz' then begin

		gas_ids= id(0:npart(0)-1)
		idx= where(gas_ids EQ idtofollow)
		if idx(0) LE 0 or n_elements(idx) NE 1 then begin
			print, "Hey, what's up with idx?  printing"
			print, idx
			return, [-1,-1,-1]
		endif

		x= xgas(idx)-c(0)
		y= ygas(idx)-c(1)
		z= zgas(idx)-c(2)
		return, [x,y,z]
	endif
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




