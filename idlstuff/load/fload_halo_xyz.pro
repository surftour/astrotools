function fload_halo_xyz, dummy, center=center


    COMMON GalaxyHeader
    COMMON HaloData
    COMMON Center

    if keyword_set(center) then begin
	c= center
    endif else begin
	c= com
    endelse

    if npart(1) GT 0 then begin
	if dummy EQ 'x' then return, xhalo-c(0)
	if dummy EQ 'y' then return, yhalo-c(1)
	if dummy EQ 'z' then return, zhalo-c(2)

	if dummy EQ 'phi' then begin
           phi = atan(yhalo-c(1),xhalo-c(0))
           return, phi
	endif

	if dummy EQ 'theta' then begin
           theta= atan(sqrt((xhalo-c(0))*(xhalo-c(0)) + (yhalo-c(1))*(yhalo-c(1))),(zhalo-c(2)))
           return, theta
	endif

        if dummy EQ 'rxy' then begin
           rxy= sqrt((xhalo-c(0))*(xhalo-c(0)) + (yhalo-c(1))*(yhalo-c(1)))
           return, rxy
        endif

	if dummy EQ 'r' then begin
           r= sqrt((xhalo-c(0))*(xhalo-c(0)) + (yhalo-c(1))*(yhalo-c(1)) + (zhalo-c(2))*(zhalo-c(2)))
           return, r
	endif

    endif

end


