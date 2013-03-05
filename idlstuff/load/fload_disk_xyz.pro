function fload_disk_xyz, dummy, center=center


    COMMON GalaxyHeader
    COMMON DiskData
    COMMON Center

    if keyword_set(center) then begin
	c= center
    endif else begin
	c= com
    endelse

    if npart(2) GT 0 then begin
	if dummy EQ 'x' then return, xdisk-c(0)
	if dummy EQ 'y' then return, ydisk-c(1)
	if dummy EQ 'z' then return, zdisk-c(2)

	if dummy EQ 'phi' then begin
           phi = atan(ydisk-c(1),xdisk-c(0))
           return, phi
	endif

	if dummy EQ 'theta' then begin
           theta= atan(sqrt((xdisk-c(0))*(xdisk-c(0)) + (ydisk-c(1))*(ydisk-c(1))),(zdisk-c(2)))
           return, theta
	endif

        if dummy EQ 'rxy' then begin
           rxy= sqrt((xdisk-c(0))*(xdisk-c(0)) + (ydisk-c(1))*(ydisk-c(1)))
           return, rxy
        endif

	if dummy EQ 'r' then begin
           r= sqrt((xdisk-c(0))*(xdisk-c(0)) + (ydisk-c(1))*(ydisk-c(1)) + (zdisk-c(2))*(zdisk-c(2)))
           return, r
	endif

    endif

end


