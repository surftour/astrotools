function fload_bulge_xyz, dummy, center=center


    COMMON GalaxyHeader
    COMMON BulgeData
    COMMON Center


    if keyword_set(center) then begin
	c= center
    endif else begin
	c= com
    endelse

    if npart(3) GT 0 then begin
	if dummy EQ 'x' then return, xbulge-c(0)
	if dummy EQ 'y' then return, ybulge-c(1)
	if dummy EQ 'z' then return, zbulge-c(2)

	if dummy EQ 'phi' then begin
	   phi = atan(ybulge-c(1),xbulge-c(0))
	   return, phi
	endif

	if dummy EQ 'theta' then begin
	   theta= atan(sqrt((xbulge-c(0))*(xbulge-c(0)) + (ybulge-c(1))*(ybulge-c(1))),(zbulge-c(2)))
	   return, theta
	endif

	if dummy EQ 'r' then begin
	   r= sqrt((xbulge-c(0))*(xbulge-c(0)) + (ybulge-c(1))*(ybulge-c(1)) + (zbulge-c(2))*(zbulge-c(2)))
	   return, r
	endif

	if dummy EQ 'rxy' then begin
           rxy= sqrt((xbulge-c(0))*(xbulge-c(0)) + (ybulge-c(1))*(ybulge-c(1)))
           return, rxy
	endif

    endif

end


