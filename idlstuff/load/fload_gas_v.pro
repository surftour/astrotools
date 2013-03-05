Function fload_gas_v, dummy, center=center, comvel=comvel, idtofollow=idtofollow

    COMMON GalaxyHeader
    COMMON GasData
    COMMON OtherData
    COMMON Center


    ; ----------------------------------
    if keyword_set(center) then begin
        c= center
    endif else begin
        c= com
    endelse


    if not keyword_set(comvel) then begin
        comvel= [0.0, 0.0, 0.0]
    endif


    ; -------------------------
    if npart(0) GT 0 then begin
	allvxs = [vxgas] - comvel(0)
	allvys = [vygas] - comvel(1)
	allvzs = [vzgas] - comvel(2)
    endif else begin
	return, [0]
    endelse




    if dummy EQ 'mag' then return, sqrt(allvxs*allvxs + allvys*allvys + allvzs*allvzs)
    if dummy EQ 'v2' then return, allvxs*allvxs + allvys*allvys + allvzs*allvzs
    if dummy EQ 'x' then return, allvxs
    if dummy EQ 'y' then return, allvys
    if dummy EQ 'z' then return, allvzs

    xx= xgas-c(0)
    yy= ygas-c(1)
    zz= zgas-c(2)
    xydist= sqrt(xx*xx + yy*yy)
    xyzdist= sqrt(xx*xx + yy*yy + zz*zz)

    if dummy EQ 'phi' then begin
	phi = atan(ygas-c(1),xgas-c(0))
	v_phi= -allvxs*sin(phi) + allvys*cos(phi)
	;v_phi = (allvxs*yy - allvys*xx)/xydist
	return, v_phi
    endif

    if dummy EQ 'theta' then begin
        phi= atan(ygas-c(1),xgas-c(0))
        theta= atan(sqrt((xgas-c(0))*(xgas-c(0)) + (ygas-c(1))*(ygas-c(1))),(zgas-c(2)))
        v_theta = cos(phi)*cos(theta)*allvxs + cos(theta)*sin(phi)*allvys + cos(theta)*allvzs
	;v_theta = (-allvxs*xx*zz - allvys*yy*zz + allvzs*xydist*xydist)/sqrt(xx*xx*zz*zz + yy*yy*zz*zz + (xydist*xydist)^2)
        return, v_theta
    endif

    if dummy EQ 'r' then begin
	phi= atan(ygas-c(1),xgas-c(0))
	theta= atan(sqrt((xgas-c(0))*(xgas-c(0)) + (ygas-c(1))*(ygas-c(1))),(zgas-c(2)))
	v_r = cos(phi)*sin(theta)*allvxs + sin(theta)*sin(phi)*allvys + cos(theta)*allvzs
	;v_r = (allvxs*xx + allvys*yy + allvzs*zz)/xyzdist
	return, v_r
    endif


    if keyword_set(idtofollow) then begin

	gas_ids= id(0:npart(0)-1)
	idx= where(gas_ids EQ idtofollow)
	if idx(0) EQ -1 or n_elements(idx) NE 1 then begin
		print, "Hey, what's up with idx?  printing"
		print, idx
		return, [0,0,0]
	endif

	vx= allvxs(idx)
	vy= allvys(idx)
	vz= allvzs(idx)
	x= xgas(idx)-c(0)
	y= ygas(idx)-c(1)
	z= zgas(idx)-c(2)


	if dummy EQ '0_id' then return, sqrt(vx*vx + vy*vy + vz*vz)
	if dummy EQ 'r_id' then begin
		phi= atan(y,x)
		theta= atan(sqrt(x*x + y*y),z)
		return, cos(phi)*sin(theta)*vx + sin(theta)*sin(phi)*vy + cos(theta)*vz
	endif

	return, [vx,vy,vz]

    endif

end




