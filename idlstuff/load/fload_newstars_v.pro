Function fload_newstars_v, dummy, center=center, vcom=vcom, idtofollow=idtofollow

    COMMON GalaxyHeader
    COMMON NewStarData
    COMMON OtherData
    COMMON Center


    ; ----------------------------------
    if keyword_set(center) then begin
        c= center
    endif else begin
        c= com
    endelse


    if not keyword_set(vcom) then begin
        vcom= [0.0, 0.0, 0.0]
    endif


    ; -------------------------
    if npart(4) GT 0 then begin
	allvxs = [vxstars] - vcom(0)
	allvys = [vystars] - vcom(1)
	allvzs = [vzstars] - vcom(2)
    endif else begin
	allvxs = [0]
	allvys = [0]
	allvzs = [0]
	return, [0,0,0]
    endelse




    if dummy EQ 'mag' then return, sqrt(allvxs*allvxs + allvys*allvys + allvzs*allvzs)
    if dummy EQ 'x' then return, allvxs
    if dummy EQ 'y' then return, allvys
    if dummy EQ 'z' then return, allvzs

    if dummy EQ 'phi' then begin
	phi = atan(ystars-c(1),xstars-c(0))
	v_phi= -vxstars*sin(phi) + vystars*cos(phi)
	return, v_phi
    endif

    if dummy EQ 'theta' then begin
        phi= atan(ystars-c(1),xstars-c(0))
        theta= atan(sqrt((xstars-c(0))*(xstars-c(0)) + (ystars-c(1))*(ystars-c(1))),(zstars-c(2)))
        v_theta = cos(phi)*cos(theta)*allvxs + cos(theta)*sin(phi)*allvys + cos(theta)*allvzs
        return, v_theta
    endif

    if dummy EQ 'r' then begin
	phi= atan(ystars-c(1),xstars-c(0))
	theta= atan(sqrt((xstars-c(0))*(xstars-c(0)) + (ystars-c(1))*(ystars-c(1))),(zstars-c(2)))
	v_r = cos(phi)*sin(theta)*allvxs + sin(theta)*sin(phi)*allvys + cos(theta)*allvzs
	return, v_r
    endif


    if keyword_set(idtofollow) then begin

        sid= npart(0)+npart(1)+npart(2)+npart(3)
        eid= npart(0)+npart(1)+npart(2)+npart(3)+npart(4)-1
        stars_ids= id(sid:eid)
	idx= where(stars_ids EQ idtofollow)
	if idx(0) EQ -1 or n_elements(idx) NE 1 then begin
		print, "Hey, what's up with idx?  printing"
		print, idx
		return, [0,0,0]
	endif

	vx= allvxs(idx)
	vy= allvys(idx)
	vz= allvzs(idx)
	x= xstars(idx)-c(0)
	y= ystars(idx)-c(1)
	z= zstars(idx)-c(2)


	if dummy EQ '0_id' then return, sqrt(vx*vx + vy*vy + vz*vz)
	if dummy EQ 'r_id' then begin
		phi= atan(y,x)
		theta= atan(sqrt(x*x + y*y),z)
		return, cos(phi)*sin(theta)*vx + sin(theta)*sin(phi)*vy + cos(theta)*vz
	endif

    endif

end




