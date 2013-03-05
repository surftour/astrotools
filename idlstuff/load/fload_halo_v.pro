Function fload_halo_v, dummy, center=center, comvel=comvel


    COMMON GalaxyHeader
    COMMON HaloData
    COMMON Center

    if keyword_set(center) then c= center else c= com

    if keyword_set(comvel) then cvel= comvel else cvel=[0,0,0]

    allvxs= vxhalo-cvel(0)
    allvys= vyhalo-cvel(1)
    allvzs= vzhalo-cvel(2)

    if dummy EQ 'x' then return, allvxs
    if dummy EQ 'y' then return, allvys
    if dummy EQ 'z' then return, allvzs

    xx= xhalo-c(0)
    yy= yhalo-c(1)
    zz= zhalo-c(2)
    xydist= sqrt(xx*xx + yy*yy)
    xyzdist= sqrt(xx*xx + yy*yy + zz*zz)

;print, "using cvel=", cvel
;print, "using    c=",c

    if dummy eq 'tot' then begin
        v_tot = sqrt(allvxs*allvxs + allvys*allvys + allvzs*allvzs)
        return, v_tot
    endif

    if dummy EQ 'phi' then begin
        ;phi = atan(yy,xx)
        ;v_phi= -allvxs*sin(phi) + allvys*cos(phi)
        v_phi = (allvxs*yy - allvys*xx)/xydist
        return, v_phi
    endif

    if dummy EQ 'theta' then begin
        ;phi= atan(yy,xx)
        ;theta= atan(sqrt(xx*xx + yy*yy),zz)
        ;v_theta = cos(phi)*cos(theta)*allvxs + cos(theta)*sin(phi)*allvys + cos(theta)*allvzs
        v_theta = (-allvxs*xx*zz - allvys*yy*zz + allvzs*xydist*xydist)/sqrt(xx*xx*zz*zz + yy*yy*zz*zz + (xydist*xydist)^2)
        return, v_theta
    endif


    if dummy EQ 'r' then begin
        ;print, "using angles"
        ;phi= atan(yy,xx)
        ;theta= atan(sqrt(xx*xx + yy*yy),zz)
        ;v_r = cos(phi)*sin(theta)*allvxs + sin(theta)*sin(phi)*allvys + cos(theta)*allvzs
        ;print, "using coordinates"
        v_r = (allvxs*xx + allvys*yy + allvzs*zz)/xyzdist
        return, v_r
    endif


    if dummy eq 'tan' then begin
        v_tot = allvxs*allvxs + allvys*allvys + allvzs*allvzs
        v_r = (allvxs*xx + allvys*yy + allvzs*zz)/xyzdist
        v_tan = sqrt(v_tot - v_r*v_r)
        return, v_tan
    endif




end


