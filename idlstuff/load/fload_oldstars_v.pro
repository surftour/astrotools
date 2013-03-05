Function fload_oldstars_v, dummy, center=center, comvel=comvel, $
                idtofollow=idtofollow


    COMMON GalaxyHeader
    COMMON DiskData
    COMMON BulgeData
    COMMON Center


    if keyword_set(center) then c= center else c= com

    if keyword_set(comvel) then cvel= comvel else cvel=[0,0,0]


    if (npart(2)+npart(3)) eq 0 then begin
             allvxs = [0]
             allvys = [0]
             allvzs = [0]
             allxs = [0]
             allys = [0]
             allzs = [0]
    endif


        ; is there a disk?
        ; -----------------
        if npart(2) GT 0 then begin
             allvxs = [vxdisk]
             allvys = [vydisk]
             allvzs = [vzdisk]
             allxs = [xdisk]
             allys = [ydisk]
             allzs = [zdisk]
        endif

        ; is there a bulge?
        ; -----------------
        if npart(3) GT 0 then begin
	 if n_elements(allvxs) gt 0 then begin
          allvxs = [allvxs, vxbulge]
          allvys = [allvys, vybulge]
          allvzs = [allvzs, vzbulge]
	  allxs = [allxs, xbulge]
          allys = [allys, ybulge]
          allzs = [allzs, zbulge]
	 endif else begin
          allvxs = [vxbulge]
          allvys = [vybulge]
          allvzs = [vzbulge]
          allxs = [xbulge]
          allys = [ybulge]
          allzs = [zbulge]
	 endelse
        endif

    allvxs= allvxs-cvel(0)
    allvys= allvys-cvel(1)
    allvzs= allvzs-cvel(2)


        ; load star ID's
        ; ---------------
        if keyword_set(idtofollow) then begin
          ;startids=npart(0)+npart(1)
          ;allstarids= id(startids:startids+npart(2)+npart(3)+npart(4)-1)
          oldstarids= fload_oldstars_id(1)

          idx= where(oldstarids EQ idtofollow)
          if idx(0) eq -1 then begin
                print, "Hey, what's up with IDX?"
                print, "are you giving me bogus information?"
                return, [-1,-1,-1]
          endif

          vx= allvxs(idx)
          vy= allvys(idx)
          vz= allvzs(idx)

          return, [vx,vy,vz]
        endif



    if dummy EQ 'x' then return, allvxs
    if dummy EQ 'y' then return, allvys
    if dummy EQ 'z' then return, allvzs


    xx= allxs-c(0)
    yy= allys-c(1)
    zz= allzs-c(2)
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


