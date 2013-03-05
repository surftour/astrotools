function fload_newstars_xyz, dummy, idtofollow=idtofollow, center=center, xy=xy, xz=xz, yz=yz


    COMMON GalaxyHeader
    COMMON NewStarData 
    COMMON OtherData
    COMMON Center


    if keyword_set(center) then begin
        c= center
    endif else begin
        c=com
    endelse
    print, "fload_newstars_xyz using com= ", c



if npart(4) GT 0 then begin
	if dummy EQ 'x' then return, xstars-c(0)
	if dummy EQ 'y' then return, ystars-c(1)
	if dummy EQ 'z' then return, zstars-c(2)

    if dummy EQ 'phi' then begin
        phi = atan(ystars-c(1),xstars-c(0))
        return, phi
    endif

    if dummy EQ 'theta' then begin
        theta= atan(sqrt((xstars-c(0))*(xstars-c(0)) + (ystars-c(1))*(ystars-c(1))),(zstars-c(2)))
        return, theta
    endif

    if dummy EQ 'rxy' then begin
        rxy= sqrt((xstars-c(0))*(xstars-c(0)) + (ystars-c(1))*(ystars-c(1)))
        return, rxy
    endif

    if dummy EQ 'r' then begin
        r= sqrt((xstars-c(0))*(xstars-c(0)) + (ystars-c(1))*(ystars-c(1)) + (zstars-c(2))*(zstars-c(2)))
        return, r
    endif

    ; grab a specific id numbers xyz
    if keyword_set(idtofollow) then begin
        if dummy EQ 'xyz' then begin

        	sid= npart(0)+npart(1)+npart(2)+npart(3)
        	eid= npart(0)+npart(1)+npart(2)+npart(3)+npart(4)-1
        	stars_ids= id(sid:eid)
                idx= where(stars_ids EQ idtofollow)
                if idx(0) LE 0 or n_elements(idx) NE 1 then begin
                        print, "Hey, what's up with idx?  printing"
                        print, idx
                        return, [-1,-1,-1]
                endif

                x= xstars(idx)-c(0)
                y= ystars(idx)-c(1)
                z= zstars(idx)-c(2)
                return, [x,y,z]
        endif
    endif


        if dummy eq 'reff' then begin
                ; prepare projected r and mass
                if keyword_set(xy) then begin
                        gx= xstars-c(0)
                        gy= ystars-c(1)
                        rxy= sqrt(gx*gx + gy*gy)
                endif

                if keyword_set(xz) then begin
                        gx= xstars-c(0)
                        gz= zstars-c(2)
                        rxy= sqrt(gx*gx + gz*gz)
                endif

                if keyword_set(yz) then begin
                        gy= ystars-c(1)
                        gz= zstars-c(2)
                        rxy= sqrt(gy*gy + gz*gz)
                endif

                ms = mstars
                mtot = total(ms)
                hm = 0.5*mtot
                nidx = n_elements(ms)

                ; ok, preparations done
                sorta= sort(rxy)
                r= rxy(sorta)
                m= ms(sorta)

                ; find, effective radius
                n_guess = long(nidx/2.0)
                n_iterations= 0

                repeat begin
                    m_guess = total(m[0:n_guess])
                    dm = m_guess-hm
                    n_guess=long(n_guess - nidx*dm/mtot)
                    n_iterations= n_iterations+1
                endrep until ((abs(dm/mtot) lt 0.0001) or (n_iterations gt 100))

                r_eff = r[n_guess]

                return, r_eff

        endif


endif

end


