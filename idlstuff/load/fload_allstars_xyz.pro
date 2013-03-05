Function fload_allstars_xyz, dummy, center=center, xy=xy, xz=xz, yz=yz, $
				idtofollow=idtofollow

    COMMON GalaxyHeader
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON OtherData
    COMMON Center

    if dummy NE '-1' then begin

	if (npart(2)+npart(3)+npart(4)) eq 0 then begin
	   allxs=[0]
	   allys=[0]
	   allzs=[0]
	endif

	if keyword_set(center) then begin
		c= center
	endif else begin
		c= com
	endelse
	print, "fload_allstars_xyz using com= ", c

; -------------------------------------------------------------------
; now recalculate allxs, allys and allzs for returning correct x,y,z
; -------------------------------------------------------------------

        ; is there a disk?
        ; -----------------
        if npart(2) GT 0 then begin
	    allxs = [xdisk]
	    allys = [ydisk]
	    allzs = [zdisk]
	    allmass = [mdisk]
        endif

        ; is there a bulge?
        ; -----------------
        if npart(3) GT 0 then begin
	  if n_elements(allxs) gt 0 then begin
            allxs = [allxs, xbulge]
            allys = [allys, ybulge]
	    allzs = [allzs, zbulge]
	    allmass = [allmass, mbulge]
	  endif else begin
            allxs = [xbulge]
            allys = [ybulge]
            allzs = [zbulge]
            allmass = [mbulge]
	  endelse
        endif

        ; is there stars?
        ; -----------------
        if npart(4) GT 0 then begin
          if n_elements(allxs) gt 0 then begin
            allxs = [allxs, xstars]
            allys = [allys, ystars]
            allzs = [allzs, zstars]
            allmass = [allmass, mstars]
          endif else begin
            allxs = [xstars]
            allys = [ystars]
            allzs = [zstars]
            allmass = [mstars]
          endelse
        endif


	; load star ID's
	; ---------------
	if keyword_set(idtofollow) then begin
	  ;startids=npart(0)+npart(1)
	  ;allstarids= id(startids:startids+npart(2)+npart(3)+npart(4)-1)
	  allstarids= fload_allstars_ids(1)

          idx= where(allstarids EQ idtofollow)
          if idx(0) eq -1 then begin
                print, "Hey, what's up with IDX?"
                print, "are you giving me bogus information?"
                return, [-1,-1,-1]
          endif 

          x= allxs(idx)-c(0)
          y= allys(idx)-c(1)
          z= allzs(idx)-c(2)

	  return, [x,y,z]
	endif


	if dummy eq 'x' then return, allxs-c(0)
	if dummy eq 'y' then return, allys-c(1)
	if dummy eq 'z' then return, allzs-c(2)
	if dummy eq 'r' then begin
		r2= (allxs-c(0))*(allxs-c(0)) + (allys-c(1))*(allys-c(1)) + (allzs-c(2))*(allzs-c(2))
		r= sqrt(r2)
		return, r
	endif


	if dummy eq 'rxy' then begin
                r2= (allxs-c(0))*(allxs-c(0)) + (allys-c(1))*(allys-c(1))
                r= sqrt(r2)
                return, r
        endif


	if dummy eq 'reff' then begin
        	; prepare projected r and mass
        	if keyword_set(xy) then begin
        	        gx= allxs-c(0)
        	        gy= allys-c(1)
        	        rxy= sqrt(gx*gx + gy*gy)
        	endif

        	if keyword_set(xz) then begin
        	        gx= allxs-c(0)
        	        gz= allzs-c(2)
        	        rxy= sqrt(gx*gx + gz*gz)
        	endif

        	if keyword_set(yz) then begin
        	        gy= allys-c(1)
        	        gz= allzs-c(2)
        	        rxy= sqrt(gy*gy + gz*gz)
        	endif

        	ms = allmass
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


