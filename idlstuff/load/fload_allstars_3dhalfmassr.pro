Function fload_allstars_3dhalfmassr, dummy, center=center

    COMMON GalaxyHeader
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON OtherData
    COMMON Center

    if dummy NE '-1' then begin

	if (npart(2)+npart(3)+npart(4)) eq 0 then return, 0.0

	if keyword_set(center) then begin
		c= center
	endif else begin
		c= com
	endelse

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


	; ---------------
	r2= (allxs-c(0))*(allxs-c(0)) + (allys-c(1))*(allys-c(1)) + (allzs-c(2))*(allzs-c(2))
	r= sqrt(r2)

        	ms = allmass
        	mtot = total(ms)
        	hm = 0.5*mtot
        	nidx = n_elements(ms)

        	; ok, preparations done
        	sorta= sort(r)
        	r= r(sorta)
        	m= ms(sorta)

        	; find, effective radius
        	n_guess = long(nidx/2.0)
		n_iterations= 0

        	repeat begin
        	    m_guess = total(m[0:n_guess])
        	    dm = m_guess-hm
        	    ;n_guess=long(n_guess - nidx*dm/mtot)
        	    n_guess=long(n_guess - nidx*0.5*dm/mtot)
		    n_iterations= n_iterations+1

		    last_n_guess= n_guess
		    second_to_last_n_guess= last_n_guess

        	endrep until ((abs(dm/mtot) lt 0.001) or (n_iterations gt 100))

        	r_eff = r[n_guess]

        	return, r_eff


    endif

end


