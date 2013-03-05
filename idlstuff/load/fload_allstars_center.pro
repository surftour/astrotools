Function fload_allstars_center, dummy


    COMMON GalaxyHeader
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData


    if dummy EQ 1 then begin
	com = [0.0,0.0,0.0]
	lastcom = [0.0,0.0,0.0]
	dcom = 100.0

        ; is there a disk
        ; ------------------------
        if npart(2) GT 0 then begin
          allmasses = [mdisk]
          allxs = [xdisk]
          allys = [ydisk]
          allzs = [zdisk]
        endif

	; is there a bulge
	; ------------------------
	if npart(3) GT 0 then begin
	  if n_elements(allmasses) gt 1 then begin
	    allmasses = [allmasses, mbulge]
	    allxs = [allxs, xbulge]
	    allys = [allys, ybulge]
	    allzs = [allzs, zbulge]
	  endif else begin
            allmasses = [mbulge]
            allxs = [xbulge]
            allys = [ybulge]
            allzs = [zbulge]
          endelse
	endif

	; is there new stars
	; ----------------------
	if npart(4) GT 0 then begin
	  if n_elements(allmasses) gt 1 then begin
	    allmasses = [allmasses, mstars]
	    allxs = [allxs, xstars]
	    allys = [allys, ystars]
	    allzs = [allzs, zstars]
	  endif else begin
            allmasses = [mbulge]
            allxs = [xbulge]
            allys = [ybulge]
            allzs = [zbulge]
          endelse

	endif

	allrs = sqrt(allxs*allxs + allys*allys + allzs*allzs)
	rmax = max(allrs)
	
	while (((dcom GT 0.01) AND (n_elements(allmasses) GT 1000)) OR (rmax GT 20.0)) do begin
	
	; calculate distance from com and select less than some radius
	   allrs = sqrt((allxs-com(0))*(allxs-com(0)) + (allys-com(1))*(allys-com(1)) + (allzs-com(2))*(allzs-com(2)))
	   idx = where(allrs LE rmax)
	   allmasses = allmasses(idx)
	   allxs = allxs(idx)
	   allys = allys(idx)
	   allzs = allzs(idx)

	; save last com and find new center of mass
	   lastcom = com
	   com(0) = total(allmasses*allxs)/total(allmasses)
	   com(1) = total(allmasses*allys)/total(allmasses)
	   com(2) = total(allmasses*allzs)/total(allmasses)

	; shift coordinates so com is at 0,0,0
	; no - don't want to do this because we loose our coordinate system
	;   allxs = allxs - com(0)
	;   allys = allys - com(1)
	;   allzs = allzs - com(2)

	; see what error is from 0
	   comdiff = com-lastcom
	   dcom = sqrt(total(comdiff*comdiff))

;	   print, rmax, n_elements(allmasses), " com = ", com(0),com(1),com(2), " dcom= ",dcom

	; cut search circle in half
	   rmax = rmax / 2.0

	endwhile

	print, "stellar (all) center: ",com

	return, com
    endif

end


