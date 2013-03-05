Function fload_center, dummy


    COMMON GalaxyHeader
    COMMON HaloData
    COMMON GasData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON BlackHoleData


    if dummy EQ 1 then begin
	com = [0.0,0.0,0.0]
	lastcom = [0.0,0.0,0.0]
	dcom = 100.0

        ; we know there is a halo
        ; -------------------------
	if npart(0) GT 0 then begin
          allmasses = [mgas]
          allxs = [xgas]
          allys = [ygas]
          allzs = [zgas]
	endif

        ; halo?
        ; ------
	if npart(1) GT 0 then begin
	  if n_elements(allmasses) gt 0 then begin
		allmasses = [mgas, mhalo]
		allxs = [xgas, xhalo]
		allys = [ygas, yhalo]
		allzs = [zgas, zhalo]
	  endif else begin
		allmasses = [mhalo]
		allxs = [xhalo]
		allys = [yhalo]
		allzs = [zhalo]
	  endelse
	endif

	; is there a disk?
	; ----------------
	if npart(2) GT 0 then begin
	  allmasses = [allmasses, mdisk]
	  allxs = [allxs, xdisk]
	  allys = [allys, ydisk]
	  allzs = [allzs, zdisk]
	endif

        ; is there a bulge?
        ; -----------------
        if npart(3) GT 0 then begin
          allmasses = [allmasses, mbulge]
          allxs = [allxs, xbulge]
          allys = [allys, ybulge]
          allzs = [allzs, zbulge]
        endif

        ; is there stars?
        ; -----------------
        if npart(4) GT 0 then begin
          allmasses = [allmasses, mstars]
          allxs = [allxs, xstars]
          allys = [allys, ystars]
          allzs = [allzs, zstars]
        endif

	; is there a black hole?
	; ----------------------
        if npart(5) GT 0 then begin
          allmasses = [allmasses, mbh]
          allxs = [allxs, xbh]
          allys = [allys, ybh]
          allzs = [allzs, zbh]
        endif

	allrs = sqrt(allxs*allxs + allys*allys + allzs*allzs)
	rmax = max(allrs)
	
	orign= n_elements(allmasses)

	iteration= 0

	while (((dcom GT 0.01) AND (n_elements(allmasses) GT 1000)) OR (rmax GT 20.0)) do begin
;print, "i= ", iteration, dcom, n_elements(allmasses), rmax
	
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

	   ;if iteration eq 0 then print, "       rmax   n(allmasses) "
	   ;print, rmax, n_elements(allmasses), "    com = ", com(0),com(1),com(2), " dcom= ",dcom

	; cut search circle in half
	   rmax = rmax / 1.5

	   iteration= iteration+1
	   if iteration gt 100 then break

	endwhile

;	print, "center= ", com," (n=",orign,")"
;
        print, FORMAT= '("center= ",F10.5,"  ",F10.5,"  ",F10.5,"   n= ",F10.1)', com, orign

	return, com
    endif

end


