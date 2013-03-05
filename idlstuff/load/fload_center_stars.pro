Function fload_center_stars, dummy


    COMMON GalaxyHeader
    COMMON HaloData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON BlackHoleData


    if dummy EQ 1 then begin
	com = [0.0,0.0,0.0]
	lastcom = [0.0,0.0,0.0]
	dcom = 100.0

	; is there a disk?
	; ----------------
	if npart(2) GT 0 then begin
		allmasses = [mdisk]
		allxs = [xdisk]
		allys = [ydisk]
		allzs = [zdisk]
	endif

        ; is there a bulge?
        ; -----------------
        if npart(3) GT 0 then begin
	  if (npart(0)+npart(1)) GT 0 then begin
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

	; ------------------------------------------------
	; from this point forward, we're following felix
	; ------------------------------------------------
	rmax = 1000.0
	minnumber= 100
	
	com(0)= total(allmasses*allxs,/double)/total(allmasses,/double)
	com(1)= total(allmasses*allys,/double)/total(allmasses,/double)
	com(2)= total(allmasses*allzs,/double)/total(allmasses,/double)

	allrs = (allxs-com(0))*(allxs-com(0)) + (allys-com(1))*(allys-com(1)) + (allzs-com(2))*(allzs-com(2))
	
	N= n_elements(allmasses)
	count= N

	while (count gt minnumber) do begin

		rmax= 0.8 * rmax

		select = where(allrs lt rmax, count)

		if select(0) eq -1 then break

		oldcenter= com
		com(0) = total(allxs(select)*allmasses(select),/double)/total(allmasses(select),/double)
		com(1) = total(allys(select)*allmasses(select),/double)/total(allmasses(select),/double)
		com(2) = total(allzs(select)*allmasses(select),/double)/total(allmasses(select),/double)

		allrs = (allxs-com(0))*(allxs-com(0)) + (allys-com(1))*(allys-com(1)) + (allzs-com(2))*(allzs-com(2))

	endwhile
	
	;print, "center= ", com," (baryons only)"

        print, FORMAT= '("center= ",F10.5,"  ",F10.5,"  ",F10.5,"   ---- baryons only")', com

	return, com
    endif

end


