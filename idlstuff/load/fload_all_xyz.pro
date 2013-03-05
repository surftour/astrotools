Function fload_all_xyz, dummy, center=center

; --------------------------------------------------------
; Note: this is essentially fload_center, with a minor
;       modification, to save allrs and return with
;       center subtracted
; --------------------------------------------------------

    COMMON GalaxyHeader
    COMMON HaloData
    COMMON GasData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON Center


    if keyword_set(center) then begin
	c= center
    endif else begin
	c= com
    endelse


        ; is there a gas?
        ; -----------------
        if npart(0) GT 0 then begin
          allxs = [xgas]
          allys = [ygas]
          allzs = [zgas]
        endif

        ; is there a halo?
        ; -----------------
        if npart(1) GT 0 then begin
	  if npart(0) GT 0 then begin
             allxs = [allxs, xhalo]
             allys = [allys, yhalo]
             allzs = [allzs, zhalo]
	  endif else begin
	     allxs = [xhalo]
	     allys = [yhalo]
	     allzs = [zhalo]
	  endelse
        endif

        ; is there a disk?
        ; -----------------
        if npart(2) GT 0 then begin
          allxs = [allxs, xdisk]
          allys = [allys, ydisk]
          allzs = [allzs, zdisk]
        endif

        ; is there a bulge?
        ; -----------------
        if npart(3) GT 0 then begin
          allxs = [allxs, xbulge]
          allys = [allys, ybulge]
	  allzs = [allzs, zbulge]
        endif

        ; is there stars?
        ; -----------------
        if npart(4) GT 0 then begin
          allxs = [allxs, xstars]
          allys = [allys, ystars]
	  allzs = [allzs, zstars]
        endif


        if dummy eq 'x' then return, allxs-c(0)
        if dummy eq 'y' then return, allys-c(1)
        if dummy eq 'z' then return, allzs-c(2)

	if dummy eq 'r' then begin
	   r2= (allxs-c(0))*(allxs-c(0)) + (allys-c(1))*(allys-c(1)) + (allzs-c(2))*(allzs-c(2))
	   r= sqrt(r2)
	   return, r
	endif

end


