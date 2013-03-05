Function fload_all_comvel, dummy, center=center, justcenter=justcenter, rfac=rfac


    COMMON GalaxyHeader
    COMMON HaloData
    COMMON GasData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON OtherData
    COMMON Center


    if keyword_set(center) then c= center else c= com


    if dummy GE 0 then begin

	; gas?
	; -------------------------
	if npart(0) GT 0 then begin
	  allxs = [xgas]
	  allys = [ygas]
	  allzs = [zgas]
	  allmasses = [mgas]
	  allvxs = [vxgas]
	  allvys = [vygas]
	  allvzs = [vzgas]
	endif

        ; halo?
        ; -------------------------
        if npart(1) GT 0 then begin
	  if npart(0) GT 0 then begin
		allxs = [allxs, xhalo]
		allys = [allys, yhalo]
		allzs = [allzs, zhalo]
		allmasses = [allmasses, mhalo]
		allvxs = [allvxs, vxhalo]
		allvys = [allvys, vyhalo]
		allvzs = [allvzs, vzhalo]
	  endif else begin
		allxs = [xhalo]
		allys = [yhalo]
		allzs = [zhalo]
		allmasses = [mhalo]
		allvxs = [vxhalo]
		allvys = [vyhalo]
		allvzs = [vzhalo]
	  endelse
        endif

        ; disk?
        ; -------------------------
        if npart(2) GT 0 then begin
          allxs = [allxs, xdisk]
          allys = [allys, ydisk]
          allzs = [allzs, zdisk]
          allmasses = [allmasses, mdisk]
          allvxs = [allvxs, vxdisk]
          allvys = [allvys, vydisk]
          allvzs = [allvzs, vzdisk]
        endif

	; is there a bulge?
	; -----------------
	if npart(3) GT 0 then begin
	  allxs = [allxs, xbulge]
	  allys = [allys, ybulge]
	  allzs = [allzs, zbulge]
	  allmasses = [allmasses, mbulge]
	  allvxs = [allvxs, vxbulge]
	  allvys = [allvys, vybulge]
	  allvzs = [allvzs, vzbulge]
	endif

	; is there stars?
	; -----------------
        if npart(4) GT 0 then begin
          allxs = [allxs, xstars]
          allys = [allys, ystars]
          allzs = [allzs, zstars]
          allmasses = [allmasses, mstars]
          allvxs = [allvxs, vxstars]
          allvys = [allvys, vystars]
          allvzs = [allvzs, vzstars]
        endif

print, "N= ", n_elements(allmasses)

	if keyword_set(justcenter) then begin
	   allxs= allxs-c(0)
	   allys= allys-c(1)
	   allzs= allzs-c(2)

	   radius= sqrt(allxs*allxs + allys*allys + allzs*allzs)

	   idx= where(radius le justcenter)
	   allvxs= allvxs(idx)
	   allvys= allvys(idx)
	   allvzs= allvzs(idx)
	   allmasses= allmasses(idx)

	   print, "justcenter= ", justcenter
	   print, "N (within radius)= ", n_elements(allmasses)

	endif

	comvel= fltarr(3)

        idx1= where(finite(allvxs) eq 0)
        idx2= where(finite(allvys) eq 0)
        idx3= where(finite(allvzs) eq 0)
        idx4= where(finite(allmasses) eq 0)

        idx1thru4= [idx1, idx2, idx3, idx4]
        idx= where(idx1thru4 gt 0)
        if idx(0) ne -1 then begin
           idx1thru4= idx1thru4(where(idx1thru4 gt 0))
           print, "fixing"
           print, "idx1thru4= ", idx1thru4
           print, "allvxs(idx1thru4)= ", allvxs(idx1thru4)
           print, "allvys(idx1thru4)= ", allvys(idx1thru4)
           print, "allvzs(idx1thru4)= ", allvzs(idx1thru4)
           print, "allmasses(idx1thru4)= ", allmasses(idx1thru4)
           allvxs(idx1thru4)= 0.0
           allvys(idx1thru4)= 0.0
           allvzs(idx1thru4)= 0.0
           allmasses(idx1thru4)= 0.0
        endif

        idx= where(abs(allvxs) gt 1.0e+6) & if idx(0) ne -1 then allvxs(idx)= 0.0 & if idx(0) ne -1 then help, "idx= ", idx
        idx= where(abs(allvys) gt 1.0e+6) & if idx(0) ne -1 then allvys(idx)= 0.0 & if idx(0) ne -1 then help, "idx= ", idx
        idx= where(abs(allvzs) gt 1.0e+6) & if idx(0) ne -1 then allvzs(idx)= 0.0 & if idx(0) ne -1 then help, "idx= ", idx
        idx= where(abs(allmasses) gt 1.0e+6) & if idx(0) ne -1 then allmasses(idx)= 0.0 & if idx(0) ne -1 then help, "idx= ", idx


	comvel[0]= total(allvxs*allmasses)/total(allmasses)
	comvel[1]= total(allvys*allmasses)/total(allmasses)
	comvel[2]= total(allvzs*allmasses)/total(allmasses)

	print, "total com velocity",comvel
	return, comvel

    endif

end


