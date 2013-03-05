pro fload_all_data, allmasses, allxs, allys, allzs, allvxs, allvys, allvzs


    COMMON GalaxyHeader
    COMMON HaloData
    COMMON GasData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON ParentID
    COMMON OtherData
    COMMON BlackHoleData
    COMMON Center


        ; gas!
        ; -------------------------
	if npart(0) GT 0 then begin
          allmasses = [mgas]
	  allxs= [xgas]
	  allys= [ygas]
	  allzs= [zgas]
          allvxs = [vxgas]
          allvys = [vygas]
          allvzs = [vzgas]
	endif

        ; is there a halo?
        ; -----------------
        if npart(1) GT 0 then begin
	  if npart(0) EQ 0 then begin
            allmasses = [mhalo]
            allxs = [xhalo]
            allys = [yhalo]
            allzs = [zhalo]
            allvxs = [vxhalo]
            allvys = [vyhalo]
            allvzs = [vzhalo]
	  endif else begin
            allmasses = [allmasses, mhalo]
	    allxs = [allxs, xhalo]
	    allys = [allys, yhalo]
	    allzs = [allzs, zhalo]
            allvxs = [allvxs, vxhalo]
            allvys = [allvys, vyhalo]
            allvzs = [allvzs, vzhalo]
	  endelse
        endif

        ; is there a disk?
        ; -----------------
        if npart(2) GT 0 then begin
          allmasses = [allmasses, mdisk]
	  allxs = [allxs, xdisk]
	  allys = [allys, ydisk]
	  allzs = [allzs, zdisk]
          allvxs = [allvxs, vxdisk]
          allvys = [allvys, vydisk]
          allvzs = [allvzs, vzdisk]
        endif

        ; is there a bulge?
        ; -----------------
        if npart(3) GT 0 then begin
          allmasses = [allmasses, mbulge]
	  allxs = [allxs, xbulge]
	  allys = [allys, ybulge]
	  allzs = [allzs, zbulge]
          allvxs = [allvxs, vxbulge]
          allvys = [allvys, vybulge]
          allvzs = [allvzs, vzbulge]
        endif

        ; is there stars?
        ; -----------------
        if npart(4) GT 0 then begin
          allmasses = [allmasses, mstars]
	  allxs = [allxs, xstars]
	  allys = [allys, ystars]
	  allzs = [allzs, zstars]
          allvxs = [allvxs, vxstars]
          allvys = [allvys, vystars]
          allvzs = [allvzs, vzstars]
        endif



end


