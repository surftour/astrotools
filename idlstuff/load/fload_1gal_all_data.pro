pro fload_1gal_all_data, startid, numpart, allmasses, allxs, allys, allzs, allvxs, allvys, allvzs


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
	;
	; take this out because of our new star forming
	; scheme this makes some of these id's not in
	; the id range we search for in the following step.
	; instead we'll deal with this later.
	;
	;
	; (put this back in, and we'll deal with the
	;  whacky ID's below)


        ; just get the particles which have id's for galaxy 1
        ; -------------------------------------------------------
	if keyword_set(startid) then begin

	    ;allids= id(0:nparts_minus_ns-1)
	    allids= id

	    ; need to correct for whacky ID's
	    whacky_ids= where((allids lt 0) or (allids gt 5.0e7))
	    num= long(1)
	    if whacky_ids(0) ne -1 then allids(whacky_ids)= allids(whacky_ids) - ishft(num,31)


	    idx = where((allids ge long(startid)) and (allids lt long(startid)+long(numpart)))

            allmasses = allmasses(idx)
            allxs = allxs(idx)
            allys = allys(idx)
            allzs = allzs(idx)
	    allvxs = allvxs(idx)
	    allvys = allvys(idx)
	    allvzs = allvzs(idx)
	endif


        ; deal with stars
        ; -------------------------
        ;if npart(4) GT 0 then begin
;
        ;    if flag_stargens eq 1 and flag_parentid eq 1 then begin
        ;        allids=parentid
        ;    endif else begin
;		nsstartid= npart(0)+npart(1)+npart(2)+npart(3)
;		Ntot= npart(0)+npart(1)+npart(2)+npart(3)+npart(4)
;                allids=id(nparts_minus_ns:Ntot-1)
;            endelse 
;
;	    idx= where((allids ge long(startid)) and (allids lt long(startid)+long(numpart)))
;
;	   if idx(0) ne -1 then begin
;		allmasses = [allmasses, mstars(idx)]
;		allxs = [allxs, xstars(idx)]
;		allys = [allys, ystars(idx)]
;		allzs = [allzs, zstars(idx)]
;		allvxs = [allvxs, vxstars(idx)]
;		allvys = [allvys, vystars(idx)]
;		allvzs = [allvzs, vzstars(idx)]
;	   endif
;        endif




end


