Function fload_1gal_mass, startid, numpart


    COMMON GalaxyHeader
    COMMON HaloData
    COMMON GasData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON OtherData


    if numpart GT 0 then begin

	; gas?
	; -------------------------
	if npart(0) gt 0 then begin
	  allmasses = [mgas]
	endif

        ; halo?
        ; -------------------------
        if npart(1) gt 0 then begin
          if npart(0) gt 0 then allmasses = [allmasses, mhalo] else allmasses= [mhalo]
        endif

        ; disk?
        ; -------------------------
        if npart(2) gt 0 then begin
          allmasses = [allmasses, mdisk]
        endif

	; is there a bulge?
	; -----------------
	if npart(3) GT 0 then begin
	  allmasses = [allmasses, mbulge]
	endif

	; is there stars?
	; -----------------
        if npart(4) GT 0 then begin
          allmasses = [allmasses, mstars]
        endif

	; just get the particles which have id's for galaxy 1
	; -------------------------------------------------------
	idx = where((id ge long(startid)) and (id lt long(startid)+long(numpart)))
	allmasses = allmasses(idx)


	tot_mass= total(allmasses)
	;print, "gal1 mass:", tot_mass
	return, tot_mass

    endif

end


