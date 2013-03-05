function fload_allstars_z, dummy

    COMMON GalaxyHeader
    COMMON SfrData
    COMMON OldStarData

    if dummy NE -1 then begin

	if npart(2) GT 0 then zmets=[diskmetals]
        if npart(3) GT 0 then begin
		if n_elements(zmets) gt 0 then zmets=[zmets,bulgemetals] else zmets=[bulgemetals]
	endif
        if npart(4) GT 0 then begin
                if n_elements(zmets) gt 0 then zmets=[zmets,smetals] else zmets=[smetals]
        endif

	return, zmets
    endif

end



