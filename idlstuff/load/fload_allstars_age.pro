function fload_allstars_age, dummy

    COMMON GalaxyHeader
    COMMON SfrData
    COMMON OldStarData

    if dummy NE -1 then begin

	if npart(2) GT 0 then allage=[diskage]
        if npart(3) GT 0 then begin
		if n_elements(allage) gt 0 then allage=[allage,bulgeage] else allage=[bulgeage]
	endif
        if npart(4) GT 0 then begin
                if n_elements(allage) gt 0 then allage=[allage,stellage] else allage=[stellage]
        endif

	return, allage
    endif

end



