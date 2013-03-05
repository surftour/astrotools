function fload_allstars_mass, dummy

    COMMON GalaxyHeader
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData

    if dummy NE -1 then begin

	if (npart(2)+npart(3)+npart(4)) eq 0 then m=[0]

	if npart(2) GT 0 then m=[mdisk]
        if npart(3) GT 0 then begin
		if n_elements(m) gt 0 then m=[m,mbulge] else m=[mbulge]
	endif
        if npart(4) GT 0 then begin
                if n_elements(m) gt 0 then m=[m,mstars] else m=[mstars]
        endif
help, mdisk, mbulge, mstars, m

	return, m
    endif

end



