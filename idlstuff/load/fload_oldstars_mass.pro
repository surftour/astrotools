function fload_oldstars_mass, dummy

    COMMON GalaxyHeader
    COMMON DiskData
    COMMON BulgeData

    if dummy NE -1 then begin

	if npart(2) GT 0 then m=[mdisk]
        if npart(3) GT 0 then begin
		if n_elements(m) gt 0 then m=[m,mbulge] else m=[mbulge]
	endif

	return, m
    endif

end



