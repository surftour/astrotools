function fload_baryon_mass, dummy

    COMMON GalaxyHeader
    COMMON GasData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData

    if dummy NE -1 then begin

	if npart(0) GT 0 then m=[mgas]
	if npart(2) GT 0 then begin
		if npart(0) GT 0 then m=[m,mdisk] else m=[mdisk]
	endif
        if npart(3) GT 0 then begin
		if n_elements(m) gt 0 then m=[m,mbulge] else m=[mbulge]
	endif
        if npart(4) GT 0 then m=[m,mstars]

	return, m
    endif

end



