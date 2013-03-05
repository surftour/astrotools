function fload_1gal_allstars_mass, startid, numpart

    COMMON GalaxyHeader
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData

    if numpart gt 0 then begin

	if (npart(2)+npart(3)+npart(4)) eq 0 then m=[0]

	startid= long(startid)
	numpart= long(numpart)

	; first grab id's for disk and bulge particles
	gid = fload_allstars_ids(1)

	; then grab appropriate ones for galaxy 1
	idx = where((gid ge startid) and (gid lt startid+numpart))

	if npart(2) GT 0 then m=[mdisk]
        if npart(3) GT 0 then begin
		if n_elements(m) gt 0 then m=[m,mbulge] else m=[mbulge]
	endif
        if npart(4) GT 0 then begin
                if n_elements(m) gt 0 then m=[m,mstars] else m=[mstars]
        endif

	return, m(idx)
    endif

end



