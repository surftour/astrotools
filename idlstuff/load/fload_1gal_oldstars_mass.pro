function fload_1gal_oldstars_mass, startid, numpart

    COMMON GalaxyHeader
    COMMON DiskData
    COMMON BulgeData

    if numpart gt 0 then begin

	m= fload_oldstars_mass(1)
	gid= fload_oldstars_id(1)

	; then grab appropriate ones for galaxy 1
	idx = where((gid ge startid) and (gid lt startid+numpart))

	return, m(idx)
    endif

end



