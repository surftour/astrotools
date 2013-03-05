Function fload_1gal_dm_comvel, startid, numpart, center=center, rfact=rfact

    if not keyword_set(rfact) then rfact= 0.05


    COMMON GalaxyHeader
    COMMON HaloData
    COMMON ParentID
    COMMON OtherData



    if numpart GT 0 then begin

        ; is there a halo?
        ; -----------------
        if npart(1) GT 0 then begin
            allmasses = [mhalo]
            allxs = [xhalo]
            allys = [yhalo]
            allzs = [zhalo]
            allvxs = [vxhalo]
            allvys = [vyhalo]
            allvzs = [vzhalo]
        endif else return, [0,0,0]


        ; just get the particles which have id's for galaxy 1
        ; -------------------------------------------------------
        if keyword_set(startid) then begin
            hids = id(npart(0):npart(0)+npart(1)-1)
            idx = where((hids ge long(startid)) and (hids lt long(startid)+long(numpart)))
            allmasses = allmasses(idx)
            allxs = allxs(idx)
            allys = allys(idx)
            allzs = allzs(idx)
            allvxs = allvxs(idx)
            allvys = allvys(idx)
            allvzs = allvzs(idx)
        endif


	if keyword_set(center) then begin
		allrs= sqrt((allxs-center(0))*(allxs-center(0)) + (allys-center(1))*(allys-center(1)) + (allzs-center(2))*(allzs-center(2)))
		idx= where(allrs LE (rfact*max(allrs)))

		allmasses= allmasses(idx)
		allvxs= allvxs(idx)
		allvys= allvys(idx)
		allvzs= allvzs(idx)
	endif

	comvel= fltarr(3)

	comvel[0]= total(allvxs*allmasses)/total(allmasses)
	comvel[1]= total(allvys*allmasses)/total(allmasses)
	comvel[2]= total(allvzs*allmasses)/total(allmasses)

	print, "gal1 dm com velocity   ",comvel, " (n= ",n_elements(allmasses),")"
	return, comvel

    endif

end


