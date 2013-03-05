Function fload_1gal_comvel, startid, numpart, center=center, rfact=rfact

    if not keyword_set(rfact) then rfact= 0.05

    if numpart GT 0 then begin

	fload_1gal_all_data, startid, numpart, allmasses, allxs, allys, allzs, allvxs, allvys, allvzs

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

	print, "gal1 com velocity   ",comvel, " (n= ",n_elements(allmasses),")"
	return, comvel

    endif

end


