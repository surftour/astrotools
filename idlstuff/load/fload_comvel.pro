Function fload_comvel, dummy, center=center, rfact=rfact

    if not keyword_set(rfact) then rfact= 0.05

    if dummy NE -23 then begin

	fload_all_data, allmasses, allxs, allys, allzs, allvxs, allvys, allvzs

	if keyword_set(center) then begin
		allrs= sqrt((allxs-center(0))*(allxs-center(0)) + (allys-center(1))*(allys-center(1)) + (allzs-center(2))*(allzs-center(2)))

		maxrs= max(allrs)
		print, "max(allrs)= ", maxrs
		searchrs= rfact*max(allrs)
		if rfact le 0.1 then begin
			if searchrs gt 10.0 then searchrs= 10.0
		endif
		idx= where(allrs LE searchrs)

		allmasses= allmasses(idx)
		allvxs= allvxs(idx)
		allvys= allvys(idx)
		allvzs= allvzs(idx)
	endif

	comvel= fltarr(3)

        idx1= where(finite(allvxs) eq 0)
        idx2= where(finite(allvys) eq 0)
        idx3= where(finite(allvzs) eq 0)
        idx4= where(finite(allmasses) eq 0)

        idx1thru4= [idx1, idx2, idx3, idx4]
        idx= where(idx1thru4 gt 0)
        if idx(0) ne -1 then begin
	   idx1thru4= idx1thru4(where(idx1thru4 gt 0))
	   print, "fixing"
           print, "idx1thru4= ", idx1thru4
           print, "allvxs(idx1thru4)= ", allvxs(idx1thru4)
           print, "allvys(idx1thru4)= ", allvys(idx1thru4)
           print, "allvzs(idx1thru4)= ", allvzs(idx1thru4)
           print, "allmasses(idx1thru4)= ", allmasses(idx1thru4)
           allvxs(idx1thru4)= 0.0
           allvys(idx1thru4)= 0.0
           allvzs(idx1thru4)= 0.0
           allmasses(idx1thru4)= 0.0
        endif

        idx= where(abs(allvxs) gt 1.0e+6) & if idx(0) ne -1 then allvxs(idx)= 0.0 & if idx(0) ne -1 then help, "idx= ", idx
        idx= where(abs(allvys) gt 1.0e+6) & if idx(0) ne -1 then allvys(idx)= 0.0 & if idx(0) ne -1 then help, "idx= ", idx
        idx= where(abs(allvzs) gt 1.0e+6) & if idx(0) ne -1 then allvzs(idx)= 0.0 & if idx(0) ne -1 then help, "idx= ", idx
        idx= where(abs(allmasses) gt 1.0e+6) & if idx(0) ne -1 then allmasses(idx)= 0.0 & if idx(0) ne -1 then help, "idx= ", idx

	comvel[0]= total(allvxs*allmasses)/total(allmasses)
	comvel[1]= total(allvys*allmasses)/total(allmasses)
	comvel[2]= total(allvzs*allmasses)/total(allmasses)

	print, "com velocity   ",comvel, " (n= ",n_elements(allmasses),")"
	return, comvel

    endif

end


