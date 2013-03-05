Function fload_1gal_center, startid, numpart
;-------------------------------------------
;  WARNING:
;    Sometimes the algorithm works better
; if only one component of the center is
; searched for.  i.e. sometimes the baryonic
; center is more stable than the dm one.
;
;-------------------------------------------


    if numpart GT 0 then begin
	com = [0.0,0.0,0.0]
	lastcom = [0.0,0.0,0.0]
	dcom = 100.0


	fload_1gal_all_data, startid, numpart, allmasses, allxs, allys, allzs, allvxs, allvys, allvzs

	allrs = sqrt(allxs*allxs + allys*allys + allzs*allzs)
	rmax = max(allrs)
	
	while (((dcom GT 0.01) AND (n_elements(allmasses) GT 1000)) OR (rmax GT 20.0)) do begin
	
	; calculate distance from com and select less than some radius
	   allrs = sqrt((allxs-com(0))*(allxs-com(0)) + (allys-com(1))*(allys-com(1)) + (allzs-com(2))*(allzs-com(2)))
	   idx = where(allrs LE rmax)
	   allmasses = allmasses(idx)
	   allxs = allxs(idx)
	   allys = allys(idx)
	   allzs = allzs(idx)

	; save last com and find new center of mass
	   lastcom = com
	   com(0) = total(allmasses*allxs)/total(allmasses)
	   com(1) = total(allmasses*allys)/total(allmasses)
	   com(2) = total(allmasses*allzs)/total(allmasses)

	; see what error is from 0
	   comdiff = com-lastcom
	   dcom = sqrt(total(comdiff*comdiff))

	; cut search circle in half
	   rmax = rmax / 2.0
	endwhile

	print, "gal1 center         ",com
	return, com
    endif

end


