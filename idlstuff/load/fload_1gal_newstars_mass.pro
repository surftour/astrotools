function fload_1gal_newstars_mass, startid, numpart, $
					useparentinfo=useparentinfo

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON NewStarData
    COMMON Center
    COMMON ParentID


    if npart(4) le 0 then return, [0]


    nsid= fload_newstars_id(1)

    if keyword_set(useparentinfo) then begin
	; now instead of grabbing new star id number's
	;  we'll grab their parent id's
	if flag_stargens eq 1 and flag_parentid eq 1 then begin
		nsid= [long(parentid)]
	endif else begin
		nsstartid= npart(0)+npart(1)+npart(2)+npart(3)
		N= total(npart)
		nsid= id(nsstartid:N-1)
	endelse
    endif


    idx = where((nsid ge long(startid)) and (nsid lt long(startid)+long(numpart)))


    allmasses= mstars(idx)

    return, allmasses


end


