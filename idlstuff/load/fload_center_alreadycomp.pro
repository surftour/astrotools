function fload_center_alreadycomp, time, denorbh=denorbh, bhlbl=bhlbl


    COMMON Center


    if keyword_set(denorbh) and keyword_set(bhlbl) then begin
	fname=fload_frun(1)+"/centers_bh_"+bhlbl+".txt"
	if denorbh eq "den" then fname=fload_frun(1)+"/centers_den_"+bhlbl+".txt"

	read_file_centers, ftime, fcenter, filename=fname

	idx=where(ftime ge (time-1.0e-6))
	if idx(0) eq -1 then begin
		print, " "
		print, " PROBLEM: can't find time in centers file"
		print, " "
		return, com
	endif

	return, fcenter[*,idx(0)]

    endif

    return, com

end


