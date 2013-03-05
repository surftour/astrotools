function fload_all_mass, dummy


    COMMON GalaxyHeader
    COMMON GasData
    COMMON HaloData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData


    if dummy NE -276 then begin
	if npart(0) GT 0 then allmasses= [mgas]
	if npart(1) GT 0 then begin
	   if npart(0) EQ 0 then begin
		allmasses= [mhalo]
	   endif else begin
		allmasses= [allmasses, mhalo]
	   endelse
	endif
	if npart(2) GT 0 then allmasses= [allmasses, mdisk]
	if npart(3) GT 0 then allmasses= [allmasses, mbulge]
	if npart(4) GT 0 then allmasses= [allmasses, mstars]
	return, allmasses
    endif

end




