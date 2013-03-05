function fload_1gal_disk_massinr, radius, startid, numpart, $
				total_d=total_d, $
				center=center

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON DiskData

    if numpart eq 0 then return, [0]

    if npart(2) eq 0 then return, [0]


    if not keyword_set(center) then begin
	c=fload_1gal_center(startid,numpart)
    endif else begin
	c=center
    endelse

    ; do this to only get disk particles
    startingid= npart(0)+npart(1)-1
    endid= startingid+npart(2)
    did = id(startingid:endid)
    idx = where((did ge long(startid)) and (did lt long(startid)+long(numpart)))

    g1x = xdisk(idx)
    g1y = ydisk(idx)
    g1z = zdisk(idx)
    g1m = mdisk(idx)

    ; return the total dark matter 
    ;if keyword_set(total_d) then begin
	total_d= total(g1m)
    ;endif

    r= sqrt((g1x-c(0))*(g1x-c(0)) + (g1y-c(1))*(g1y-c(1)) + (g1z-c(2))*(g1z-c(2)))

    idx= where(r lt radius)
    if idx(0) GE 0 then m_inr= g1m(idx) else m_inr= 0.0


    tm= total(m_inr)

    print, n_elements(idx), " part. < ",radius, " kpc/h  mass= ",tm,"   ,(Total:",total_d,")"
    return, tm

end


