function fload_1gal_bulge_massinr, radius, startid, numpart, $
				total_b=total_b, $
				center=center

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON BulgeData

    if numpart eq 0 then return, [0]

    if npart(3) eq 0 then return, [0]


    if not keyword_set(center) then begin
	c=fload_1gal_center(startid,numpart)
    endif else begin
	c=center
    endelse

    ; do this to only get bulge particles
    startingid= npart(0)+npart(1)+npart(2)-1
    endid= startingid+npart(3)
    did = id(startingid:endid)
    idx = where((did ge long(startid)) and (did lt long(startid)+long(numpart)))

    if idx(0) eq -1 then return, [0]

    g1x = xbulge(idx)
    g1y = ybulge(idx)
    g1z = zbulge(idx)
    g1m = mbulge(idx)

    ; return the total dark matter 
    ;if keyword_set(total_b) then begin
	total_b= total(g1m)
    ;endif

    r= sqrt((g1x-c(0))*(g1x-c(0)) + (g1y-c(1))*(g1y-c(1)) + (g1z-c(2))*(g1z-c(2)))

    idx= where(r lt radius)
    if idx(0) GE 0 then m_inr= g1m(idx) else m_inr= 0.0


    tm= total(m_inr)

    print, n_elements(idx), " part. < ",radius, " kpc/h  mass= ",tm,"   ,(Total:",total_b,")"
    return, tm

end


