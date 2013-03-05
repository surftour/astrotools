function fload_1gal_halo_massinr, radius, startid, numpart, $
				total_dm=total_dm, $
				center=center

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON HaloData

    if numpart eq 0 then return, [0]

    if npart(1) eq 0 then return, [0]


    if not keyword_set(center) then begin
	c=fload_1gal_center(startid,numpart)
    endif else begin
	c=center
    endelse


    ; do this to only get disk particles
    if npart(0) gt 0 then begin
	startingid= npart(0)-1
        endid= startingid+npart(1)
    endif else begin
        startingid= 0
        endid= npart(1)-1
    endelse
    hid = id(startingid:endid)
    idx = where((hid ge long(startid)) and (hid lt long(startid)+long(numpart)))

    g1x = xhalo(idx)
    g1y = yhalo(idx)
    g1z = zhalo(idx)
    g1m = mhalo(idx)

    ; return the total dark matter 
    if keyword_set(total_dm) then begin
	total_dm= total(g1m)
    endif

    r= sqrt((g1x-c(0))*(g1x-c(0)) + (g1y-c(1))*(g1y-c(1)) + (g1z-c(2))*(g1z-c(2)))

    idx= where(r lt radius)
    if idx(0) GE 0 then m_inr= g1m(idx) else m_inr= 0.0


    tm= total(m_inr)

    print, n_elements(idx), " part. < ",radius, " kpc/h  mass= ",tm,"   ,(Total:",total_dm,")"
    return, tm

end


