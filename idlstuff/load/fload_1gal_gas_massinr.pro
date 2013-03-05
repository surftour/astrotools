function fload_1gal_gas_massinr, radius, startid, numpart, gas_startid, gas_numpart, $
				total_mgas=total_mgas, mfs_inr=mfs_inr, total_mfs=total_mfs, $
				center=center

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON GasData
    COMMON SfrData

    if gas_numpart eq 0 then return, [0]


    if not keyword_set(center) then begin
	c=fload_1gal_center(startid,numpart)
    endif else begin
	c=center
    endelse

    ; i did this to only get gas particles, and not new star particles
    gid = id(0:npart(0)-1)
    idx = where((gid ge long(gas_startid)) and (gid lt long(gas_startid)+long(gas_numpart)))

    g1xgas = xgas(idx)
    g1ygas = ygas(idx)
    g1zgas = zgas(idx)
    g1mgas = mgas(idx)
    g1mfs = mfs(idx)

    ; return the total gas (gas)
    if keyword_set(total_mgas) then begin
	total_mgas= total(g1mgas)-total(g1mfs)
    endif

    rgas= sqrt((g1xgas-c(0))*(g1xgas-c(0)) + (g1ygas-c(1))*(g1ygas-c(1)) + (g1zgas-c(2))*(g1zgas-c(2)))

    idx= where(rgas lt radius)
    if idx(0) GE 0 then mgas_inr= g1mgas(idx)-g1mfs(idx) else mgas_inr= 0.0


    tmgas= total(mgas_inr)

    print, n_elements(idx), " part. < ",radius, " kpc/h  mass= ",tmgas,"   ,(Total:",total_mgas,")"
    return, tmgas

end


