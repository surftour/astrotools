function fload_halo_kinetic_e, dummy, comv=comv, total=total

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON HaloData

    if keyword_set(comv) then begin
	c=comv
    endif else begin
	c=[0,0,0]
    endelse

    if npart(1) GT 0 then begin
        ; there must be gas!
        ; -------------------------
        allmasses = [mhalo]
        allvxs = [vxhalo]
        allvys = [vyhalo]
        allvzs = [vzhalo]


	v2= (allvxs-c(0))*(allvxs-c(0)) + (allvys-c(1))*(allvys-c(1)) + (allvzs-c(2))*(allvzs-c(2))

	kinetic_energy= 0.5*allmasses*v2
	total_kinetic_energy= total(kinetic_energy)

	return, kinetic_energy
    endif else begin
	return, [0]
    endelse

end


