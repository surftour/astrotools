function fload_all_kinetic_e, dummy, comv=comv, total=total

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON HaloData
    COMMON GasData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData

    if keyword_set(comv) then begin
	c=comv
    endif else begin
	c=[0,0,0]
    endelse

        ; is there gas!
        ; -------------------------
	if npart(0) GT 0 then begin
          allmasses = [mgas]
          allvxs = [vxgas]
          allvys = [vygas]
          allvzs = [vzgas]
	endif

        ; is there a halo?
        ; -----------------
        if npart(1) GT 0 then begin
	  if npart(0) EQ 0 then begin    ; require either gas or halo
            allmasses = [mhalo]
            allvxs = [vxhalo]
            allvys = [vyhalo]
            allvzs = [vzhalo]
	  endif else begin
            allmasses = [allmasses, mhalo]
            allvxs = [allvxs, vxhalo]
            allvys = [allvys, vyhalo]
            allvzs = [allvzs, vzhalo]
	  endelse
        endif

        ; is there a disk?
        ; -----------------
        if npart(2) GT 0 then begin
          allmasses = [allmasses, mdisk]
          allvxs = [allvxs, vxdisk]
          allvys = [allvys, vydisk]
          allvzs = [allvzs, vzdisk]
        endif

        ; is there a bulge?
        ; -----------------
        if npart(3) GT 0 then begin
          allmasses = [allmasses, mbulge]
          allvxs = [allvxs, vxbulge]
          allvys = [allvys, vybulge]
          allvzs = [allvzs, vzbulge]
        endif

        ; is there stars?
        ; -----------------
        if npart(4) GT 0 then begin
          allmasses = [allmasses, mstars]
          allvxs = [allvxs, vxstars]
          allvys = [allvys, vystars]
          allvzs = [allvzs, vzstars]
        endif


    v2= (allvxs-c(0))*(allvxs-c(0)) + (allvys-c(1))*(allvys-c(1)) + (allvzs-c(2))*(allvzs-c(2))

    kinetic_energy= 0.5*allmasses*v2
    total_kinetic_energy= total(kinetic_energy)

    print, n_elements(allmasses), " particles have total kinetic energy      ", total_kinetic_energy

    if keyword_set(total) then return, total_kinetic_energy

    ; else return it all
    return, kinetic_energy

end


