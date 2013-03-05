function fload_1gal_kinetic_e, startid, numpart, comv=comv

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

        ; there must be gas!
        ; -------------------------
        allmasses = [mgas]
        allvxs = [vxgas]
        allvys = [vygas]
        allvzs = [vzgas]

        ; is there a halo?
        ; -----------------
        if npart(1) GT 0 then begin
          allmasses = [allmasses, mhalo]
          allvxs = [allvxs, vxhalo]
          allvys = [allvys, vyhalo]
          allvzs = [allvzs, vzhalo]
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


    ; just get specific id's
    ; -------------------------
    idx = where((id ge long(startid)) and (id lt long(startid)+long(numpart)))

    mass = allmasses(idx)
    vx = allvxs(idx)
    vy = allvys(idx)
    vz = allvzs(idx)

    v2= (vx-c(0))*(vx-c(0)) + (vy-c(1))*(vy-c(1)) + (vz-c(2))*(vz-c(2))

    kinetic_energy= 0.5*mass*v2
    total_kinetic_energy= total(kinetic_energy)

    print, n_elements(idx), " particles have total kinetic energy", total_kinetic_energy
    return, kinetic_energy

end


