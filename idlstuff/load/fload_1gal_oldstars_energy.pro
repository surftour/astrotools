function fload_1gal_oldstars_energy, startid, numpart, comv=comv, total=total, $
                                        kinetic=kinetic, $
                                        potential=potential, $
					specific=specific


    COMMON GalaxyHeader
    COMMON OtherData
    COMMON DiskData
    COMMON BulgeData
    COMMON PotData

    if keyword_set(comv) then begin
        c=comv
    endif else begin
        c=[0,0,0]
    endelse 

    if npart(2)+npart(3) le 0 then begin
        print, " "
        print, " no old star particles"
        print, " "
        return, [0]
    endif


        ; is there stuff
        ; -------------------------
        if npart(2) gt 0 then begin
           allmasses = [mdisk]
           allvxs = [vxdisk]
           allvys = [vydisk]
           allvzs = [vzdisk]
           allps = [pdisk]
        endif

        if npart(3) gt 0 then begin
          if n_elements(allmasses) gt 0 then begin
             allmasses = [allmasses,mbulge]
             allvxs = [allvxs,vxbulge]
             allvys = [allvys,vybulge]
             allvzs = [allvzs,vzbulge]
             allps = [allps,pbulge]
          endif else begin
             allmasses = [mbulge]
             allvxs = [vxbulge]
             allvys = [vybulge]
             allvzs = [vzbulge]
             allps = [pbulge]
          endelse
        endif

    ; first grab id's for all baryons
    ;if (npart(2)+npart(3) gt 0) then begin
    ;    gid= [id(npart(0)+npart(1):N-1-npart(4))]
    ;endif else begin
    ;    print, "PROBLEM"
    ;    return, [0]
    ;endelse

    gid= fload_oldstars_id(1)


    ; then grab appropriate ones for galaxy 1
    idx = where((gid ge startid) and (gid lt startid+numpart))


    allmasses= allmasses(idx)
    allvxs= allvxs(idx)
    allvys= allvys(idx)
    allvzs= allvzs(idx)
    allps= allps(idx)


    v2= (allvxs-c(0))*(allvxs-c(0)) + (allvys-c(1))*(allvys-c(1)) + (allvzs-c(2))*(allvzs-c(2))

    kinetic_energy= 0.5*allmasses*v2

    if keyword_set(specific) then kinetic_energy= kinetic_energy / allmasses
    if keyword_set(kinetic) then return, kinetic_energy

    potential_energy= 0.5*allps*allmasses

    if keyword_set(specific) then potential_energy= potential_energy / allmasses
    if keyword_set(potential) then return, potential_energy


    total_energy= potential_energy+kinetic_energy

    return, total_energy

end





