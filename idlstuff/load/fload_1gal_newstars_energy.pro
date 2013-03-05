function fload_1gal_newstars_energy, startid, numpart, comv=comv, total=total, $
                                        kinetic=kinetic, $
                                        potential=potential, $
					specific=specific

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON NewStarData
    COMMON PotData

    if keyword_set(comv) then begin
        c=comv
    endif else begin
        c=[0,0,0]
    endelse

    if npart(4) le 0 then begin
        print, " "
        print, " no new star particles"
        print, " "
        return, [0]
    endif


        ; is there stuff
        ; -------------------------
        allmasses = [mstars]
        allvxs = [vxstars]
        allvys = [vystars]
        allvzs = [vzstars]


    nsid= fload_newstars_id(1)

    if keyword_set(useparentinfo) then begin
        ; now instead of grabbing new star id number's
        ;  we'll grab their parent id's
        if flag_stargens eq 1 and flag_parentid eq 1 then begin
                nsid= [parentid]
        endif else begin
                nsstartid= npart(0)+npart(1)+npart(2)+npart(3)
                N= total(npart)
                nsid= id(nsstartid:N-1)
        endelse
    endif


    idx = where((nsid ge long(startid)) and (nsid lt long(startid)+long(numpart)))

    allmasses= allmasses(idx)
    allvxs= allvxs(idx)
    allvys= allvys(idx)
    allvzs= allvzs(idx)


        v2= (allvxs-c(0))*(allvxs-c(0)) + (allvys-c(1))*(allvys-c(1)) + (allvzs-c(2))*(allvzs-c(2))

        kinetic_energy= 0.5*allmasses*v2

    if keyword_set(specific) then kinetic_energy= kinetic_energy / allmasses
    if keyword_set(kinetic) then return, kinetic_energy

    potential_energy= 0.5*pstars*mstars

    if keyword_set(specific) then potential_energy= potential_energy / allmasses
    if keyword_set(potential) then return, potential_energy


    total_energy= potential_energy+kinetic_energy

    return, total_energy

end


