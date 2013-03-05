function fload_newstars_energy, dummy, comv=comv, total=total, $
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


