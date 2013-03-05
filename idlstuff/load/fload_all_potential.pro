Function fload_all_potential, dummy, total=total

    if not keyword_set(dummy) then dummy= 1

    COMMON GalaxyHeader
    COMMON HaloData
    COMMON GasData
    COMMON DiskData
    COMMON BulgeData
    COMMON NewStarData
    COMMON PotData

    if dummy NE -278 and flag_snaphaspot then begin
        if npart(0) GT 0 then begin allpotential= [pgas] & allmasses= [mgas] & endif
        if npart(1) GT 0 then begin
           if npart(0) EQ 0 then begin        ; require either gas or halo
                allpotential= [phalo] & allmasses= [mhalo]
           endif else begin
                allpotential= [allpotential, phalo] & allmasses= [allmasses, mhalo]
           endelse
        endif
        if npart(2) GT 0 then begin allpotential= [allpotential, pdisk] & allmasses= [allmasses, mdisk] & endif
        if npart(3) GT 0 then begin allpotential= [allpotential, pbulge] & allmasses= [allmasses, mbulge] & endif
        if npart(4) GT 0 then begin allpotential= [allpotential, pstars] & allmasses= [allmasses, mstars] & endif


	total_potential_energy= 0.5 * total(allpotential * allmasses)

	print, n_elements(allpotential), " particles have total potential energy    ", total_potential_energy

	if keyword_set(total) then return, total_potential_energy

        return, allpotential
    endif


end


