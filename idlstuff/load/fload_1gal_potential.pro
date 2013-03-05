Function fload_1gal_potential, startid, numpart

    COMMON GalaxyHeader
    COMMON PotData
    COMMON OtherData

    if numpart GT 0 then begin

        ; we know there are these
        ; -------------------------
        allps = [pgas, phalo, pdisk]

        ; is there a bulge?
        ; -----------------
        if npart(3) GT 0 then allps = [allps, pbulge]

        ; is there stars?
        ; -----------------
        if npart(4) GT 0 then allps = [allps, pstars]

	
        ; just get the particles which have id's for galaxy 1
        ; -------------------------------------------------------
        idx = where((id ge startid) and (id lt startid+numpart))
        allps = allps(idx)


; ----------------------------------------------------------------
; ----------------------------------------------------------------

	return, allps

    endif

end


