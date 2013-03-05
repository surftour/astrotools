function fload_gas_id, dummy

    COMMON GalaxyHeader
    COMMON OtherData


    if dummy GT 0 then begin

	if npart(0) eq 0 then return, 0

	gids= id(0:npart(0)-1)

	; -------------------------------------------------
	;  Old Method
	;
	;
        ; recover same ID as original gas particle
        ;print, "HEY !!!!!!  WARNING:  turned off new star ID correction"
        ;whacky_ids= where((gids lt 0) or (gids gt 5.0e7))
        ;num= long(1)
        ;if whacky_ids(0) ne -1 then gids(whacky_ids)= gids(whacky_ids) - ishft(num,31)

        ;return, gids
	; -------------------------------------------------

        ns_generations= 0 * gids

        ;limit_lower= 0L
        ;limit_upper= long(5.0e7)

        ; binary stuff
        ; IDL> print, binary(ishft(num,0))
        ; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
        ; IDL> print, binary(ishft(num,31))
        ; 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0


        ; find the number of generations (should never get to 10 bits, i.e., 1024 generations)
        for nbit= 1, 10L do begin
                one_exist= where(ishft(gids, -32+nbit) eq 1)
                if one_exist(0) eq -1 then break
        endfor
        nbit = nbit -1    ; take out the failed attempt that triggered the break
        print, "nbit= ", nbit
        print, "expected number of generations= ", max(ishft(gids, -32+nbit))


        ns_generations= ishft(gids, -32+nbit)


        ; recover same ID as original gas particle ( essentially
        ;    undoing the bit flipping done in sfr_eff.c )
        ;
        ; note: sfr_eff.c affects the first x bits, where x is the miminum
        ;       number of bits needed to keep track of the number of generations,
        ;       i.e., 10 generations needs 4 bits, 8 generations needs 3 bits, etc.
        ;
        ;

        good_ids= gids - ishft(ns_generations, 32-nbit)

        ; could check that they're good if needed
        ;good_ids= where((gids ge limit_lower) and (gids le limit_upper))


        return, good_ids

    endif

end

