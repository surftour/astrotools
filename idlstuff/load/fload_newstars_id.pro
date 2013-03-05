function fload_newstars_id, dummy, ns_generations=ns_generations


    COMMON GalaxyHeader
    COMMON OtherData


    if dummy NE -1 then begin
        if npart(4) EQ 0 then return, 0

        sid= npart(0)+npart(1)+npart(2)+npart(3)
        eid= npart(0)+npart(1)+npart(2)+npart(3)+npart(4)-1
        stars_ids= id(sid:eid)

	ns_generations= 0 * stars_ids

	;limit_lower= 0L
	;limit_upper= long(5.0e7)

	; binary stuff
	; num= 1L
	; IDL> print, binary(ishft(num,0))
	; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
	; IDL> print, binary(ishft(num,31))
	; 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0


	; find the number of generations (should never get to 10 bits, i.e., 1024 generations)
	for nbit= 1, 10L do begin
		one_exist= where(ishft(stars_ids, -32+nbit) eq 1)
		if one_exist(0) eq -1 then break
	endfor
	nbit = nbit -1    ; take out the failed attempt that triggered the break
	print, "nbit= ", nbit
	print, "expected number of generations= ", max(ishft(stars_ids, -32+nbit))
	

	ns_generations= ishft(stars_ids, -32+nbit)


	; recover same ID as original gas particle ( essentially
	;    undoing the bit flipping done in sfr_eff.c )
	;
	; note: sfr_eff.c affects the first x bits, where x is the miminum
	;       number of bits needed to keep track of the number of generations,
	;       i.e., 10 generations needs 4 bits, 8 generations needs 3 bits, etc.
	;
	;

	good_ids= stars_ids - ishft(ns_generations, 32-nbit)

	; could check that they're good if needed
	;good_ids= where((stars_ids ge limit_lower) and (stars_ids le limit_upper))

    endif

        return, good_ids


end


