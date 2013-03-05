function fload_gas_potential, dummy

    COMMON GalaxyHeader
    COMMON PotData


    if dummy EQ 1 and flag_snaphaspot then begin
	return, pgas
    endif else begin
	return, 0
    endelse

end



