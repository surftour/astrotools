function fload_halo_potential, dummy

    COMMON GalaxyHeader
    COMMON PotData

    if dummy EQ 1 and flag_snaphaspot then begin
	return, phalo
    endif else begin
	return, 0
    endelse

end



