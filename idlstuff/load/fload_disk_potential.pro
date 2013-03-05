function fload_disk_potential, dummy

    COMMON GalaxyHeader
    COMMON PotData


    if dummy EQ 1 and flag_snaphaspot then begin
	return, pdisk
    endif else begin
	return, 0
    endelse

end



