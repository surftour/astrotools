function fload_bulge_mass, dummy

    COMMON GalaxyHeader
    COMMON BulgeData

    if dummy EQ 1 then begin
	if npart(3) EQ 0 then return, 0
	return, mbulge
    endif

end


