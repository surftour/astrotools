function fload_1gal_halo_mass, startid, numpart

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON HaloData

    if numpart eq 0 then return, [0]

    ;gid = id(0:npart(0)-1)
    gid = fload_halo_id(1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))

    return, mhalo(idx)

end


