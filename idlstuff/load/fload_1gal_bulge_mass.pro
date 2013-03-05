function fload_1gal_bulge_mass, dummy, startid, numpart

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON BulgeData

    if numpart eq 0 then return, [0]
    if npart(3) le 0 then return, [0]

    ;gid = id(0:npart(0)-1)
    gid = fload_bulge_id(1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))

    return, mbulge(idx)

end


