function fload_1gal_gas_sfr, dummy, startid, numpart

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON SfrData

    if numpart eq 0 then return, [0]


    ;gid = id(0:npart(0)-1)
    gid = fload_gas_id(1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))

    return, sfr(idx)

end


