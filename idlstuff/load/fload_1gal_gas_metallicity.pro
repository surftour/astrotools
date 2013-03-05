function fload_1gal_gas_metallicity, dummy, startid, numpart

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON SfrData

    if numpart eq 0 then return, [0]
    if flag_metals le 0 then return, [0]


    gid = id(0:npart(0)-1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))

    return, gmetals(idx)

end


