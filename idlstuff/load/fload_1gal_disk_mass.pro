function fload_1gal_disk_mass, dummy, startid, numpart

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON DiskData

    if numpart eq 0 then return, [0]

    gid = fload_disk_id(1)
    idx = where((gid ge long(startid)) and (gid lt long(startid)+long(numpart)))

    return, mdisk(idx)

end


