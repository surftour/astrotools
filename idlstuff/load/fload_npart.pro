function fload_npart, dummy


    COMMON GalaxyHeader


    if dummy EQ 0 then return, npart(0)
    if dummy EQ 1 then return, npart(1)
    if dummy EQ 2 then return, npart(2)
    if dummy EQ 3 then return, npart(3)
    if dummy EQ 4 then return, npart(4)
    if dummy EQ 5 then return, npart(5)

    if dummy GT 5 then return, npart

end


