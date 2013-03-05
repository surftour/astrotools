function fload_gas_kinetic_e, dummy, comv=comv, total=total, idtofollow=idtofollow

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON GasData

    if keyword_set(comv) then begin
	c=comv
    endif else begin
	c=[0,0,0]
    endelse


        ; there must be gas!
        ; -------------------------
        allmasses = [mgas]
        allvxs = [vxgas]
        allvys = [vygas]
        allvzs = [vzgas]


    if dummy eq 1 then begin
	v2= (allvxs-c(0))*(allvxs-c(0)) + (allvys-c(1))*(allvys-c(1)) + (allvzs-c(2))*(allvzs-c(2))

	kinetic_energy= 0.5*allmasses*v2
	total_kinetic_energy= total(kinetic_energy)

	;print, n_elements(allmasses), " particles have total kinetic energy      ", total_kinetic_energy
	return, kinetic_energy
    endif



    if keyword_set(idtofollow) then begin

        gas_ids= id(0:npart(0)-1)
        idx= where(gas_ids EQ idtofollow)
        if idx(0) EQ -1 or n_elements(idx) NE 1 then begin
                print, "Hey, what's up with idx?  printing"
                print, idx
                return, [0,0,0]
        endif

        vx= allvxs(idx)
        vy= allvys(idx)
        vz= allvzs(idx)
	m= allmasses(idx)

	v2= vx*vx + vy*y + vz*vz
        if dummy EQ '99' then return, 0.5*m*v2

    endif


end


