function fload_gas_energy, junk, comv=comv, idtofollow=idtofollow, $
				kinetic=kinetic, $
				potential=potential, $
				radial_kinetic=radial_kinetic, $
				thermal=thermal, $
				specific=specific

    COMMON GalaxyHeader
    COMMON OtherData
    COMMON GasData
    COMMON FeedbackData
    COMMON PotData

    if keyword_set(comv) then begin
	c=comv
    endif else begin
	c=[0,0,0]
    endelse

    if npart(0) le 0 then begin
	print, " "
	print, " no gas particles"
	print, " "
	return, [0]
    endif

        ; there must be gas!
        ; -------------------------
        allmasses = [mgas]
        allvxs = [vxgas]
        allvys = [vygas]
        allvzs = [vzgas]
	allus = [u]


    v2= (allvxs-c(0))*(allvxs-c(0)) + (allvys-c(1))*(allvys-c(1)) + (allvzs-c(2))*(allvzs-c(2))
    kinetic_energy= 0.5*allmasses*v2


    if keyword_set(specific) then kinetic_energy= kinetic_energy/allmasses
    if keyword_set(kinetic) then return, kinetic_energy


    v_r= fload_gas_v('r', comvel=c)
    kinetic_energy_radial= 0.5*allmasses*v_r*v_r

    if keyword_set(specific) then kinetic_energy_radial= kinetic_energy_radial/allmasses
    if keyword_set(kinetic_radial) then return, kinetic_energy_radial


    potential_energy= 0.5*pgas*mgas

    if keyword_set(specific) then potential_energy= potential_energy/allmasses
    if keyword_set(potential) then return, potential_energy

    thermal_energy= allus*allmasses

    if keyword_set(specific) then thermal_energy= thermal_energy/allmasses
    if keyword_set(thermal) then return, thermal_energy


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

	v2= vx*vx + vy*vy + vz*vz

	ke= 0.5*m*v2
	pe= m*pgas(idx)
	th= u(idx)*m
	;th_q= u(idx)*m + tpu(idx)*m

	etot= ke + pe + th
	
	if dummy EQ 'all' then begin
		alles= [etot,ke,pe,th,th_q]
		return, alles
	endif

    endif


    e= potential_energy+kinetic_energy+thermal_energy
    return, e

end


