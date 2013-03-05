function fload_blackhole_accrate, frun, Ttime, idtofollow=idtofollow

    COMMON GalaxyHeader
    COMMON BlackHoleData
    COMMON OtherData

    Nbh= npart(5)

    if Ttime lt 0 or Nbh le 0 then return, [0]

    bhids= fload_blackhole_id(1)

    bhaccrate= fltarr(Nbh)

    ; get bolometric luminosity from the accretion rate
    ; --------------------------------------------------
    fload_blackhole_details, frun, ids=bhdids, time=bhdtime, bhmass=bhdmass, mdot=bhdmdot, local_rho=lr, local_c=lc

    for i=1,Nbh do begin
	idx= where(bhdids eq bhids[i-1])
	this_time= bhdtime(idx)
	this_bhmass= bhdmass(idx)
	this_mdot= 10.2247*bhdmdot(idx)   ; 10.2247 converts from gadget units to m_solar/yr

	sorted_time= this_time(sort(this_time))
	sorted_bhmass= this_bhmass(sort(this_time))
	sorted_mdot= this_mdot(sort(this_time))

	idx= where(sorted_time ge Ttime)
	if idx(0) ne -1 then begin
		; instantaneous
		current_mdot= sorted_mdot(idx(0))
		current_bhmass= sorted_bhmass(idx(0))

		; averaged
		lowidx= idx(0)
		if n_elements(idx) gt 5 then highidx= idx(4) else highidx= lowidx
		avg_c_bhmass= mean(sorted_bhmass(lowidx:highidx))
		avg_c_mdot= mean(sorted_mdot(lowidx:highidx))

		; cgs units, cm, sec, g
		L_solar= 3.826d+33                               ; in ergs per sec
		cc= 2.9979d+10                                 ; speed of light in cm/sec
		;convert_sunyr_gsec= 6.446d+27                  ; convert msun/yr -> g/sec
		convert_sunyr_gsec= 6.30428d+25                  ; convert msun/yr -> g/sec

		;boloL_ergss= 0.1*current_mdot*cc*cc*convert_sunyr_gsec
		BHAccEfficiency= 0.1
		boloL_ergss= BHAccEfficiency*avg_c_mdot*cc*cc*convert_sunyr_gsec
		boloL_sun= boloL_ergss / L_solar

		print, "Time= ", Ttime, "  time= ", sorted_time(idx(0))
		print, "BH Id= ", bhids[i-1]
		print, " Mass= ", current_bhmass;,   "    Averge of 5= ", mean(sorted_bhmass(idx(0):idx(4)))
		print, "  averaged mass= ", avg_c_bhmass
		print, " Mdot= ", current_mdot;,   "    Averge of 5= ", mean(sorted_mdot(idx(0):idx(4)))
		print, "  averaged mdot= ", avg_c_mdot
		print, "BoloLum (L_sun) ", boloL_sun
		print, "BoloLum (erg/s) ", boloL_ergss

		bhaccrate[i-1]= avg_c_mdot
	endif else begin
		print, "Problem with BH ID= ",bhids[i-1]
		print, "  can't match the two times"
		print, "  sorted_time",sorted_time(n_elements(sorted_time)-1)
		print, "  Ttime=",Ttime
		print, " "
		bhaccrate[i-1]= 0.0
	endelse
    endfor


    ; default is bolometric
    return, bhaccrate


end






