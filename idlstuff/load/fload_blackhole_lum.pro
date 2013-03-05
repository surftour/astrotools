function fload_blackhole_lum, frun, Ttime, idtofollow=idtofollow, $
				bolometric=bolometric, $
				Bband=Bband, $
				softX=softX, hardX=hardX

    COMMON GalaxyHeader
    COMMON BlackHoleData
    COMMON OtherData

    Nbh= npart(5)

    if Ttime lt 0 or Nbh le 0 then return, [0]

    bhids= fload_blackhole_id(1)

    bhbololum= fltarr(Nbh)

    ; get bolometric luminosity from the accretion rate
    ; --------------------------------------------------
    fload_blackhole_details, frun, ids=bhdids, time=bhdtime, bhmass=bhdmass, mdot=bhdmdot, local_rho=lr, local_c=lc

    ; don't have blackhole_details file
    if bhdtime(0) lt 0.0 then begin
	; old method just assumed something
	;bhbololum(*)= 10.0
    	;return, bhbololum

	; used blackholes.txt, if possible
	;
	; need to .run bh_multi
	open_blackhole_txt, frun, bhtime, bh_num, bh_mass, bh_mdot_gu, bh_mdot_sunyr, $
                                bh_totalmass, bh_mdot_edd
	Nbh= 1
	bhdids= 0.0*bhtime + bhids[0]
	bhdtime= bhtime
	bhdmass= bh_mass  / 1.0d+10   ; the above multiplied by this
	bhdmdot= bh_mdot_gu
	
    endif

    for i=1,Nbh do begin
	idx= where(bhdids eq bhids[i-1])
	this_time= bhdtime(idx)
	this_bhmass= bhdmass(idx)
	this_mdot= 10.2247*bhdmdot(idx)   ; 10.2247 converts from gadget units to m_solar/yr

	sorted_time= this_time(sort(this_time))
	sorted_bhmass= this_bhmass(sort(this_time))
	sorted_mdot= this_mdot(sort(this_time))

	idx= where(sorted_time ge Ttime)

	; trap for the case of snapshot time and last timestep time
	; are close but timestep is slightly more
	if idx(0) eq -1 then begin
		lastidx= n_elements(sorted_time)-1
		tdiff= double(Ttime) - double(sorted_time(lastidx))
		if tdiff lt 0.0001 then idx(0)= lastidx
	endif

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
		boloL_ergss= 0.1*avg_c_mdot*cc*cc*convert_sunyr_gsec
		boloL_sun= boloL_ergss / L_solar

		print, "Time= ", Ttime, "  time= ", sorted_time(idx(0))
		print, "BH Id= ", bhids[i-1]
		print, " Mass= ", current_bhmass;,   "    Averge of 5= ", mean(sorted_bhmass(idx(0):idx(4)))
		print, "  averaged mass= ", avg_c_bhmass
		print, " Mdot= ", current_mdot;,   "    Averge of 5= ", mean(sorted_mdot(idx(0):idx(4)))
		print, "  averaged mdot= ", avg_c_mdot
		print, "BoloLum (L_sun) ", boloL_sun
		print, "BoloLum (erg/s) ", boloL_ergss

		bhbololum[i-1]= boloL_sun
	endif else begin
		print, "Problem with BH ID= ",bhids[i-1]
		print, "can't match sorted_time with T=",time
		bhbololum[i-1]= 0.0
	endelse
    endfor


    ; default is bolometric
    return, bhbololum


end






