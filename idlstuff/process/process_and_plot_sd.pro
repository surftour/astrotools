;================================================================================
;
;
;
pro process_and_plot_sd, frun, snapnum, proftype, $
			xmin, xmax, bins, $
			x_is_devac=x_is_devac, $
			x_is_log=x_is_log, $
			linecolor=linecolor, $
			includeerr=includeerr, $
			fitdevac=fitdevac, $
			fitsersic=fitsersic, $
			frommanyproj=frommanyproj, $
			h=h, $
			thispsym=thispsym, thisthick=thisthick, thislinest=thislinest, $
			special_x=special_x, $
			special_y=special_y, $
			special_z=special_z, $
			special_m=special_m


	if snapnum lt 0 then begin
		snapnum= fload_frun_nsnaps(frun) - 2
	endif


	if frun ne "loaded" then begin
		;ok=fload_snapshot_bh(frun,snapnum)
		ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap)
		;ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap,/arepo)
	endif


	; -----------------------------------------------

	if not keyword_set(proftype) then proftype= "allstars"

        ; total (all star) profile
        ; --------------------------
	if proftype eq "allstars" then begin
        	x=fload_allstars_xyz('x')
        	y=fload_allstars_xyz('y')
        	z=fload_allstars_xyz('z')
        	m=fload_allstars_mass(1)
	endif

        ; gas profile
        ; --------------
        if proftype eq "gas" then begin
		if fload_npart(0) le 0 then return
		x=fload_gas_xyz('x')
		y=fload_gas_xyz('y')
		z=fload_gas_xyz('z')
		m=fload_gas_mass(1)
        endif


        ; disk profile
        ; --------------
        if proftype eq "disk" then begin
		if fload_npart(2) le 0 then return
                x=fload_disk_xyz('x')
                y=fload_disk_xyz('y')
                z=fload_disk_xyz('z')
                m=fload_disk_mass(1)
        endif


        ; bulge profile
        ; --------------
        if proftype eq "bulge" then begin
		if fload_npart(3) le 0 then return 
                x=fload_bulge_xyz('x')
                y=fload_bulge_xyz('y')
                z=fload_bulge_xyz('z')
                m=fload_bulge_mass(1)
        endif


        ; new stars profile
        ; -------------------
        if proftype eq "newstars" then begin
		if fload_npart(4) le 0 then return
                x=fload_newstars_xyz('x')
                y=fload_newstars_xyz('y')
                z=fload_newstars_xyz('z')
                m=fload_newstars_mass(1)
	endif



	; special cases
	; ---------------
	if proftype eq "special" then begin
		x= special_x
		y= special_y
		z= special_z
		m= special_m
	endif



	if keyword_set(h) then begin
		x= x/h
		y= y/h
		z= z/h
		m= m/h
	endif




	; compute (appropriate) profile
	; ------------------------------
	r_s = fltarr(bins)
	mass_sd_avg = fltarr(bins)
	mass_sd_1sig = fltarr(bins)

	if not keyword_set(frommanyproj) then begin
		a= sqrt(x*x + y*y)
                c= m

		if keyword_set(x_is_devac) then a=a^(0.25)
		if keyword_set(x_is_log) then a=alog10(a)

		process_prof_sd, a, c, bins, xmax, xmin, r_s, mass_sd_avg, $
					sd_1sig=mass_sd_1sig, $
					x_is_devac=x_is_devac, $
					x_is_log=x_is_log
	endif else begin
		process_prof_frommanyprojections, x, y, z, m, $
					xmin, xmax, bins, $
					x_is_devac=x_is_devac, $
					x_is_log=x_is_log, $
					mass_sd_avg, mass_sd_1sig, r_s
	endelse

;
; info
;
;print, "total mass= ", total(c)
;idx=where(a lt 5.5)
;print, "mass inside 5.5= ", total(c(idx)), "  (", total(c(idx))/total(c)," )"
;idx=where(a lt 25.0)
;print, "mass inside 25= ", total(c(idx)), "  (", total(c(idx))/total(c)," )"
;idx=where(a lt 50.0)
;print, "mass inside 50= ", total(c(idx)), "  (", total(c(idx))/total(c)," )"


	; take out zeros
	; ----------------
	idx= where(mass_sd_avg gt 0)
	if idx(0) ne -1 then begin
		mass_sd_avg= mass_sd_avg(idx)
		mass_sd_1sig= mass_sd_1sig(idx)
		r_s= r_s(idx)
	endif


	;
	; make log
	;    &
	; standard conversion to m_solar pc-2
	;-------------------------------------
	mass_sd_p1sig= alog10(mass_sd_avg+mass_sd_1sig) + 4.0
	mass_sd_m1sig= alog10(mass_sd_avg-mass_sd_1sig) + 4.0

	;mass_sd_p1sig= alog10(mass_sd_avg)+alog10(mass_sd_1sig) + 4.0
        ;mass_sd_m1sig= alog10(mass_sd_avg)-alog10(mass_sd_1sig) + 4.0
	mass_sd= alog10(mass_sd_avg) + 4.0

	help, mass_sd, mass_sd_p1sig, mass_sd_m1sig
	print, "mass_sd max/min= ", max(mass_sd), min(mass_sd)

	
	; now printing
	; -------------
	if not keyword_set(thislinest) then thislinest= 0
	if not keyword_set(thisthick) then thisthick= 4.0
	if not keyword_set(thispsym) then thispsym= 3
	;if linecolor eq 200 then begin thislinest= 2 & thisthick=1.0 & thispsym= 3 & endif


	; plot the total
	oplot, r_s, mass_sd, psym=-thispsym, linestyle=thislinest, color= linecolor, thick=thisthick

	if keyword_set(includeerr) then begin
		oplot, r_s, mass_sd_p1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
		oplot, r_s, mass_sd_m1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
	endif


	; try to fit de Vac. to it
	; -------------------------
	if keyword_set(fitdevac) then begin

		if keyword_set(x_is_devac) then r_tofit= r_s^(4.0)
		if keyword_set(x_is_log) then r_tofit= 10^r_s

		sd_tofit= 10^(mass_sd)
		weight= mass_sd_1sig
		;weight= 1.0e4*abs(mass_sd_1sig/mass_sd)

		fitandoverplot_devacprofile, r_tofit, sd_tofit, $
					weight, /ylogaxis, $
					x_is_devac=x_is_devac, $
					x_is_log=x_is_log

	endif



        ; try to fit Sersic Profile to it
        ; --------------------------------
        if keyword_set(fitsersic) then begin

		r_tofit= r_s
                if keyword_set(x_is_devac) then r_tofit= r_s^(4.0)
		if keyword_set(x_is_log) then r_tofit= 10^r_s

                sd_tofit= 10^(mass_sd)
                weight= mass_sd_1sig
                ;weight= 1.0e4*abs(mass_sd_1sig/mass_sd)

		;idx=where(r_tofit gt 1.0)
		;r_tofit= r_tofit(idx)
		;sd_tofit= sd_tofit(idx)
		;weight= weight(idx)
		minfitradius= 0.5

                fitandoverplot_sersicprofile, r_tofit, sd_tofit, weight, $
						/ylogaxis, $
						minfitradius=minfitradius, $
						x_is_devac=x_is_devac, $
                                        	x_is_log=x_is_log


        endif




end





;================================================================================
;================================================================================




