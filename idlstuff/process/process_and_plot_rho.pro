;================================================================================
;
;
;
pro process_and_plot_rho, frun, snapnum, proftype, $
			xmin, xmax, bins, $
			x_is_log=x_is_log, $
			linecolor=linecolor, $
			includeerr=includeerr, $
			fitnfw=fitnfw, $
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

	if not keyword_set(proftype) then proftype= "dm"

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


        ; halo profile
        ; --------------
        if proftype eq "halo" or proftype eq "dm" then begin
                if fload_npart(1) le 0 then return
                x=fload_halo_xyz('x')
                y=fload_halo_xyz('y')
                z=fload_halo_xyz('z')
                m=fload_halo_mass(1)
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
	rho= fltarr(bins)
	rho_1sig = fltarr(bins)

	radius= sqrt(x*x + y*y + z*z)

	if keyword_set(x_is_log) then radius= alog10(radius)

	process_prof_rho, radius, m, bins, xmax, xmin, r_s, rho, rho_1sig, $
					x_is_log=x_is_log


	; take out zeros
	; ----------------
	idx= where(rho gt 0)
	if idx(0) ne -1 then begin
		rho= rho(idx)
		rho_1sig= rho_1sig(idx)
		r_s= r_s(idx)
	endif


	;
	; make log
	;    &
	; standard conversion to m_solar pc-3
	;-------------------------------------
	rho_p1sig= alog10(rho+rho_1sig) + 1.0
	rho_m1sig= alog10(rho-rho_1sig) + 1.0

	;rho_p1sig= alog10(rho)+alog10(rho_1sig) + 1.0
        ;rho_m1sig= alog10(rho)-alog10(rho_1sig) + 1.0
	rho= alog10(rho) + 1.0

	help, rho, rho_p1sig, rho_m1sig
	print, "rho max/min= ", max(rho), min(rho)

	
	; now printing
	; -------------
	if not keyword_set(thislinest) then thislinest= 0
	if not keyword_set(thisthick) then thisthick= 4.0
	if not keyword_set(thispsym) then thispsym= 3
	;if linecolor eq 200 then begin thislinest= 2 & thisthick=1.0 & thispsym= 3 & endif


	; plot the total
	oplot, r_s, rho, psym=-thispsym, linestyle=thislinest, color= linecolor, thick=thisthick

	if keyword_set(includeerr) then begin
		oplot, r_s, rho_p1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
		oplot, r_s, rho_m1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
	endif


	; try to fit NFW profile
	; -------------------------
	if keyword_set(fitnfw) then begin

		if keyword_set(x_is_log) then r_tofit= 10^r_s

		sd_tofit= 10^(rho)
		weight= rho_1sig
		;weight= 1.0e4*abs(rho_1sig/rho)

		fitandoverplot_nfw, r_tofit, sd_tofit, $
					weight, /ylogaxis, $
					xmax=xmax, xmin=xmin, $
					x_is_log=x_is_log

	endif




end





;================================================================================
;================================================================================




