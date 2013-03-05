;================================================================================
;
;
;
pro process_and_plot_3dvel, frun, snapnum, plotdetail, $
			xmin, xmax, bins, $
			x_is_log=x_is_log, $
			linecolor=linecolor, $
			includeerr=includeerr, $
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

	galcomp= strmid(plotdetail,0,strpos(plotdetail,'/'))
	print, "galcomp= ", galcomp
	printit= strmid(plotdetail,strpos(plotdetail,'/')+1)
	print, "printit= ", printit


	; -----------------------------------------------

	if not keyword_set(galcomp) then galcomp= "dm"

        ; total (all star) profile
        ; --------------------------
	if galcomp eq "allstars" then begin
        	x=fload_allstars_xyz('x')
        	y=fload_allstars_xyz('y')
        	z=fload_allstars_xyz('z')
        	m=fload_allstars_mass(1)
        	vx=fload_allstars_v('x')
        	vy=fload_allstars_v('y')
        	vz=fload_allstars_v('z')
        	vr=fload_allstars_v('r')
        	vtan=fload_allstars_v('tan')
        	vtheta=fload_allstars_v('theta')
	endif

        ; gas profile
        ; --------------
        if galcomp eq "gas" then begin
		if fload_npart(0) le 0 then return
		x=fload_gas_xyz('x')
		y=fload_gas_xyz('y')
		z=fload_gas_xyz('z')
		m=fload_gas_mass(1)
        endif


        ; halo profile
        ; --------------
        if galcomp eq "halo" or galcomp eq "dm" then begin
                if fload_npart(1) le 0 then return
                x=fload_halo_xyz('x')
                y=fload_halo_xyz('y')
                z=fload_halo_xyz('z')
                m=fload_halo_mass(1)
        	vx=fload_halo_v('x')
        	vy=fload_halo_v('y')
        	vz=fload_halo_v('z')
        	vr=fload_halo_v('r')
        	vtan=fload_halo_v('tan')
        	vtheta=fload_halo_v('theta')
        endif


        ; disk profile
        ; --------------
        if galcomp eq "disk" then begin
		if fload_npart(2) le 0 then return
                x=fload_disk_xyz('x')
                y=fload_disk_xyz('y')
                z=fload_disk_xyz('z')
                m=fload_disk_mass(1)
                vx=fload_disk_v('x')
                vy=fload_disk_v('y')
                vz=fload_disk_v('z')
                vr=fload_disk_v('r')
                vtan=fload_disk_v('tan')
                vtheta=fload_disk_v('theta')
        endif


        ; bulge profile
        ; --------------
        if galcomp eq "bulge" then begin
		if fload_npart(3) le 0 then return 
                x=fload_bulge_xyz('x')
                y=fload_bulge_xyz('y')
                z=fload_bulge_xyz('z')
                m=fload_bulge_mass(1)
                vx=fload_bulge_v('x')
                vy=fload_bulge_v('y')
                vz=fload_bulge_v('z')
                vr=fload_bulge_v('r')
                vtan=fload_bulge_v('tan')
                vtheta=fload_bulge_v('theta')
        endif


        ; new stars profile
        ; -------------------
        if galcomp eq "newstars" then begin
		if fload_npart(4) le 0 then return
                x=fload_newstars_xyz('x')
                y=fload_newstars_xyz('y')
                z=fload_newstars_xyz('z')
                m=fload_newstars_mass(1)
	endif



	; special cases
	; ---------------
	if galcomp eq "special" then begin
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
	avg_s= fltarr(bins)
	avg_1sig = fltarr(bins)
	avg_disp= fltarr(bins)

	radius= sqrt(x*x + y*y + z*z)

	if keyword_set(x_is_log) then radius= alog10(radius)

	; what quantity are we averaging?
	avgq= m
	if printit eq "sigr" then avgq= vr
	;if printit eq "sigtan" then avgq= vtan
	if printit eq "sigtheta" then avgq= vtheta
	if printit eq "beta" then avgq= vr

	;
	; determine the radia average of avgq
	;
	process_prof_avg, radius, avgq, bins, xmax, xmin, r_s, avg_s, avg_1sig, avg_disp=avg_disp ;, /xaxis_min_sn


	if printit eq "beta" then begin
		tempa= avg_s & tempb= avg_1sig & tempc= avg_disp
		;avgq= vtan
		avgq= vtheta
		process_prof_avg, radius, avgq, bins, xmax, xmin, r_s, avg_s, avg_1sig, avg_disp=avg_disp ;, /xaxis_min_sn

		idx= where(tempc gt 0.0)
		if idx(0) ne -1 then begin
			r_s= r_s(idx)
			avg_s= 1.0 - ((avg_disp(idx) * avg_disp(idx)) / (tempc(idx) * tempc(idx)))
		endif
	endif



	; take out zeros (don't necessarily want to do this with averages)
	; ----------------
	;idx= where(avg_s gt 0)
	;if idx(0) ne -1 then begin
	;	avg_s= avg_s(idx)
	;	avg_1sig= avg_1sig(idx)
	;	avg_disp= avg_disp(idx)
	;	r_s= r_s(idx)
	;endif



	; determine +/- 1 sigma
	;-------------------------------------
	avg_p1sig= avg_s+avg_1sig
	avg_m1sig= avg_s-avg_1sig

	help, avg_s, avg_p1sig, avg_m1sig
	print, "avg_s    max/min= ", max(avg_s), min(avg_s)
	print, "avg_1sig max/min= ", max(avg_1sig), min(avg_1sig)
	print, "avg_disp max/min= ", max(avg_disp), min(avg_disp)



	

	; now printing
	; -------------

	if not keyword_set(thislinest) then thislinest= 0
	if not keyword_set(thisthick) then thisthick= 4.0
	if not keyword_set(thispsym) then thispsym= 3
	;if linecolor eq 200 then begin thislinest= 2 & thisthick=1.0 & thispsym= 3 & endif


        if printit eq "sigr" then begin 
		; overplot the expected sigma_r for a H90 profile
		rs= 10^(r_s)
		rplusa= rs + 69.03
		rdiva= rs/69.03
		sig2r= 43007.1 * 1088.46 / 12. / 69.03 * $
			((12.*rs * rplusa^3.)/((69.03)^4.) * alog(rplusa/rs) - $
				(rs/rplusa) * (25. + 52.*rdiva + 42.*rdiva*rdiva + 12.*rdiva*rdiva*rdiva))
		oplot, r_s, sqrt(sig2r), psym=-3, linestyle= 2, color= 0, thick=8.0

		; show the mean, in addition to the dispersion (which is shown below)
		;oplot, r_s, avg_s, psym=-3, linestyle=1, color= linecolor, thick=(thisthick/2.0)
                avg_s= avg_disp
        endif

        if printit eq "sigtheta" then begin
		; show the mean, in addition to the dispersion (which is shown below)
		;oplot, r_s, avg_s, psym=-3, linestyle=1, color= linecolor, thick=(thisthick/2.0)
                avg_s= avg_disp
        endif



	; plot the total
	oplot, r_s, avg_s, psym=-thispsym, linestyle=thislinest, color= linecolor, thick=thisthick

	if keyword_set(includeerr) then begin
		oplot, r_s, avg_p1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
		oplot, r_s, avg_m1sig, psym=-3, linestyle=1, color= linecolor, thick=2.0
	endif




end





;================================================================================
;================================================================================




