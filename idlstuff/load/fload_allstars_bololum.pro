function fload_allstars_bololum, junk


        ; grab luminosities
        ; ----------------------
        ndisk= fload_npart(2)                           ; disk particles
        nbulge= fload_npart(3)                          ; bulge particles
        nstars= fload_npart(4)
        npart= long(ndisk) + long(nbulge) + long(nstars)
        print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars
        TTime= float(fload_time(1))
        N= npart
        m= 1.0e+10*fload_allstars_mass(1)
        age=fload_allstars_age(1)
        age=float(TTime-age)
        zmets=fload_allstars_z(1)


        ; get the luminosities
        ;  - in units of solar luminosities
        print, "load luminosities"
	print, "N= ",N

	if(N gt 0) then begin

	        NewStarLuminosities= fltarr(14,N)      ;luminosities
	        MassInfo= fltarr(2,N)      ;luminosities

	
	
	;  new version of colors code (see below)
	;------------------------------
	;model= 0   ; salpeter
	model= 1   ; chabrier
	silent= 0
	;silent= 1
	spawn, 'echo $TJHOME', result
	homedir= strcompress(result,/remove_all)
	;libfile= homedir+'/Tools/C-Routines_for_IDL/ComputeColors/colors'   ; older version
	libfile= homedir+'/Tools/C-Routines_for_IDL/colors/colors'
	S = CALL_EXTERNAL(libfile, $
	        'main', $
	        long(N), $
	        age, $           ; age lin gyr
	        zmets, $             ; metallicity in solar units
	        NewStarLuminosities, $ ; luminosities in abs magnitudes
	        MassInfo, $
	        long(model), $         ; model==0 -> Salpeter, model==1 Chabrier
	        long(silent), $        ; produce output - or not
	        /F_VALUE)

		BolMag= NewStarLuminosities(0,*)
		usdssMag= NewStarLuminosities(1,*)
		gsdssMag= NewStarLuminosities(2,*)
		rsdssMag= NewStarLuminosities(3,*)
		isdssMag= NewStarLuminosities(4,*)
		zsdssMag= NewStarLuminosities(5,*)
		UMag= NewStarLuminosities(6,*)
		BMag= NewStarLuminosities(7,*)
		VMag= NewStarLuminosities(8,*)
		RMag= NewStarLuminosities(9,*)
		IMag= NewStarLuminosities(10,*)
		JMag= NewStarLuminosities(11,*)
		HMag= NewStarLuminosities(12,*)
		KMag= NewStarLuminosities(13,*)

	endif else begin
        	print,'No stars, no luminosities (or magnitudes, in this case)!'
        	return, -1
	endelse


	; -----------------------
	;  Convert to Luminosity
	; -----------------------
	if keyword_set(notsolar) then begin
        	Lum_Bol= (10^(-0.4*BolMag)) * m
        	Lum_U= (10^(-0.4*UMag)) * m
        	Lum_B= (10^(-0.4*BMag)) * m
        	Lum_V= (10^(-0.4*VMag)) * m
        	Lum_R= (10^(-0.4*RMag)) * m
        	Lum_I= (10^(-0.4*IMag)) * m
        	Lum_J= (10^(-0.4*JMag)) * m
        	Lum_H= (10^(-0.4*HMag)) * m
        	Lum_K= (10^(-0.4*KMag)) * m
        	Lum_sdss_u= (10^(-0.4*usdssMag)) * m
        	Lum_sdss_g= (10^(-0.4*gsdssMag)) * m
        	Lum_sdss_r= (10^(-0.4*rsdssMag)) * m
        	Lum_sdss_i= (10^(-0.4*isdssMag)) * m
        	Lum_sdss_z= (10^(-0.4*zsdssMag)) * m
	endif else begin
        	;Lum_Bol= (10^(-0.4*(BolMag-4.75))) * m
        	;Lum_U= (10^(-0.4*(UMag-5.60))) * m
        	;Lum_B= (10^(-0.4*(BMag-5.51))) * m
        	;Lum_V= (10^(-0.4*(VMag-4.84))) * m
        	;Lum_R= (10^(-0.4*(RMag-4.48))) * m
        	;Lum_I= (10^(-0.4*(IMag-4.13))) * m
        	;Lum_J= (10^(-0.4*(JMag-3.70))) * m
        	;Lum_H= (10^(-0.4*(HMag-3.37))) * m
        	;Lum_K= (10^(-0.4*(KMag-3.33))) * m
        	;Lum_sdss_u= (10^(-0.4*(usdssMag-6.28))) * m
        	;Lum_sdss_g= (10^(-0.4*(gsdssMag-4.95))) * m
        	;Lum_sdss_r= (10^(-0.4*(rsdssMag-4.45))) * m
        	;Lum_sdss_i= (10^(-0.4*(isdssMag-4.35))) * m
        	;Lum_sdss_z= (10^(-0.4*(zsdssMag-4.36))) * m

        	; brant's new values - see below
        	Lum_Bol= (10^(-0.4*(BolMag-4.74))) * m
        	Lum_U= (10^(-0.4*(UMag-5.56))) * m
        	Lum_B= (10^(-0.4*(BMag-5.45))) * m
        	Lum_V= (10^(-0.4*(VMag-4.80))) * m
        	Lum_R= (10^(-0.4*(RMag-4.46))) * m
        	Lum_I= (10^(-0.4*(IMag-4.10))) * m
        	Lum_J= (10^(-0.4*(JMag-3.66))) * m
        	Lum_H= (10^(-0.4*(HMag-3.32))) * m
        	Lum_K= (10^(-0.4*(KMag-3.28))) * m
        	Lum_sdss_u= (10^(-0.4*(usdssMag-6.75))) * m
        	Lum_sdss_g= (10^(-0.4*(gsdssMag-5.33))) * m
        	Lum_sdss_r= (10^(-0.4*(rsdssMag-4.67))) * m
        	Lum_sdss_i= (10^(-0.4*(isdssMag-4.48))) * m
        	Lum_sdss_z= (10^(-0.4*(zsdssMag-4.42))) * m
	endelse


	;======================

        ; trap for any NaN's
        ;idx=where(finite(Lum_B) eq 0)
        idx=where(finite(Lum_K) eq 0)
        if idx(0) ne -1 then begin
                Lum_Bol(idx)= 100.0 & Lum_U(idx)= 100.0 & Lum_B(idx)= 100.0 & Lum_V(idx)= 100.0 & Lum_R(idx)= 100.0
                Lum_I(idx)= 100.0 & Lum_J(idx)= 100.0 & Lum_H(idx)= 100.0 & Lum_K(idx)= 100.0
                Lum_sdss_u(idx)= 100.0 & Lum_sdss_g(idx)= 100.0 & Lum_sdss_r(idx)= 100.0
                Lum_sdss_i(idx)= 100.0 & Lum_sdss_z(idx)= 100.0
        print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
        endif

        print, "Bolometric luminosity= ", total(Lum_Bol)


	return, Lum_Bol

end

