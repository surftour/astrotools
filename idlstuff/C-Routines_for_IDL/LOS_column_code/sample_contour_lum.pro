pro contour_lum, Band, frun, snapnum, xlen, sendto, $
                        center=center, $
						use_calc_center=use_calc_center, $
                        colorbar=colorbar, $
                        crude=crude, $
                        filename=filename, $
                        fitstoo=fitstoo, $
                        loadedsnap=loadedsnap, $
                        msg=msg, $
						nodust=nodust, $
                        nolabels=nolabels, $
                        particlesonly=particlesonly, $
                        pixels=pixels, $
                        pubstyle=pubstyle, $
                        plotpts=plotpts, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        startid=startid, numpart=numpart, $
                        thumbnail=thumbnail, $
                        track_to_draw= track_to_draw, $
                        xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
						showbhs=showbhs, $
                        xz=xz, yz=yz

if not keyword_set(xlen) then xlen=100.0
if not keyword_set(sendto) then sendto='x'
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(Band) then begin
   print, "  "
   print, "contour_lum, Band, frun, snapnum, xlen, sendto, /xz, /yz, /colorbar, "
   print, "              filename=filename, /thumbnail, /fitstoo, /crude, /loadedsnap, "
   print, "              xthickness=xthickness, ythickness=ythickness, zthickness=zthickness"
   print, "              /nolabels, center=center, /plotpts/, track_to_draw= track_to_draw"
   print, "              /use_calc_center"
   print, "  "
   print, " -> see contour_gas for information on available options"
   print, "  "
   print, "  "
   print, "  "
   return
endif



; --------------------------------
;  Center this on something other
;  than 0,0,0
; --------------------------------
if not keyword_set(center) then center=[0.0,0.0,0.0]


set_maxden=10.0
set_dynrng=1.0e+5

; --------------------------------
;  Load variables for smoothing
; --------------------------------
    if not keyword_set(loadedsnap) then begin
      ;if (fload_snapshot(frun, snapnum)) then begin
      if (fload_snapshot_bh(frun, snapnum)) then begin
        print, "PROBLEM: opening file"
        return
      endif
    endif

    if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)

    ndisk= fload_npart(2)				; disk particles
    nbulge= fload_npart(3)				; bulge particles
    nstars= fload_npart(4)
    npart= long(ndisk) + long(nbulge) + long(nstars)
    print, "Ntot,Ndisk,Nbulge,Nstars= ",npart,ndisk,nbulge,nstars

    x= fload_allstars_xyz('x',center=[0,0,0])
    y= fload_allstars_xyz('y',center=[0,0,0])
    z= fload_allstars_xyz('z',center=[0,0,0])
    ;m= fload_allstars_mass(1)

    ; if we ever want to seperate out just one galaxy
    ;if keyword_set(startid) and keyword_set(numpart) then begin
	;x= fload_1gal_allstars_xyz('x',startid,numpart,center=[0,0,0])
	;y= fload_1gal_allstars_xyz('y',startid,numpart,center=[0,0,0])
	;z= fload_1gal_allstars_xyz('z',startid,numpart,center=[0,0,0])
	;m= fload_1gal_allstars_xyz(startid,numpart)
    ;endif

    ;hsml= 0.0*x + 0.1    ; assume a standard smoothing size
    hsml= 0.5*fload_allstars_hsml(1)

    TTime= float(fload_time(1))

    N= npart
    m= 1.0e+10*fload_allstars_mass(1)
    age=fload_allstars_age(1)
    age=float(TTime-age)
    zmets=fload_allstars_z(1)


; ---------------------------------
;  Add BH Luminosity
; ---------------------------------
include_BH_lum= 1
;include_BH_lum= 0
if include_BH_lum eq 1 and fload_npart(5) gt 0 then begin
	; get BH lum
	; -------------
	Nbh= fload_npart(5)
	bhids= fload_blackhole_id(1)
	bhbololum= fload_blackhole_lum(frun,Ttime,/bolometric)
	; Add BH to other lists
	; -----------------------
	x= [fload_blackhole_xyz('x',center=[0,0,0]), x]
	y= [fload_blackhole_xyz('y',center=[0,0,0]), y]
	z= [fload_blackhole_xyz('z',center=[0,0,0]), z]

	for i=1,fload_npart(5) do hsml= [0.1, hsml]

	set_bh_to_zero= 0

endif else begin
	set_bh_to_zero= 1
endelse


ProjectThisLum=fltarr(Nstars+Nbh)+1.0e+6


; ---------------------------------
;  Attenuate the Light
; ---------------------------------

;attenuate_light= 0
;attenuate_light= 1
;if keyword_set(attenuate_light) then begin
if not keyword_set(nodust) then begin


        Ngas= fload_npart(0)
        Nstars= fload_npart(2)+fload_npart(3)+fload_npart(4)
        Nbh= fload_npart(5)
	if set_bh_to_zero eq 1 then Nbh= 0

	; need it in radians
        if keyword_set(rotate_theta) then theta = (!PI / 180.0) * rotate_theta else theta= 0.0
        if keyword_set(rotate_phi) then phi = (!PI / 180.0) * rotate_phi else phi= 0.0


        ; fields to pass
        Coord = fltarr(10,Ngas+Nstars+Nbh)

        ; gas
        Coord(0,0:Ngas-1) = fload_gas_xyz('x', center=[0,0,0])
        Coord(1,0:Ngas-1) = fload_gas_xyz('y', center=[0,0,0])
        Coord(2,0:Ngas-1) = fload_gas_xyz('z', center=[0,0,0])

        Coord(3,0:Ngas-1) = fload_gas_u(1)
        Coord(4,0:Ngas-1) = fload_gas_rho(1)
        Coord(5,0:Ngas-1) = fload_gas_hsml(1)
        Coord(6,0:Ngas-1) = fload_gas_numh(1)
        Coord(7,0:Ngas-1) = fload_gas_nume(1)
        Coord(8,0:Ngas-1) = fload_gas_metallicity(1)
        Coord(9,0:Ngas-1) = fload_gas_mass(1)
        ;;print, fload_gas_mass(1)

        ; stars
        Coord(0,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('x', center=[0,0,0])
        Coord(1,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('y', center=[0,0,0])
        Coord(2,Ngas:Nstars+Ngas-1) = fload_allstars_xyz('z', center=[0,0,0])

        ; black holes
		if Nbh gt 0 then begin
        	Coord(0,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('x', center=[0,0,0])
        	Coord(1,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('y', center=[0,0,0])
        	Coord(2,Ngas+Nstars:Nstars+Ngas+Nbh-1) = fload_blackhole_xyz('z', center=[0,0,0])
        endif
      

        los_NH= fltarr(Nstars+Nbh)
        los_Z= fltarr(Nstars+Nbh)

        print, "PASSING: "
        print, "N_gas= ", Ngas
        print, "N_stars= ", Nstars
        print, "N_bh= ", Nbh
        print, "theta= ", theta
        print, "phi= ", phi
        help, Coord

	S = CALL_EXTERNAL('./getnh.so', $
                'getnh', $
                Ngas, $ 
                Nstars, $
                Nbh, $ 
                theta, $
                phi, $ 
                Coord, $
                los_NH, $
                los_Z)

	; trap for really low NH values
	;   and zero metallicity (make it small instead)
	idx= where(los_NH lt 1.0e+10)
	if idx(0) ne -1 then los_NH(idx)= 1.0e+10
	idx= where(los_Z le 0.0)
	if idx(0) ne -1 then los_Z(idx)= 1.0e-5


	;  OK, now we've got the N_h values,
	; what is the attenuated luminosity

	print, "Now calling attenuation program."

	; pass in anstroms, but it gets converted
	; to microns in ism_absorption
	;Band_Lambda= float(Band_Lambda)
	;print, "Band_Lambda= ",Band_Lambda

        AttLums= fltarr(Nstars+Nbh)      ;attenuated luminosities
	PreLums= float(ProjectThisLum)

	AttLums = PreLums * EXP(-cross_section(6.818e14)*los_NH)

print, "total light= ", total(ProjectThisLum), total(PreLums)
print, "attenuated light= ", total(AttLums)
print, "    abs. f= ", total(AttLums)/total(ProjectThisLum)

print, " "
print, " Stars "
print, " ------"
print, "att. light = ", total(AttLums[Nbh:Nbh+Nstars-1])
print, "     abs. f= ", total(AttLums[Nbh:Nbh+Nstars-1])/total(ProjectThisLum[Nbh:Nbh+Nstars-1])

if Nbh gt 0 then begin
	print, " "
	print, " BH  "
	print, " ------"
	print, "att. light = ", total(AttLums[0:Nbh-1])
	print, "     abs. f= ", total(AttLums[0:Nbh-1])/total(ProjectThisLum[0:Nbh-1])
	print, "     los_NH= ", los_NH[0:Nbh-1]
	print, "     los_Z =  ", los_Z[0:Nbh-1]
	print, " "
endif else begin
	print, " NO Black Hole !!!"
endelse


; actually project the Attenuated Light
; ----------------------------------------
ProjectThisLum= AttLums



endif



; ----------------------------------------------
;  Call our generic contour plotting procedure
; ----------------------------------------------

;set_maxden=10.0
;set_maxden=1.0e+10
set_maxden=1.0e+9
set_dynrng=1.0e+5


contour_makeplot, x, y, z, ProjectThisLum, hsml, xlen, sendto, xz=xz, yz=yz, $
                        filename=filename, thumbnail=thumbnail, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
                        pixels=pixels, zthickness=zthickness, $
                        crude=crude, center=center, msg=msg, $
			plotpts=plotpts, particlesonly=particlesonly, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			nolabels=nolabels, pubstyle=pubstyle, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			track_to_draw=track_to_draw, showbhs=showbhs





; -------------
;  Done
; -------------





end


