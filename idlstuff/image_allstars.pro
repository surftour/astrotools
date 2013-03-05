pro image_allstars, frun, snapnum, xlen, $
                        center=center, $
			use_calc_center=use_calc_center, $
                        crude=crude, $
                        filename=filename, $
                        loadedsnap=loadedsnap, $
                        msg=msg, $
			old_tj_snap=old_tj_snap, $
                        pixels=pixels, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        startid=startid, numpart=numpart, $
                        xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
                        xz=xz, yz=yz

if not keyword_set(xlen) then xlen=100.0
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(frun) then begin
   print, "  "
   print, "  image_allstars:  check the file  "
   print, "  "
   print, "  "
   return
endif



; --------------------------------
;  Center this on something other
;  than 0,0,0
; --------------------------------
if not keyword_set(center) then center=[0.0,0.0,0.0]

if not keyword_set(msg) then msg='all stars'

set_maxden=5.0
;set_dynrng=5.0e+3
set_dynrng=5.0e+7

; --------------------------------
;  Load variables for smoothing
; --------------------------------
    if not keyword_set(loadedsnap) then begin
      if keyword_set(old_tj_snap) then begin
        ok= fload_snapshot(frun, snapnum)
      endif else begin
        ;ok= fload_snapshot_bh(frun, snapnum)
        ok= fload_snapshot_bh(frun, snapnum, /nopot_in_snap)
      endelse
      if ok eq 1 then begin
        print, "PROBLEM: opening file:", frun, snapnum
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
    m= fload_allstars_mass(1)

    if keyword_set(startid) and keyword_set(numpart) then begin
	x= fload_1gal_allstars_xyz('x',startid,numpart,center=[0,0,0])
	y= fload_1gal_allstars_xyz('y',startid,numpart,center=[0,0,0])
	z= fload_1gal_allstars_xyz('z',startid,numpart,center=[0,0,0])
	m= fload_1gal_allstars_mass(1,startid,numpart)
    endif

    ; version 1
    ;hsml= 0.1+0.0*m
    ;print, "WARNING: multiplying the stellar hsml by 2.0"

    ; version 2
    ;r=fload_allstars_xyz('r')
    ;hsml= (r*1.0)
    ;idx= where(hsml lt 0.4)
    ;if idx(0) ne -1 then hsml(idx)= 0.4

    ; version 3
    if not keyword_set(particlesonly) then hsml= fload_allstars_hsml(1) else hsml= 0.2+0.0*m



; ----------------------------------------------
;  Make Image
; ----------------------------------------------


img_makeimage, x, y, z, m, hsml, xlen, 'jpg', xz=xz, yz=yz, $
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


