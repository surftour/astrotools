pro jy_contour_uv, frun, snapnum, xlen, sendto, xz=xz, yz=yz, $
			filename=filename, thumbnail=thumbnail, fitstoo=fitstoo, $
			xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
			colorbar=colorbar, crude=crude, loadedsnap=loadedsnap, $
			pubstyle=pubstyle, old_tj_snap=old_tj_snap, $
			msg=msg, particlesonly=particlesonly, nolabels=nolabels, $
			showbhs=showbhs, use_calc_center=use_calc_center, $
			center=center, plotpts=plotpts, track_to_draw= track_to_draw,$
                        rotate_phi=rotate_phi,rotate_theta=rotate_theta

if not keyword_set(xlen) then xlen=100.0
if not keyword_set(sendto) then sendto='x'
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(frun) then begin
   print, "  "
   print, "contour_uv, frun, snapnum, xlen, sendto, /xz, /yz,"
   print, "                  filename=filename, /thumbnail, /fitstoo,/colorbar,"
   print, "                  xthickness=xthickness, ythickness=ythickness, zthickness=zthickness,"
   print, "                  "
   print, "  "
   print, "  Pixel Map of newstars component in xy projection."
   print, "  Can do contour, but need to turn it on within!"
   print, "  "
   print, "  "
   print, "  "
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


if not keyword_set(msg) then msg=' '



; --------------------------------
;  Load variables for smoothing
; --------------------------------
if not keyword_set(loadedsnap) then begin
      if keyword_set(old_tj_snap) then begin
        ok= fload_snapshot(frun, snapnum,/nopot)
      endif else begin
        ok= fload_snapshot_bh(frun, snapnum,/nopot)
      endelse
      if ok eq 1 then begin
        print, "PROBLEM: opening file"
        return
    endif
endif



    if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)


    ; new stellar particles
    ; -------------------------------------------
    N= 0L
    if fload_npart(4) gt 0 then begin

	xs= fload_allstars_xyz('x',center=[0,0,0])
	ys= fload_allstars_xyz('y',center=[0,0,0])
	zs= fload_allstars_xyz('z',center=[0,0,0])
	ms= fload_allstars_mass(1)
        zs= fload_allstars_z(1)

        ; what time is it?
        Ttime= float(fload_time(1))
        nsage= float(Ttime-fload_allstars_age(1))

	hsmls= fload_allstars_hsml(1)
	;hsmls= 0.1 + 0.0*xs

    endif 

    x = xs
    y = ys
    z = zs
    hsml = hsmls

    uvlum = jy_get_uv('/home/tcox/Tools/idlstuff/jy_uvtable_bc03_salp.fits',nsage,abs(zs),ms)
    ms = uvlum
    set_maxden = 1.0e8
    set_dynrng=1.0e5

; ----------------------------------------------
;  Call our generic contour plotting procedure
; ----------------------------------------------

;
; this is the same as mine only he turned off the overplotting 
; of young stars
;
;jy_contour_makeplot, x, y, z, ms, hsml, float(xlen), sendto, xz=xz, yz=yz, $
;

contour_makeplot, x, y, z, ms, hsml, float(xlen), sendto, xz=xz, yz=yz, $
                        filename=filename, thumbnail=thumbnail, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
                        pixels=pixels, zthickness=zthickness, $
                        crude=crude, center=center, msg=msg, plotpts=plotpts, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        nolabels=nolabels, pubstyle=pubstyle, particlesonly=particlesonly, $
                        showbhs=showbhs, set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        track_to_draw=track_to_draw



; -------------
;  Done
; -------------





end


