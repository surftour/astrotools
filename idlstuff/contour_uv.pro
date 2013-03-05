pro contour_uv, frun, snapnum, xlen, sendto, xz=xz, yz=yz, $
			filename=filename, thumbnail=thumbnail, fitstoo=fitstoo, $
			xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
			colorbar=colorbar, crude=crude, loadedsnap=loadedsnap, $
			pubstyle=pubstyle, old_tj_snap=old_tj_snap, $
			msg=msg, particlesonly=particlesonly, nolabels=nolabels, $
			showbhs=showbhs, use_calc_center=use_calc_center, $
			center=center, plotpts=plotpts, track_to_draw= track_to_draw

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
        ok= fload_snapshot(frun, snapnum)
      endif else begin
        ok= fload_snapshot_bh(frun, snapnum)
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

	xs= fload_newstars_xyz('x',center=[0,0,0])
	ys= fload_newstars_xyz('y',center=[0,0,0])
	zs= fload_newstars_xyz('z',center=[0,0,0])
	ms= fload_newstars_mass(1)

        ; what time is it?
        Ttime= float(fload_time(1))
        nsage= float(Ttime-fload_newstars_age(1))

	;hsml= fload_newstars_hsml(1)
	hsmls= 0.1 + 0.0*x

	idx=where(nsage lt 0.07)
	if idx(0) ne -1 then begin
	    N= long(n_elements(idx))
	    x= xs(idx)
	    y= ys(idx)
	    z= zs(idx)
	    m= ms(idx)
	    hsml= hsmls(idx)
	endif

    endif 


    ; sf gas
    ; ---------------------------------------
    sfr= fload_gas_sfr(1)
    xg= fload_gas_xyz('x',center=[0,0,0])
    yg= fload_gas_xyz('y',center=[0,0,0])
    zg= fload_gas_xyz('z',center=[0,0,0])
    mg= fload_gas_mass(1)
    hsmlg= fload_gas_hsml(1)

    idx= where(sfr gt 0)
    if idx(0) ne -1 then begin

	if N gt 0 then begin
	     N= long(N + n_elements(idx))
	     x= [x, xg(idx)]
	     y= [y, yg(idx)]
	     z= [z, zg(idx)]
	     m= [m, mg(idx)]
	     hsml= [hsml, hsmlg(idx)]
	endif else begin
	     N= long(n_elements(idx))
	     x= xg(idx)
	     y= yg(idx)
	     z= zg(idx)
	     m= mg(idx)
	     hsml= hsmlg(idx)
	endelse

    endif



; ----------------------------------------------
;  Call our generic contour plotting procedure
; ----------------------------------------------



contour_makeplot, x, y, z, m, hsml, float(xlen), sendto, xz=xz, yz=yz, $
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


