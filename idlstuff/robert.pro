pro robert, frun, snapnum, xlen, sendto, $
			center=center, $
			colorbar=colorbar, $
			crude=crude, $
			filename=filename, $
			fitstoo=fitstoo, $
			loadedsnap=loadedsnap, $
			msg=msg, $
			nolabels=nolabels, $
			old_tj_snap=old_tj_snap, $
			particlesonly=particlesonly, $
			pixels=pixels, $
			pubstyle=pubstyle, $
			plotpts=plotpts, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			showbhs=showbhs, $
			startid=startid, numpart=numpart, $
			thumbnail=thumbnail, $
			track_to_draw= track_to_draw, $
			use_calc_center=use_calc_center, $
			xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
			xz=xz, yz=yz

if not keyword_set(xlen) then xlen=100.0
if not keyword_set(sendto) then sendto='ps'
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(filename) then filename='contour.eps'
if not keyword_set(frun) then begin
   print, "  "
   print, "  "
   print, "  "
   print, "  This used to be contour_gas.pro "
   print, "  "
   print, "  "
   print, "  "
   print, "-------------------------------------------"
   print, "   Pixel Map of Density in Projection"
   print, "-------------------------------------------"
   print, "robert, frun, snapnum, xlen, sendto, /xz, /yz,"
   print, "              filename=filename, /thumbnail, /fitstoo, /colorbar, "
   print, "              pixels=pixels, zthickness=zthickness,"
   print, "              xthickness=xthickness, ythickness=ythickness,"
   print, "              /loadedsnap, /crude, center=center"
   print, "              /nolabels, /showbhs"
   print, "  "
   print, " keyword [default]       - description "
   print, " crude [off]             - make a crude pixel map rather than c-program one"
   print, " center [0,0,0]          - zoom in on center"
   print, " colorbar [off]          - prints color bar and scale on right side"
   print, " contour [off]           - Can do, but need to turn it on within"
   print, " filename [contour.eps]  - same file as"
   print, " fitstoo [off]           - write a fits file"
   print, " loadedsnap [off]        - snap is already load it, routine does not"
   print, " msg [gas]               - print a message on figure"
   print, " nolabels [off]          - suppress labels, so that the figure fills the plotted area"
   print, " numpart [all]           - with startid, can select to print certain particles"
   print, " old_tj_snap [off]       - data is in the old 'TJ' snapshot format"
   print, " particlesonly [off]     - plot the particles rather than the density map"
   print, " plotpts [off]           - plot the actual data points, can be prohibitive if there are many"
   print, " pixels [480]            - "
   print, " pubstyle [off]          - omit run name and msg printing"
   print, " rotate_phi [0]          - rotate the image by phi, about z-axis"
   print, " rotate_theta [0]        - rotate the image by theta, about y-axis"
   print, " sendto [ps]             - ps, x, jpg (not working), gif, z(png)"
   print, " showbhs [off]           - put an x where the black hole is"
   print, " snapnum [100]           - snapshot number, reads snapshot_<snapnum>"
   print, " startid [all]           - with numpart, can select to print certain particles"
   print, " thumbnail [off]         - creates small image, only works for sendto='gif'"
   print, " track_to_draw [off]     - will add the path as found in centers.txt"
   print, " use_calc_center [off]   - uses the center automatically calculated when a snap is opened"
   print, " xlen [100.0]            - length of each side is 2x(xlen)"
   print, " x|y|zthickness [off]    - takes slice of this thickness centered at center"
   print, " xz [off]                - project along xz, default is xy"
   print, " yz [off]                - project along yz, default is xy"
   print, "  "
   print, "  "
   return
endif





; --------------------------------
;  Center this on something other
;  than 0,0,0
; --------------------------------
if not keyword_set(center) then center=[0.0,0.0,0.0]

if not keyword_set(msg) then msg='gas'

; standard
; ---------
set_maxden= 75.0
set_dynrng= 3.0
;set_maxden= 3.0e+3
;set_dynrng= 8.0

; zoom into center
; -----------------
;set_maxden= 1.0
;set_dynrng= 1.0e+4



; --------------------------------
;  Load variables for smoothing
; --------------------------------
    if not keyword_set(loadedsnap) then begin
      if keyword_set(old_tj_snap) then begin
	ok= fload_snapshot(frun, snapnum)
      endif else begin
	;ok= fload_snapshot_bh(frun, snapnum)
	;ok= fload_snapshot_bh(frun, snapnum,/nopot_in_snap)
	;ok= fload_snapshot_bh(frun, snapnum,/ics)
	ok= robert_load(frun, snapnum,/ics)
      endelse
      if ok eq 1 then begin
        print, "PROBLEM: opening file"
        return
      endif
    endif
    npart= fload_npart(0)                        ; return gas particles
    N= long(npart)



    ; default - Fe
    x= fload_gas_xyz('x',center=[0,0,0])
    y= fload_gas_xyz('y',center=[0,0,0])
    z= fload_gas_xyz('z',center=[0,0,0])
    m= fload_gas_mass(1)
    hsml= fload_gas_hsml(1)

    ; default - Si
    ;x= fload_halo_xyz('x',center=[0,0,0])
    ;y= fload_halo_xyz('y',center=[0,0,0])
    ;z= fload_halo_xyz('z',center=[0,0,0])
    ;m= fload_gas_mass(1)
    ;hsml= fload_gas_hsml(1)

    m= x*0.0 + 1.
    hsml= x*0.0 + 0.05

    x= x/1.0e4
    y= y/1.0e4
    z= z/1.0e4
    ;x= x/1.0e8
    ;y= y/1.0e8
    ;z= z/1.0e8
    ;m= m/1.0e22

    ; 1 galaxy
    if keyword_set(numpart) and keyword_set(startid) then begin
	x= fload_1gal_gas_xyz('x',startid,numpart,center=[0,0,0])
	y= fload_1gal_gas_xyz('y',startid,numpart,center=[0,0,0])
	z= fload_1gal_gas_xyz('z',startid,numpart,center=[0,0,0])
	m= fload_1gal_gas_mass(startid,numpart)
	hsml= fload_1gal_gas_hsml(startid,numpart)
    endif


    if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)


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


