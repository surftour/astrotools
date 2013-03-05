pro contour_stars_plus_gas, frun, snapnum, xlen, sendto, $
                        center=center, $
			use_calc_center=use_calc_center, $
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
                        xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
                        xz=xz, yz=yz

if not keyword_set(xlen) then xlen=100.0
if not keyword_set(sendto) then sendto='ps'
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(frun) then begin
   print, "  "
   print, "contour_stars_plus_gas, frun, snapnum, xlen, sendto, /xz, /yz, /colorbar, "
   print, "              filename=filename, /thumbnail, /fitstoo, /crude, /loadedsnap, "
   print, "              xthickness=xthickness, ythickness=ythickness, zthickness=zthickness"
   print, "              /nolabels, center=center, /plotpts/, /track_to_draw, "
   print, "              /use_calc_center, /old_tj_snap"
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
    hsml= 0.2+0.0*m
    ;print, "WARNING: multiplying the stellar hsml by 2.0"

    ; version 2
    ;r=fload_allstars_xyz('r')
    ;hsml= (r*1.0)
    ;idx= where(hsml lt 0.4)
    ;if idx(0) ne -1 then hsml(idx)= 0.4

    ; version 3
    ;if not keyword_set(particlesonly) then hsml= fload_allstars_hsml(1) else hsml= 0.2+0.0*m






; ----------------------------------------------
;  Rather then the generic contour plotting procedures,
;    we'll do things by hand.
; ----------------------------------------------


        contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, $
                        filename=filename, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, $
                        pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage, HalfMassSB=HalfMassSB

        Pic= NxNImage
        PIXELS= n_elements(Pic[0,*])


    xmin= -xlen+center[0]
    xmax=  xlen+center[0]
    ymin= -xlen+center[1]
    ymax=  xlen+center[1]




        ; -----------------------------
        ;   Send it to postscript
        ; -----------------------------

        initialize_plotinfo, 1

        setup_plot_stuff, sendto, filename=filename, colortable= 4
        ;setup_plot_stuff, sendto, filename=filename, colortable= 3
        ;setup_plot_stuff, sendto, filename=filename, colortable= 1
        ;setup_plot_stuff, sendto, filename=filename, colortable= 0
        ;setup_plot_stuff, sendto, filename=filename, colortable= 4, imgxsize=14
        ;setup_plot_stuff, sendto, filename=filename, colortable= 0, imgxsize=14


        if keyword_set(nolabels) then begin
            x0= 0.01
            x1= 0.99
            y0= 0.01
            y1= 0.99
        endif else begin
            x0= 0.20
            if keyword_set(colorbar) then x1= 0.85 else x1= 0.97
            y0= 0.15
            y1= 0.97
        endelse


        !p.position=[x0,y0,x1,y1]
        !p.ticklen=0.03



        ;  place the image down
        ; ----------------------
        tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


        ; creates axes and plot style
        ; ---------------------------
        if keyword_set(nolabels) then begin
                plot,[0],[0], psym=3, xstyle=1, ystyle=1, color= 0, /noerase, $
                      xrange=[xmin,xmax], yrange=[ymin,ymax], $
                      xcharsize=0.01, ycharsize=0.01, $
                      xthick=4.0, ythick=4.0, /normal, /nodata, $
                      xtickformat='(a1)', ytickformat='(a1)'
        endif else begin
                plot,[0],[0], psym=3, xstyle=1, ystyle=1, color= 0, /noerase, $
                      xrange=[xmin,xmax], yrange=[ymin,ymax], $
                      xcharsize=1.50, ycharsize=1.50, $
                      xthick=4.0, ythick=4.0, /normal, /nodata, $
                      charthick=3.0, $
                      xtitle=xtit, ytitle=ytit
        endelse 







        ; --------------------
        ; done, close this up
        ; --------------------
        if sendto eq 'ps' then device, /close






; -------------
;  Done
; -------------





end


