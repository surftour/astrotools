pro contour_separategals, frun, snapnum, xlen, sendto, $
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
   print, "contour_separategals, frun, snapnum, xlen, sendto, /xz, /yz, /colorbar, "
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




startid_1= 1
numpart_1= 200001

startid_2= 200002
numpart_2= 200001

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




    ;
    ; galaxy 1
    ;---------------------------------
    x= fload_1gal_allstars_xyz('x',startid_1,numpart_1,center=[0,0,0])
    y= fload_1gal_allstars_xyz('y',startid_1,numpart_1,center=[0,0,0])
    z= fload_1gal_allstars_xyz('z',startid_1,numpart_1,center=[0,0,0])
    m= fload_1gal_allstars_mass(startid_1,numpart_1)
    hsml= 0.2+0.0*m

    contour_makepic, x, y, z, m, xlen, hsml=hsml, xz=xz, yz=yz, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage, HalfMassSB=HalfMassSB

    GalaxyOne= NxNImage






    ;
    ; galaxy 2
    ;---------------------------------
    x= fload_1gal_allstars_xyz('x',startid_2,numpart_2,center=[0,0,0])
    y= fload_1gal_allstars_xyz('y',startid_2,numpart_2,center=[0,0,0])
    z= fload_1gal_allstars_xyz('z',startid_2,numpart_2,center=[0,0,0])
    m= fload_1gal_allstars_mass(startid_2,numpart_2)
    hsml= 0.2+0.0*m

    contour_makepic, x, y, z, m, xlen, hsml=hsml, xz=xz, yz=yz, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage, HalfMassSB=HalfMassSB

    GalaxyTwo= NxNImage







;------------------

; Now Combine the two images





;Pic= GalaxyOne
;Pic= GalaxyTwo

pixels= (size(GalaxyOne))[1]

Pic = fltarr(3,pixels,pixels)    ; red, green, blue
Pic(0,*,*) = byte(GalaxyOne)
Pic(1,*,*) = byte(GalaxyTwo)
Pic(2,*,*) = 0


WRITE_JPEG, 'test.jpg', Pic, TRUE=1, quality= 95

;--------------------------------
;   Send it to postscript
;
;   Most of this is taken
;   straight from contour_makeplot
;
;
; -----------------------------

	;initialize_plotinfo, 1
	;setup_plot_stuff, sendto, filename=filename, colortable= 1

	set_plot, 'ps'
        device, filename= filename, /encapsulated,/color,bits_per_pixel=8
        device, SET_CHARACTER_SIZE=[200,300], xsize=12, ysize=12
	;loadct, 5


	x0= 0.20
	if keyword_set(colorbar) then x1= 0.85 else x1= 0.97
	y0= 0.15
	y1= 0.97


	!p.position=[x0,y0,x1,y1]
	!p.ticklen=0.03


	; default is xy axis
	; ------------------
	;xtit="!8x!3 (kpc/h)"
	;ytit="!8y!3 (kpc/h)"
	xtit="x (kpc)"
	ytit="y (kpc)"

	xmin= -xlen+center[0]
	xmax=  xlen+center[0]
	ymin= -xlen+center[1]
	ymax=  xlen+center[1]


	;  place the image down
	; ----------------------
	tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal, true= 1


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




        ; -----------------------------------
        ;  Puts std information in upper left
        ; -----------------------------------
        if not keyword_set(nolabels) then begin
           if keyword_set(pubstyle) then begin
                xyouts, 0.23, 0.88, fload_timelbl(1,2), /normal, size= 1.7, charthick=3.0, color= 0
           endif else begin
                ; plot other information
                ; ------------------------
                xyouts, 0.22, 0.88, fload_fid(1), /normal, size= 1.7, charthick=3.0, color= 0
                xyouts, 0.22, 0.83, fload_timelbl(1,3), /normal, size= 1.7, charthick=3.0, color= 0

                if keyword_set(msg) then begin
                   xyouts, 0.65, 0.88, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
                endif
           endelse
        endif else begin
             if keyword_set(pubstyle) then begin
                ;xyouts, 0.07, 0.90, fload_timelbl(1,2,/noteq), /normal, size= 1.7, charthick=3.0, color= 0
                ;xyouts, 0.07, 0.09, fload_timelbl(1,2,/noteq), /normal, size= 1.3, charthick=3.0, color= 0
                xyouts, 0.07, 0.09, fload_timelbl(0.7,2,/noteq), /normal, size= 1.3, charthick=3.0, color= 0
                sidelbl='side= '+strcompress(string(long(2*xlen)), /remove_all)+' kpc'
                ;xyouts, 0.07, 0.05, sidelbl, /normal, size= 1.3, charthick=3.0, color= 0

                ;;xyouts, 0.07, 0.83, 'side view', /normal, size= 1.7, charthick=3.0, color= 0
                ;;xyouts, 0.07, 0.83, fload_fid('firsttwo'), /normal, size= 1.7, charthick=3.0, color= 0

                ;xyouts, 0.07, 0.83, '!6no black hole', /normal, size= 1.7, charthick=3.0, color= 0

             endif else begin
                ; plot other information
                ; ------------------------
                xyouts, 0.07, 0.88, fload_fid(1), /normal, size= 1.7, charthick=3.0, color= 0
                xyouts, 0.07, 0.82, fload_timelbl(1,4), /normal, size= 1.7, charthick=3.0, color= 0
             endelse

             if keyword_set(msg) then begin
                ;xyouts, 0.07, 0.91, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
                ;xyouts, 0.07, 0.83, msg, /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
             endif
        endelse





        ; --------------------
        ; done, close this up
        ; --------------------
        if sendto eq 'ps' then device, /close 



; -------------
;  Done
; -------------



end


