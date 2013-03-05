pro produce_all_images, junk


        ; Make Gas and All star images with the center on the primary galaxy
        ; ----------------------------------
	;make_image, 100, filename='gassd1.eps'
	;make_image, 110, filename='gassd2.eps'
	;make_image, 130, filename='gassd3.eps'
	;make_image, 200, filename='gassd4.eps'
	;make_image, 300, filename='gassd5.eps'
	;make_image, 330, filename='gassd6.eps'
	;make_image, 360, filename='gassd7.eps'
	;make_image, 380, filename='gassd8.eps'
	;make_image, 400, filename='gassd9.eps'
	;make_image, 500, filename='gassd10.eps'
	;make_image, 600, filename='gassd11.eps'

	make_image, 600, filename='gassd12.eps', /xz



        ; Make All star images with the center on the primary galaxy
        ; ----------------------------------
        ;make_image, 100, filename='stellarsd1.eps', /diskimage
        ;make_image, 110, filename='stellarsd2.eps', /diskimage
        ;make_image, 130, filename='stellarsd3.eps', /diskimage
        ;make_image, 200, filename='stellarsd4.eps', /diskimage
        ;make_image, 300, filename='stellarsd5.eps', /diskimage
        ;make_image, 330, filename='stellarsd6.eps', /diskimage
        ;make_image, 360, filename='stellarsd7.eps', /diskimage
        ;make_image, 380, filename='stellarsd8.eps', /diskimage
        ;make_image, 400, filename='stellarsd9.eps', /diskimage
        ;make_image, 500, filename='stellarsd10.eps', /diskimage
        ;make_image, 600, filename='stellarsd11.eps', /diskimage

        make_image, 600, filename='stellarsd12.eps', /xz, /diskimage




	; Make just satellite image
	; ---------------------------
        ;make_1gal_image, 100, filename='satgsd1.eps'
        ;make_1gal_image, 110, filename='satgsd2.eps'
        ;make_1gal_image, 130, filename='satgsd3.eps'
        ;make_1gal_image, 200, filename='satgsd4.eps'
        ;make_1gal_image, 280, filename='satgsd5.eps'
        ;make_1gal_image, 300, filename='satgsd6.eps'
        ;make_1gal_image, 360, filename='satgsd7.eps'
        ;make_1gal_image, 380, filename='satgsd8.eps'
        ;make_1gal_image, 400, filename='satgsd9.eps'
        ;make_1gal_image, 420, filename='satgsd10.eps'
        ;make_1gal_image, 500, filename='satgsd11.eps'

        ;make_1gal_image, 500, filename='satgsd12.eps', /xz
	

end






pro make_image, snapnum, filename=filename, xz=xz, diskimage=diskimage

        frun= 'execute/G3G1-u3'
        xlen= 60.0
        sendto= 'ps'

        ; Individual Galaxy Images (note that these ONLY show each galaxy)
        ; --------------------------
        ; G3
        startid= 1
        numpart= 240000
        gstartid= 1
        gnumpart= 50000

        rotate_theta= 30.0
        rotate_phi= 0.0

        ok=fload_snapshot(frun,snapnum)
        cen1= fload_1gal_center(startid,numpart)

        ; get center info
        ; --------------- 
        read_centerfiles, frun, time=time, csep=csep
        idx= where(time le fload_time(1))
        if idx(0) ne -1 then begin
                track_to_draw= transpose(csep[*,idx])
        endif

        if snapnum gt 380 then track_to_draw= 0

        if not keyword_set(diskimage) then begin
                if keyword_set(xz) then begin 
                  contour_gas, frun, snapnum, xlen, sendto, /loadedsnap, center=cen1, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
                                                /nolabels, /pubstyle, filename=filename, track_to_draw= track_to_draw, /xz, msg='side view'
                endif else begin
                  contour_gas, frun, snapnum, xlen, sendto, /loadedsnap, center=cen1, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
                                                /nolabels, /pubstyle, filename=filename, track_to_draw= track_to_draw, msg=' '
                endelse
        endif else begin
                if keyword_set(xz) then begin
                  contour_allstars, frun, snapnum, xlen, sendto, /loadedsnap, center=cen1, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
                                                /nolabels, /pubstyle, filename=filename, track_to_draw= track_to_draw, /xz, msg='side view'
                endif else begin
                  contour_allstars, frun, snapnum, xlen, sendto, /loadedsnap, center=cen1, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
                                                /nolabels, /pubstyle, filename=filename, track_to_draw= track_to_draw, msg=' '
                endelse
        endelse


end












pro make_images, junk

        frun= 'execute/G3il-u1a'
        xlen= 50.0
        sendto= 'ps'
        snapnum= 1200
        ;snapnum= 999

        rotate_theta= 0.0
        rotate_phi= 0.0


        contour_gas, frun, snapnum, xlen, sendto, /center_it, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
                                               /nolabels, /pubstyle, filename='gasiso.eps', /xz, msg=' '
        contour_allstars, frun, snapnum, xlen, sendto, /center_it, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
                                               /nolabels, /pubstyle, filename='starsiso.eps', /xz, msg=' '


	;frun= 'execute/G3G1-u3'
        ;xlen= 50.0
        ;sendto= 'ps'
	;snapnum= 100

        ;rotate_theta= 30.0
        ;rotate_phi= 0.0


	;contour_gas, frun, snapnum, xlen, sendto, /center_it, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
	;					/nolabels, /pubstyle, filename='gas31.eps', /xz, msg=' '
	;contour_allstars, frun, snapnum, xlen, sendto, /center_it, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
	;					/nolabels, /pubstyle, filename='stars31.eps', /xz, msg=' '

	;frun= 'execute/G3G2-u3-rerun'
        ;xlen= 50.0
        ;sendto= 'ps'
	;snapnum= 75

        ;rotate_theta= 30.0
        ;rotate_phi= 0.0


	;contour_gas, frun, snapnum, xlen, sendto, /center_it, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
	;					/nolabels, /pubstyle, filename='gas32.eps', /xz, msg=' '
	;contour_allstars, frun, snapnum, xlen, sendto, /center_it, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
	;					/nolabels, /pubstyle, filename='stars32.eps', /xz, msg=' '


	;frun= 'execute/G3G3b-u1-cont'
        ;xlen= 50.0
        ;sendto= 'ps'
	;snapnum= 600

        ;rotate_theta= 30.0
        ;rotate_phi= 0.0


	;contour_gas, frun, snapnum, xlen, sendto, /center_it, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
	;					/nolabels, /pubstyle, filename='gas33.eps', /xz, msg=' '
	;contour_allstars, frun, snapnum, xlen, sendto, /center_it, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
	;					/nolabels, /pubstyle, filename='stars33.eps', /xz, msg=' '



end








pro make_1gal_image, snapnum, filename=filename, xz=xz, diskimage=diskimage

        frun= 'execute/G3G1-u3'
        xlen= 60.0
        sendto= 'ps'

        ; Individual Galaxy Images (note that these ONLY show each galaxy)
        ; --------------------------
        ; G3
        startid= 1
        numpart= 240000
        gstartid= 1
        gnumpart= 50000

	;G3
	satstartid= 240001
	satnumpart= 51000

        rotate_theta= 30.0
        rotate_phi= 0.0

        ok=fload_snapshot(frun,snapnum)
        cen1= fload_1gal_center(startid,numpart)

        ; get center info
        ; ---------------
        read_centerfiles, frun, time=time, csep=csep
        idx= where(time le fload_time(1))
        if idx(0) ne -1 then begin
                track_to_draw= transpose(csep[*,idx])
        endif

        if snapnum gt 380 then track_to_draw= 0

        if not keyword_set(diskimage) then begin
                if keyword_set(xz) then begin
                  contour_gas, frun, snapnum, xlen, sendto, /loadedsnap, center=cen1, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
						startid= satstartid, numpart=satnumpart, $
                                                /nolabels, /pubstyle, filename=filename, track_to_draw= track_to_draw, /xz
                endif else begin
                  contour_gas, frun, snapnum, xlen, sendto, /loadedsnap, center=cen1, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
						startid= satstartid, numpart=satnumpart, $
                                                /nolabels, /pubstyle, filename=filename, track_to_draw= track_to_draw
                endelse
        endif else begin
                if keyword_set(xz) then begin
                  contour_allstars, frun, snapnum, xlen, sendto, /loadedsnap, center=cen1, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
						startid= satstartid, numpart=satnumpart, $
                                                /nolabels, /pubstyle, filename=filename, track_to_draw= track_to_draw, /xz
                endif else begin
                  contour_allstars, frun, snapnum, xlen, sendto, /loadedsnap, center=cen1, rotate_theta= rotate_theta, rotate_phi= rotate_phi, $
						startid= satstartid, numpart=satnumpart, $
                                                /nolabels, /pubstyle, filename=filename, track_to_draw= track_to_draw
                endelse
        endelse


end









pro read_centerfiles, filename, time=time, csep=csep

filename=filename+"/centers.txt"

; read this with the following
openr, 1, filename
junk= ''
lines= 0
readf, 1, junk     ; galaxy 1
print, junk
readf, 1, lines
data= fltarr(4,lines)
readf, 1, data
time= data[0,*]
gal1_xyz= data[1:3,*]
readf, 1, junk     ; galaxy 1
print, junk
readf, 1, lines
data= fltarr(4,lines)
readf, 1, data
gal2_xyz= data[1:3,*]
close, 1


csep= fltarr(3,lines)

; wrt gal 1
csep[0,*] = gal2_xyz[0,*] - gal1_xyz[0,*]
csep[1,*] = gal2_xyz[1,*] - gal1_xyz[1,*]
csep[2,*] = gal2_xyz[2,*] - gal1_xyz[2,*]
; wrt gal 2
;csep[0,*] = gal1_xyz[0,*] - gal2_xyz[0,*]
;csep[1,*] = gal1_xyz[1,*] - gal2_xyz[1,*]
;csep[2,*] = gal1_xyz[2,*] - gal2_xyz[2,*]


end
























;===========================================================================
;===========================================================================




pro produce_iso_images, junk


	snapnum= 150

	frun= 'execute/G0i-u1a'
	make_iso_images, frun, snapnum


        frun= 'execute/G1i-u1a'
        make_iso_images, frun, snapnum


        frun= 'execute/G2im-u1a'
        make_iso_images, frun, snapnum


        frun= 'execute/G3il-u1a'
        make_iso_images, frun, snapnum



end






pro make_iso_images, frun, snapnum, filename=filename, xz=xz, diskimage=diskimage

        xlen= 30.0
        sendto= 'ps'

        ok=fload_snapshot(frun,snapnum)

	runid= fload_fid(1)

        ;contour_gas, frun, snapnum, xlen, sendto, msg='gas (side view)', /loadedsnap, /nolabels, /pubstyle, filename=runid+'gxz.eps', /xz
	;contour_gas, frun, snapnum, xlen, sendto, msg='gas', /loadedsnap, /nolabels, /pubstyle, filename=runid+'gxy.eps'

	;contour_allstars, frun, snapnum, xlen, sendto, msg='stars (side view)', /loadedsnap, /nolabels, /pubstyle, filename=runid+'sxz.eps', /xz
	contour_allstars, frun, snapnum, xlen, sendto, msg='stars', /loadedsnap, /nolabels, /pubstyle, filename=runid+'sxy.eps'


	;make_fancy_xz_image, 1


end






pro make_fancy_xz_image, junk

		xlen= 30.0
		runid= fload_fid(1)


		; grab stars xz image
		; --------------------
                x= fload_allstars_xyz('x',center=[0,0,0])
                y= fload_allstars_xyz('y',center=[0,0,0])
                z= fload_allstars_xyz('z',center=[0,0,0])
                m= fload_allstars_mass(1)

                contour_makepic, x, y, z, m, xlen, center=[0,0,0], NxNImage=NxNImage, /xz

		; assumes 480 x 480 images
                StarsImage= NxNImage[*,121:360]


		; grab gas xz image
		; --------------------
                x= fload_gas_xyz('x',center=[0,0,0])
                y= fload_gas_xyz('y',center=[0,0,0])
                z= fload_gas_xyz('z',center=[0,0,0])
                m= fload_gas_mass(1) - fload_gas_mfs(1)

                contour_makepic, x, y, z, m, xlen, center=[0,0,0], NxNImage=NxNImage, /xz

                GasImage= NxNImage[*,121:360]




	; now plot the thing
	; -------------------
            initialize_plotinfo, 1
            setup_plot_stuff, 'ps', filename=runid+'_xz.eps', colortable= 4


	; first one
	; ---------
        x0= 0.01
        x1= 0.99
        y0= 0.01
        y1= 0.50
        !p.position=[x0,y0,x1,y1]

            ;  place the image down
            ; ----------------------
            tv, StarsImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

            ; creates axes and plot style
            ; ---------------------------
            plot,[0],[0], psym=3,  $
                      xrange=[-xlen,xlen], $
                      yrange=[-0.5*xlen,0.5*xlen], $
                      /NOERASE, $
                      color= 0, $
                      xcharsize=0.01, ycharsize=0.01, $
                      xthick=4.0, ythick=4.0, $
                      xstyle=1, ystyle=1, $
                      xtickformat='(a1)', ytickformat='(a1)', $
                      /normal, $
                      /nodata



        ; first one
        ; ---------
        x0= 0.01
        x1= 0.99
        y0= 0.50
        y1= 0.99
        !p.position=[x0,y0,x1,y1]

            ;  place the image down
            ; ----------------------
            tv, GasImage, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

            ; creates axes and plot style
            ; ---------------------------
            plot,[0],[0], psym=3,  $
                      xrange=[-xlen,xlen], $
                      yrange=[-0.5*xlen,0.5*xlen], $
                      /NOERASE, $
                      color= 0, $
                      xcharsize=0.01, ycharsize=0.01, $
                      xthick=4.0, ythick=4.0, $
                      xstyle=1, ystyle=1, $
                      xtickformat='(a1)', ytickformat='(a1)', $
                      /normal, $
                      /nodata




	; info
	; ------
	;xyouts, 0.07, 0.90, fload_timelbl(1,2,/noteq), /normal, size= 1.7, charthick=3.0, color= 0
	xyouts, 0.07, 0.92, 'gas', /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white
	xyouts, 0.07, 0.42, 'stars', /normal, size= 1.7, charthick=3.0, color= 0   ; 0=black, 1=white


	device, /close


end








;============================================================================












;============================================================================




;random scripts for mino merger paper,
;or talk
pro random_shit, junk


contour_gas, "execute/G0i-u1a", 200, 40.0, 'ps', filename='g0.eps', /nolabels, /pubstyle
contour_gas, "execute/G1i-u1a", 200, 40.0, 'ps', filename='g1.eps', /nolabels, /pubstyle
contour_gas, "execute/G2im-u1a", 200, 40.0, 'ps', filename='g2.eps', /nolabels, /pubstyle
contour_gas, "execute/G3il-u1a", 200, 40.0, 'ps', filename='g3.eps', /nolabels, /pubstyle


end




