
;===============================================================================




pro do3, junk


frun= "/n/home/tcox/data1/lg_m33/SbSbm33_7/"

imagename= 'image4'
xlen= 280.0
pixels= 1200L




;--------------------------------------
;  Make images
;--------------------------------------


make_one_image, frun, 70, pixels, xlen, "image1_temp.jpg", "+1.8 Gyr"
;
;convert -antialias -pointsize 60 -fill black -draw 'text  40, 86 "+1.5 Gyr"'  -quality 95   image1_temp.jpg  image1.jpg
;convert -antialias -pointsize 30 -fill black -draw 'text 370,820 "Milky Way"' -quality 95   image1.jpg       image1_a.jpg
;convert -antialias -pointsize 30 -fill black -draw 'text 720,380 "M31"'       -quality 95   image1_a.jpg     image1_b.jpg
;convert -antialias -pointsize 30 -fill black -draw 'text 890,630 "M33"'       -quality 95   image1_b.jpg     image1_c.jpg
;

make_one_image, frun, 168, pixels, xlen, "image2_temp.jpg", "+2.5 Gyr"
;
;convert -antialias -pointsize 60 -fill black -draw 'text  40, 86 "+2.5 Gyr"'  -quality 95   image2_temp.jpg  image2_a.jpg
;convert -antialias -pointsize 30 -fill black -draw 'text 890,680 "Milky Way"' -quality 95   image2_a.jpg     image2_b.jpg
;convert -antialias -pointsize 30 -fill black -draw 'text 370,500 "M31"'       -quality 95   image2_b.jpg     image2_c.jpg
;convert -antialias -pointsize 30 -fill black -draw 'text 260,850 "M33"'       -quality 95   image2_c.jpg     image2.jpg
;

make_one_image, frun, 700, pixels, xlen, "image3_temp.jpg", "+7.5 Gyr"
;
;convert -antialias -pointsize 60 -fill black -draw 'text  40, 86 "+7.5 Gyr"'  -quality 95   image3_temp.jpg  image3_a.jpg
;convert -antialias -pointsize 30 -fill black -draw 'text 920,500 "Milkomeda"' -quality 95   image3_a.jpg     image3_b.jpg
;convert -antialias -pointsize 30 -fill black -draw 'text 320,650 "M33"'       -quality 95   image3_b.jpg     image3.jpg
;


end




;------------------------------------------------





pro make_one_image, frun, snapnum, pixels, xlen, imagename, msg


;-------------------------------------


; get plotting stuff ready
;SET_PLOT, 'ps'
set_plot, 'z'
dxx= 120 & dyy= 40
device, set_resolution= [dxx, dyy]

LOADCT, 1
v1=[0,255] 
v2=[0,255] 
v3=[0,255]
tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
TVLCT, blue0, blue1, blue2, /GET

LOADCT, 3
v1=[0,255]
v2=[0,255]
v3=[0,255]
tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
TVLCT, red0, red1, red2, /GET

LOADCT, 4
;LOADCT, 13
v1=[0,255]
v2=[0,255]
v3=[0,255]
tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
TVLCT, c0, c1, c2, /GET



	if snapnum ge 1000 then begin
		ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap,/skip_center,/do_four)
	endif else begin
		ok=fload_snapshot_bh(frun,snapnum,/nopot_in_snap,/skip_center)
	endelse

        ; get allstar image
        ; --------------------
        x= fload_allstars_xyz('x')
        y= fload_allstars_xyz('y')
        z= fload_allstars_xyz('z')
        m= fload_allstars_mass(1)
        ;hsml= fload_allstars_hsml(1)
        hsml= 0.2 + m*0.0

	set_maxden= 1.0e+0
	set_dynrng= 1.0e+6

	; top view
        movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, pixels=pixels, center=[0.0, 0.0, 0.0],$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        StellarPic= NxNImage


	;   write jpg files
	; --------------------
	image3d = BYTARR(3, pixels, pixels)

	image3d(0, *, *)= [c0(StellarPic)]
	image3d(1, *, *)= [c1(StellarPic)]
	image3d(2, *, *)= [c2(StellarPic)]


	;WRITE_JPEG, 'movietempimage.jpg', image3d, TRUE=1, quality= 95
	WRITE_JPEG, imagename, image3d, TRUE=1, quality= 95

	; add msg to the figure
	;cmd = 'convert -font arial -antialias -pointsize 30 -fill black -draw '
	;cmd = cmd + "'"
	;cmd = cmd + 'text 20,46 '
	;cmd = cmd + """
	;cmd = cmd + msg
	;cmd = cmd + """
	;cmd = cmd + "'"
	;cmd = cmd + ' -quality 95 '
	;cmd = cmd + "  movietempimage.jpg  " + imagename

	;print, cmd
	;spawn, cmd


end






;===============================================================================






;===============================================================================



