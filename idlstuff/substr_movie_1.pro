;======================================================================
;
;
;
;======================================================================
pro mov_2, junk


frun= 'data/substr/t_23_0'
;frun= 'data/substr/t_23_1'
;frun= 'data/substr/t_23_2'
imagedir= frun+'/movieimages'

imagename= 'image'

startsnap= 0
;startsnap= 42
;endsnap= 10
;endsnap= 44
endsnap= 98
;endsnap= 140
;endsnap= 201
;endsnap= 202
;endsnap= 300

;xlen= 120.0

pixels= 640L


startid= 2150001L
numpart= 22748            ; std, _0
;numpart= 2274            ; _1
;numpart= 227480          ; _2

tempimagename= 'temp_'+imagename


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






;--------------------------------------
;  Make images
;--------------------------------------


for i=startsnap, endsnap do begin


	ok=fload_snapshot_bh(frun,i,/nopot_in_snap,/skip_center)


        ; get XY image
        ; --------------------
        x= fload_1gal_halo_xyz('x',startid,numpart)
        y= fload_1gal_halo_xyz('y',startid,numpart)
        z= fload_1gal_halo_xyz('z',startid,numpart)
        m= fload_1gal_halo_mass(startid,numpart)
        ;hsml= fload_allstars_hsml(1)
        hsml= 0.2 + m*0.0

	set_maxden= 1.0e+0
	set_dynrng= 1.0e+6


	; -------------------

        movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, pixels=pixels, center=center,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        XYPic= NxNImage


	; -------------------


        movie_makepic, y, z, x, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, pixels=pixels, center=center,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        YZPic= NxNImage


        ; -------------------


	ilbl= '000'+strcompress(string(i),/remove_all)
	ilbl= strmid(ilbl,strlen(ilbl)-3,3)


	;   write jpg files
	; --------------------

	; gas and stars image
	image3d = BYTARR(3, 2*pixels, pixels)
	image3d(0, *, *)= [c0(XYPic), c0(YZPic)]
	image3d(1, *, *)= [c1(XYPic), c1(YZPic)]
	image3d(2, *, *)= [c2(XYPic), c2(YZPic)]


	; black line to separate the two panels
	image3d[*,pixels:pixels+1,*]= 0


        ;filename= imagedir+'/'+tempimagename+'_'+ilbl+'.jpg'
        filename= imagedir+'/temp.jpg'
	WRITE_JPEG, filename, image3d, TRUE=1, quality= 95



	; add text onto the figure
	; method #2 - good
	;
	; basically uses the convert command (could also use mogrify)
	;
	;convert -font arial -antialias -pointsize 30 -fill white \                                               
        ;     -draw 'text 20,50 "z = 2.0"'  \                                                                 
        ;     mypic.jpg  annotated_pic.jpg        
	;

	fout= imagedir+'/'+imagename+'_'+ilbl+'.jpg'
	timelbl=fload_timelbl(0.7,2,/noteq)

	cmd = 'convert -font arial -antialias -pointsize 30 -fill black -draw '
	cmd = cmd + "'"
	cmd = cmd + 'text 20,46 '
	cmd = cmd + """
	cmd = cmd + timelbl
	cmd = cmd + """
	cmd = cmd + "'"
	cmd = cmd + ' -quality 95 '
	cmd = cmd + filename + " " + fout

	print, cmd
	spawn, cmd


endfor


; -------------
;  Done
; -------------


end





;===============================================================================



