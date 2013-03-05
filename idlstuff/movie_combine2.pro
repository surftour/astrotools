;======================================================================
;
;
;   Script to combine two images into one.
;
;
;
;======================================================================
pro movie_combine_2, junk


outdir= 'movie_090/'
;outdir= 'movie_150/'
spawn, "mkdir "+outdir
outimage= outdir+'img'

;img1= 'data/ds/vc3vc3e_no/movie_v/pic_gas_'
;img2= 'data/ds/vc3vc3e_2/movie_v/pic_gas_'
;img1= 'data/ds/vc3vc3h_no/movie_v/pic_gas_'
;img2= 'data/ds/vc3vc3h_2/movie_v/pic_gas_'
;img1= 'data/minor/Sbfg0.4Scfg0.4_090/movie_v2/pic_stars_v2_'
;img2= 'data/minor/Sbfg0.4Scfg0.4_090/movie_v2/pic_gas_'
img1= 'data/minor/Sbfg0.4Scfg0.4_150/movie_v2/pic_stars_v2_'
img2= 'data/minor/Sbfg0.4Scfg0.4_150/movie_v2/pic_gas_'

startsnap= 0
;startsnap= 42
;startsnap= 202
;endsnap= 2
;endsnap= 10
;endsnap= 44
;endsnap= 140
;endsnap= 201
;endsnap= 202
;endsnap= 300
endsnap= 347
;endsnap= 399


; get plotting stuff ready
;SET_PLOT, 'ps'
set_plot, 'z'
dxx= 120 & dyy= 40
device, set_resolution= [dxx, dyy]

LOADCT, 0
v1=[0,255]
v2=[0,255]
v3=[0,255]
tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
TVLCT, bw0, bw1, bw, /GET

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


;



;--------------------------------------
;  Make images
;--------------------------------------


for i=startsnap, endsnap do begin


	; -------------------


	ilbl= '0000'+strcompress(string(i),/remove_all)
	ilbl= strmid(ilbl,strlen(ilbl)-4,4)



	;   read jpg file 1
	; --------------------

	; original image
	fimg1= img1+ilbl+'.jpg'


	;
	; put time on the image
	;
	;ok=fload_snapshot_bh(frun,i,/header_only)
	;timelbl=fload_timelbl(0.7,2,/noteq)
	;cmd = 'convert -font Courier-Bold -antialias -pointsize 30 -fill white -draw '
	;cmd = cmd + "'"
	;cmd = cmd + 'text 20,46 '
	;cmd = cmd + """
	;cmd = cmd + timelbl
	;cmd = cmd + """
	;cmd = cmd + "'"
	;cmd = cmd + ' -quality 95 '
	;cmd = cmd + fimg1 + " image1.jpg"
	;print, cmd
	;spawn, cmd


	;
	; put "no BH" on there
	;
        ;cmd = "convert -font Courier-Bold -antialias -pointsize 50 -fill white -draw 'text 400,400 "
        cmd = "convert -font Courier-Bold -antialias -pointsize 20 -fill black -draw 'text 400,400 "
        cmd = cmd + """
        ;cmd = cmd + 'No BH'
        cmd = cmd + 'stars'
        cmd = cmd + """
        cmd = cmd + "'"
        cmd = cmd + ' -quality 95 '
        cmd = cmd + fimg1 + " image1.jpg"
        print, cmd
        spawn, cmd


	READ_JPEG, "image1.jpg", image1, TRUE=1




	;   read jpg file 2
	; --------------------

	; original image
	fimg2= img2+ilbl+'.jpg'

        ;cmd = "convert -font Courier-Bold -antialias -pointsize 50 -fill white -draw 'text 480,400 "
        cmd = "convert -font Courier-Bold -antialias -pointsize 20 -fill black -draw 'text 480,400 "
        cmd = cmd + """
        ;cmd = cmd + 'BH'
        cmd = cmd + 'gas'
        cmd = cmd + """
        cmd = cmd + "'"
        cmd = cmd + ' -quality 95 '
        cmd = cmd + fimg2 + " image2.jpg"
        print, cmd
        spawn, cmd


	READ_JPEG, "image2.jpg", image2, TRUE=1




	;   combine images
	; -------------------
	xpixels= (size(image1))[2]
	ypixels= (size(image1))[3]
	image_final = BYTARR(3, 2*xpixels, ypixels)
	;image_final(0, *, *)= [image1[0,*,*], image2[0,*,*]]
	;image_final(1, *, *)= [image1[1,*,*], image2[1,*,*]]
	;image_final(2, *, *)= [image1[2,*,*], image2[2,*,*]]
	image_final[*,0:xpixels-1,*]=image1
	image_final[*,xpixels:2*xpixels-1,*]=image2


	; white line to separate the two panels
	;image_final[*,xpixels:xpixels+1,*]= 255
	image_final[*,xpixels:xpixels+1,*]= 0



	; write jpeg
	outfile= outimage+ilbl+'.jpg'
	WRITE_JPEG, outfile, image_final, TRUE=1, quality= 95



endfor


; -------------
;  Done
; -------------


end







;======================================================================
;
;
;   Script to combine four images into one.
;
;
;
;======================================================================
pro movie_combine_4, junk, thumbnail=thumbnail


;frun= "data1/ring/v3v2b_1x_a"
frun= junk


outdir= frun+'/movie_4/'
spawn, "mkdir "+outdir
outimage= outdir+'img_'

;
;  -----------------
;  |       |       |
;  |  1    |  2    |
;  |       |       |
;  -----------------
;  |       |       |
;  |  3    |  4    |
;  |       |       |
;  -----------------
;
img1= frun+'/movie_images/img_stars_top_'
img2= frun+'/movie_images/img_gas_top_'
img3= frun+'/movie_images/img_stars_side_'
img4= frun+'/movie_images/img_gas_side_'

startsnap= 0
;startsnap= 42
;startsnap= 202

spawn, "/bin/ls "+frun+"/snapshot_* | wc ",result
endsnap=long(result[0])-1

;endsnap= 2
;endsnap= 10
;endsnap= 44
;endsnap= 140
;endsnap= 201
;endsnap= 202
;endsnap= 300
;endsnap= 347
;endsnap= 399





;--------------------------------------
;  Make images
;--------------------------------------


for i=startsnap, endsnap do begin


	; -------------------


	ilbl= '0000'+strcompress(string(i),/remove_all)
	ilbl= strmid(ilbl,strlen(ilbl)-4,4)


	if keyword_set(thumbnail) then begin
		cmd= 'cp '+img1+ilbl+'.jpg  image11.jpg'
		spawn, cmd
		cmd= 'cp '+img2+ilbl+'.jpg  image22.jpg'
		spawn, cmd
		cmd= 'cp '+img3+ilbl+'.jpg  image33.jpg'
		spawn, cmd
		cmd= 'cp '+img4+ilbl+'.jpg  image44.jpg'
		spawn, cmd
		goto, combineimages
	endif

	;   read jpg file 1
	; --------------------

	; original image
	fimg1= img1+ilbl+'.jpg'

	if file_size(fimg1) le 0 then continue
	img_add_text_to_jpg, fimg1, "image1.jpg", "face-on", 40, 100
	img_add_text_to_jpg, "image1.jpg", "image11.jpg", "stars", 40, 40




	;   read jpg file 2
	; --------------------

	; original image
	fimg2= img2+ilbl+'.jpg'

	if file_size(fimg2) le 0 then continue
	img_add_text_to_jpg, fimg2, "image2.jpg", "face-on", 40, 100
	img_add_text_to_jpg, "image2.jpg", "image22.jpg", "gas", 40, 40




	;   read jpg file 3
	; --------------------

	; original image
	fimg3= img3+ilbl+'.jpg'

	if file_size(fimg3) le 0 then continue
	img_add_text_to_jpg, fimg3, "image3.jpg", "edge-on", 40, 100
	img_add_text_to_jpg, "image3.jpg", "image33.jpg", "stars", 40, 40




	;   read jpg file 2
	; --------------------

	; original image
	fimg4= img4+ilbl+'.jpg'

	if file_size(fimg4) le 0 then continue
	img_add_text_to_jpg, fimg4, "image4.jpg", "edge-on", 40, 100
	img_add_text_to_jpg, "image4.jpg", "image44.jpg", "gas", 40, 40





combineimages:
	; -------------------
	; -------------------
	;   combine images
	; -------------------
	; -------------------

	READ_JPEG, "image11.jpg", image1, TRUE=1
	READ_JPEG, "image22.jpg", image2, TRUE=1
	READ_JPEG, "image33.jpg", image3, TRUE=1
	READ_JPEG, "image44.jpg", image4, TRUE=1

	xpixels= (size(image1))[2]
	ypixels= (size(image1))[3]
	image_final = BYTARR(3, 2*xpixels, 2*ypixels)
	;image_final(0, *, *)= [image1[0,*,*], image2[0,*,*]]
	;image_final(1, *, *)= [image1[1,*,*], image2[1,*,*]]
	;image_final(2, *, *)= [image1[2,*,*], image2[2,*,*]]
	pix0= 0
	pix1= xpixels-1
	pix2= xpixels
	pix3= 2*xpixels-1
	image_final[*,0:xpixels-1,ypixels:2*ypixels-1]=image1
	image_final[*,xpixels:2*xpixels-1,ypixels:2*ypixels-1]=image2
	image_final[*,0:xpixels-1,0:ypixels-1]=image3
	image_final[*,xpixels:2*xpixels-1,0:ypixels-1]=image4


	; vertical white line
	;image_final[*,xpixels:xpixels+1,*]= 255
	image_final[*,xpixels:xpixels+1,*]= 0

	; horizontal white line
	;image_final[*,*,ypixels:ypixels+1]= 255
	image_final[*,*,ypixels:ypixels+1]= 0


	; write jpeg
	outfile= outimage+ilbl+'.jpg'
	WRITE_JPEG, outfile, image_final, TRUE=1, quality= 95



endfor


; -------------
;  Done
; -------------


end


