;======================================================================
;
;
;   General scripts to make movie images.  You still will need to
;   compile all the images into a movie using either ffmpeg, QuickTimePro
;   or some equivalent program.
;
;   and you might need to combine them into multiple panels too.
;
;
;
;======================================================================



pro make_movie_panels, frun


frun="data1/pg/pg32_g"

if not keyword_set(frun) then begin
	print, " "
	print, " make_movie, frun"
	print, " "
	return
endif

;
; make directory for images
;
imagedir= frun+'/movie_images_bw'
spawn, "mkdir "+imagedir
imagedir= frun+'/movie_images_ct'
spawn, "mkdir "+imagedir



;
; get snapshot range
;
startsnap= 0
;startsnap= 40
;startsnap= 750
;startsnap= 815

;spawn, "/bin/ls "+frun+"/snapshot_* | wc ",result
;endsnap=long(result[0])-1
;endsnap=10
;endsnap=800
;endsnap=820
endsnap=1250
;endsnap=3400


;
; other parameter choices
;

rotateby= 0.0

;xlen= 350.0
;xlen= 280.0
;xlen= 210.0
;xlen= 140.0
;xlen= 105.0
;xlen= 75.0
xlen= 70.0
;xlen= 52.5
;xlen= 40.0
;xlen= 20.0

;pixels= 640L
;pixels= 480L
;pixels= 420L
;pixels= 320L


;
; standard screen resolutions
;   1280 x  960
;   1920 x 1200  (this is my mac)
;
;
; standard projector resolutions
;
;    800 x   600
;   1024 x   768
;   1400 x  1050
;
;

xpixels= 1024L
ypixels= 768L

;xpixels= 512L
;ypixels= 384L

center= [0.0, 0.0, 0.0]  ; get's reset later

; get center from file
read_center, cfile_time, cfile_center, filename=frun+"/centers_den_550001.txt"


;-------------------------------------


print, "frun= ", frun
print, "imagedir= ", imagedir
print, "startsnap= ", startsnap
print, "endsnap= ", endsnap
print, "xlen= ", xlen
print, "xpixels= ", xpixels
print, "ypixels= ", ypixels
print, "center= ", center


;-------------------------------------


img_load_colormaps, 1





;--------------------------------------
;  Make images
;--------------------------------------


for i=startsnap, endsnap do begin


	if i ge 1000 then begin
		ok=fload_snapshot_bh(frun,i,/nopot_in_snap,/skip_center,/do_four)
	endif else begin
		ok=fload_snapshot_bh(frun,i,/nopot_in_snap,/skip_center)
	endelse

	if ok lt 0 then begin 
		print, " ============================================== "
		print, "  Something is wrong with the file, continuing  "
		print, " ============================================== "
		continue
	endif

        ; -------------------

        ilbl= '0000'+strcompress(string(i),/remove_all)
        ilbl= strmid(ilbl,strlen(ilbl)-4,4)

        ; -------------------


	; center the image
	; -----------------

	; default
	center= [0.0, 0.0, 0.0]  ; get's reset later


	center_bh= 0
	;center_bh= 1
	if keyword_set(center_bh) then begin
		bhid= fload_blackhole_id(1)
        	bhid= long(min(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
        	center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

        	print, "Blackhole ID: ", bhid
        	print, "Blackhole center: ", center_bh

		center= center_bh
	endif


	; set center from file
	time= fload_time(1)
	idx= where(cfile_time ge (time-1.0e-6))
	if idx(0) ne -1 then center= cfile_center(*,idx(0))



	ngas= fload_npart(0)
	if ngas eq 0 then goto, dostars


dostars:

        ; get allstar image
        ; --------------------
        x= fload_allstars_xyz('x')
        y= fload_allstars_xyz('y')
        z= fload_allstars_xyz('z')
	;movie_process_rotation, x, y, z, theta, phi, x_new, y_new, z_new
	x_new= cos(rotateby)*x - sin(rotateby)*z
	y_new= y
	z_new= sin(rotateby)*x + cos(rotateby)*z
        m= fload_allstars_mass(1)
        ;hsml= fload_allstars_hsml(1)
        hsml= 0.2 + m*0.0


	; top view

	;
	; color table
	;
	set_maxden= 1.0e+0
	set_dynrng= 1.0e+5
        movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, xpixels=xpixels, ypixels=ypixels, center=center,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        StellarPicTop= rotate(NxNImage,6)
	image3d= img_map_to_color(StellarPicTop, "colortable", xpixels, ypixels, /iinvert)
	;WRITE_JPEG, 'data1/pg/pg32_g/movie_images_ct/img_stars_top_'+ilbl+'.jpg', image3d, TRUE=1, quality= 95
	WRITE_JPEG, 'temp.jpg', image3d, TRUE=1, quality= 95

        timelbl=fload_timelbl(0.7,2,/noteq)

        cmd = 'convert -font arial -antialias -pointsize 30 -fill white -draw '
        cmd = cmd + "'"
        cmd = cmd + 'text 20,46 '
        cmd = cmd + """
        cmd = cmd + timelbl
        cmd = cmd + " Gyr ""
        cmd = cmd + "'"
        cmd = cmd + ' -quality 95 '
        cmd = cmd + " temp.jpg data1/pg/pg32_g/movie_images_ct/img_stars_top_"+ilbl+".jpg"

        print, cmd
        spawn, cmd


        ;
        ; b/w
        ;
        set_maxden= 1.0e+0 / 50.
        set_dynrng= 1.0e+5 * 10.
        movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, xpixels=xpixels, ypixels=ypixels, center=center,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage, /floatingbckgrnd

        StellarPicTop= rotate(NxNImage,6)
	image3d= img_map_to_color(StellarPicTop, "bw", xpixels, ypixels)
	;WRITE_JPEG, 'data1/pg/pg32_g/movie_images_bw/img_stars_top_'+ilbl+'.jpg', image3d, TRUE=1, quality= 95
	WRITE_JPEG, 'temp.jpg', image3d, TRUE=1, quality= 95


        timelbl=fload_timelbl(0.7,2,/noteq)

        cmd = 'convert -font arial -antialias -pointsize 30 -fill black -draw '
        cmd = cmd + "'"
        cmd = cmd + 'text 20,46 '
        cmd = cmd + """
        cmd = cmd + timelbl
        cmd = cmd + " Gyr ""
        cmd = cmd + "'"
        cmd = cmd + ' -quality 95 '
        cmd = cmd + " temp.jpg data1/pg/pg32_g/movie_images_bw/img_stars_top_"+ilbl+".jpg"

        print, cmd
        spawn, cmd




	openr, 1, imagedir+"/stop", error=err
	close, 1
	if (err eq 0) then begin
	 cmd= '/bin/rm -f '+imagedir+'/stop'
	 spawn, cmd
	 print, " "
	 print, " STOP FILE FOUND - stopping "
	 print, " "
	 return
	endif


endfor


; -------------
;  Done
; -------------


end







;===============================================================================





function img_map_to_color, image, colormap, xpixels, ypixels, iinvert=iinvert

COMMON Colormaps

image3d = BYTARR(3, xpixels, ypixels)


if keyword_set(iinvert) then begin
	print, "inverting image"
	print, "min/max= ", min(image), max(image)
	image= 256-image   ; typically my stuff does it to 254, not 256
	; background is white=1 (now 256-1=255)
	idx= where(image eq 255)
	if idx(0) ne -1 then image(idx)= 0
endif


if colormap eq "blue" then goto, blue
if colormap eq "std" then goto, std
if colormap eq "red" then goto, red
if colormap eq "bw" then goto, bw


;----------------------
; color plane
colplane:
for i=0, xpixels-1 do begin
   for j=0, ypixels-1 do begin
      ;image3d(0,i,j)= ColPlane(image(i,j),imagesecondary(i,j),0)
      ;image3d(1,i,j)= ColPlane(image(i,j),imagesecondary(i,j),1)
      ;image3d(2,i,j)= ColPlane(image(i,j),imagesecondary(i,j),2)
      image3d(0,i,j)= ColPlane(image(i,j),image(i,j),0)
      image3d(1,i,j)= ColPlane(image(i,j),image(i,j),1)
      image3d(2,i,j)= ColPlane(image(i,j),image(i,j),2)
      ;image3d(0,i,j)= ColPlane(image(i,j),0,0)
      ;image3d(1,i,j)= ColPlane(image(i,j),0,1)
      ;image3d(2,i,j)= ColPlane(image(i,j),0,2)
   endfor
endfor
return, image3d


;----------------------
; blue
blue:
image3d(0, *, *)= [blue0(image)]
image3d(1, *, *)= [blue1(image)]
image3d(2, *, *)= [blue2(image)]
return, image3d


;----------------------
; std
std:
image3d(0, *, *)= [c0(image)]
image3d(1, *, *)= [c1(image)]
image3d(2, *, *)= [c2(image)]
return, image3d


;----------------------
; red
red:
image3d(0, *, *)= [red0(image)]
image3d(1, *, *)= [red1(image)]
image3d(2, *, *)= [red2(image)]
return, image3d


;----------------------
; bw
bw:
image3d(0, *, *)= [bw0(image)]
image3d(1, *, *)= [bw1(image)]
image3d(2, *, *)= [bw2(image)]
return, image3d




end






;===============================================================================



pro img_load_colormaps, junk

COMMON Colormaps, ColPlane, blue0, blue1, blue2, red0, red1, red2, c0, c1, c2, bw0, bw1, bw2

; get plotting stuff ready
;SET_PLOT, 'ps'
set_plot, 'z'
dxx= 120 & dyy= 40
device, set_resolution= [dxx, dyy]

LOADCT, 0
;v1=[0,255]
;v2=[0,255]
;v3=[0,255]
;tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
TVLCT, bw0, bw1, bw2, /GET

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



cols=256
ColPlane=fltarr(cols,cols,3)

;openr,1,"colplane2.dat"
openr,1,"/n/home/tcox/idlstuff/image/colplane_flip2.dat"
readu,1,ColPlane
close,1



end





;==========================================================================
;==========================================================================
;==========================================================================
;==========================================================================




pro test, junk

xlen= 70.0
xpixels= 1024L
ypixels= 768L

center= [0.0, 0.0, 0.0]  ; get's reset later

img_load_colormaps, 1

ok=fload_snapshot_bh("data1/pg/pg32_g/",817,/nopot_in_snap,/skip_center)

x= fload_allstars_xyz('x')
y= fload_allstars_xyz('y')
z= fload_allstars_xyz('z')
m= fload_allstars_mass(1)
hsml= 0.2 + m*0.0

set_maxden= 1.0e+0
set_dynrng= 1.0e+5


; test 0 = the standard
; ----------------------
movie_makepic, x, y, z, m, xlen, hsml=hsml, xpixels=xpixels, ypixels=ypixels, center=center,$
                        set_maxden=set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage

StellarPicTop= transpose(NxNImage)
image3d= img_map_to_color(StellarPicTop, "colortable", xpixels, ypixels, /iinvert)
WRITE_JPEG, 'data1/pg/pg32_g/movie_images/test_std.jpg', image3d, TRUE=1, quality= 95

; test 4
; ----------------------
;StellarPicTop= transpose(NxNImage)
;image3d= img_map_to_color(StellarPicTop, "bw", xpixels, ypixels, /iinvert)
;WRITE_JPEG, 'data1/pg/pg32_g/movie_images/test_4.jpg', image3d, TRUE=1, quality= 95

; test 5
; ----------------------
movie_makepic, x, y, z, m, xlen, hsml=hsml, xpixels=xpixels, ypixels=ypixels, center=center,$
                        set_maxden=set_maxden/50., set_dynrng=set_dynrng*10., NxNImage=NxNImage, /floatingbckgrnd

;StellarPicTop= transpose(NxNImage)
StellarPicTop= rotate(NxNImage,1)
image3d= img_map_to_color(StellarPicTop, "bw", xpixels, ypixels)
WRITE_JPEG, 'data1/pg/pg32_g/movie_images/test_5.jpg', image3d, TRUE=1, quality= 95

StellarPicTop= rotate(NxNImage,3)
image3d= img_map_to_color(StellarPicTop, "bw", xpixels, ypixels)
WRITE_JPEG, 'data1/pg/pg32_g/movie_images/test_6.jpg', image3d, TRUE=1, quality= 95

StellarPicTop= rotate(NxNImage,4)
image3d= img_map_to_color(StellarPicTop, "bw", xpixels, ypixels)
WRITE_JPEG, 'data1/pg/pg32_g/movie_images/test_7.jpg', image3d, TRUE=1, quality= 95

StellarPicTop= rotate(NxNImage,6)
image3d= img_map_to_color(StellarPicTop, "bw", xpixels, ypixels)
WRITE_JPEG, 'data1/pg/pg32_g/movie_images/test_8.jpg', image3d, TRUE=1, quality= 95


; test 1
; ----------------------
;movie_makepic, x, y, z, m, xlen, hsml=hsml*5.0, xpixels=xpixels, ypixels=ypixels, center=center,$
;                        set_maxden=set_maxden/10., set_dynrng=set_dynrng/100., NxNImage=NxNImage
;
;StellarPicTop= transpose(NxNImage)
;image3d= img_map_to_color(StellarPicTop, "colortable", xpixels, ypixels, /iinvert)
;WRITE_JPEG, 'data1/pg/pg32_g/movie_images/test_1.jpg', image3d, TRUE=1, quality= 95


; test 2
; ----------------------
;movie_makepic, x, y, z, m, xlen, hsml=hsml*5.0, xpixels=xpixels, ypixels=ypixels, center=center,$
;                        set_maxden=set_maxden, set_dynrng=set_dynrng, NxNImage=NxNImage
;
;StellarPicTop= transpose(NxNImage)
;image3d= img_map_to_color(StellarPicTop, "colortable", xpixels, ypixels, /iinvert)
;WRITE_JPEG, 'data1/pg/pg32_g/movie_images/test_2.jpg', image3d, TRUE=1, quality= 95


; test 3
; ----------------------
;movie_makepic, x, y, z, m, xlen, hsml=hsml*7.5, xpixels=xpixels, ypixels=ypixels, center=center,$
;                        set_maxden=set_maxden, set_dynrng=set_dynrng*10., NxNImage=NxNImage
;
;StellarPicTop= transpose(NxNImage)
;image3d= img_map_to_color(StellarPicTop, "colortable", xpixels, ypixels, /iinvert)
;WRITE_JPEG, 'data1/pg/pg32_g/movie_images/test_3.jpg', image3d, TRUE=1, quality= 95

end










;===========================================================================


pro do_movie_zoom, junk

frun= "data1/pg/pg32_g"

;
; make directory for images
;
imagedir= frun+'/movie_images_zoom'
spawn, "mkdir "+imagedir


;
; original is 1024 x 768
;
; Note: original is different depending on whether
;       it has been "converted" or not, i.e., whether
;       convert was used to put text on there.  This
;       generates the duplication on some of these
;       calls.
;

backgrndclr= 0
;backgrndclr= 255
;READ_JPEG, 'data1/pg/pg32_g/movie_images_bw/img_stars_top_0816.jpg', image3d, TRUE=1
;READ_JPEG, 'data1/pg/pg32_g/movie_images_bw/img_stars_top_0817.jpg', image3d, TRUE=1
READ_JPEG, 'data1/pg/pg32_g/movie_images_ct/img_stars_top_0817.jpg', image3d, TRUE=1
help, image3d



i0= 0
i1= 20      ; frames before & after
i2= 50      ; number of frames for zoom 
i3= 125     ; frames for comparison image

for i=i0,i1 do begin
        ilbl= '0000'+strcompress(string(i),/remove_all)
        ilbl= strmid(ilbl,strlen(ilbl)-4,4)
	WRITE_JPEG, 'data1/pg/pg32_g/movie_images_zoom/img_stars_zoom_'+ilbl+'.jpg', image3d, TRUE=1, quality= 95
	;WRITE_JPEG, 'data1/pg/pg32_g/movie_images_zoom/img_stars_zoom_'+ilbl+'.jpg', image3d
        ilbl= '0000'+strcompress(string(i1+i2+i3+i2+i),/remove_all)
        ilbl= strmid(ilbl,strlen(ilbl)-4,4)
	WRITE_JPEG, 'data1/pg/pg32_g/movie_images_zoom/img_stars_zoom_'+ilbl+'.jpg', image3d, TRUE=1, quality= 95
	;WRITE_JPEG, 'data1/pg/pg32_g/movie_images_zoom/img_stars_zoom_'+ilbl+'.jpg', image3d
endfor

for i=i1,i1+i2 do begin

	di= i-i1
	newx= long(1024 * (100.0 - di)/100.0)
	newy= long(768  * (100.0 - di)/100.0)
	smimg= congrid(image3d, 3, newx, newy)
	;smimg= congrid(image3d, newx, newy)

	; ct
	newimage= bytarr(3,1024,768)
	newimage[*,*,*]= backgrndclr
	newimage[*,0:newx-1,(768-newy)/2:(768-newy)/2+newy-1]= smimg
	; bw
	;newimage= bytarr(1024,768)
	;newimage[*,*]= backgrndclr
	;newimage[0:newx-1,(768-newy)/2:(768-newy)/2+newy-1]= smimg

        ilbl= '0000'+strcompress(string(i),/remove_all)
        ilbl= strmid(ilbl,strlen(ilbl)-4,4)
	WRITE_JPEG, 'data1/pg/pg32_g/movie_images_zoom/img_stars_zoom_'+ilbl+'.jpg', newimage, TRUE=1, quality= 95
	;WRITE_JPEG, 'data1/pg/pg32_g/movie_images_zoom/img_stars_zoom_'+ilbl+'.jpg', newimage

        ilbl= '0000'+strcompress(string(i3+i2+i2+i1+i1+1-i),/remove_all)
        ilbl= strmid(ilbl,strlen(ilbl)-4,4)
	WRITE_JPEG, 'data1/pg/pg32_g/movie_images_zoom/img_stars_zoom_'+ilbl+'.jpg', newimage, TRUE=1, quality= 95
	;WRITE_JPEG, 'data1/pg/pg32_g/movie_images_zoom/img_stars_zoom_'+ilbl+'.jpg', newimage

endfor

;
; this image is 646 x 614
;
READ_JPEG, 'data1/pg/sdssbinaryqso.jpg', sdssbinaryqso, TRUE=1
qsoimg= congrid(sdssbinaryqso,3,512,487)

for i= i1+i2,i1+i2+i3 do begin
	
	; ct
        newimage= bytarr(3,1024,768)
	newimage[*,*,*]= backgrndclr
        newimage[*,0:newx-1,(newy-1)/2:(newy-1)/2 + newy-1]= smimg
	newimage[*,512:1023,140:626]= qsoimg
	; bw
        ;newimage= bytarr(1024,768)
	;newimage[*,*]= backgrndclr
        ;newimage[0:newx-1,(newy-1)/2:(newy-1)/2 + newy-1]= smimg
	;newimage[512:1023,140:626]= qsoimg(0,*,*)

        ilbl= '0000'+strcompress(string(i),/remove_all)
        ilbl= strmid(ilbl,strlen(ilbl)-4,4)
	;WRITE_JPEG, 'data1/pg/pg32_g/movie_images_zoom/img_stars_zoom_'+ilbl+'.jpg', newimage, TRUE=1, quality= 95
	WRITE_JPEG, 'temp1.jpg', newimage, TRUE=1, quality= 95
	;WRITE_JPEG, 'temp1.jpg', newimage

        ;cmd = "convert -font Courier-Bold -antialias -pointsize 30 -fill black -draw 'text 150,180 "
        cmd = "convert -font Courier-Bold -antialias -pointsize 30 -fill white -draw 'text 150,180 "
        cmd = cmd + """
        cmd = cmd + 'Simulation'
        cmd = cmd + """
        cmd = cmd + "'"
        cmd = cmd + ' -quality 95 '
        cmd = cmd + " temp1.jpg temp2.jpg"
        print, cmd
        spawn, cmd

        ;cmd = "convert -font Courier-Bold -antialias -pointsize 30 -fill black -draw 'text 620,120 "
        cmd = "convert -font Courier-Bold -antialias -pointsize 30 -fill white -draw 'text 620,120 "
        cmd = cmd + """
        cmd = cmd + 'SDSS J1254+0846'
        cmd = cmd + """
        cmd = cmd + "'"
        cmd = cmd + ' -quality 95 '
        cmd = cmd + " temp2.jpg data1/pg/pg32_g/movie_images_zoom/img_stars_zoom_"+ilbl+".jpg"
        print, cmd
        spawn, cmd



endfor


end






