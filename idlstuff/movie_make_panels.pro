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



pro make_movie_panels, frun, thumbnail=thumbnail, $
			xlen=xlen

if not keyword_set(frun) then begin
	print, " "
	print, " make_movie, frun"
	print, " "
	return
endif

;
; make directory for images
;
imagedir= frun+'/movie_images'
spawn, "mkdir "+imagedir



;
; get snapshot range
;
startsnap= 0
;startsnap= 40
;startsnap= 120
;startsnap= 750
;startsnap= 815

spawn, "/bin/ls "+frun+"/snapshot_* | wc ",result
endsnap=long(result[0])-1
;endsnap=5
;endsnap=125
;endsnap=800
;endsnap=820
;endsnap=1250
;endsnap=3400


;
; other parameter choices
;

rotateby= 0.0

;
; side = 2 * xlen (in kpc/h)
;
if not keyword_set(xlen) then begin
;xlen= 350.0
;xlen= 280.0
;xlen= 210.0
;xlen= 140.0
;xlen= 105.0
;xlen= 75.0
;xlen= 70.0
;xlen= 52.5
xlen= 42.0
;xlen= 40.0
;xlen= 35.0
;xlen= 20.0
endif

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


if keyword_set(thumbnail) then begin
	xpixels= 256L
	ypixels= 192L
endif

center= [0.0, 0.0, 0.0]  ; gets reset later

; get center from file
;read_center, cfile_time, cfile_center, filename=frun+"/centers_den_550001.txt"


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
	center= [0.0, 0.0, 0.0]  ; gets reset later

	; manual
	;if i eq 1 then center= [-100.0, 0.0, 0.0]
	;if i eq 2 then center= [-132.0, 0.0, 0.0]

	; find angular momentum vector
	do_angular_j= 0
	;do_angular_j= 1
	if do_angular_j eq 1 then begin
	 if i eq startsnap then begin
	   com=fload_center(1)
	   r= fload_allstars_xyz('r', center=com)
	   comvel= fload_all_comvel(1, center=com, justcenter=10.0)
	   jx= fload_allstars_j(21, center=com, vcom=comvel)   ; specific j_x
	   jy= fload_allstars_j(22, center=com, vcom=comvel)   ; specific j_y
	   jz= fload_allstars_j(23, center=com, vcom=comvel)   ; specific j_z

	   ; angular momentum within some radius 
	   radius= 5.0
	   idx= where(r lt radius)
	   if idx(0) eq -1 then stop
	   jtot_x= total(jx(idx))
	   jtot_y= total(jy(idx))
	   jtot_z= total(jz(idx))
	   jtot= sqrt(jtot_x*jtot_x + jtot_y*jtot_y + jtot_z*jtot_z)
	   nx= jtot_x/jtot
	   ny= jtot_y/jtot
	   nz= jtot_z/jtot

	   n_hat= [nx, ny, nz]
	   n_hat_tot= sqrt(nx*nx + ny*ny + nz*nz)

	   theta= acos(nz)*180./!PI
	   phi= atan(ny,nx)*180./!PI

	   print, "angular momentum vector= ", nx, ny, nz
	   print, "              n_hat_tot= ", n_hat_tot
	   print, "                  theta= ", theta
	   print, "                    phi= ", phi
	 endif
	endif



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
	;idx= where(cfile_time ge (time-1.0e-6))
	;if idx(0) ne -1 then center= cfile_center(*,idx(0))



	ngas= fload_npart(0)
	if ngas eq 0 then goto, dostars


dogas:

	; get gas image
	; ---------------
	x= fload_gas_xyz('x')
	y= fload_gas_xyz('y')
	z= fload_gas_xyz('z')
	;movie_process_rotation, x, y, z, theta, phi, x_new, y_new, z_new
	;x_new= x * n_hat[0]
	;y_new= y * n_hat[0]
	;z_new= z * n_hat[0]
	x_new= cos(rotateby)*x - sin(rotateby)*z
	y_new= y
	z_new= sin(rotateby)*x + cos(rotateby)*z
	m= fload_gas_mass(1)
	hsml= fload_gas_hsml(1)

	set_maxden= 1.0e+0
	set_dynrng= 1.0e+7


	; top view
	movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, x_new, y_new, z_new, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, xpixels=xpixels, ypixels=ypixels, center=center, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			NxNImage=NxNImage, NxNTempImage=NxNTempImage

	GasPicTop= transpose(NxNImage)
	TempPicTop= transpose(NxNTempImage)

        ;   write jpg files
        ; --------------------
	image3d= img_map_to_color(GasPicTop, "blue", xpixels, ypixels)
        WRITE_JPEG, imagedir+'/img_gas_top_'+ilbl+'.jpg', image3d, TRUE=1, quality= 95



	; side view
	;theta= rotateby * (!PI / 180.0)
	;x_new= x*cos(theta) - z*sin(theta)
	;y_new= x*sin(theta) + z*cos(theta)
	;z_new= -y
	;movie_makepic, x_new, z_new, -y_new, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, z_new, x_new, -y_new, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, y, z, x, m, xlen, xz=xz, yz=yz, $
	center_side= [center(2), center(0), center(1)]
	;movie_makepic, z, x, y, m, xlen, xz=xz, yz=yz, $
	movie_makepic, z, y, x, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, xpixels=xpixels, ypixels=ypixels, center=center_side, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			NxNImage=NxNImage, NxNTempImage=NxNTempImage

	GasPicSide= transpose(NxNImage)
	TempPicSide= transpose(NxNTempImage)

        ;   write jpg files
        ; --------------------
	image3d= img_map_to_color(GasPicSide, "blue", xpixels, ypixels)
        WRITE_JPEG, imagedir+'/img_gas_side_'+ilbl+'.jpg', image3d, TRUE=1, quality= 95


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

	set_maxden= 1.0e+0
	set_dynrng= 1.0e+5


	; top view
        movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
        ;movie_makepic, x_new, y_new, z_new, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, xpixels=xpixels, ypixels=ypixels, center=center,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        StellarPicTop= transpose(NxNImage)

        ;   write jpg files
        ; --------------------
	image3d= img_map_to_color(StellarPicTop, "std", xpixels, ypixels)
	;image3d= img_map_to_color(StellarPicTop, "colortable", xpixels, ypixels)
        WRITE_JPEG, imagedir+'/img_stars_top_'+ilbl+'.jpg', image3d, TRUE=1, quality= 95

	; side view
	;theta= rotateby * (!PI / 180.0)
	;x_new= x*cos(theta) - z*sin(theta)
	;y_new= x*sin(theta) + z*cos(theta)
	;z_new= -y
	;movie_makepic, x_new, z_new, -y_new, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, z_new, x_new, -y_new, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, x, z, -y, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, y, z, x, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, z, x, y, m, xlen, xz=xz, yz=yz, $    ; original
	movie_makepic, z, y, x, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, xpixels=xpixels, ypixels=ypixels, center=center_side,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        StellarPicSide= transpose(NxNImage)

        ;   write jpg files
        ; --------------------
	image3d= img_map_to_color(StellarPicSide, "std", xpixels, ypixels)
        WRITE_JPEG, imagedir+'/img_stars_side_'+ilbl+'.jpg', image3d, TRUE=1, quality= 95




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
openr,1,"/n/home03/tcox/idlstuff/image/colplane_flip2.dat"
readu,1,ColPlane
close,1



end





;==========================================================================
;==========================================================================
;==========================================================================
;==========================================================================




