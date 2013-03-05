;======================================================================
;
;
;   General scripts to make movie images.  You still will need to
;   compile all the images into a movie using either ffmpeg, QuickTimePro
;   or some equivalent program.
;
;
;
;======================================================================



pro make_movie_4, frun

if not keyword_set(frun) then begin
	print, " "
	print, " make_movie, frun"
	print, " "
	return
endif

;
; make directory for images
;
imagedir= frun+'/movie4'
spawn, "mkdir "+imagedir



;
; which temp file should we use?
; 
tempi= 0
repeat begin
 tempilbl= strcompress(string(tempi),/remove_all)
 tempimagename= "temp"+tempilbl+".jpg"
 openr, 1, tempimagename, error=err
 close, 1
 print, "tempi= ", tempi
 print, "tempilbl= ", tempilbl
 print, "tempimagename= ", tempimagename
 print, "err= ", err
 print, "--------------"
 tempi= tempi + 1
endrep until (err ne 0)



;
; get snapshot range
;
startsnap= 0

spawn, "/bin/ls "+frun+"/snapshot_* | wc ",result
endsnap=long(result[0])-1



;
; other parameter choices
;

imagename= 'image4'
rotateby= 0.0

;xlen= 350.0
;xlen= 280.0
;xlen= 210.0
xlen= 140.0
;xlen= 75.0
;xlen= 70.0
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


xpixels= 512L
ypixels= 384L

center= [0.0, 0.0, 0.0]  ; get's reset later



;-------------------------------------


print, "frun= ", frun
print, "imagedir= ", imagedir
print, "imagename= ", imagename
print, "tempimagename= ", tempimagename
print, "startsnap= ", startsnap
print, "endsnap= ", endsnap
print, "xlen= ", xlen
print, "xpixels= ", xpixels
print, "ypixels= ", ypixels
print, "center= ", center


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






;--------------------------------------
;  Make images
;--------------------------------------


for i=startsnap, endsnap do begin


	if i ge 1000 then begin
		ok=fload_snapshot_bh(frun,i,/nopot_in_snap,/skip_center,/do_four)
	endif else begin
		ok=fload_snapshot_bh(frun,i,/nopot_in_snap,/skip_center)
	endelse


	; center the image
	; -----------------

	; default
	center= [0.0, 0.0, 0.0]  ; get's reset later

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

	set_maxden= 1.0e+1
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


	; side view
	;theta= rotateby * (!PI / 180.0)
	;x_new= x*cos(theta) - z*sin(theta)
	;y_new= x*sin(theta) + z*cos(theta)
	;z_new= -y
	;movie_makepic, x_new, z_new, -y_new, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, z_new, x_new, -y_new, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, y, z, x, m, xlen, xz=xz, yz=yz, $
	center_side= [center(2), center(0), center(1)]
	movie_makepic, z, x, y, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, xpixels=xpixels, ypixels=ypixels, center=center_side, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			NxNImage=NxNImage, NxNTempImage=NxNTempImage

	GasPicSide= transpose(NxNImage)
	TempPicSide= transpose(NxNTempImage)



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


	; side view
	;theta= rotateby * (!PI / 180.0)
	;x_new= x*cos(theta) - z*sin(theta)
	;y_new= x*sin(theta) + z*cos(theta)
	;z_new= -y
	;movie_makepic, x_new, z_new, -y_new, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, z_new, x_new, -y_new, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, x, z, -y, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, y, z, x, m, xlen, xz=xz, yz=yz, $
	movie_makepic, z, x, y, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, xpixels=xpixels, ypixels=ypixels, center=center_side,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        StellarPicSide= transpose(NxNImage)



	; -------------------


	ilbl= '0000'+strcompress(string(i),/remove_all)
	ilbl= strmid(ilbl,strlen(ilbl)-4,4)



	;   write jpg files
	; --------------------

	image3d = BYTARR(3, 2*xpixels, 2*ypixels)

	image3d(0, *, 0:ypixels-1)= [blue0(GasPicSide), c0(StellarPicSide)]
	image3d(1, *, 0:ypixels-1)= [blue1(GasPicSide), c1(StellarPicSide)]
	image3d(2, *, 0:ypixels-1)= [blue2(GasPicSide), c2(StellarPicSide)]

	image3d(0, *, ypixels:2*ypixels-1)= [blue0(GasPicTop), c0(StellarPicTop)]
	image3d(1, *, ypixels:2*ypixels-1)= [blue1(GasPicTop), c1(StellarPicTop)]
	image3d(2, *, ypixels:2*ypixels-1)= [blue2(GasPicTop), c2(StellarPicTop)]


	; black line to separate the different panels
	image3d[*,xpixels:xpixels+1,*]= 0    ; verticle line
	image3d[*,*,ypixels:ypixels+1]= 0    ; horizontal line


        ;tempimagename= imagedir+'/'+tempimagename+'_'+ilbl+'.jpg'
        ;tempimagename= 'temp.jpg'
        ;tempimagename= 'temp1.jpg'
        ;tempimagename= 'temp2.jpg'
	WRITE_JPEG, tempimagename, image3d, TRUE=1, quality= 95



	; add text onto the figure
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
	cmd = cmd + tempimagename + " " + fout

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



