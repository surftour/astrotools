;======================================================================
;
;
;   General scripts to make movie images.  You still will need to
;   compile all the images into a movie using either ffmpeg, QuickTimePro
;   or some equivalent program.
;
;
;
;
;   first script make_movie_2 does a two panel gas/stars comparison
;
;   second script make_movie_4 does a four panel gas/stars 
;   comparison from two viewing angles.
;
;
;
;======================================================================
pro make_movie_2, junk


;frun= '/raid4/tcox/localgroup/v5'
;frun= '/raid4/tcox/bs/b2e_igm3'
frun= 'data/altsf/vc3vc3e_N2'
;frun= '/raid4/tcox/z3/b6e'
;frun= '/raid4/tcox/minor/min_30'
;frun= '/raid4/tcox/vc3vc3e_2'
;frun= '/raid4/tcox/vc3vc3e_no'
;frun= '/raid4/tcox/vc3vc3i_rem4'
imagedir= frun+'/movieimages'
spawn, "mkdir "+imagedir

imagename= 'image'
;imagename= 'gs_image'
;imagename= 'gt_image'

;startsnap= 0
startsnap= 42
;endsnap= 10
endsnap= 44
;endsnap= 140
;endsnap= 201
;endsnap= 202
;endsnap= 300

;imagename= 'testimage'
;imagename= 'test2image'
;imagename= 'test3image'
;startsnap=0
;startsnap= 50
;startsnap= 56
;startsnap= 73
;startsnap= 187
;endsnap= 1
;endsnap= 74
;endsnap=57
;endsnap= 188



;xlen= 200.0
;xlen= 120.0
;xlen= 75.0
;xlen= 70.0
xlen= 40.0
;xlen= 5.0

pixels= 640L


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


	; get gas image
	; ---------------
	x= fload_gas_xyz('x')
	y= fload_gas_xyz('y')
	z= fload_gas_xyz('z')
	m= fload_gas_mass(1)
	hsml= fload_gas_hsml(1)

	set_maxden= 1.0e-2
	set_dynrng= 1.0e+6

	movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, pixels=pixels, center=center, $
			idtofollow=idtofollow, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			NxNImage=NxNImage, NxNTempImage=NxNTempImage

	GasPic= NxNImage
	TempPic= NxNTempImage



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

        ;movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, /crude, $
        movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, pixels=pixels, center=center,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        StellarPic= NxNImage


	; -------------------


	ilbl= '000'+strcompress(string(i),/remove_all)
	ilbl= strmid(ilbl,strlen(ilbl)-3,3)



	;   write jpg files
	; --------------------

	; gas only image
	;s = SIZE(GasPic)
	;image3d = BYTARR(3, s(1), s(2))
	;image3d(0, *, *) = blue0(GasPic)
	;image3d(1, *, *) = blue1(GasPic)
	;image3d(2, *, *) = blue2(GasPic)
	;WRITE_JPEG, imagedir+'/image2_'+ilbl+'.jpg', image3d, TRUE=1, quality= 95

	; gas and gas temp image
	;image3d = BYTARR(3, 2*pixels, pixels)
	;image3d(0, *, *)= [blue0(GasPic), red0(TempPic)]
	;image3d(1, *, *)= [blue1(GasPic), red1(TempPic)]
	;image3d(2, *, *)= [blue2(GasPic), red2(TempPic)]

	; gas and stars image
	image3d = BYTARR(3, 2*pixels, pixels)
	image3d(0, *, *)= [blue0(GasPic), c0(StellarPic)]
	image3d(1, *, *)= [blue1(GasPic), c1(StellarPic)]
	image3d(2, *, *)= [blue2(GasPic), c2(StellarPic)]


	; black line to separate the two panels
	image3d[*,pixels:pixels+1,*]= 0


	; add the time to the image
	; method #1 - bad
	;
	;timelbl=fload_timelbl(0.7,2,/noteq)
	;;im = pixfont(timelbl, xsize=dxx, ysize=dyy, align=1)
	;xyouts, 0.99, 0.1, timelbl, charthick=5.0, align=1, charsize=3.0
	;im=tvrd()
	;for ii=0,(dxx-1) do begin
	;  for jj=0,(dyy-1) do begin
	;	if im[ii,jj] gt 0 then begin
	;		;print, "ii, jj, old= ", ii, jj, image3d[*,ii,jj]
	;		;image3d[*,ii+20,pixels-dyy-15+jj]= 255     ; white numbers
	;		image3d[*,ii+20,pixels-dyy-15+jj]= 0       ; black numbers
	;	endif
	;  endfor
	;endfor
	;erase

        filename= imagedir+'/'+tempimagename+'_'+ilbl+'.jpg'


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


	; ----------------------------------
	; write a eps file
	; ----------------------------------
	;write_eps_file, image3d, imagedir+'/f_'+ilbl+'.eps', xlen




	; ----------------------------------
	; use the z buffer
	; ----------------------------------
	;write_zbuff_to_jpg, image3D, imagedir+'/z_'+ilbl+'.jpg', xlen, pixels





        ;   write png file
        ; --------------------
        ;DoublePic= [GasPic, StellarPic]
        ;DoublePic[pixels:pixels+1,*]= 0

        ;filename= imagedir+'/'+imagename+'_'+ilbl+'.png'
        ;write_png, filename, DoublePic, blue0, blue1, blue2



endfor


; -------------
;  Done
; -------------


end





;===============================================================================




pro make_movie_4, junk


;frun= '/raid4/tcox/localgroup/v5'
;frun= '/raid4/tcox/bs/b2e_igm3' & rotateby= 90.0 * (!PI / 180.)
;frun= '/raid4/tcox/minor/min_30' & rotateby= -60.0
;frun= '/raid4/tcox/minor/min_150' & rotateby= 150.0 * (!PI / 180.)
;frun= '/raid4/tcox/gtest' & rotateby= 90.0 * (!PI / 180.)
;frun= '/raid4/tcox/gurtina' & rotateby= 90.0 * (!PI / 180.)
;frun= '/raid4/tcox/z3/b6e'
;frun= 'data/z3/b4f'
;frun= "/n/home/tcox/data/kecktests/Sb2x1Sb2x1_e_32/"
;frun= "/n/home/tcox/data/kecktests/Sb10xSb10xhs_e_256/"
;frun= "/n/home/tcox/data1/lg_m33/Sbm33_1/"
;frun= "/n/home/tcox/data1/lg_m33/Sbm33_2/"
;frun= "/n/home/tcox/data1/lg_m33/Sbm33_3/"
;frun= "/n/home/tcox/data1/lg_m33/SbSbm33_1/"
;frun= "/n/home/tcox/data1/lg_m33/SbSbm33_2/"
;frun= "/n/home/tcox/data1/lg_m33/SbSbm33_3/"
;frun= "/n/home/tcox/data1/lg_m33/SbSbm33_4/"
;frun= "/n/home/tcox/data1/lg_m33/SbSbm33_5/"
;frun= "/n/home/tcox/data1/lg_m33/SbSbm33_6/"
;frun= "/n/home/tcox/data1/lg_m33/SbSbm33_7/"
;frun= "/n/home/tcox/data1/remergers/iso_b3e/"
;frun= "/n/home/tcox/data1/remergers/iso_sph/"
;frun= "/n/home/tcox/data1/remergers/iso_v4/"
;frun= "/n/home/tcox/data1/remergers/b3ev2/"
;frun= "/n/home/tcox/data1/remergers/b3ev2v2/"
;frun= "/n/home/tcox/data1/remergers/b3ev2sph/"
;frun= "/n/home/tcox/data1/remergers/b3ev3/"
;frun= "/n/home/tcox/data1/remergers/b3ev3v2/"
;frun= "/n/home/tcox/data1/remergers/b3ev3v3/"
;frun= "/n/home/tcox/data1/remergers/b3ev4/"
;frun= "/n/home/tcox/data1/remergers/b3esph/"
;frun= "/n/home/tcox/data1/remergers/b3esphv2/"
;frun= "/n/home/tcox/data1/remergers/b3esphv3/"
;frun= "/n/home/tcox/data1/remergers/b3eb3e2/"

if not keyword_set(frun) then frun= junk

;
; make directory for images
;
imagedir= frun+'/movie4'
spawn, "mkdir "+imagedir


imagename= 'image4'
rotateby= 0.0


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

startsnap= 0
;startsnap= 1
;startsnap= 72
;startsnap= 240
;startsnap= 500
;startsnap= 541
;startsnap= 551
;startsnap= 572
;startsnap= 599
;startsnap= 1495

;endsnap= 2
;endsnap= 5
;endsnap= 10
;endsnap= 55
;endsnap= 74
;endsnap= 140
;endsnap= 224
;endsnap= 266
;endsnap= 300
;endsnap= 302
;endsnap= 360
;endsnap= 400
;endsnap= 500
;endsnap= 545
;endsnap= 553
;endsnap= 620
;endsnap= 699
;endsnap= 700
;endsnap= 999
;endsnap= 1200
;endsnap= 1500
spawn, "/bin/ls "+frun+"/snapshot_* | wc ",result
endsnap=long(result[0])-1



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
pixels= 320L

center= [0.0, 0.0, 0.0]  ; get's reset later



;-------------------------------------


print, "frun= ", frun
print, "imagedir= ", imagedir
print, "imagename= ", imagename
print, "tempimagename= ", tempimagename
print, "startsnap= ", startsnap
print, "endsnap= ", endsnap
print, "xlen= ", xlen
print, "pixels= ", pixels
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

	set_maxden= 1.0e-0
	set_dynrng= 1.0e+7


	; top view
	movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
	;movie_makepic, x_new, y_new, z_new, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, pixels=pixels, center=center, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			NxNImage=NxNImage, NxNTempImage=NxNTempImage

	GasPicTop= NxNImage
	TempPicTop= NxNTempImage


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
			hsml=hsml, pixels=pixels, center=center_side, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			NxNImage=NxNImage, NxNTempImage=NxNTempImage

	GasPicSide= NxNImage
	TempPicSide= NxNTempImage



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
	set_dynrng= 1.0e+6


	; top view
        movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
        ;movie_makepic, x_new, y_new, z_new, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, pixels=pixels, center=center,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        StellarPicTop= NxNImage


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
                        hsml=hsml, pixels=pixels, center=center_side,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        StellarPicSide= NxNImage



	; -------------------


	ilbl= '0000'+strcompress(string(i),/remove_all)
	ilbl= strmid(ilbl,strlen(ilbl)-4,4)



	;   write jpg files
	; --------------------

	image3d = BYTARR(3, 2*pixels, 2*pixels)

	image3d(0, *, 0:pixels-1)= [blue0(GasPicSide), c0(StellarPicSide)]
	image3d(1, *, 0:pixels-1)= [blue1(GasPicSide), c1(StellarPicSide)]
	image3d(2, *, 0:pixels-1)= [blue2(GasPicSide), c2(StellarPicSide)]

	image3d(0, *, pixels:2*pixels-1)= [blue0(GasPicTop), c0(StellarPicTop)]
	image3d(1, *, pixels:2*pixels-1)= [blue1(GasPicTop), c1(StellarPicTop)]
	image3d(2, *, pixels:2*pixels-1)= [blue2(GasPicTop), c2(StellarPicTop)]


	; black line to separate the different panels
	image3d[*,pixels:pixels+1,*]= 0    ; verticle line
	image3d[*,*,pixels:pixels+1]= 0    ; horizontal line


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




pro write_eps_file, image, fname, xlen



set_plot, 'ps'
device, filename= fname, /encapsulated,/color,bits_per_pixel=8
device, SET_CHARACTER_SIZE=[200,300], xsize=24, ysize=12



x0= 0.01
y0= 0.01

xhalf= 0.50

x1=0.99
y1=0.99

tv, image, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal, true=1
;tv, image, true=1


xmin= -xlen
xmax=  xlen
ymin= -xlen
ymax=  xlen


!p.ticklen=0.01
!p.position=[x0,y0,xhalf,y1]
plot,[0],[0], psym=3, xstyle=1, ystyle=1, color= 0, /noerase, $
                      xrange=[xmin,xmax], yrange=[ymin,ymax], $
                      xcharsize=0.01, ycharsize=0.01, $
                      xthick=4.0, ythick=4.0, /normal, /nodata, $
                      xtickformat='(a1)', ytickformat='(a1)'

!p.position=[xhalf,y0,x1,y1]
plot,[0],[0], psym=3, xstyle=1, ystyle=1, color= 0, /noerase, $
                      xrange=[xmin,xmax], yrange=[ymin,ymax], $
                      xcharsize=0.01, ycharsize=0.01, $
                      xthick=4.0, ythick=4.0, /normal, /nodata, $
                      xtickformat='(a1)', ytickformat='(a1)'

xyouts, 0.06, 0.92, fload_timelbl(1,2,/noteq), /normal, size= 1.3, charthick=3.0, color= 0



device, /close



end






;===============================================================================





pro write_zbuff_to_jpg, image, fname, xlen, pixels

        ; ----------------------------------
        ; use the z buffer
        ; ----------------------------------
        ; this is a nice way to combine just the pixel image
        ; and words or these types of things
        ;

	set_plot, 'z'
	;device, set_resolution= [2*pixels,pixels], set_character_size= [12,18]
	device, set_resolution= [2*pixels,pixels]
	erase


	x0= 0.01
	y0= 0.01

	xhalf= 0.50

	x1=0.99
	y1=0.99

	tv, image, true=1


	xmin= -xlen
	xmax=  xlen
	ymin= -xlen
	ymax=  xlen



	;!p.ticklen=0.01
	;!p.position=[x0,y0,xhalf,y1]
	;plot,[0],[0], psym=3, xstyle=1, ystyle=1, color= 0, /noerase, $
        ;              xrange=[xmin,xmax], yrange=[ymin,ymax], $
        ;              xcharsize=0.01, ycharsize=0.01, $
        ;              xthick=4.0, ythick=4.0, /normal, /nodata, $
        ;              xtickformat='(a1)', ytickformat='(a1)'

	;!p.position=[xhalf,y0,x1,y1]
	;plot,[0],[0], psym=3, xstyle=1, ystyle=1, color= 0, /noerase, $
        ;              xrange=[xmin,xmax], yrange=[ymin,ymax], $
        ;              xcharsize=0.01, ycharsize=0.01, $
        ;              xthick=4.0, ythick=4.0, /normal, /nodata, $
        ;              xtickformat='(a1)', ytickformat='(a1)'

	;xyouts, 0.06, 0.92, fload_timelbl(1,2,/noteq), /normal, size= 1.3, charthick=3.0, color= 0


        WRITE_JPEG, fname, tvrd(true=1), true=1, quality=100

end






;===============================================================================






pro movie_process_rotation, x_orig, y_orig, z_orig, rotate_theta, rotate_phi, x_new, y_new, z_new

        print, "*** process_rotation ***"
        print, "rotate: by (theta,phi)=", rotate_theta,rotate_phi

        ; transform to radians
        theta= !PI*rotate_theta/180.0
        phi= !PI*rotate_phi/180.0

        x= x_orig
        y= y_orig
        z= z_orig

        ; rotate the coordinates
        ; ------------------------
        ; first around z axis (phi)
        ; then around y axis (theta
        x_new= x*(cos(theta)*cos(phi)) + y*(cos(theta)*sin(phi)) + z*sin(theta)
        y_new= -x*sin(phi) + y*cos(phi)
        z_new= -x*(sin(theta)*cos(phi)) - y*(sin(theta)*sin(phi)) + z*cos(theta)


end

                                  
