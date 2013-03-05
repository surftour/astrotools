;======================================================================
;
;
;   new script based upon movie_makeit, but we instead want to
;   zoom into one of the centers.
;
;
;
;======================================================================
pro make_movie_2, junk


frun= '/raid4/tcox/vc3vc3e_2'
;frun= '/raid4/tcox/vc3vc3e_no'
imagedir= frun+'/movieimages'

imagename= 'image'
;startsnap= 0
startsnap= 36
;endsnap= 2
endsnap= 37
;endsnap= 107

xlen= 5.0

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

; grab one galaxy's center

; warning, need to load time_centers for this to work
read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers.txt"
cen_x= cp_cen2[0,*]
cen_y= cp_cen2[1,*]
cen_z= cp_cen2[2,*]






;--------------------------------------
;  Make images
;--------------------------------------


for i=startsnap, endsnap do begin


	ok=fload_snapshot_bh(frun,i,/skip_center)

	thiscenter= [cen_x[i], cen_y[i], cen_z[i]]

	; get gas image
	; ---------------
	x= fload_gas_xyz('x', center=thiscenter)
	y= fload_gas_xyz('y', center=thiscenter)
	z= fload_gas_xyz('z', center=thiscenter)
	m= fload_gas_mass(1)
	hsml= fload_gas_hsml(1)

	set_maxden= 1.0e-2
	set_dynrng= 1.0e+5

	movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, pixels=pixels, center=center, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			NxNImage=NxNImage, NxNTempImage=NxNTempImage

	GasPic= NxNImage
	TempPic= NxNTempImage



        ; get allstar image
        ; --------------------
        x= fload_allstars_xyz('x', center=thiscenter)
        y= fload_allstars_xyz('y', center=thiscenter)
        z= fload_allstars_xyz('z', center=thiscenter)
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


frun= '/raid4/tcox/vc3vc3e_2'
rotateby= 30.0 * (!PI / 180.)
imagedir= frun+'/movieimages'

imagename= 'image'
startsnap= 0
endsnap= 2
;endsnap= 107

xlen= 5.0 


pixels= 480L


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

; grab one galaxy's center

; warning, need to load time_centers for this to work
read_centerpositions, cp_time, cp_cen1, cp_cen2, filename=frun+"/centers.txt"
cen_x= cp_cen2[0,*]
cen_y= cp_cen2[1,*]
cen_z= cp_cen2[2,*]







;--------------------------------------
;  Make images
;--------------------------------------


for i=startsnap, endsnap do begin


	ok=fload_snapshot_bh(frun,i,/skip_center)

	thiscenter= [cen_x[i], cen_y[i], cen_z[i]]


	calc_j= 0
	;calc_j= 1
	if calc_j eq 1 then begin
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



	; get gas image
	; ---------------
	x= fload_gas_xyz('x', center= thiscenter)
	y= fload_gas_xyz('y', center= thiscenter)
	z= fload_gas_xyz('z', center= thiscenter)
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
	set_dynrng= 1.0e+6


	; top view
	;movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
	movie_makepic, x_new, y_new, z_new, m, xlen, xz=xz, yz=yz, $
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
	movie_makepic, z_new, x_new, -y_new, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, pixels=pixels, center=center, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			NxNImage=NxNImage, NxNTempImage=NxNTempImage

	GasPicSide= NxNImage
	TempPicSide= NxNTempImage



        ; get allstar image
        ; --------------------
        x= fload_allstars_xyz('x', center= thiscenter)
        y= fload_allstars_xyz('y', center= thiscenter)
        z= fload_allstars_xyz('z', center= thiscenter)
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
        ;movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
        movie_makepic, x_new, y_new, z_new, m, xlen, xz=xz, yz=yz, $
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
	movie_makepic, z_new, x_new, -y_new, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, pixels=pixels, center=center,$
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage

        StellarPicSide= NxNImage



	; -------------------


	ilbl= '000'+strcompress(string(i),/remove_all)
	ilbl= strmid(ilbl,strlen(ilbl)-3,3)



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


        filename= imagedir+'/'+tempimagename+'_'+ilbl+'.jpg'
	WRITE_JPEG, filename, image3d, TRUE=1, quality= 95



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
	cmd = cmd + filename + " " + fout

	print, cmd
	spawn, cmd



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

                                  
