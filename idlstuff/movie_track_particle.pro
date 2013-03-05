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
pro make_movie, junk


;frun= '/raid4/tcox/localgroup/v5'
;frun= '/raid4/tcox/bs/b2e_igm3'
frun= 'data/altsf/vc3vc3e_N2'
idtofollow= 160000L
;frun= '/raid4/tcox/z3/b6e'
;frun= '/raid4/tcox/minor/min_30'
;frun= '/raid4/tcox/vc3vc3e_2'
;frun= '/raid4/tcox/vc3vc3e_no'
;frun= '/raid4/tcox/vc3vc3i_rem4'
imagedir= frun+'/movieimages'

imagename= 'image'
;imagename= 'gs_image'
;imagename= 'gt_image'

startsnap= 0
;startsnap= 42
;endsnap= 2
;endsnap= 10
;endsnap= 44
;endsnap= 140
endsnap= 201
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


rotateby= 0.0

;xlen= 200.0
;xlen= 120.0
;xlen= 75.0
xlen= 70.0
;xlen= 40.0
;xlen= 5.0

pixels= 640L


tempimagename= 'temp_'+imagename


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






;--------------------------------------
;  Make images
;--------------------------------------


for i=startsnap, endsnap do begin


	ok=fload_snapshot_bh(frun,i,/nopot_in_snap,/skip_center)
	; skip_center effectively sets the center to [0,0,0]


	bhid= fload_blackhole_id(1)
        ;bhid= bhid[0]
        ;bhid= bhid[1]
        ;bhid= 200001L
        ;bhid= 280002L   ; used for z3/b4e
        ;bhid= 400002L   ; used for ds/vc3vc3e_2
        bhid= long(max(bhid))   ; uses minimum bhid - then you don't need to know the BH ID
        center_bh= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid)

        print, "Blackhole ID: ", bhid
        print, "Blackhole center: ", center_bh
	xpix= long(pixels * (0.5 + center_bh[0] /(2. * xlen)))
	ypix= long(pixels * (0.5 + center_bh[1] /(2. * xlen)))
	zpix= long(pixels * (0.5 + center_bh[2] /(2. * xlen)))
	print, " BH pix= ", xpix, ypix, zpix




        ; get allstar image
        ; --------------------
        x= fload_allstars_xyz('x')
        y= fload_allstars_xyz('y')
        z= fload_allstars_xyz('z')
        ;movie_process_rotation, x, y, z, theta, phi, x_new, y_new, z_new
        ;x_new= cos(rotateby)*x - sin(rotateby)*z
        ;y_new= y
        ;z_new= sin(rotateby)*x + cos(rotateby)*z
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



	;-----------------------------------------------

	if idtofollow  gt 0 then begin
		;
		; used id_earthr.pro to produce idlist
		;
		;idlist_fmfile= fload_id_list(frun+'/earth_idlist.txt')
		idlist_fmfile= fload_id_list('data/ds/vc3vc3h/earth_idlist.txt')

		print, "id's in file= ",n_elements(idlist_fmfile)
		sidx=sort(idlist_fmfile)
		idlist_fmfile= idlist_fmfile(sidx)


		; find x,y of earth stars
		gid= fload_allstars_ids(1)
		;
			;--------------------------------------
			; straight from "oplot_contour_add_idlist"
			;--------------------------------------
                        idx= intarr(n_elements(gid))
                        lstid= -1
                        duplicates= 0
                        for ii=0L,n_elements(idlist_fmfile)-1 do begin
                            if idlist_fmfile[ii] eq lstid then duplicates= duplicates+1
                            inlst= where(gid eq idlist_fmfile[ii])
                            if inlst(0) ne -1 then begin
                                idx(inlst)= 1
                            endif else begin
                                ;print, "couldn't find id=",idlist_fmfile[i]
                            endelse
                            lstid= idlist_fmfile[ii]
                        endfor
                        midx= where(idx eq 1)
                        print, "matching gid's= ",n_elements(midx)
                        print, "duplicate id's= ",duplicates
                        ;if midx(0) ne -1 then oplot, x(midx), y(midx), psym=7, color= 150, symsize= 0.2
                        ;if idx(0) ne -1 then oplot, x(midx), y(midx), psym=3, color= 150
			;--------------------------------------

		; now, ex(midx), etc. are the x,y,z of earths
		;
		for ii=0, n_elements(midx)-1 do begin
			xpix= long(pixels * (0.5 + x(midx[ii]) /(2. * xlen)))
			ypix= long(pixels * (0.5 + y(midx[ii]) /(2. * xlen)))
			zpix= long(pixels * (0.5 + z(midx[ii]) /(2. * xlen)))
			
			; remember that some of the xyz stuff is flipped
			if pixok(xpix,ypix,pixels) then StellarPicTop[ypix,xpix]=0
			if pixok(xpix+1,ypix,pixels) then StellarPicTop[ypix,xpix+1]=0
			if pixok(xpix+2,ypix,pixels) then StellarPicTop[ypix,xpix+2]=0
			if pixok(xpix-1,ypix,pixels) then StellarPicTop[ypix,xpix-1]=0
			if pixok(xpix-2,ypix,pixels) then StellarPicTop[ypix,xpix-2]=0
			if pixok(xpix,ypix+1,pixels) then StellarPicTop[ypix+2,xpix]=0
			if pixok(xpix,ypix+2,pixels) then StellarPicTop[ypix+1,xpix]=0
			if pixok(xpix,ypix-1,pixels) then StellarPicTop[ypix-1,xpix]=0
			if pixok(xpix,ypix-2,pixels) then StellarPicTop[ypix-2,xpix]=0

			; remember that some of the xyz stuff is flipped
			if pixok(zpix,xpix,pixels) then StellarPicSide[xpix,zpix]=0
			if pixok(zpix+1,xpix,pixels) then StellarPicSide[xpix,zpix+1]=0
			if pixok(zpix+2,xpix,pixels) then StellarPicSide[xpix,zpix+2]=0
			if pixok(zpix-1,xpix,pixels) then StellarPicSide[xpix,zpix-1]=0
			if pixok(zpix-2,xpix,pixels) then StellarPicSide[xpix,zpix-2]=0
			if pixok(zpix,xpix+1,pixels) then StellarPicSide[xpix+2,zpix]=0
			if pixok(zpix,xpix+2,pixels) then StellarPicSide[xpix+1,zpix]=0
			if pixok(zpix,xpix-1,pixels) then StellarPicSide[xpix-1,zpix]=0
			if pixok(zpix,xpix-2,pixels) then StellarPicSide[xpix-2,zpix]=0

		endfor

	endif


	; -------------------


	ilbl= '000'+strcompress(string(i),/remove_all)
	ilbl= strmid(ilbl,strlen(ilbl)-3,3)



	;   write jpg files
	; --------------------

	; gas and stars image
	image3d = BYTARR(3, 2*pixels, pixels)
	image3d(0, *, *)= [c0(StellarPicTop), c0(StellarPicSide)]
	image3d(1, *, *)= [c1(StellarPicTop), c1(StellarPicSide)]
	image3d(2, *, *)= [c2(StellarPicTop), c2(StellarPicSide)]


	; black line to separate the two panels
	image3d[*,pixels:pixels+1,*]= 0


        ;filename= imagedir+'/'+tempimagename+'_'+ilbl+'.jpg'
        filename= imagedir+'/temp.jpg'


	WRITE_JPEG, filename, image3d, TRUE=1, quality= 95



	; add text onto the figure
	;
	;convert -font Courier-Bold -antialias -pointsize 30 -fill white \                                               
        ;     -draw 'text 20,50 "z = 2.0"'  \                                                                 
        ;     mypic.jpg  annotated_pic.jpg        
	;

	fout= imagedir+'/'+imagename+'_'+ilbl+'.jpg'
	timelbl=fload_timelbl(0.7,2,/noteq)

	cmd = 'convert -font Courier-Bold -antialias -pointsize 30 -fill black -draw '
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






function pixok, xp, yp, pixnum


	if xp lt 0 then return, 0
	if yp lt 0 then return, 0
	if xp gt pixnum then return, 0
	if yp gt pixnum then return, 0

	return, 1

end



