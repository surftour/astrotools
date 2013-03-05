pro img_makeimage, x, y, z, m, hsml, xlen, sendto, xz=xz, yz=yz, $
			filename=filename, thumbnail=thumbnail, fitstoo=fitstoo, $
			xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
			pixels=pixels, zthickness=zthickness, $
			crude=crude, center=center, msg=msg, plotpts=plotpts, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			pubstyle=pubstyle, particlesonly=particlesonly, $
			showbhs=showbhs, set_maxden=set_maxden, set_dynrng=set_dynrng, $
			track_to_draw=track_to_draw, nolabels=nolabels, $
			drawbox=drawbox, xslit=xslit, yslit=yslit


;
;  This is our general contour plotting procedure, note that this will
; make matters less confusing since we can simply plot everything from here now.
;
;
;
;
;
;
;
;
;


if not keyword_set(xlen) then xlen=100.0
if not keyword_set(x) then begin
   print, "  "
   print, "PROBLEM: contour_makeplot"
   print, "  "
   print, "  "
   print, "  "
   return
endif


if keyword_set(center) then begin
	orig_center= center
endif else begin
	orig_center=[0,0,0]
	center=[0,0,0]
endelse




;--------------------------------------
;  Process Raw Data - if needed
;--------------------------------------


; rotate
; ---------
if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin

        print, "rotate: moving to center", center
        x=x-center[0]
        y=y-center[1]
        z=z-center[2]

        print, "rotate: set center = [0,0,0]"
        center= [0,0,0]

        process_rotation, x, y, z, rotate_theta, rotate_phi, x_new, y_new, z_new

        x= x_new
        y= y_new
        z= z_new
endif





    N= long(n_elements(x))

    ; cut in z direction
    ; -------------------
    if keyword_set(zthickness) then begin
        idx= where((z LE center[2]+zthickness) and (z GE center[2]-zthickness))
        if idx(0) ne -1 then begin
                N= long(n_elements(idx))
                x= x(idx)
                y= y(idx)
                z= z(idx)
                m= m(idx)
                hsml= hsml(idx)
        endif
    endif


    ; cut in x direction
    ; -------------------
    if keyword_set(xthickness) then begin
        idx= where((x LE center[0]+xthickness) and (x GE center[0]-xthickness))
        if idx(0) ne -1 then begin
                N= long(n_elements(idx))
                x= x(idx)
                y= y(idx)
                z= z(idx)
                m= m(idx)
        endif
    endif

    ; cut in y direction
    ; -------------------
    if keyword_set(ythickness) then begin
        idx= where((y LE center[1]+ythickness) and (y GE center[1]-ythickness))
        if idx(0) ne -1 then begin
                N= long(n_elements(idx))
                x= x(idx)
                y= y(idx)
                z= z(idx)
                m= m(idx)
        endif
    endif




    if keyword_set(center) then begin

        print, "moving to center", center
        x=x-center[0]
        y=y-center[1]
        z=z-center[2]

        print, "set center = [0,0,0]"
        center= [0,0,0]

    endif





    ;--------------------------------------
    ;  Make image
    ;--------------------------------------


    img_makepic_raw, x, y, z, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, $
                        pixels=pixels, $
                        crude=crude, $
                        NxNImage=NxNImage, $
                        NxNTempImage=NxNTempImage





   ;------------------------------
   ;   Scale Image
   ;------------------------------


   ;  Option 1 
   ;--------------
   ;do_opt1= 0
   do_opt1= 1
   if do_opt1 eq 1 then begin
	;ResultW = transpose(NxNImage)
	ResultW = NxNImage
	ResultW /= total(ResultW)

	print, "ResultW  max/min= ", max(ResultW), min(ResultW)

	xs= pixels
	ys= pixels

	Map = ResultW

	;;; Now need to map the result onto a color table

	ma = max(Map)
	;mi= ma /10000000
	mi= ma /100000
	;mi= ma /10000
	;mi= ma /6000
	;mi= ma /4000
	;mi= ma /1000
	;mi= ma /600

	print, mi, ma

	;; now do clipping in Map

	ind=where(Map lt mi)
	if ind(0) ne -1 then Map(ind)=mi
	ind=where(Map gt ma)
	if ind(0) ne -1 then Map(ind)=ma

	;; now map the log of the field linearly onto the color table

	cols= 255
   
	ImageDens= byte((alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-2) +2)
   
	;   tv, imagedens

gasfrac_max= 1.0
gasfrac_min= 0.0

   Mgas =total(fload_gas_mass(1))
   Mstars= total(fload_allstars_mass(1))
   

   gasfrac =  mgas/(mgas+mstars)
   print, "gas fraction=", gasfrac
   gasfrac =  0.8
   print, "gas fraction=", gasfrac


	imageTemp = fltarr(pixels, pixels)
	Imagetemp(*) = (gasfrac_max-gasfrac)/(gasfrac_max-gasfrac_min) * 180

	imageTemp = fix(ImageTemp)


   endif



   ;------------------------------
   ;   Map to Color Table 
   ;------------------------------

   ;
   ; Volker's Pink/White/Blue Colors
   ; ----------------------------------
   ;do_volker= 0
   do_volker= 1
   if do_volker eq 1 then begin
        cols=256
        ColPlane=fltarr(cols,cols,3)

        ;openr,1,"colplane2.dat"
        openr,1,"/home2/tcox/Tools/idlstuff/image/colplane_flip2.dat"
        readu,1,ColPlane
        close,1


        JpegPic=bytarr(pixels,pixels,3)

        for i=0, pixels-1 do begin 
              for j=0, pixels-1 do begin
                 ;JPegPic(i,j,0)= ColPlane(imageDens(i,j),imageTemp(i,j),0)
                 ;JPegPic(i,j,1)= ColPlane(imageDens(i,j),imageTemp(i,j),1)
                 ;JPegPic(i,j,2)= ColPlane(imageDens(i,j),imageTemp(i,j),2)
                 JPegPic(i,j,0)= ColPlane(imageDens(i,j),imageDens(i,j),0)
                 JPegPic(i,j,1)= ColPlane(imageDens(i,j),imageDens(i,j),1)
                 JPegPic(i,j,2)= ColPlane(imageDens(i,j),imageDens(i,j),2)
              endfor
        endfor
   endif


   ;
   ; Black/White, or one of the IDL presets
   ; -----------------------------------------
   do_opt2= 0
   ;do_opt2= 1
   if do_opt2 eq 1 then begin

	;imageDens= 255-imageDens

        LOADCT, 0     ; black/white
        ;LOADCT, 1     ; blue/white
        ;LOADCT, 3     ; red temperature
        ;LOADCT, 4     ; rainbow colors
        v1=[0,255]
        v2=[0,255]
        v3=[0,255]
        tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
        TVLCT, c0, c1, c2, /GET

        JpegPic=bytarr(pixels,pixels,3)

        for i=0, pixels-1 do begin 
              for j=0, pixels-1 do begin
                 JPegPic(i,j,0)= c0(imageDens(i,j))
                 JPegPic(i,j,1)= c1(imageDens(i,j))
                 JPegPic(i,j,2)= c2(imageDens(i,j))
              endfor
        endfor
   endif


   ; -------------
   ;  Write Image
   ; -------------


   ; load background
   read_jpeg, "/home2/tcox/Images/Astro/fields_cosmos.jpg", background_image_full
   bkgrnd_image= transpose(background_image_full[*, 1024:2047, 0:1023])




   write_jpeg, filename , Jpegpic, true=3, quality=90




   ;time= fload_time(1)
   ;cmd = 'mogrify -font arial -antialias -pointsize 58 -fill white -draw '
   ;cmd = cmd + "'"
   ;cmd = cmd + 'text 40,72 '
   ;cmd = cmd + """
   ;cmd = cmd + 'T = ' + string(time/0.7,format='(F4.2)')+' Gyr'
   ;cmd = cmd + """
   ;cmd = cmd + "'"
   ;cmd = cmd + ' -quality 95 '
   ;cmd = cmd + fout
   ;print, cmd
   ;spawn, cmd







   ; -------------
   ;  Done
   ; -------------



end


