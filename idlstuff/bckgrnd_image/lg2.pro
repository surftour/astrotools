
Object_File='/home2/tcox/Tools/C-Routines_for_IDL/ComputeHsmlAndProject/HsmlAndProject.so'

;mw_bhid= 500001L
mw_bhid= 817955L

frun = "/raid4/tcox/localgroup/bhires/" 
gasfrac_max= 1.0
gasfrac_min= 0.0

Snap=   "snapshot"

pic_base= "lg_"

;ehist2, ' ', snapnum=150, filename=frun+'/earthist_150.eps', msg='Today'
;ehist2, ' ', snapnum=154, filename=frun+'/earthist_154.eps', msg='+0.15 Gyr'
;ehist2, ' ', snapnum=197, filename=frun+'/earthist_197.eps', msg='+1.8 Gyr, first passage'
;ehist2, ' ', snapnum=209, filename=frun+'/earthist_209.eps', msg='+2.2 Gyr, post-first passage'
;ehist2, ' ', snapnum=279, filename=frun+'/earthist_279.eps', msg='+3.5 Gyr, second passage'
;ehist2, ' ', snapnum=291, filename=frun+'/earthist_291.eps', msg='+4.0 Gyr, post-second passage'
;ehist2, ' ', snapnum=317, filename=frun+'/earthist_317.eps', msg='+4.5 Gyr, post-merger'

;snaparr= [150, 154, 197, 209, 279, 291, 317]
snaparr= [154, 197, 209, 279, 291, 317]
;snaparr= [209]
n_pics= n_elements(snaparr)

PIX = 1024
Len = 80.0

bw= 200.0


xpixels= long(PIX)              ; pixels in x-direction
ypixels= long(PIX)              ; pixels in y-direction
;window, xsize=XPixels,ysize=Ypixels,retain=2
seed=42L




;
;
;
;
;
;---------------------------------------------------------
for rep = 0,n_pics-1 do begin
   
   
   ok= fload_snapshot_bh(frun,snaparr[rep])


   Mgas =total(fload_gas_mass(1))
   Mstars= total(fload_allstars_mass(1))
   
   ; get blackhole id's
   bhid= fload_blackhole_id(1)

   mwcenter= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=mw_bhid)
   print, "using mw bh as center: ", mwcenter


;;;;;;;;;;;;; start  X-Y- projection of density

   x= fload_allstars_xyz('x',center=mwcenter)
   y= fload_allstars_xyz('y',center=mwcenter)
   z= fload_allstars_xyz('z',center=mwcenter)
   N=   long(n_elements(x))

   ; where is the earth
   ;star_to_follow= 484400L      ; for v7
   ;earth_pos= fload_allstars_xyz('dummy',center=mwcenter,idtofollow=star_to_follow)
   ;print, "earth xyz= ", earth_pos
   ;r= sqrt(earth_pos[0]*earth_pos[0] + earth_pos[1]*earth_pos[1])
   ;theta= acos(earth_pos[2]/r)
   ;phi= atan(earth_pos[1], earth_pos[0]) + !PI   ; (now -x is 0 deg, -y is 90, +x is 180 deg)
   ;print, "earth, r, theta, phi= ", r, theta*180./!PI, phi*180./!PI


   Xyz=  fltarr(3,N)            ; generate coordinate array for particles
   ;Xyz(0,*)= x
   ;Xyz(1,*)= y
   ;Xyz(2,*)= z
   Xyz(0,*)= y
   Xyz(1,*)= z
   Xyz(2,*)= x


   Hsml = fltarr(N)    
   ;Hsml = fload_allstars_hsml(1)
   ;Hsml(*) = 0.2
   Hsml(*) = 0.01

   Weight = fltarr(N)           ; generate space for mass array
   Weight(0:*)= fload_allstars_mass(1)   ; fill in mass values
   quantity = fltarr(N)       
   quantity(*)= 1.0             ; fill in mass values
   Axis1=   long(0)             ; select horizontal axis
   Axis2=   long(1)             ; select vertical axis
   Axis3=   long(2)             ; select projection axis
   
   xmin= float(-len)            ; left edge of map
   xmax= float( len)            ; right edge of map 
   ymin= float(-len)            ; lower
   ymax= float( len)            ; upper edge
   Zmin= float(-1.0e6)          ; lower range of particles that are included
   Zmax= float(1.0e6)           ; upper range


   xpixels= long(PIX)           ; pixels in x-direction
   ypixels= long(PIX)           ; pixels in y-direction

   Hmax=    float( (xmax-xmin)/xpixels * 50) ; maximum used smoothing length
   BoxSize  = double(0.0)

   DesNgb = 32L                 ; for those particles that have no zero smoothing length

   ResultW = fltarr(YPixels, XPixels) ; projected weights (normally projected mass map)
   ResultQ = fltarr(YPixels, XPixels) ; projected `Weight' weighted quantity



   gasfrac =  mgas/(mgas+mstars)
   print, "gas fraction=", gasfrac
   gasfrac =  0.8
   print, "gas fraction=", gasfrac



   ; will calculate hsml and project
   ;
   ;S = CALL_EXTERNAL(Object_File, $
   ;                  'findHsmlAndProject', $
   ;                  N, $
   ;                  Xyz, $
   ;                  Hsml, $
   ;                  Weight, $
   ;                  Quantity, $
   ;                  Xmin, $
   ;                  Xmax, $
   ;                  Ymin, $
   ;                  Ymax, $
   ;                  Zmin, $
   ;                  Zmax, $
   ;                  XPixels, $
   ;                  YPixels, $
   ;                  DesNgb, $
   ;                  Axis1, $
   ;                  Axis2, $
   ;                  Axis3, $
   ;                  Hmax,  $
   ;                  BoxSize, $
   ;                  ResultW, $
   ;                  ResultQ)



   ; good if already know Hsml
   ;-------------------------------

   ;bw= 10.0

   pos2D= fltarr(2,N)
   pos2D(0,*)= Xyz(0,*) * bw/(bw-Xyz(2,*))
   pos2D(1,*)= Xyz(1,*) * bw/(bw-Xyz(2,*))

   Hsml= Hsml * bw/(bw-Xyz(2,*))
   Weight= Weight * bw/(bw-Xyz(2,*))

   ind =where(Hsml gt 20.0)
   if ind(0) ne -1 then begin
       ;Weight(ind) = 0
       ;Hsml(ind) = 0.01
   endif

   ind =where(Xyz(2,*) ge bw)
   if ind(0) ne -1 then begin
       Weight(ind) = 0
       Hsml(ind) = 0.001
   endif

   quantity= fltarr(N)
   quantity(*)= 1.0
   ;if fload_npart(0) gt 0 then quantity(*)= fload_gas_temperature(1) else quantity(*)= Masses


   S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so', $
                     'slice', $
                     N, $
                     pos2D, $
   		     Hsml, $
   		     Weight, $
   		     Quantity, $
                     Xmin, $
   		     Xmax, $
   		     Ymin, $
   		     Ymax, $
                     Xpixels, $
   		     Ypixels, $
                     ResultW, $
   		     ResultQ)

   ;-------------------------------


   ResultW = transpose(ResultW)
   ResultW /= total(ResultW)


   print, "ResultW  max/min= ", max(ResultW), min(ResultW)


   R = 100

   xs= xpixels
   ys= ypixels


   Map = ResultW

;;; Now need to map the result onto a color table


   ma = max(Map)
   mi= ma /10000000
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

   imageTemp = fltarr(XPixels, ypixels)
   Imagetemp(*) = (gasfrac_max-gasfrac)/(gasfrac_max-gasfrac_min) * 180 


   imageTemp = fix(ImageTemp)



   ;
   ; OPTION 1
   ; ------------
   ;do_volker= 0
   do_volker= 1
   if do_volker eq 1 then begin
	cols=256
	ColPlane=fltarr(cols,cols,3)

	openr,1,"colplane2.dat"
	readu,1,ColPlane
	close,1


	JpegPic=bytarr(XPixels,YPixels,3)

	for i=0, XPixels-1 do begin
	      for j=0, YPixels-1 do begin
	         JPegPic(i,j,0)= ColPlane(imageDens(i,j),imageTemp(i,j),0)
	         JPegPic(i,j,1)= ColPlane(imageDens(i,j),imageTemp(i,j),1)
	         JPegPic(i,j,2)= ColPlane(imageDens(i,j),imageTemp(i,j),2)
	      endfor
	endfor
   endif



   ;
   ; OPTION 1
   ; ------------
   do_other= 0
   ;do_other= 1
   if do_other eq 1 then begin
	LOADCT, 1
	v1=[0,255]
	v2=[0,255]
	v3=[0,255]
	tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white (respectively)
	TVLCT, blue0, blue1, blue2, /GET

	JpegPic=bytarr(XPixels,YPixels,3)

	for i=0, XPixels-1 do begin
	      for j=0, YPixels-1 do begin
	         JPegPic(i,j,0)= blue0(imageDens(i,j))
	         JPegPic(i,j,1)= blue1(imageDens(i,j))
	         JPegPic(i,j,2)= blue2(imageDens(i,j))
	      endfor
	endfor
   endif




   ; load background
   read_jpeg, "/home2/tcox/Images/Astro/fields_cosmos.jpg", background_image_full
   bkgrnd_image= transpose(background_image_full[*, 1024:2047, 0:1023])

   ; v0
   ;JpegPic += 0.6*bkgrnd_image
   ;idx= where(JpegPic gt max(ColPlane))
   ;if idx(0) ne -1 then JpegPic(idx)= max(ColPlane)

   ; v1
   ;for i=0, XPixels-1 do begin
   ;   for j=0, YPixels-1 do begin
   ;      if bkgrnd_image(i,j,0) gt JPegPic(i,j,0)= ColPlane(imageDens(i,j),imageTemp(i,j),0)
   ;      JPegPic(i,j,1)= ColPlane(imageDens(i,j),imageTemp(i,j),1)
   ;      JPegPic(i,j,2)= ColPlane(imageDens(i,j),imageTemp(i,j),2)
   ;   endfor
   ;endfor


   ;tv, JPegPic, true=3,0

   exts='000'
   ;exts=exts+strcompress(string(rep),/remove_all)
   exts=exts+strcompress(string(snaparr[rep]),/remove_all)
   exts=strmid(exts,strlen(exts)-3,3)


   ;fout= "pic_hr_"+exts+".jpg"
   fout= pic_base+exts+".jpg"
   fout=strcompress(fout,/remove_all)

   write_jpeg, fout , Jpegpic, true=3, quality=90


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


endfor

end





 










 
