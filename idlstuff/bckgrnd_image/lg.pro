
Object_File='/home2/tcox/Tools/C-Routines_for_IDL/ComputeHsmlAndProject/HsmlAndProject.so'

frun = "/raid4/tcox/localgroup/bhires/" 
gasfrac_max= 1.0
gasfrac_min= 0.0

Snap=   "snapshot"

pic_base= "lg_"


snapnum= 209
n_pics= 1


PIX = 1024
Len = 55.0


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
   
   
   ok= fload_snapshot_bh(frun,snapnum)


   Mgas =total(fload_gas_mass(1))
   Mstars= total(fload_allstars_mass(1))
   
   ; get blackhole id's
   bhid= fload_blackhole_id(1)
   Nholes= n_elements(bhid)


   center= [0.0, 0.0, 0.0]
   if (rep ge 5) then center= fload_center_alreadycomp(1)




;;;;;;;;;;;;; start  X-Y- projection of density

   x= fload_allstars_xyz('x',center=center)
   y= fload_allstars_xyz('y',center=center)
   z= fload_allstars_xyz('z',center=center)
   N=   long(n_elements(x))

   Xyz=  fltarr(3,N)            ; generate coordinate array for particles
   Xyz(0,*)= x
   Xyz(1,*)= y
   Xyz(2,*)= z

   Hsml = fltarr(N)    

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



   S = CALL_EXTERNAL(Object_File, $
                     'findHsmlAndProject', $
                     N, $
                     Xyz, $
                     Hsml, $
                     Weight, $
                     Quantity, $
                     Xmin, $
                     Xmax, $
                     Ymin, $
                     Ymax, $
                     Zmin, $
                     Zmax, $
                     XPixels, $
                     YPixels, $
                     DesNgb, $
                     Axis1, $
                     Axis2, $
                     Axis3, $
                     Hmax,  $
                     BoxSize, $
                     ResultW, $
                     ResultQ)



   ; good if already know Hsml
   ;-------------------------------
   ;os2D= fltarr(2,N)
   ;os2D(0,*)= Xyz(0,*)
   ;os2D(1,*)= Xyz(1,*)

   ;quantity= fltarr(N)
   ;if fload_npart(0) gt 0 then quantity(*)= fload_gas_temperature(1) else quantity(*)= Masses


   ; = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so', $
   ;                  'slice', $
   ;                  N, $
   ;                  pos2D, $
   ;		      Hsml, $
   ;		      Weight, $
   ;		      Quantity, $
   ;                  Xmin, $
   ;		      Xmax, $
   ;		      Ymin, $
   ;		      Ymax, $
   ;                  Xpixels, $
   ;		      Ypixels, $
   ;                  ResultW, $
   ;		      ResultQ)

   ;-------------------------------


   ResultW = transpose(ResultW)
   ResultW /= total(ResultW)


   print, "ResultW  max/min= ", max(ResultW), min(ResultW)


   R = 100

   xs= xpixels
   ys= ypixels


   pic = fltarr(Xs, Ys)


   ;if (rep eq 14) or (rep eq 9) or (rep eq 13) then begin
   ;if (rep eq 4) or (rep eq 5) then begin
   if (rep eq 5) then begin

      for bhnr=0, Nholes-1 do begin

         N = 200

         bh_pos= fload_blackhole_xyz('xyz',center=center,idtofollow=bhid[bhnr])
         xc = (bh_pos[0]-xmin)/(xmax-xmin)*xs
         yc=  (bh_pos[1]-ymin)/(ymax-ymin)*ys

         print, "blackhole position: ", xc, yc


         Nspoke = 80


         spoke = randomn(seed, Nspoke)

         Radius= 5


         for row = 0, ys-1 do begin
            ;print, row
            for col = 0, Xs-1 do begin

               u =  double(col - xc) 
               v =  double(row - yc) 
               l = sqrt (u * u + v * v)/2

               t = (atan (u, v) / (2 * !PI) + .51) * Nspoke 
               i = floor (t) 
               t -= i 
               i = i mod Nspoke 

               w1 = spoke[i] * (1 - t) + spoke[(i + 1)  mod nspoke] * t 
               
               w1 = w1 * w1 

               if l lt radius then begin
                  w1 = max(spoke^2)*(1-  l/ (radius))  + l/ (radius) * w1
               endif

               w = 1.0 / (l + 2*radius)

               c =  w1 * w^3.8

               pic(col, row) += c
            endfor
         endfor
      endfor


      ;if rep eq 14 then begin
      ;if rep eq 4 then begin
      if rep eq -1 then begin

         Pic = pic /total(pic) * total(resultW)*2.0

         resultW *= 0.8

         imageTemp2 = ((+alog10(Pic)-alog10(max(resultw)/100) +3  )/4.0)*75.0

         ind=where(imageTemp2 lt 0) 
         if ind(0) ne -1 then ImageTemp2(ind)=0
         ind=where(ImageTemp2 gt 75.0)
         if ind(0) ne -1 then Imagetemp2(ind)=75.0
      endif


      ;if rep eq 13 then begin
      if rep eq 5 then begin

         Pic = pic /total(pic) * total(resultW)/1.5

         imageTemp2 = ((+alog10(Pic)-alog10(max(resultw)/200) +3  )/4.0)*50.0
         
         ind=where(imageTemp2 lt 0) 
         if ind(0) ne -1 then ImageTemp2(ind)=0
         ind=where(ImageTemp2 gt 50.0)
         if ind(0) ne -1 then Imagetemp2(ind)=50.0
      endif


      if rep eq 9 then begin

         Pic = pic /total(pic) * total(resultW)/4

         imageTemp2 = ((+alog10(Pic)-alog10(max(resultw)/800) +3  )/4.0)*50.0
         
         ind=where(imageTemp2 lt 0) 
         if ind(0) ne -1 then ImageTemp2(ind)=0
         ind=where(ImageTemp2 gt 50.0)
         if ind(0) ne -1 then Imagetemp2(ind)=50.0
      endif

   endif

   Map = ResultW  + Pic

;;; Now need to map the result onto a color table


   if rep eq 0 then begin
      ma = max(Map)
      ;mi= ma /10000
      ;mi= ma /6000
      ;mi= ma /4000
      mi= ma /1000
      ;mi= ma /600
   endif else begin
      ma*= 1.12884
      mi*= 1.02
   endelse


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

   ;if (rep eq 14) or (rep eq 9) or (rep eq 13) then begin
   ;if (rep eq 4) or (rep eq 5) then begin
   if (rep eq 5) then begin
      Imagetemp += imageTemp2
   endif


   imageTemp = fix(ImageTemp)

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



   ; load background
   read_jpeg, "/home2/tcox/Images/Astro/fields_cosmos.jpg", background_image_full
   bkgrnd_image= transpose(background_image_full[*, 1024:2047, 0:1023])

   ; v0
   JpegPic += 0.6*bkgrnd_image
   idx= where(JpegPic gt max(ColPlane))
   if idx(0) ne -1 then JpegPic(idx)= max(ColPlane)

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
   exts=exts+strcompress(string(rep),/remove_all)
   exts=strmid(exts,strlen(exts)-3,3)


   ;fout= "pic_hr_"+exts+".jpg"
   fout= pic_base+exts+".jpg"
   fout=strcompress(fout,/remove_all)

   write_jpeg, fout , Jpegpic, true=3, quality=90


   time= fload_time(1)

   cmd = 'mogrify -font arial -antialias -pointsize 58 -fill white -draw '
   cmd = cmd + "'"
   cmd = cmd + 'text 40,72 '
   cmd = cmd + """
   cmd = cmd + 'T = ' + string(time/0.7,format='(F4.2)')+' Gyr'
   cmd = cmd + """
   cmd = cmd + "'"


   cmd = cmd + ' -quality 95 '
   cmd = cmd + fout

   print, cmd
   spawn, cmd


endfor

end





 










 
