
Object_File='/home2/tcox/Tools/C-Routines_for_IDL/ComputeHsmlAndProject/HsmlAndProject.so'

frun = "/raid4/tcox/bs/b4e/" 
gasfrac_max= 0.8
gasfrac_min= 0.0

Snap=   "snapshot"

pic_base= "pic_stars"


get_time_list, Frun, snap, NumList, TiList

time_begin= 0.15
time_end= 1.35
n_pics= 20

time_arr= [0.1, 0.3, 0.5, 0.8, 1.1, 1.2, 1.8]
n_pics= n_elements(time_arr)

PIX = 1024
Len = 70.0


xpixels= long(PIX)              ; pixels in x-direction
ypixels= long(PIX)              ; pixels in y-direction
;window, xsize=XPixels,ysize=Ypixels,retain=2
seed=42L

for rep = 0,n_pics-1 do begin
   
   
   ti = time_begin + (time_end-time_begin)/(n_pics-1)*rep

   if (size(time_arr))[0] gt 0 then ti= time_arr[rep]

   mm = abs(tilist -ti)
   ind = sort(mm)
   num = numlist(ind(0))


   UnitLength_in_cm =        3.085678d21 ;  1.0 kpc 
   UnitMass_in_g    =        1.989d43 ;  1.0e10 solar masses 
   UnitVelocity_in_cm_per_s = 1d5 ;  1 km/sec 
   UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s 
   UnitDensity_in_cgs= UnitMass_in_g/ UnitLength_in_cm^3
   UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ UnitTime_in_s^2
   UnitEnergy_in_cgs= UnitMass_in_g * UnitLength_in_cm^2 / UnitTime_in_s^2
   GRAVITY   = 6.672d-8 
   BOLTZMANN = 1.3806d-16
   PROTONMASS = 1.6726e-24
   G=GRAVITY/ UnitLength_in_cm^3 * UnitMass_in_g * UnitTime_in_s^2
   Xh=0.76D                     ; mass fraction of hydrogen
   HubbleParam= 0.65
   gamma= 5.0/3

   exts='000'
   exts=exts+strcompress(string(num),/remove_all)
   exts=strmid(exts,strlen(exts)-3,3)
   f=frun + snap+ "_"+exts
   f=strcompress(f,/remove_all)

   npart=lonarr(6)	
   massarr=dblarr(6)
   time=0.0D
   redshift=0.0D
   flag_sfr=0L
   flag_feedback=0L
   npartall=lonarr(6)	
   bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4
   la=intarr(bytesleft/2)
   
   print, "opening: ", f
   openr,1,f,/f77_unformatted   ;, /swap_endian
   readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,la
   print
   print,num, time,redshift
   print
   N=total(npart)
   pos=fltarr(3,N)
   id=lonarr(N)
   readu,1,pos
   readu,1                      ;,vel
   readu,1 ,id
   ind=where((npart gt 0) and (massarr eq 0))
   if ind(0) ne -1 then begin
      Nm= total(npart(ind))
      mass=fltarr(Nm)	
      readu,1,mass
   endif

   NGas=npart(0)
   NHalo=npart(1)
   NDisk=npart(2)
   NBulge=npart(3)
   NStars=npart(4)
   NHoles=npart(5)

   Nskip = NGas + NHalo + NDisk +  NBulge +  NStars



   u=fltarr(Ngas)
   readu,1,u
   rho=fltarr(Ngas)
   readu,1,rho
   Nelec=fltarr(Ngas)
   readu,1 ,Nelec
   readu,1                      ;,NH0

   Hsml=fltarr(Ngas)
   readu,1, Hsml
   Sfr= fltarr(Ngas)
   readu,1, Sfr

   close,1


   MeanWeight= 4.0/(3*Xh+1+4*Xh*Nelec) * PROTONMASS
   Temp = MeanWeight/BOLTZMANN * (gamma-1) * U * UnitEnergy_in_cgs/ UnitMass_in_g


   xc = Pos(0,Nskip)            ; position of BH
   yc = Pos(1,Nskip) 
   zc = Pos(2,Nskip)


   if Ngas gt 0 then begin
      xgas=fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas) & mgas=fltarr(Ngas)
      xgas(*)=pos(0,0:Ngas-1)
      ygas(*)=pos(1,0:Ngas-1)
      zgas(*)=pos(2,0:Ngas-1)

      if massarr(0) eq 0 then begin
         mgas(*)=mass(0:Ngas-1)	
      endif else begin
         mgas(*)= massarr(0)
      endelse
   endif

   ind =where((abs(xgas) lt len) and (abs(ygas) lt len) and (abs(zgas) lt len))
   Mgas =total(mgas(ind))
   


   if Nstars gt 0 then begin
      xstars=fltarr(Nstars) &  ystars=fltarr(Nstars)  & zstars=fltarr(Nstars)  & mstars=fltarr(Nstars)
      xstars(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
      ystars(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
      zstars(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
      if massarr(4) eq 0 then begin
         skip=0L
         for t=0,2 do begin
            if (npart(t) gt 0) and (massarr(t) eq 0) then begin
               skip=skip + npart(t)
            endif
         endfor
         mstars(*)=mass(0+skip:Nstars-1+skip)
      endif else begin
         mstars(*)= massarr(4)
      endelse
   endif

   if Ndisk gt 0 then begin
      xdisk=fltarr(NDisk) &  ydisk=fltarr(NDisk) &  zdisk=fltarr(NDisk) & mdisk=fltarr(NDisk)
      xdisk(*)=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
      ydisk(*)=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
      zdisk(*)=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
      if massarr(2) eq 0 then begin
         skip=0L
         for t=0,1 do begin
            if (npart(t) gt 0) and (massarr(t) eq 0) then begin
               skip=skip + npart(t)
            endif
         endfor
         mdisk(*)=mass(0+skip:Ndisk-1+skip)
      endif else begin
         mdisk(*)= massarr(2)
      endelse
   endif




;;;;;;;;;;;;;


   center= [0.0, 0.0, 0.0]
   if (rep eq 6) then center= [0.0, -10.0, 0.0]




;;;;;;;;;;;;; start  X-Y- projection of density

   N=   long(Nstars+Ndisk)      ; number of particles

   Xyz=  fltarr(3,N)            ; generate coordinate array for particles
   Xyz(0,0:nstars-1)= xstars(*) ; fill in the x-coordinates
   Xyz(1,0:nstars-1)= ystars(*) ; ...         y-coordinates
   Xyz(2,0:nstars-1)= zstars(*) ; ...         y-coordinates
   Xyz(0,Nstars:*)= xdisk(*)    ; fill in the x-coordinates
   Xyz(1,Nstars:*)= ydisk(*)    ; ...         y-coordinates
   Xyz(2,Nstars:*)= zdisk(*)    ; ...         y-coordinates

   Xyz(0,*)= Xyz(0,*) - center(0)
   Xyz(1,*)= Xyz(1,*) - center(1)
   Xyz(2,*)= Xyz(2,*) - center(2)

   Hsml = fltarr(N)    

   Weight = fltarr(N)           ; generate space for mass array
   Weight(0:nstars-1)= mstars(*) ; fill in mass values
   Weight(nstars:*)= mdisk(*)   ; fill in mass values
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




   ind =where((abs(xyz(0,*)) lt len) and (abs(xyz(1,*)) lt len) and (abs(xyz(2,*)) lt len))
   Mstars =total(Weight(ind))

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

         xc = (Pos(0,Nhalo+Ndisk+Ngas+NBulge+Nstars+bhnr)-xmin)/(xmax-xmin)*xs
         yc=  (Pos(1,Nhalo+Ndisk+Ngas+NBulge+Nstars+bhnr)-ymin)/(ymax-ymin)*ys

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
      mi= ma /4000
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


   ;tv, JPegPic, true=3,0

   exts='000'
   exts=exts+strcompress(string(rep),/remove_all)
   exts=strmid(exts,strlen(exts)-3,3)


   ;fout= "pic_hr_"+exts+".jpg"
   fout= pic_base+exts+".jpg"
   fout=strcompress(fout,/remove_all)

   write_jpeg, fout , Jpegpic, true=3, quality=90



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





 










 
