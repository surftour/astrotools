
;Object_File='C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so'
;Object_File='/n/home/tcox/C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so'
Object_File='/n/home/tcox/C-Routines_for_IDL/ComputeHsmlAndProject/HsmlAndProject.so'

add_time_lbl= 0

;Frun =  "data/"
;Frun =  "/raid4/tcox/vc3vc3e/"
;Frun = "/raid5/yxli/Runs/merger_tree/run4/Sims/1gpc_zoom1/test_2/z_8.50_0/"
;Frun = "/n/home/tcox/data/ds/vc3vc3e_2/"
Frun = "/n/home/tcox/data/ds/vc3vc3e_no/" ; & add_time_lbl= 1
;Frun = "/n/home/tcox/data/ds/vc3vc3h_2/"
;Frun = "/n/home/tcox/data/ds/vc3vc3h_no/" ; & add_time_lbl= 1

;PicsDir = "vc3vc3e_2/"
;PicsDir = "vc3vc3e_no/"
PicsDir = Frun+"/movie_v/"

Snap = "snapshot"
;Snap = "z_8.50_0"

spawn, "mkdir "+PicsDir

get_time_list, Frun, Snap, NumList, TiList

snap_min = 0
;snap_max = 5
;snap_max = 29
;snap_max = 913
snap_max= 107

;NumPix = 10
NumPix = 400

PIX=   480

LenBase =  70.0
Len =  LenBase

gasfrac_max= 0.8
gasfrac_min= 0.0


for picnr=0, NumPix-1 do begin

; testing purposes
;if picnr gt 4 then goto, skipit


    exts='0000'
    exts=exts+strcompress(string(picnr),/remove_all)
    exts=strmid(exts,strlen(exts)-4,4)
    f=PicsDir +"/pic_stars_"+exts+".jpg"
    f=strcompress(f,/remove_all)
    picfname= f


    openr,1,picfname, error=err ;, /swap_endian
    if err eq 0 then begin
        close,1
        goto,skipit
    endif

    ; add this to help speed things up
    if picnr lt 281 then begin
	print, "skipping ", picfname
	goto, skipit
    endif

    openw,1,picfname,/f77_unformatted
    writeu,1,Len
    close,1

    print
    print, "PicNr=", Picnr
    print

    Phi = 0.0

    Time = TiList(snap_min) + picnr*(TiList(snap_max)-TiList(snap_min))/float(NumPix-1)

    timelbl= strcompress(string(time/0.7),/remove_all)
    timelbl= strmid(timelbl,0,5)+" Gyr"

    print,'time = ',time
    print,'frun = ',frun
    print,'snap = ',snap
    interpolate_stars_data, Time, Frun, Snap, NumList, TiList, N, Pos, Vel, Mass, Id, gasfrac




;;;;;;;;;;;;; start  X-Y- projection of density

    N=   long(N)

    Pos2d=  fltarr(2,N)      ; generate coordinate array for particles

    
    phi=!pi/2.0
    
    x = Pos(0,*)  
    z = Pos(1,*)*cos(phi) - Pos(2,*)*sin(phi)    ; basically flips z to -z
    y = Pos(1,*)*sin(phi) + Pos(2,*)*cos(phi)

    bw = 65.0

    x = bw/(bw - z)*x
    y = bw/(bw - z)*y
    Mass = (bw/(bw - z))*Mass


    Pos2d(0,*)= x        ; fill in the x-coordinates
    Pos2d(1,*)= y        ; ...         y-coordinates

    Mass = float(Mass)

    ind =where(z ge bw)
    if ind(0) ne -1 then begin
       Mass(ind) = 0
    endif

 

    ;
    ;  Good if we know Hsml
    ;-------------------------------------------------------------
    ;quantity = fltarr(N)        ; generate space for mass array
    ;quantity(*)= Temp(*)        ; fill in temperature
    ;
    ;xmin= float(-len*4.0/3)     ; left edge of map
    ;xmax= float( len*4.0/3)     ; right edge of map 
    ;ymin= float(-len)           ; lower
    ;ymax= float( len)           ; upper edge
    ;
    ;xpixels= long(PIX)*4/3      ; pixels in x-direction
    ;ypixels= long(PIX)          ; pixels in y-direction
    ;
    ;ValueXY= fltarr(ypixels,xpixels)      ; holds the result
    ;ValueXY_Temp= fltarr(ypixels,xpixels) ; holds the result
    ;
    ;S = CALL_EXTERNAL(Object_File, $
    ;                  'slice', $
    ;                  N, $
    ;                  pos2D, hsml, mass, quantity, $
    ;                  xmin,xmax,ymin,ymax,$
    ;                  Xpixels,Ypixels,$
    ;                  ValueXY, ValueXY_Temp)
    ;
    ;ValueXY = transpose(ValueXY)
    ;ValueXY_Temp= transpose(ValueXY_Temp)
    ;-------------------------------------------------------------


   Xyz=  fltarr(3,N)            ; generate coordinate array for particles
   Xyz(0,*)= x
   Xyz(1,*)= y
   Xyz(2,*)= z
   
   Hsml = fltarr(N)

   Weight = fltarr(N)           ; generate space for mass array
   Weight(0:*)= Mass   ; fill in mass values
   quantity = fltarr(N)       
   quantity(*)= 1.0             ; fill in mass values
   Axis1=   long(0)             ; select horizontal axis
   Axis2=   long(1)             ; select vertical axis
   Axis3=   long(2)             ; select projection axis
   
   ;xmin= float(-len)            ; left edge of map
   ;xmax= float( len)            ; right edge of map 
   xmin= float(-len*4.0/3)            ; left edge of map
   xmax= float( len*4.0/3)            ; right edge of map 
   ymin= float(-len)            ; lower
   ymax= float( len)            ; upper edge
   Zmin= float(-1.0e6)          ; lower range of particles that are included
   Zmax= float(1.0e6)           ; upper range
   

   xpixels= long(PIX*4/3)           ; pixels in x-direction
   ypixels= long(PIX)           ; pixels in y-direction
   
   Hmax=    float( (xmax-xmin)/xpixels * 50) ; maximum used smoothing length
   BoxSize  = double(0.0)

   DesNgb = 32L                 ; for those particles that have no zero smoothing length

   ResultW = fltarr(YPixels, XPixels) ; projected weights (normally projected mass map)
   ResultQ = fltarr(YPixels, XPixels) ; projected `Weight' weighted quantity
   
   
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


   ResultW = transpose(ResultW)
   ;ResultW /= total(ResultW)

   print, "ResultW  max/min= ", max(ResultW), min(ResultW)



	ValueXY= ResultW

;----------------------------------------------------------------

; put in "spikes" of BH emision



;----------------------------------------------------------------

;;; Now need to map the result onto a color table

    Map= ValueXY

;;    ma = max(Map)
    ma = 10*3.0e-5 / (16.0/len)^2  
    mi=  ma /8000

    print, "max= ", ma, "min=", mi, max(map)



;; now do clipping in Map

    ind=where(Map lt mi) 
    if ind(0) ne -1 then Map(ind)=mi
    ind=where(Map gt ma)
    if ind(0) ne -1 then Map(ind)=ma

;; now map the log of the field linearly onto the color table

    cols= 255

    ImageDens= byte((alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-2) +2)

    

   imageTemp = fltarr(XPixels, ypixels)
   Imagetemp(*) = (gasfrac_max-gasfrac)/(gasfrac_max-gasfrac_min) * 180
   imageTemp = fix(ImageTemp)

;; Load two-dimensional color table    

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

;;;    window, xsize=PIX*4/3, ysize=PIX
;;;    tv, JPegPic, true=3,0

    if add_time_lbl eq 1 then begin
	temppicfname= PicsDir + "/temp.jpg"
	write_jpeg, temppicfname, JPegPic, true=3, quality=95

	; add time onto the figure
	;
	;
	cmd = 'convert -font Courier-Bold -antialias -pointsize 30 -fill white -draw '
	cmd = cmd + "'"
	cmd = cmd + 'text 25,51 '
	cmd = cmd + """
	cmd = cmd + timelbl
	cmd = cmd + """
	cmd = cmd + "'"
	cmd = cmd + ' -quality 95 '
	cmd = cmd + temppicfname + " " + picfname

	print, cmd
	spawn, cmd
    endif else begin
	write_jpeg, picfname, JPegPic, true=3, quality=95
    endelse

skipit:

endfor

end





 










 
