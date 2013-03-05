
;Object_File='C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so'
Object_File='/n/home/tcox/C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so'

add_time_lbl= 0

;Frun =  "data/"
;Frun =  "/raid4/tcox/vc3vc3e/"
;Frun = "/raid5/yxli/Runs/merger_tree/run4/Sims/1gpc_zoom1/test_2/z_8.50_0/"
;Frun = "/n/home/tcox/data/ds/vc3vc3e_2/"
;Frun = "/n/home/tcox/data/ds/vc3vc3e_no/" & add_time_lbl= 1
;Frun = "/n/home/tcox/data/ds/vc3vc3h_2/"
;Frun = "/n/home/tcox/data/ds/vc3vc3h_no/" & add_time_lbl= 1
;Frun = "/n/home/tcox/data/bs/b3f/"
;Frun = "/n/home/tcox/data/minor/Sbfg0.4Scfg0.4_090/"
;Frun = "/n/home/tcox/data/minor/Sbfg0.4Scfg0.4_150/"
;Frun = "/n/home/tcox/data/minor/Sbfg0.4Scfg0.4_150/"
Frun= "/n/home/tcox/data/kecktests/Sb2x1Sb2x1_e_32"

;PicsDir = "vc3vc3e_2/"
;PicsDir = "vc3vc3e_no/"
;PicsDir = Frun+"/movie_v2/"
PicsDir = Frun+"/movie_xy/"
;PicsDir = Frun+"/movie_xz/"

Snap = "snapshot"
;Snap = "z_8.50_0"

spawn, "mkdir "+PicsDir

get_time_list, Frun, Snap, NumList, TiList

snap_min = 0
;snap_min= 113
snap_max = 3
;snap_max = 5
;snap_max = 29
;snap_max= 30
;snap_max= 116
;snap_max= 615
;snap_max = 913

NumPix = 9
;NumPix = 10
;NumPix = 300
;NumPix = 348

PIX=   480

LenBase =  50.0

Len =  LenBase


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



for picnr=0, NumPix-1 do begin

    exts='0000'
    exts=exts+strcompress(string(picnr),/remove_all)
    exts=strmid(exts,strlen(exts)-4,4)
    f=PicsDir +"/pic_gas_"+exts+".jpg"
    f=strcompress(f,/remove_all)
    picfname= f


    openr,1,picfname, error=err ;, /swap_endian
    if err eq 0 then begin
        close,1
        goto,skipit
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
    interpolate_gas_data_v2, Time, Frun, Snap, NumList, TiList, N, Pos, Vel, Hsml, Mass, Temp



;;;;;;;;;;;;; do a rotation - if you want to, of course

; 090 rotation angles
;theta_deg= -130.5
;phi_deg= -21.0
; 150 rotation angles
theta_deg= -141.3
phi_deg= -104.0 + 90.0

x= Pos(0,*)
y= Pos(1,*)
z= Pos(2,*)

        ; rotate the coordinates
        ; ------------------------
	theta= theta_deg * !PI / 180.0
	phi= phi_deg * !PI / 180.0

        ; first around z axis (phi)
        ; then around y axis (theta)
        x_new= x*(cos(theta)*cos(phi)) + y*(cos(theta)*sin(phi)) + z*sin(theta)
        y_new= -x*sin(phi) + y*cos(phi)
        z_new= -x*(sin(theta)*cos(phi)) - y*(sin(theta)*sin(phi)) + z*cos(theta)

Pos(0,*)= x_new
Pos(1,*)= y_new
Pos(2,*)= z_new

;;;;;;;;;;;;; start  X-Y- projection of density

    N=   long(N)

    Pos2d=  fltarr(2,N)      ; generate coordinate array for particles

    
    phi=!pi/2.0
    
    x = Pos(0,*)  
    z = Pos(1,*)*cos(phi) - Pos(2,*)*sin(phi)
    y = Pos(1,*)*sin(phi) + Pos(2,*)*cos(phi)


    ; if doing xz projection, then swap some things
    ;temp= y
    ;y= z
    ;z= temp


    bw = 65.0

    x = bw/(bw - z)*x
    y = bw/(bw - z)*y
    hsml = bw/(bw - z)*hsml
    Mass = (bw/(bw - z))*Mass


    Pos2d(0,*)= x        ; fill in the x-coordinates
    Pos2d(1,*)= y        ; ...         y-coordinates

    Hsml= float(Hsml)

    Mass = float(Mass)

    ind =where(Hsml gt 20.0)
    if ind(0) ne -1 then begin
       Mass(ind) = 0
       Hsml(ind) = 0.01
    endif

    ind =where(z ge bw)
    if ind(0) ne -1 then begin
       Mass(ind) = 0
       Hsml(ind) = 0.001
    endif

 

    quantity = fltarr(N)        ; generate space for mass array
    quantity(*)= Temp(*)        ; fill in temperature

    xmin= float(-len*4.0/3)     ; left edge of map
    xmax= float( len*4.0/3)     ; right edge of map 
    ymin= float(-len)           ; lower
    ymax= float( len)           ; upper edge

    xpixels= long(PIX)*4/3      ; pixels in x-direction
    ypixels= long(PIX)          ; pixels in y-direction

    ValueXY= fltarr(ypixels,xpixels)      ; holds the result
    ValueXY_Temp= fltarr(ypixels,xpixels) ; holds the result

    S = CALL_EXTERNAL(Object_File, $
                      'slice', $
                      N, $
                      pos2D, hsml, mass, quantity, $
                      xmin,xmax,ymin,ymax,$
                      Xpixels,Ypixels,$
                      ValueXY, ValueXY_Temp)

    ValueXY = transpose(ValueXY)
    ValueXY_Temp= transpose(ValueXY_Temp)



;;; Now need to map the result onto a color table

    Map= ValueXY

;;    ma = max(Map)
    ma = 50*3.0e-5 / (16.0/len)^2  
    mi=  ma /10000

    print, "max= ", ma, "min=", mi, max(map)



;; now do clipping in Map

    ind=where(Map lt mi) 
    if ind(0) ne -1 then Map(ind)=mi
    ind=where(Map gt ma)
    if ind(0) ne -1 then Map(ind)=ma

;; now map the log of the field linearly onto the color table

    cols= 255

    ImageDens= byte((alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-2) +2)

    
    ind=where(ValueXY eq 0)
    if ind(0) ne -1 then ValueXy(ind)=1
    TempMap= ValueXY_Temp / ValueXY ; define mass-weighted temperature map

    Map= TempMap

    ma= 3.0e6                 ; select maximum value (brightest color)
    mi= 8.0e3                 ; select minimum value


;; now do clipping in Map

    ind=where(Map lt mi) 
    if ind(0) ne -1 then Map(ind)=mi
    ind=where(Map gt ma)
    if ind(0) ne -1 then Map(ind)=ma

;; now map the log of the field linearly onto the color table

    cols= 255

    ImageTemp= byte((alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-2) +2)


;; Load two-dimensional color table    

    cols=256
    ;ColPlane=fltarr(cols,cols,3)

    ;openr,1,"colplane2.dat"
    ;readu,1,ColPlane
    ;close,1

    invertcolors= 1
    if invertcolors eq 1 then begin
        ImageDens=256-ImageDens                          ; invert color table
        idx= where(ImageDens EQ 254)               ; set background to white
    endif else begin
        idx= where(ImageDens EQ 2)
    endelse

    ImageDens(idx)= 1    ; white

    JpegPic=bytarr(XPixels,YPixels,3)

    for i=0, XPixels-1 do begin
        for j=0, YPixels-1 do begin
            ;JPegPic(i,j,0)= ColPlane(imageDens(i,j),imageTemp(i,j),0)
            ;JPegPic(i,j,1)= ColPlane(imageDens(i,j),imageTemp(i,j),1)
            ;JPegPic(i,j,2)= ColPlane(imageDens(i,j),imageTemp(i,j),2)
            JPegPic(i,j,0)= blue0(imageDens(i,j))
            JPegPic(i,j,1)= blue1(imageDens(i,j))
            JPegPic(i,j,2)= blue2(imageDens(i,j))
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





 










 
