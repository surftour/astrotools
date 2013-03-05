; since we don't have a contour_makepic equivalent
; for the velocity map, we'll generate one
; -----------------------------------------------------
pro contour_makevelpic, x, y, z, vx, vy, vz, m, xlen, xz=xz, yz=yz, $
                        filename=filename, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, $
                        pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        center=center, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
			velocitymap=velocitymap, dispersionmap=dispersionmap, $
                        NxNImage=NxNImage
        



;XYPIXELS= 500L
XYPIXELS= 480L

if keyword_set(pixels) then XYPIXELS= long(pixels)

NxNImage= fltarr(XYPIXELS, XYPIXELS)
NxNImage= NxNImage*0.0 + 1.0


if not keyword_set(xlen) then xlen=100.0
if not keyword_set(sendto) then sendto='x'
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(x) then begin
   print, "  "
   print, "PROBLEM: contour_makepic"
   print, "  "
   print, "  "
   print, "  "
   return
endif



;--------------------------------------
;  Set up Parameters
;--------------------------------------





    ; Smoothing Parameters
    ; ---------------------
    ;DesNgb=30L
    ;DesNgb=1000L  ; make this really large so it always uses Hmax
    ;DesNgb=500L  ; make this really large so it always uses Hmax
    DesNgb=200L  ; make this really large so it always uses Hmax
    ;DesNgb=100L  ; make this really large so it always uses Hmax
    Hmax= 5.0


    ; Plot ranges
    ; ------------------------------------
    MaxDensXY= 300.0
    DynRange= 0.0

    if keyword_set(set_maxden) then MaxDensXY= set_maxden
    if keyword_set(set_dynrng) then DynRange= set_dynrng

    ma=MaxDensXY
    mi=MaxDensXY-DynRange



;----------------------------------
;----------------------------------

if not keyword_set(center) then center= [0,0,0]

if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin

        print, "rotate: moving to center to",center
        x2=x-center[0]
        y2=y-center[1]
        z2=z-center[2]

        print, "rotate: set center = [0,0,0]"
        center= [0,0,0]

        process_rotation, x2, y2, z2, rotate_theta, rotate_phi, x_new, y_new, z_new

        ;x= x_new
        ;y= y_new
        ;z= z_new

        process_rotation, vx, vy, vz, rotate_theta, rotate_phi, vx_new, vy_new, vz_new

        ;vx= vx_new
        ;vy= vy_new
        ;vz= vz_new

endif else begin
	x_new= x
	y_new= y
	z_new= z
	vx_new= vx
	vy_new= vy
	vz_new= vz
endelse



    N= long(n_elements(x))
    
    ; cut in z direction
    ; -------------------
    if keyword_set(zthickness) then begin
        idx= where((z_new LE center[2]+zthickness) and (z_new GE center[2]-zthickness))
        if idx(0) ne -1 then begin
                N= long(n_elements(idx))
                x_new= x_new(idx)
                y_new= y_new(idx)
                z_new= z_new(idx)
                vx_new= vx_new(idx)
                vy_new= vy_new(idx)
                vz_new= vz_new(idx)
                m= m(idx)
        endif
    endif

    ; cut in x direction
    ; -------------------
    if keyword_set(xthickness) then begin
        idx= where((x_new LE center[0]+xthickness) and (x_new GE center[0]-xthickness))
        if idx(0) ne -1 then begin
                N= long(n_elements(idx))
                x_new= x_new(idx)
                y_new= y_new(idx)
                z_new= z_new(idx)
                vx_new= vx_new(idx)
                vy_new= vy_new(idx)
                vz_new= vz_new(idx)
                m= m(idx)
        endif
    endif

    ; cut in y direction
    ; -------------------
    if keyword_set(ythickness) then begin
        idx= where((y_new LE center[1]+ythickness) and (y_new GE center[1]-ythickness))
        if idx(0) ne -1 then begin
                N= long(n_elements(idx))
                x_new= x_new(idx)
                y_new= y_new(idx)
                z_new= z_new(idx)
                vx_new= vx_new(idx)
                vy_new= vy_new(idx)
                vz_new= vz_new(idx)
                m= m(idx)
        endif
    endif



    print, "contour_makevelpic: changing center to= ", center


    ; OK we got all field set, send it away
    ; -------------------------------------
    Coord=fltarr(3,N)
    Velocity=fltarr(N)                 ; velocity out of the plane

    Coord(0,*)=x_new
    Coord(1,*)=y_new
    Coord(2,*)=z_new
    Velocity(*)=vz_new

    xmin= -xlen+center[0]
    xmax=  xlen+center[0]
    ymin= -xlen+center[1]
    ymax=  xlen+center[1]

    xpixels=XYPIXELS
    ypixels=XYPIXELS

    ; default is xy axis
    ; ------------------
    Axis1=0L    ; select x-axis
    Axis2=1L    ; select y-axis

    if keyword_set(xz) then begin
        Axis1=0L    ; select x-axis
        Axis2=2L    ; select z-axis
        ymin= -xlen+center[2]
        ymax=  xlen+center[2]

        Velocity(*)= vy_new
    endif

    if keyword_set(yz) then begin
        Axis1=1L    ; select y-axis
        Axis2=2L    ; select z-axis
        xmin= -xlen+center[1]
        xmax=  xlen+center[1]
        ymin= -xlen+center[2]
        ymax=  xlen+center[2]

        Velocity(*)= vz_new
    endif





    ; -----------------------------------
    ; Call external c program to smooth
    ; -----------------------------------
    Map=fltarr(xpixels,ypixels)
    if keyword_set(velocitymap) then begin
	spawn, 'echo $HOME', result
	homedir= strcompress(result,/remove_all)
	libfile= homedir+'/Tools/C-Routines_for_IDL/SmoothedAvg/smoothedavg.so'
;	libfile= homedir+'/Tools/C-Routines_for_IDL/NearestValue/nearesval.so'
	S = CALL_EXTERNAL(libfile, $
                      'smoothedavg', $
                      N, $
                      Coord,Velocity,$
                      xmin,xmax,ymin,ymax,$
                      Xpixels,Ypixels,$
                      DesNgb, $
                      Axis1, $
                      Axis2, $
                      Hmax, $
                      Map  )

; --------------
;  write fits
; --------------

    if keyword_set(fitstoo) then begin

        header = strarr((size(Map))[0] + 7) + string(' ',format='(a80)')      ;Create empty array
        header[0] = 'END' + string(replicate(32b,77))

        sxaddpar, header, 'AUTHOR', 'T.J. Cox',' Written by contour_makevelpic:  '+ systime()

        s= size(Map)
        stype = s[s[0]+1]              ;Type of data    
        case stype of
                0:      message,'ERROR: Input data array is undefined'
                1:      bitpix = 8
                2:      bitpix = 16
                3:      bitpix = 32
                4:      bitpix = -32
                5:      bitpix = -64
                6:      message,'Complex types not allowed as FITS primary arrays'
                7:      bitpix = 8
                12:      bitpix = 16
                13:      bitpix = 32
                14:      bitpix = 64
        else:   message,'ERROR: Illegal Image Datatype'
        endcase

        sxaddpar, header, 'BITPIX', -32, ' Number of bits per data pixel'      ; manually set this, see mkhdr for how to automatically do this
        sxaddpar, header, 'NAXIS', s[0],' Number of data axes'       ;# of dimensions

        if ( s[0] GT 0 ) then begin
           for i = 1, s[0] do sxaddpar,header,'NAXIS' + strtrim(i,2),s[i]
        endif

        if (s[0] EQ 0) then $
            sxaddpar, header, 'EXTEND', 'T', ' FITS data may contain extensions'
        Get_date, dte                       ;Get current date as CCYY-MM-DD
        sxaddpar, header, 'DATE', dte, 'Creation UTC (CCCC-MM-DD) date of FITS header'

        sxaddpar, header, 'XLEN', xlen, 'xlen used for image generation'
        sxaddpar, header, 'SIDELEN', 2.0*xlen, 'image side length = 2.0 * xlen'

        header = header[0:s[0]+7]

        sxaddpar,header,'COMMENT ', "TJ test"


        extrainfo= 'vel'
        n2= strlen(filename)-4
        fitsfilename= strmid(filename,0,n2)+'_'+extrainfo+'.fits'
	print, "writing: ", fitsfilename
	print, "max(Map)= ", max(Map)
	print, "min(Map)= ", min(Map)
        writefits, fitsfilename, Map, header
    endif



    endif



    if keyword_set(dispersionmap) then begin
	spawn, 'echo $HOME', result
	homedir= strcompress(result,/remove_all)
	libfile= homedir+'/Tools/C-Routines_for_IDL/SmoothedDisp/smootheddisp.so'
	S = CALL_EXTERNAL(libfile, $
                      'smootheddisp', $
                      N, $
                      Coord,Velocity,$
                      xmin,xmax,ymin,ymax,$
                      Xpixels,Ypixels,$
                      DesNgb, $
                      Axis1, $
                      Axis2, $
                      Hmax, $
                      Map  )


; --------------
;  write fits
; --------------

    if keyword_set(fitstoo) then begin

        header = strarr((size(Map))[0] + 7) + string(' ',format='(a80)')      ;Create empty array
        header[0] = 'END' + string(replicate(32b,77))

        sxaddpar, header, 'AUTHOR', 'T.J. Cox',' Written by contour_makevelpic:  '+ systime()

        s= size(Map)
        stype = s[s[0]+1]              ;Type of data    
        case stype of
                0:      message,'ERROR: Input data array is undefined'
                1:      bitpix = 8
                2:      bitpix = 16
                3:      bitpix = 32
                4:      bitpix = -32
                5:      bitpix = -64
                6:      message,'Complex types not allowed as FITS primary arrays'
                7:      bitpix = 8
                12:      bitpix = 16
                13:      bitpix = 32
                14:      bitpix = 64
        else:   message,'ERROR: Illegal Image Datatype'
        endcase

        sxaddpar, header, 'BITPIX', -32, ' Number of bits per data pixel'      ; manually set this, see mkhdr for how to automatically do this
        sxaddpar, header, 'NAXIS', s[0],' Number of data axes'       ;# of dimensions

        if ( s[0] GT 0 ) then begin
           for i = 1, s[0] do sxaddpar,header,'NAXIS' + strtrim(i,2),s[i]
        endif

        if (s[0] EQ 0) then $
            sxaddpar, header, 'EXTEND', 'T', ' FITS data may contain extensions'
        Get_date, dte                       ;Get current date as CCYY-MM-DD
        sxaddpar, header, 'DATE', dte, 'Creation UTC (CCCC-MM-DD) date of FITS header'

        sxaddpar, header, 'XLEN', xlen, 'xlen used for image generation'
        sxaddpar, header, 'SIDELEN', 2.0*xlen, 'image side length = 2.0 * xlen'

        header = header[0:s[0]+7]

        sxaddpar,header,'COMMENT ', "TJ test"


        extrainfo= 'dis'
        n2= strlen(filename)-4
        fitsfilename= strmid(filename,0,n2)+'_'+extrainfo+'.fits'
	print, "writing: ", fitsfilename
	print, "max(Map)= ", max(Map)
	print, "min(Map)= ", min(Map)
        writefits, fitsfilename, Map, header
    endif




    endif








; -----------------------------------------------------
; now map the field to a color scale
; -----------------------------------------------------


    ; print ranges
    print, "Map     max: ", max(Map), "   min: ", min(Map)
    print, "Clipping at   ma= ", ma, " mi= ", mi


    ; do some clipping

    ind=where(Map lt mi)
    if ind(0) ne -1 then Map(ind)=mi
    ind=where(Map gt ma)
    if ind(0) ne -1 then Map(ind)=ma


    Cols=254 ; number of colors

    ; ------------------
    ; do this linearly, and remembering 
    ; that we'll have neg. numbers
    ;

    ; for velocitymap
    ;if keyword_set(velocitymap) then Pic= (Map+MaxDensXY)/DynRange * (cols-3) + 2
    if keyword_set(velocitymap) then Pic= (Map+MaxDensXY)/DynRange * cols + 2

    ; for dispersionmap
    if keyword_set(dispersionmap) then Pic= (Map/MaxDensXY) * (cols-3) + 2

    ; ------------------

    ; this should never occur
    ;ind=where(Pic ge 256)
    ;if ind(0) ne -1 then Pic(ind)=255

    ; this might occur if the field was zero to begin with
    ;ind=where(Pic eq 0)
    ;if ind(0) ne -1 then Pic(ind)=2

    zero_idx=where(Pic eq 0)

    ; should be an array with values between 2 and 256 at this point
    ;invertcolors= 1
    invertcolors= 0
    if invertcolors eq 1 then begin
        Pic=258-Pic                           ; invert the color table
    endif else begin
        ;idx= where(Pic EQ 2)
    endelse

    idx= where(Pic gt 255)                ; it seems that many color tables have 256 as black
    if idx(0) ne -1 then Pic(idx)= 255

    ;if idx(0) NE -1 then Pic(idx)= 1       ; this is setting it to white
    ;if idx(0) NE -1 then Pic(idx)= 0       ; this is setting it to black

    if zero_idx(0) ne -1 then Pic(zero_idx)= 1



; -----------------------------
;   set image and we're done
; -----------------------------

    NxNImage= Pic


; -------------
;  Done 
; -------------
        

end




