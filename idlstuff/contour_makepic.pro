pro contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, $
			xthickness=xthickness, ythickness=ythickness, $
			filename=filename, fitstoo=fitstoo, $
			pixels=pixels, zthickness=zthickness, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			crude=crude, center=center, $
			NxNImage=NxNImage, set_maxden=set_maxden, set_dynrng=set_dynrng, domaxrng=domaxrng, $
			compute_halfmasscnt=compute_halfmasscnt, HalfMassSB=HalfMassSB, FindMassFactor=FindMassFactor, $
			var_already_inlog=var_already_inlog, $
			NxNTempImage=NxNTempImage, RawNxNImage=RawNxNImage


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

if not keyword_set(x) then begin
   print, "  "
   print, "PROBLEM: contour_makepic"
   print, "  "
   print, "  "
   print, "  "
   return
endif


if not keyword_set(hsml) then hsml=x*0.0 + 1.0

if not keyword_set(center) then center= [0.0, 0.0, 0.0]




;----------------------------------
;----------------------------------


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
    
    

    if keyword_set(center) then begin

        print, "moving to center", center
        x=x-center[0]
        y=y-center[1]
        z=z-center[2]

        print, "set center = [0,0,0]"
        center= [0,0,0]

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




    ; OK, we got fields, prepare to send
    ; -----------------------------------

    img_makepic_raw, x, y, z, m, xlen, xz=xz, yz=yz, $
                        hsml=hsml, $
                        pixels=pixels, $
                        crude=crude, $
                        NxNImage=NxNImage, $
                        NxNTempImage=NxNTempImage, $
			Pixel_Area=Pixel_Area, $
			Total_Mass=Total_Mass, $
			Plotted_Mass=Plotted_Mass


     RawNxNImage= NxNImage





; -----------------------------------------------------
; now find the half mass contour
; -----------------------------------------------------

HalfMassContourLevel= -1.0
if not keyword_set(FindMassFactor) then FindMassFactor= 0

;compute_halfmasscnt=0
;compute_halfmasscnt=1
if keyword_set(compute_halfmasscnt) then begin

        img_compute_halfmasscnt, NxNImage, HalfMassContourLevel, $
                                Pixel_Area=Pixel_Area, $
                                Total_Mass=Total_Mass, $
                                FindMassFactor=FindMassFactor

endif






; ----------------
;   write fits
; ----------------

    if keyword_set(fitstoo) then begin


        header = strarr((size(RawNxNImage))[0] + 7) + string(' ',format='(a80)')      ;Create empty array
        header[0] = 'END' + string(replicate(32b,77))

        sxaddpar, header, 'AUTHOR', 'T.J. Cox',' Written by contour_makevelpic:  '+ systime()

        s= size(RawNxNImage)
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
        sxaddpar, header, 'HALFMASS', HalfMassContourLevel, 'half mass contour level'

        header = header[0:s[0]+7]

        sxaddpar,header,'COMMENT ', "TJ test"


        n2= strlen(filename)-4                      ; assumed filename is *.eps
        fitsfilename= strmid(filename,0,n2)+'_sd.fits'
	print, " "
	print, "writing FITS file: "+fitsfilename
	help, RawNxNImage
	print, "max(RawNxNImage)= ", max(RawNxNImage)
	print, "min(RawNxNImage)= ", min(RawNxNImage)
	print, "xlen= ", xlen
	print, " "
        writefits, fitsfilename, RawNxNImage, header

    endif







; -----------------------------------------------------
; now map the field logarithmically to a color scale
; -----------------------------------------------------


    ; Plot ranges
    ; ------------------------------------
    ;MaxDensXY= 1.0e+3
    ;MaxDensXY= 1.0e+1
    ;MaxDensXY= 2.0
    ;MaxDensXY= 1.0
    MaxDensXY= 1.0e-1            ; (gadget units, 10^10 msun/ kpc^2 = 10^4 msun/pc^2)
    ;MaxDensXY= 5.0e-2
    ;MaxDensXY= 1.0e-2
    if keyword_set(domaxrng) then MaxDensXY= max(NxNImage)
    
    ;DynRange=1.0e2
    ;DynRange=1.0e3
    ;DynRange=5.0e3
    DynRange=1.0e4
    ;DynRange=1.0e5
    ;DynRange=1.0e6
    ;DynRange=1.0e10
    ;DynRange=50
    if keyword_set(domaxrng) then DynRange=1.0e4
    
    ma=MaxDensXY
    mi=MaxDensXY/DynRange

    if keyword_set(set_maxden) then ma= set_maxden
    if keyword_set(set_dynrng) then mi= set_maxden/set_dynrng



    ; print ranges
    print, "NxNImage     max: ", max(NxNImage), "   min: ", min(NxNImage)
    print, "Clipping at   ma= ", ma, " mi= ", mi



    ; do some clipping

    ind=where(NxNImage lt mi)
    if ind(0) ne -1 then NxNImage(ind)=mi
    ind=where(NxNImage gt ma)
    if ind(0) ne -1 then NxNImage(ind)=ma


    Cols=255 ; number of colors

    if keyword_set(var_already_inlog) then begin
	Pic=(NxNImage-mi)/(ma/mi) * (cols-3) + 2
    endif else begin
	Pic=(alog10(NxNImage)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2
    endelse

    HalfMassContourPic= (alog10(HalfMassContourLevel)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2

    ind=where(Pic ge 256)
    if ind(0) ne -1 then Pic(ind)=255

    ind=where(Pic le 0)
    if ind(0) ne -1 then Pic(ind)=1


    invertcolors= 1
    if invertcolors eq 1 then begin
	Pic=256-Pic			     ; invert color table
        HalfMassContourPic=256-HalfMassContourPic
	idx= where(Pic EQ 254)               ; set background to white
    endif else begin
	idx= where(Pic EQ 2)
    endelse

    if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
    ;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black




;--------------------
; now do the same thing for the temperature image - if we have one.

if not keyword_set(ValueXY_Temp) then ValueXY_Temp= Pic*0.0

if total(ValueXY_Temp) gt 0 then begin
    ind=where(ValueXY eq 0)
    if ind(0) ne -1 then ValueXy(ind)=1
    TempMap= ValueXY_Temp / ValueXY ; define mass-weighted temperature map
    
    ma= 1.0e6                 ; select maximum value (brightest color)
    mi= 8.0e3                 ; select minimum value

    
;; now do clipping in Map

    ind=where(TempMap lt mi)
    if ind(0) ne -1 then TempMap(ind)=mi
    ind=where(TempMap gt ma)
    if ind(0) ne -1 then TempMap(ind)=ma

;; now map the log of the field linearly onto the color table
    
    cols= 255

    TempPic= byte((alog10(TempMap)-alog10(mi))/(alog10(ma/mi)) * (cols-2) +2)
    
endif else begin
    TempPic= Pic*0.0
endelse



; -----------------------------
;   set image and we're done
; -----------------------------

    NxNImage= Pic

    HalfMassSB= HalfMassContourPic

    NxNTempImage= TempPic

; -------------
;  Done
; -------------



end


