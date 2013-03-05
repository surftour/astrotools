pro movie_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, $
			xthickness=xthickness, ythickness=ythickness, $
			filename=filename, fitstoo=fitstoo, $
			pixels=pixels, xpixels=xpixels, ypixels=ypixels, $
			zthickness=zthickness, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			crude=crude, center=center, $
			floatingbckgrnd=floatingbckgrnd, $
			NxNImage=NxNImage, set_maxden=set_maxden, set_dynrng=set_dynrng, $
			HalfMassSB=HalfMassSB, $
			var_already_inlog=var_already_inlog, $
			NxNTempImage=NxNTempImage


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


if keyword_set(pixels) then XYPIXELS= long(pixels) else XYPIXELS= 480L

if keyword_set(xpixels) then xpixels= long(xpixels) else xpixels= long(XYPIXELS)
if keyword_set(ypixels) then ypixels= long(ypixels) else ypixels= long(XYPIXELS)

NxNImage= fltarr(xpixels, ypixels)
NxNImage= NxNImage*0.0 + 1.0

NxNTempImage= fltarr(xpixels, ypixels)
NxNTempImage= NxNImage*0.0 + 1.0
ValueXY_Temp= -1


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


if not keyword_set(hsml) then hsml=x*0.0 + 1.0

if not keyword_set(center) then center= [0.0, 0.0, 0.0]



;--------------------------------------
;  Set up Parameters
;--------------------------------------




    
    ; Smoothing Parameters
    ; ---------------------
    ; good for gas
    DesNgb=32L
    Hmax= 50.0
    ; good for stars
    ;DesNgb=16L
    ;Hmax= 5.0


    ; Plot ranges
    ; ------------------------------------
    ;MaxDensXY= 1.0e+3
    ;MaxDensXY= 1.0e+1
    ;MaxDensXY= 2.0
    ;MaxDensXY= 1.0
    MaxDensXY= 1.0e-1            ; (gadget units, 10^10 msun/ kpc^2 = 10^4 msun/pc^2)
    ;MaxDensXY= 5.0e-2
    ;MaxDensXY= 1.0e-2
    
    ;DynRange=1.0e2
    ;DynRange=1.0e3
    ;DynRange=5.0e3
    DynRange=1.0e4
    ;DynRange=1.0e5
    ;`DynRange=1.0e6
    ;DynRange=1.0e10
    ;DynRange=50
    
    ma=MaxDensXY
    mi=MaxDensXY/DynRange

    if keyword_set(set_maxden) then ma= set_maxden
    if keyword_set(set_dynrng) then mi= set_maxden/set_dynrng



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


    print, "changing center to= ",center


    ; OK, we got fields, prepare to send
    ; -----------------------------------
    Coord=fltarr(3,N)
    Masses=fltarr(N)
    Hsmls=fltarr(N)
    Density=fltarr(N)

    Coord(0,*)= x
    Coord(1,*)= y
    Coord(2,*)= z
    Masses(*)= m
    Hsmls(*)= hsml
    Density(*)=1.0/hsml

    ;xmin= -xlen
    ;xmax=  xlen
    ;ymin= -xlen
    ;ymax=  xlen
    xmin= -xlen+center[0]
    xmax=  xlen+center[0]
    ymin= -xlen+center[1]
    ymax=  xlen+center[1]


    ;xpixels=XYPIXELS
    ;ypixels=XYPIXELS

    ; default is xy axis
    ; ------------------
    Axis1=0L    ; select x-axis
    Axis2=1L    ; select y-axis


    if keyword_set(xz) then begin
        Axis1=0L    ; select x-axis
        Axis2=2L    ; select z-axis
	ymin= -xlen+center[2]
	ymax=  xlen+center[2]
    endif

    if keyword_set(yz) then begin
        Axis1=1L    ; select y-axis
        Axis2=2L    ; select z-axis
	xmin= -xlen+center[1]
	xmax=  xlen+center[1]
	ymin= -xlen+center[2]
	ymax=  xlen+center[2]
    endif

    Total_Mass= total(Masses)

    ; get mass within plotted region
    ; -------------------------------
    plotmasses= Masses
    nonplotpts= where((Coord(Axis1,*) gt xmax) or (Coord(Axis1,*) lt xmin))
    if nonplotpts(0) ne -1 then plotmasses(nonplotpts)= 0

    nonplotpts= where((Coord(Axis2,*) gt ymax) or (Coord(Axis2,*) lt ymin))
    if nonplotpts(0) ne -1 then plotmasses(nonplotpts)= 0

    Plotted_mass= total(plotmasses)

    print, Plotted_mass, " out of ", Total_Mass, " (", 100.*Plotted_mass/Total_Mass,"%) will make it into the plotted area"



; --------------------------------------------------
; Actually project the particles and smooth them
; --------------------------------------------------


if not keyword_set(crude) then begin


do_old= 0
do_brant= 0
do_new= 0
do_vslicer= 1

    ; OLD VERSION OF THIS CODE !!!
    ; -----------------------------------
    ; Call external c program to smooth
    ; -----------------------------------
    if do_old eq 1 then begin
	MassMap=fltarr(xpixels,ypixels)
	S = CALL_EXTERNAL('/n/home03/tcox/C-Routines_for_IDL/ComputeSmoothedProj2/adsmooth.so', $
                      'project_and_smooth', $
                      N, $
                      Coord,Masses,$
                      xmin,xmax,ymin,ymax,$
                      Xpixels,Ypixels,$
                      DesNgb, $
                      Axis1, $
                      Axis2, $
                      Hmax, $
                      MassMap  )
    endif


    ; Brant VERSION
    ; -----------------------------------
    if do_brant eq 1 then begin
	MassMap=fltarr(xpixels,ypixels)
	S = CALL_EXTERNAL('/n/home03/tcox/C-Routines_for_IDL/ComputeSmoothedProj/smooth.so', $
                      'smooth', $
                      N, $
                      Coord,Masses,$
		      Density,Hsmls,$
                      xmin,xmax,ymin,ymax,$
                      Xpixels,Ypixels,$
                      Axis1, $
                      Axis2, $
                      MassMap  )
    endif



    ;
    ; NEW VERSION
    ; -----------------------------------
    if do_new eq 1 then begin
        MassMap=fltarr(xpixels,ypixels)
        S = CALL_EXTERNAL('/n/home03/tcox/C-Routines_for_IDL/ComputeSmoothedProj3/allnsmooth.so', $
                      'project_and_smooth', $
                      N, $
                      Coord,Masses,$
                      Hsmls,$
                      xmin,xmax,ymin,ymax,$
                      Xpixels,Ypixels,$
                      Axis1, $
                      Axis2, $
                      MassMap  )
    endif



    ;
    ; new slicer program from Volker
    ; -----------------------------------
    if do_vslicer eq 1 then begin

	pos2D= fltarr(2,N)
	pos2D(0,*)= x
	pos2D(1,*)= y

	quantity= fltarr(N)
	temp= fload_gas_temperature(1)
	if n_elements(temp) eq N then quantity(*)= temp else quantity(*)= Masses

	ValueXY= fltarr(ypixels,xpixels)      ; holds the result
	ValueXY_Temp= fltarr(ypixels,xpixels) ; holds the result

	S = CALL_EXTERNAL('/n/home03/tcox/C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so', $
                      'slice', $
                      N, $
                      pos2D, Hsmls, Masses, quantity, $
                      xmin,xmax,ymin,ymax,$
                      Xpixels,Ypixels,$
                      ValueXY, ValueXY_Temp)

	MassMap = ValueXY
	ValueXY_Temp= ValueXY_Temp

    endif

endif



if keyword_set(crude) then begin

          nbins= XYPIXELS

	  if keyword_set(yz) then begin
		x= y
		y= z
	  endif

	  if keyword_set(xz) then y= z

	  nonein= 0
          idx= where((x le xmax) and (x ge xmin))
          if idx(0) ne -1 then begin 
		x= x(idx)
		y= y(idx)
	  endif else begin
		nonein= nonein+1
	  endelse

          idx= where((y le ymax) and (y ge ymin))
          if idx(0) ne -1 then begin
		x= x(idx)
		y= y(idx)
	  endif else begin
		nonein= nonein+1
	  endelse

          x=nbins*(x-xmin)/(xmax-xmin)
          x=fix(x)

          y=nbins*(y-ymin)/(ymax-ymin)
          y=fix(y)

	  if nonein eq 0 then begin
		phist= hist_2d(x, y, max1=nbins-1, max2=nbins-1)
	  endif else begin
		phist= fltarr(nbins,nbins)
	  endelse

          ; here phist has number in each cell

	  ; units 1
	  ; ---------
          ; convert to density
          ;pixelsize= (xmax-xmin)/nbins               ; kpc
          ;pixelarea = pixelsize*pixelsize            ; kpc^2
          ;if max(m) ne min(m) then print, "WARMING: m has a range of values"
          ;imass = max(m)                           ; 1e10 msolar
          ;MassMap = phist * imass / pixelarea           ; gadget units


	  ; units 2
	  ; ---------
          ; convert to luminosity
          avg_pt_lum= mean(m)
	  MassMap= phist * avg_pt_lum


          ; method 1
	  ; ---------
          ;width= fix(0.75*Hmax/pixelsize)
          ;newMap= smooth(MassMap, width, /edge_truncate)
	  ; this actually produces a fairly nice smooth map, the only
	  ; problem is that you lose some of the high resolution center 
	  ; portions

          ; method 2
	  ; ---------
	  ; this is supposed to fix the above resolution problem, but it
	  ; didn't really do the trick
	  ;hiMapnum = 0.45*max(MassMap)
	  ;addidx= where(MassMap gt hiMapnum)
	  ;Maptoadd= 0.0*MassMap
	  ;Maptoadd(addidx) = MassMap(addidx)
	  ;width= fix(0.25*Hmax/pixelsize)
	  ;hiMap= smooth(Maptoadd, width, /edge_truncate)

          ;MassMap= newMap+hiMap
	  ;MassMap= newMap

          ; method 2
          ; kernel= [1,2,4,4,10,4,4,2,1]   ; this needs to be 2 dimensions
          ; newMap= convol(MassMap, kernel, max(kernel))


endif



; ----------------
;   write fits
; ----------------

    if keyword_set(fitstoo) then begin
        n2= strlen(filename)-4                      ; assumed filename is *.eps
        fitsfilename= strmid(filename,0,n2)+'.fits'
        writefits, fitsfilename, MassMap
    endif




; -----------------------------------------------------
; now map the field logarithmically to a color scale
; -----------------------------------------------------

HalfMassContourLevel= -1.0

if keyword_set(compute_halfmasscnt) then begin

    ; determine half-mass surface density
    Pixel_Area= 1.0*(xmax-xmin)*(ymax-ymin)/(xpixels*xpixels)   ; kpc^2
    sortedMap= MassMap(sort(MassMap))
    Nm= n_elements(sortedMap)
    print, "Pixel_Area= ", Pixel_Area
    print, "Nm= ", Nm
    print, "Total Mass= ", Total_Mass
    print, "Plotted_mass= ", Plotted_mass
    print, "Pic*Area= ",total(MassMap)*Pixel_Area,"     (should give us the plotted_mass)"
    picerr= (total(MassMap)*Pixel_Area - Plotted_mass)/Plotted_mass
    print, "          (error: ",100.0*picerr,")"

    ;skiphalfmasscontour= 1
    skiphalfmasscontour= 0
    print, " ************************************** "
    print, "  skiphalfmasscontour= ", skiphalfmasscontour
    print, " ************************************** "

    i= 0L
    repeat begin
	enclosedmass= Pixel_Area*total(sortedMap[Nm-1-i:Nm-1])
;print, "WARNING: fixing enclosedmass by hand"
;enclosedmass= 10.0e+12
	i= i+1
	if ((i gt (Nm-2)) or (Plotted_mass lt (Total_Mass*0.5))) then enclosedmass=Total_Mass
	if skiphalfmasscontour eq 1 then enclosedmass=Total_Mass
    endrep until enclosedmass ge (Total_Mass*0.5)
    if (Nm-1-i-1 ge 0) then HalfMassContourLevel= sortedMap[Nm-1-i-1] else HalfMassContourLevel=sortedMap[Nm-2]
    print, "HalfMassContourLevel= ",HalfMassContourLevel

endif


; -----------------------------------------------------
; now map the field logarithmically to a color scale
; -----------------------------------------------------


    ; print ranges
    print, "MassMap     max: ", max(MassMap), "   min: ", min(MassMap)
    print, "Clipping at   ma= ", ma, " mi= ", mi


    ; do some clipping

    ind=where(MassMap lt mi)
    if ind(0) ne -1 then MassMap(ind)=mi
    ind=where(MassMap gt ma)
    if ind(0) ne -1 then MassMap(ind)=ma


    Cols=255 ; number of colors

    if keyword_set(var_already_inlog) then begin
	Pic=(MassMap-mi)/(ma/mi) * (cols-3) + 2
    endif else begin
	Pic=(alog10(MassMap)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2
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

    ; set background
    if not keyword_set(floatingbckgrnd) then begin
	if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
	;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black
    endif




;--------------------
; now do the same thing for the temperature image - if we have one.


if total(ValueXY_Temp) gt 0 then begin
    ind=where(ValueXY eq 0)
    if ind(0) ne -1 then ValueXy(ind)=1
    TempMap= ValueXY_Temp / ValueXY ; define mass-weighted temperature map
    
    ma= 5.0e6                 ; select maximum value (brightest color)
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


