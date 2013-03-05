pro img_makepic_raw, x, y, z, m, xlen, xz=xz, yz=yz, $
			hsml=hsml, $
			pixels=pixels, $
			crude=crude, $
			NxNImage=NxNImage, $
			NxNTempImage=NxNTempImage, $
			Pixel_Area=Pixel_Area, $
			Total_Mass=Total_Mass, $
			Plotted_Mass=Plotted_Mass


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
;if keyword_set(pixels) then XYPIXELS= long(pixels) else XYPIXELS= 1024L
;if keyword_set(pixels) then XYPIXELS= long(pixels) else XYPIXELS= 4096L

NxNImage= fltarr(XYPIXELS, XYPIXELS)
NxNImage= NxNImage*0.0 + 1.0

NxNTempImage= fltarr(XYPIXELS, XYPIXELS)
NxNTempImage= NxNImage*0.0 + 1.0
ValueXY_Temp= -1


if not keyword_set(xlen) then xlen=100.0
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


    
    ; Smoothing Parameters
    ; ---------------------
    ; good for gas
    DesNgb=32L
    Hmax= 50.0
    ; good for stars
    ;DesNgb=16L
    ;Hmax= 5.0



    N= long(n_elements(x))



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


    xpixels=XYPIXELS
    ypixels=XYPIXELS

    ; default is xy axis
    ; ------------------
    Axis1= long(0)    ; select x-axis
    Axis2= long(1)    ; select y-axis
    Axis3= long(2)    ; select z-axis


    if keyword_set(xz) then begin
        Axis1=0L    ; select x-axis
        Axis2=2L    ; select z-axis
        Axis3=1L    ; select y-axis
	ymin= -xlen+center[2]
	ymax=  xlen+center[2]
    endif

    if keyword_set(yz) then begin
        Axis1=1L    ; select y-axis
        Axis2=2L    ; select z-axis
        Axis3=0L    ; select x-axis
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

    Plotted_Mass= total(plotmasses)

    print, Plotted_Mass, " out of ", Total_Mass, " (", 100.*Plotted_Mass/Total_Mass,"%) will make it into the plotted area"



; --------------------------------------------------
; Actually project the particles and smooth them
; --------------------------------------------------


if not keyword_set(crude) then begin


do_old= 0
do_brant= 0
do_new= 1
do_vslicer= 0
do_nohsml= 0


    ; OLD VERSION OF THIS CODE !!!
    ; -----------------------------------
    ; Call external c program to smooth
    ; -----------------------------------
    if do_old eq 1 then begin
	MassMap=fltarr(xpixels,ypixels)
	spawn, 'echo $TJHOME', result
	homedir= strcompress(result,/remove_all)
	libfile= homedir+'/Tools/C-Routines_for_IDL/ComputeSmoothedProj2/adsmooth.so'
	S = CALL_EXTERNAL(libfile, $
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
	spawn, 'echo $TJHOME', result
	homedir= strcompress(result,/remove_all)
	libfile= homedir+'/Tools/C-Routines_for_IDL/ComputeSmoothedProj/smooth.so'
	S = CALL_EXTERNAL(libfile, $
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
	spawn, 'echo $TJHOME', result
	homedir= strcompress(result,/remove_all)
	libfile= homedir+'/Tools/C-Routines_for_IDL/ComputeSmoothedProj3/allnsmooth.so'
	S = CALL_EXTERNAL(libfile, $
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
	if fload_npart(0) gt 0 then quantity(*)= fload_gas_temperature(1) else quantity(*)= Masses

	MassMap= fltarr(ypixels,xpixels)      ; holds the result
	ValueXY_Temp= fltarr(ypixels,xpixels) ; holds the result

	spawn, 'echo $TJHOME', result
	homedir= strcompress(result,/remove_all)
	libfile= homedir+'/Tools/C-Routines_for_IDL/Compute2dProjectionHsmlGiven/slicer.so'
	S = CALL_EXTERNAL(libfile, $
                      'slice', $
                      N, $
                      pos2D, Hsmls, Masses, quantity, $
                      xmin,xmax,ymin,ymax,$
                      Xpixels,Ypixels,$
                      MassMap, ValueXY_Temp)

    endif



    ;
    ; good if you don't have smoothing length
    ; this will calculate hsml and project
    ; -----------------------------------------
    if do_nohsml eq 1 then begin

	BoxSize= double(0.0)
	Hmax=    float( (xmax-xmin)/xpixels * 50) ; maximum used smoothing length

	quantity= fltarr(N)
	if fload_npart(0) gt 0 then quantity(*)= fload_gas_temperature(1) else quantity(*)= Masses

	MassMap= fltarr(ypixels,xpixels)      ; holds the result
	ValueXY_Temp= fltarr(ypixels,xpixels) ; holds the result

	spawn, 'echo $TJHOME', result
	homedir= strcompress(result,/remove_all)
	libfile= homedir+'/Tools/C-Routines_for_IDL/ComputeHsmlAndProject/HsmlAndProject.so'
	S = CALL_EXTERNAL(libfile, $
                     'findHsmlAndProject', $
                     N, $
                     Coord, Hsmls, Masses, $
                     quantity, $
                     Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, $
                     XPixels, YPixels, DesNgb, $
                     Axis1, Axis2, Axis3, $
                     Hmax,  $
                     BoxSize, $
                     MassMap, $
                     ValueXY_Temp)
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


    ; determine some diagnostics
    Pixel_Area= 1.0*(xmax-xmin)*(ymax-ymin)/(xpixels*xpixels)   ; kpc^2

    Nm= n_elements(MassMap)
    print, "Pixel_Area= ", Pixel_Area
    print, "Nm= ", Nm
    print, "Total Mass= ", Total_Mass
    print, "Plotted_Mass= ", Plotted_Mass
    print, "Pic*Area= ",total(MassMap)*Pixel_Area,"     (should give us the Plotted_Mass)"
    picerr= (total(MassMap)*Pixel_Area - Plotted_Mass)/Plotted_Mass
    print, "          (error: ",100.0*picerr,")"




; -----------------------------
;   set image and we're done
; -----------------------------

    NxNImage= MassMap

    NxNTempImage= ValueXY_Temp

; -------------
;  Done
; -------------



end


