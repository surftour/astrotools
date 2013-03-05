pro contour_makepic_total, x, y, z, m, xlen, xz=xz, yz=yz, $
			xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
			filename=filename, fitstoo=fitstoo, $
			pixels=pixels, zthickness=zthickness, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			crude=crude, center=center, $
			NxNImage=NxNImage, set_maxden=set_maxden, set_dynrng=set_dynrng, $
			var_already_inlog=var_already_inlog


;  This is specifically used for contour_xrays, as this totals, rather than takes the surface
;  density.
;
;
;
;
;

;XYPIXELS= 500L
;XYPIXELS= 480L
XYPIXELS= 120L

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
    DesNgb=96L
    Hmax= 4.0


    ; Plot ranges
    ; ------------------------------------
    ;MaxDensXY=  max(Map)
    ;MaxDensXY= 2.0
    ;MaxDensXY= 1.0e-1            ; (gadget units, 10^10 msun/ kpc^2 = 10^4 msun/pc^2)
    MaxDensXY= 5.0e-2
    ;MaxDensXY= 1.0e-2
    
    ;DynRange=1.0e2
    DynRange=1.0e3
    ;DynRange=5.0e3
    ;DynRange=1.0e4
    ;DynRange=1.0e6
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

	print, "rotate: by (theta,phi)=", rotate_theta,rotate_phi
	; transform to radians
	theta= !PI*rotate_theta/180.0
	phi= !PI*rotate_phi/180.0

        ; rotate the coordinates   - actually we're going to do this in contour_makepic
        ; ------------------------
        ; this transformation is 1st rotating around the y-axis such
        ; that +z -> +x (theta), and then 2nd rotating around the z-axis such
        ; that +y -> +x (phi).
        ;
        ;x_new= x*(cos(theta)*cos(phi)) - y*sin(phi) + z*(cos(phi)*sin(theta))
        ;y_new= x*(cos(theta)*sin(phi)) + y*cos(phi) + z*(sin(phi)*sin(theta))
        ;z_new= -x*sin(theta) + z*cos(theta)

        ; alternatively, this is the inverse (or not exactly) rotation of the above, 
        ;    where we rotate around z-axis (by phi) and then y-axis (by theta)
        ;x_new= x*(cos(theta)*cos(phi)) + y*(cos(theta)*sin(phi)) + z*sin(theta)
        ;y_new= -x*sin(theta) + y*cos(theta)
        ;z_new= x*(sin(theta)*cos(phi)) - y*(sin(theta)*sin(phi)) + z*cos(theta)

	; first around x axis (theta)
	; then around z axis (phi)
	;x_new= x*cos(phi) + y*(sin(phi)*cos(theta)) + z*(sin(phi)*sin(theta))
        ;y_new= -x*sin(phi) + y*(cos(phi)*cos(theta)) + z*(cos(phi)*sin(theta))
        ;z_new= -y*sin(theta) + z*cos(theta)

	; first around y axis (theta)
	; then around x axis (phi)
	x_new= x*cos(theta) + z*sin(theta)
        y_new= -x*(sin(phi)*sin(theta)) + y*cos(phi) + z*(cos(theta)*sin(phi))
        z_new= -x*(cos(phi)*sin(theta)) - y*sin(phi) + z*(cos(theta)*cos(phi))


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

    Coord(0,*)= x
    Coord(1,*)= y
    Coord(2,*)= z
    Masses(*)= m

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








; --------------------------------------------------
; Actually project the particles and smooth them
; --------------------------------------------------

if not keyword_set(crude) then begin

    ; -----------------------------------
    ; Call external c program to smooth
    ; -----------------------------------
    Map=fltarr(xpixels,ypixels)
    S = CALL_EXTERNAL('/home2/tcox/Tools/C-Routines_for_IDL/ComputeSmoothedTotal/smtotal.so', $
                      'project_smooth_and_total', $
                      N, $
                      Coord,Masses,$
                      xmin,xmax,ymin,ymax,$
                      Xpixels,Ypixels,$
                      DesNgb, $
                      Axis1, $
                      Axis2, $
                      Hmax, $
                      Map  )


endif else begin


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
          ; convert to density
          pixelsize= (xmax-xmin)/nbins               ; kpc
          pixelarea = pixelsize*pixelsize            ; kpc^2
          if max(m) ne min(m) then print, "WARMING: m has a range of values"
          imass = max(m)                           ; 1e10 msolar
          Map = phist * imass / pixelarea           ; gadget units


          ; method 1
          width= fix(0.75*Hmax/pixelsize)
          newMap= smooth(Map, width, /edge_truncate)
	  ; this actually produces a fairly nice smooth map, the only
	  ; problem is that you lose some of the high resolution center 
	  ; portions

	  ; this is supposed to fix the above resolution problem, but it
	  ; didn't really do the trick
	  ;hiMapnum = 0.45*max(Map)
	  ;addidx= where(Map gt hiMapnum)
	  ;Maptoadd= 0.0*Map
	  ;Maptoadd(addidx) = Map(addidx)
	  ;width= fix(0.25*Hmax/pixelsize)
	  ;hiMap= smooth(Maptoadd, width, /edge_truncate)

          ;Map= newMap+hiMap
	  Map= newMap

          ; method 2
          ; kernel= [1,2,4,4,10,4,4,2,1]   ; this needs to be 2 dimensions
          ; newMap= convol(Map, kernel, max(kernel))


endelse



; ----------------
;   write fits
; ----------------

    if keyword_set(fitstoo) then begin
        n2= strlen(filename)-4                      ; assumed filename is *.eps
        fitsfilename= strmid(filename,0,n2)+'.fits'
        writefits, fitsfilename, Map
    endif




; -----------------------------------------------------
; now map the field logarithmically to a color scale
; -----------------------------------------------------


    ; print ranges
    print, "Map     max: ", max(Map), "   min: ", min(Map)
    print, "Clipping at   ma= ", ma, " mi= ", mi


    ; do some clipping

    ind=where(Map lt mi)
    if ind(0) ne -1 then Map(ind)=mi
    ind=where(Map gt ma)
    if ind(0) ne -1 then Map(ind)=ma


    Cols=255 ; number of colors

    if keyword_set(var_already_inlog) then begin
	Pic=(Map-mi)/(ma/mi) * (cols-3) + 2
    endif else begin
	Pic=(alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2
    endelse

    ind=where(Pic ge 256)
    if ind(0) ne -1 then Pic(ind)=255

    ind=where(Pic eq 0)
    if ind(0) ne -1 then Pic(ind)=1


    invertcolors= 1
    if invertcolors eq 1 then begin
	Pic=256-Pic			     ; invert color table
	idx= where(Pic EQ 254)               ; set background to white
    endif else begin
	idx= where(Pic EQ 2)
    endelse

    if idx(0) ne -1 then Pic(idx)= 1       ; this is setting it to white
    ;if idx(0) ne -1 then Pic(idx)= 0       ; this sets background to black




; -----------------------------
;   set image and we're done
; -----------------------------

    NxNImage= Pic


; -------------
;  Done
; -------------



end


