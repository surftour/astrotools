pro contour_velmap, x, y, z, vx, vy, vz, m, xlen, $
			hsml=hsml, $
			pixels=pixels, zthickness=zthickness, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			center=center, $
			fitstoo=fitstoo, $
			VelxMap= VelxMap, $
			VelyMap= VelyMap



XYPIXELS= 32L

if keyword_set(pixels) then XYPIXELS= long(pixels)

NxNImage= fltarr(XYPIXELS, XYPIXELS)
NxNImage= NxNImage*0.0 + 1.0


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



;--------------------------------------
;  Set up Parameters
;--------------------------------------



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

        process_rotation, vx, vy, vz, rotate_theta, rotate_phi, vx_new, vy_new, vz_new

        vx= vx_new
        vy= vy_new
        vz= vz_new


endif
    
    
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
                vx= vx(idx)
                vy= vy(idx)
                vz= vz(idx)
        endif
    endif


    print, "changing center to= ",center


    ; OK, we got fields, prepare to send
    ; -----------------------------------
    N= long(n_elements(x))

    Coord=fltarr(3,N)
    Vel_x=fltarr(N)
    Vel_y=fltarr(N)
    Masses=fltarr(N)

    Coord(0,*)= x
    Coord(1,*)= y
    Coord(2,*)= z

    Vel_x(*)= vx
    Vel_y(*)= vy

    Masses(*)= m

    DesNgb= 32L
    Hmax= 10.0

    Axis1= 0L
    Axis2= 1L

    xmin= -xlen
    xmax=  xlen
    ymin= -xlen
    ymax=  xlen

    Xpixels= XYPIXELS
    Ypixels= XYPIXELS

    Junk= fltarr(Xpixels,Ypixels)
    VelxMap=fltarr(Xpixels,Ypixels)
    VelyMap=fltarr(Xpixels,Ypixels)

    S = CALL_EXTERNAL('Tools/C-Routines_for_IDL/VelocityField/velfield.so', $
         ;'velfield', $
	 'project_and_velfield', $
         N, $
         Coord, $
	 Vel_x, Vel_y, $
	 Masses, $
         xmin,xmax,ymin,ymax,$
         Xpixels,Ypixels,$
	 DesNgb, $
	 Axis1, $
	 Axis2, $
	 Hmax, $
         Junk, VelxMap, VelyMap) 



; ----------------
;   write fits
; ----------------

    if keyword_set(fitstoo) then begin
        n2= strlen(filename)-4                      ; assumed filename is *.eps
        fitsfilename= strmid(filename,0,n2)+'.fits'
        writefits, fitsfilename, Map
    endif



; -------------
;  Done
; -------------



end


