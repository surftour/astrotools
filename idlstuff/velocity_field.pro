pro velocity_field, frun, snapnum, xlen, sendto, $
		s_id=s_id, s_npart=s_npart, $
		xz=xz, yz=yz, zthickness=zthickness, xthickness=xthickness, ythickness=ythickness, $
		h=h, $
		filename=filename


particles_set= 1
;if keyword_set(s_id) then particles_set= 1
;if keyword_set(s_npart) then particles_set= 1


sendto= 'ps'


if not keyword_set(sendto) then sendto='ps'
if not keyword_set(snapnum) then snapnum= 0
if not keyword_set(xlen) then xlen= 20.0
if ((not keyword_set(frun)) or (not particles_set)) then begin
   print, "  "
   print, "velocity_field, frun, snapnum, xlen, sendto"
;   print, "       s_id=s_id, s_npart=s_npart"
   print, "        /xz, /yz, zthickness=zthickness, xthickness=xthickness, ythickness=ythickness,"
   print, "        filename=filename"
   print, "  "
   return
endif


	if not keyword_set(filename) then begin
            filename= 'velfield.eps'
            ans= ''
            read, ans, PROMPT='eps filename ['+filename+']:'
            if strlen(ans) GT 0 then filename= ans
	endif

        set_plot, 'ps'
        device, filename= filename, /encapsulated,/color,bits_per_pixel=8
        device, SET_CHARACTER_SIZE=[200,300], xsize=13, ysize=12
	;device, SET_FONT='Helvetica-Bold'
	;device, /helvetica, /bold, font_index=4

        ;loadct, 4    ;colors black - yellow, fully spectrum
	;loadct, 0    ;black/white
	loadct, 1    ;blue/white
	;loadct, 3    ;red/white temperature scale
        tvlct,r,g,b,/get
        v1=[0,255]
        v2=[0,255]
        v3=[0,255]
        tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white








;--------------------------------------
;  Set up Parameters
;--------------------------------------

if keyword_set(h) then begin
	h = fload_cosmology('h')
endif else begin
	h = 1.0
endelse



min_velocity= 0.0
;min_velocity= 200.0


nolabels= 1
;nolabels= 0

;draw_a_key= 1
draw_a_key= 0



;----------------------------------
;----------------------------------



; --------------------------------
;  Load variables for smoothing
; --------------------------------
;    if (fload_snapshot(frun, snapnum)) then begin
    if (fload_snapshot_bh(frun, snapnum)) then begin
	print, "PROBLEM: opening file"
	return
    endif


merger= -1

    if merger EQ -1 then begin
        npart= fload_npart(0)                        ; gas N
        N= long(npart)                               ;
    
        Coord=fltarr(3,N)
        Vel_x=fltarr(N)
        Vel_y=fltarr(N)
        Masses=fltarr(N)

        x= fload_gas_xyz('x')
        y= fload_gas_xyz('y')
        z= fload_gas_xyz('z')
        vx= fload_gas_v('x')
        vy= fload_gas_v('y')
        vz= fload_gas_v('z')
        m= fload_gas_mass(1)
    endif


    if merger EQ 0 then begin
	npart= fload_npart(99)                        ; return array
	N= long(npart(0)+npart(2)+npart(3)+npart(4))  ; baryonic N

	Coord=fltarr(3,N)
	Vel_x=fltarr(N)
	Vel_y=fltarr(N)
	Masses=fltarr(N)

	x= fload_baryon_xyz('x')
	y= fload_baryon_xyz('y')
	z= fload_baryon_xyz('z')
	vx= fload_baryon_v('x')
	vy= fload_baryon_v('y')
	vz= fload_baryon_v('z')
	m= fload_baryon_mass(1)
    endif

    if merger EQ 1 then begin
	gid= 1
        gnpart= 20000
	ngas= fload_1gal_gas_npart(gid,gnpart)
        N= long(ngas)

        Coord=fltarr(3,N)
        Vel_x=fltarr(N)
        Vel_y=fltarr(N)
        Masses=fltarr(N)

	center= fload_1gal_center(1,70000)
	comvel= fload_1gal_comvel(1,70000)
        x= fload_1gal_gas_xyz('x',gid,gnpart,center=center)
        y= fload_1gal_gas_xyz('y',gid,gnpart,center=center)
        z= fload_1gal_gas_xyz('z',gid,gnpart,center=center)
        vx= fload_1gal_gas_v('x',gid,gnpart)-comvel(0)
        vy= fload_1gal_gas_v('y',gid,gnpart)-comvel(1)
	vz= fload_1gal_gas_v('z',gid,gnpart)-comvel(2)
        m= fload_1gal_gas_mass(gid,gnpart)
;stop
    endif


    ; deal with h explicitly
    x= x/h
    y= y/h
    z= z/h
    m= m/h

    ; crop in the z direction
    if keyword_set(zthickness) then begin
	idx= where(abs(z) LE zthickness)
	if idx(0) GT -1 then begin
		N= n_elements(idx)
		x= x(idx)
		y= y(idx)
		z= z(idx)
		vx= vx(idx)
		vy= vy(idx)
		vz= vz(idx)
		m= m(idx)
	endif
    endif

    ; crop in the x direction
    if keyword_set(xthickness) then begin
        idx= where(abs(x) LE xthickness)
        if idx(0) GT -1 then begin
                N= n_elements(idx)
                x= x(idx)
                y= y(idx)
                z= z(idx)
                vx= vx(idx)
                vy= vy(idx)
		vz= vz(idx)
                m= m(idx)
        endif
    endif

    ; crop in the y direction
    if keyword_set(ythickness) then begin
        idx= where(abs(y) LE ythickness)
        if idx(0) GT -1 then begin
                N= n_elements(idx)
                x= x(idx)
                y= y(idx)
                z= z(idx)
                vx= vx(idx)
                vy= vy(idx)
		vz= vz(idx)
                m= m(idx)
        endif
    endif

    ; OK we got all field set, send it away
    ; -------------------------------------
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

    xmin= -xlen
    xmax=  xlen
    ymin= -xlen
    ymax=  xlen
    
    ;XYPIXELS= 500L
    XYPIXELS= 480L

    xpixels=XYPIXELS
    ypixels=XYPIXELS

    ; default is xy axis
    ; ------------------
    Axis1=0L    ; select x-axis
    Axis2=1L    ; select y-axis
    if h lt 1.0 then xtit="!8x!3 (kpc)" else xtit="!8x!3 (kpc/h)"
    if h lt 1.0 then xtit="!8y!3 (kpc)" else ytit="!8y!3 (kpc/h)"


    if keyword_set(xz) then begin
	Axis1=0L    ; select x-axis
	Axis2=2L    ; select z-axis
	Vel_y(*)= vz
	;if h lt 1.0 then xtit="!8x!3 (kpc)" else xtit="!8x!3 (kpc/h)"
	;if h lt 1.0 then ytit="!8z!3 (kpc)" else ytit="!8z!3 (kpc/h)"
	xtit="x (kpc)"
	ytit="z (kpc)"
    endif 

    if keyword_set(yz) then begin
	Axis1=1L    ; select y-axis
	Axis2=2L    ; select z-axis
	Vel_x(*)= vy
	Vel_y(*)= vz
        if h lt 1.0 then xtit="!8y!3 (kpc)" else xtit="!8y!3 (kpc/h)"
        if h lt 1.0 then ytit="!8z!3 (kpc)" else ytit="!8z!3 (kpc/h)"
    endif


    DesNgb=30L
    ;Hmax= 1.0 / h  
    ;Hmax= 3.0 / h  
    Hmax= 25.0 / h  



; -----------------------------------
; Call external c program to smooth
; -----------------------------------
    Map=fltarr(xpixels,ypixels)
    VelxMap=fltarr(xpixels,ypixels)
    VelyMap=fltarr(xpixels,ypixels)

    spawn, 'echo $HOME', result
    homedir= strcompress(result,/remove_all)
    libfile= homedir+'/Tools/C-Routines_for_IDL/VelocityField/velfield.so'
    S = CALL_EXTERNAL(libfile, $
                      'project_and_velfield', $
                      N, $
                      Coord,Vel_x,Vel_y,Masses,$
                      xmin,xmax,ymin,ymax,$
                      Xpixels,Ypixels,$
                      DesNgb, $
                      Axis1, $
                      Axis2, $
                      Hmax, $
                      Map, VelxMap, VelyMap  )


; -----------------------------------------------------
; now map the field logarithmically to a color scale
; -----------------------------------------------------

;    MaxDensXY=  max(Map)
    MaxDensXY= 5.0e-1            ; (gadget units, 10^10 msun/ kpc^2 = 10^4 msun/pc^2)
;    MaxDensXY= 1

    DynRange=1.0e6
    ;DynRange=1.0e5
    ;DynRange=2.0e4
    ;DynRange=1.0e3
    ;DynRange=50

;    ma=MaxDensXY/1.2
    ma=MaxDensXY
    mi=MaxDensXY/DynRange
;    mi=  min(Map)*12


    ; print ranges
    print, "Map     max: ", max(Map), "   min: ", min(Map)
    print, "Clipping at   ma= ", ma, " mi= ", mi
    
   
    ; do some clipping
  
    ind=where(Map lt mi)
    if ind(0) ne -1 then Map(ind)=mi
    ind=where(Map gt ma)
    if ind(0) ne -1 then Map(ind)=ma

	
    Cols=255 ; number of colors

    Pic=(alog10(Map)-alog10(mi))/(alog10(ma/mi)) * (cols-3) + 2

    ind=where(Pic ge 256)
    if ind(0) ne -1 then Pic(ind)=255

    ind=where(Pic eq 0)
    if ind(0) ne -1 then Pic(ind)=1
    
    Pic=256-Pic





; -----------------------------
; generate a plot
; -----------------------------


	; if we want labels and such
	; -----------------------------
	if nolabels eq 0 then begin
	    x0= 0.18
	    x1= 0.95
	    y0= 0.15
	    y1= 0.95


	    !p.position=[x0,y0,x1,y1]
	    !p.ticklen=0.03
	    !p.font= 0

	    tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

		; creates axes and plot style
		; ---------------------------
		plot,[0],[0], psym=3,  $
		      xrange=[xmin,xmax], $
		      yrange=[ymin,ymax], $
		      /NOERASE, $
		      color= 0, $
		      xstyle=1, ystyle=1, $
		      charsize=1.5, charthick= 3.0, $
		      xthick=4.0, ythick=4.0, $
		      xtitle=xtit, $
		      ytitle=ytit, /normal, $
		      /nodata

	endif 



	; no labels, please
	; -------------------
	if nolabels eq 1 then begin
	    x0= 0.02
            x1= 0.98
            y0= 0.02
            y1= 0.98


            !p.position=[x0,y0,x1,y1]
            !p.ticklen=0.03
            !p.font= 0

            tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal

                ; creates axes and plot style
                ; ---------------------------
                plot,[0],[0], psym=3,  $
                      xrange=[xmin,xmax], $
                      yrange=[ymin,ymax], $
                      /NOERASE, $
                      color= 0, $
                      xstyle=1, ystyle=1, $
		      charthick= 0.01, $
                      xcharsize=0.01, ycharsize= 0.01, $
                      xthick=4.0, ythick=4.0, $
		      xtickformat='(a1)', ytickformat='(a1)', $
		      /normal, $
                      /nodata
	endif



	; add velocity field
	; ------------------------
	;gridpts= XYPIXELS/50
	gridpts= 20
	normlen= 2.0

	pix_int= fix(XYPIXELS/(gridpts+1))
	vx= fltarr(2)
	vy= fltarr(2)
;	VNorm= max(sqrt(VelXMap*VelXMap + VelyMap*VelyMap))/normlen
;	VNorm= 20.0        ; normalize to 200  (then divide by 10, so 1 unit is 10 km/s on graph
	VNorm= 100.0

	for i=0,gridpts-1 do begin
	  for j=0,gridpts-1 do begin

		; physical location
		vx[0]= (xmax-xmin)*(i+0.5)/gridpts + xmin
		vy[0]= (ymax-ymin)*(j+0.5)/gridpts + ymin

		; add on velocity at that point
		vx[1]= vx[0] + VelxMap[(i+1)*pix_int,(j+1)*pix_int]/VNorm
		vy[1]= vy[0] + VelyMap[(i+1)*pix_int,(j+1)*pix_int]/VNorm

		vellen= sqrt((vx[1]-vx[0])*(vx[1]-vx[0]) + (vy[1]-vy[0])*(vy[1]-vy[0]))
		;if vellen gt 0.0 then begin
		if vellen gt (min_velocity/VNorm) then begin
		    ;arrow, vx[0], vy[0], vx[1], vy[1], /data, color= 1
		    arrow, vx[0], vy[0], vx[1], vy[1], /data, color= 0
		endif


	  endfor
	endfor


	; -----------------
	;  Draw a key
	; -----------------
	if draw_a_key eq 1 then begin
	   ; put key in white box
	   ;polyfill, [-55,-55,-23,-23], [-55,-33,-33,-55], color= 1
	   ;oplot, [-55,-55,-23,-23,-55], [-55,-33,-33,-55,-55], color= 0, thick=3.0
	   ;bx0=x0+0.02
	   ;bx1=x0+0.22
	   ;by0=y0+0.72
	   ;by1=y0+0.80

	   ; big images.
	   bx0=xmin+3.0
	   bx1=xmin+24.5
	   by0=ymax-21.0
	   by1=ymax-7.0
	   polyfill, [bx0,bx0,bx1,bx1], [by0,by1,by1,by0], color= 1
	   oplot, [bx0,bx0,bx1,bx1,bx0], [by0,by1,by1,by0,by0] , color= 0, thick=3.0

	   ; add a little key
	   keylen= 200.0
	   keylen= keylen/VNorm

	   ; lower left corner
	   ;x0= xmin+6.0
	   ;y0= ymin+25.0
	   ; upper left corner
	   ;x0= xmin+7.0
	   ;y0= ymax-8.0
	   ; middle left
	   ;x0= xmin+1.5
	   ;y0= ymax-12.0

	   arrow, bx0+1.0, by1-2.0, bx0+1.0+keylen, by1-2.0, /data, thick=3.0, hthick=2.0, color= 0
	   xyouts, bx0+1.0, by1-7.0, '200 km s!E-1!N', /data, size= 1.3, charthick=4.0, color= 0
	   xyouts, bx0+1.0, by1-12.0, fload_timelbl(h,2), /data, size= 1.3, charthick=4.0, color= 0
	endif




;axis,xaxis=0,xstyle=1,ystyle=1,color=1,charsize=0.001    ,xthick=4.0,ythick=4.0  
;axis,yaxis=0,xstyle=1,ystyle=1,color=1,charsize=0.001    ,xthick=4.0,ythick=4.0   
;axis,xaxis=1,xstyle=1,ystyle=1,color=1,charsize=0.001    ,xthick=4.0,ythick=4.0  
;axis,yaxis=1,xstyle=1,ystyle=1,color=1,charsize=0.001    ,xthick=4.0,ythick=4.0   

 

	; grab this and convert to GIF
	; -----------------------------
;	yesno= 'yes'
	ans= ''
;	read, ans, PROMPT='would you like to write a gif file of this? '
	if ans EQ 'yes' then begin
	        filename= 'velfield.gif'
	        ans= ''
	        read, ans, PROMPT='gif filename ['+filename+']:'
	        if strlen(ans) GT 0 then filename= ans

		image= tvrd()

		pic=byte(image)

		write_gif, filename, pic, r, g, b
	endif









; -------------
;  Done
; -------------



if (sendto EQ 'ps') then device, /close


end




