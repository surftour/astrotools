
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------
;
;     General Image of the MW-M31 System
;     -------------------------------------------
;
;------------------------------------------------------------------------
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
;------------------------------------------------------------------------

pro showit, junk


if not keyword_set(junk) then begin
	print, "  "
	print, "  showit, junk"
	print, "  "
	print, "  "
	return
endif


filename='lgshowit.eps'


boxsize= 1.4  ; in Mpc
xlen= boxsize/2.0


;  mark today's position
;
t_sep= 0.714 ; in Mpc




; ----------------------------------
;   Try this plot thingy
; ----------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=16.0, newysize=16.0





; --------------------------------
;  Make Sphere
; --------------------------------

draw_sphere= 0
;draw_sphere= 1
if draw_sphere eq 1 then begin
	; Create an empty, 3-D array: 
	;SPHERE = FLTARR(20, 20, 20) 
	SPHERE = FLTARR(8, 8, 8) 
 
	; Create the spherical dataset: 
	;FOR X=0,19 DO FOR Y=0,19 DO FOR Z=0,19 DO $ 
	;   SPHERE(X, Y, Z) = SQRT((X-10)^2 + (Y-10)^2 + (Z-10)^2) 
	FOR X=0,7 DO FOR Y=0,7 DO FOR Z=0,7 DO $ 
	   SPHERE(X, Y, Z) = SQRT((X-4)^2 + (Y-4)^2 + (Z-4)^2) 
 
	; Find the vertices and polygons for a density level of 8: 
	;SHADE_VOLUME, SPHERE, 8, V, P 
	SHADE_VOLUME, SPHERE, 3, V, P 
 
	; Set up an appropriate 3-D transformation so we can see the 
	; sphere. This step is very important: 
	;SCALE3, XRANGE=[0,20], YRANGE=[0,20], ZRANGE=[0,20] 
	SCALE3, XRANGE=[0,8], YRANGE=[0,8], ZRANGE=[0,8] 
 
	; Render the image. Note that the T3D keyword has been set so that  
	; the previously-established 3-D transformation is used: 
	image = POLYSHADE(V, P, /T3D) 
 
	; Display the image: 
	TV, image
endif







; -----------------------------------------
;
; Set the 3D coordinate space with axes.
;
; -----------------------------------------

!p.ticklen=0.0

surface, dist(5), /nodata, /save, $
	xrange=[-xlen,xlen], yrange=[-xlen,xlen], zrange=[-xlen,xlen], $
	xstyle= 1, ystyle= 1, zstyle= 1, charsize= 1.5, color= 0, $
	xtickformat='(a1)', ytickformat='(a1)', ztickformat='(a1)', $
	position=[0.02, 0.02, 0.98, 0.98, 0.02, 0.98]
	;AX = 70, AZ = 25, $





;---------------------------------
;  Print the back walls
;---------------------------------

;wall_x= [-0.1, 0.0, 0.1]
;wall_y= [ -0.1, 0.0, 0.1]
wall_x= [-0.1, 0.0, 0.1]
wall_y= [ -0.15, 0.05, 0.15]
;wall_z= [[ 0.2, 0.1, 0.1], $
;	 [ 0.3, 0.0, 0.1], $
;	 [ 0.1, 0.0, -0.1]]
;wall_z= [[ 0.5, 0.5, 0.5], $
;	 [ 0.4, 0.4, 0.4], $
;	 [ 0.3, 0.3, 0.3]]
wall_z= [[ 0.5, 0.5, 0.5], $
	 [ 0.5, 0.5, 0.5], $
	 [ 0.5, 0.5, 0.5]]

;shade_surf, wall_z, wall_x, wall_y, color= 200, /noerase, /t3d, $
;surface, wall_z, wall_x, wall_y, color= 200, /noerase, /t3d, $
	;xrange=[-xlen,xlen], yrange=[-xlen,xlen], zrange=[-xlen,xlen], $
        ;xstyle= 1, ystyle= 1, zstyle= 1, charsize= 1.5, $
        ;xtickformat='(a1)', ytickformat='(a1)', ztickformat='(a1)', $
        ;position=[0.02, 0.02, 0.98, 0.98, 0.02, 0.98]
	;AX = 70, AZ = 25, $
        ;position=[0.02, 0.02, 0.98, 0.98, 0.02, 0.98]




;---------------------------------
;  Points where M31 and MW are
;---------------------------------

; mw
plots, [0.307550], [0.179521], [0.00000], psym= 4, color= 150, symsize= 3.0, /t3d

; m31
plots, [-0.495632], [-0.289306], [0.00000], psym= 4, color= 100, symsize= 3.0, /t3d






;-----------------------------
;  Box
;-----------------------------

; front sides

plots, [-xlen,-xlen], [ xlen,-xlen], [ xlen, xlen], psym=-3, thick= 5.0, color= 0, /t3d
plots, [-xlen,-xlen], [ xlen, xlen], [ xlen,-xlen], psym=-3, thick= 5.0, color= 0, /t3d
plots, [-xlen,-xlen], [-xlen,-xlen], [ xlen,-xlen], psym=-3, thick= 5.0, color= 0, /t3d
plots, [-xlen,-xlen], [ xlen,-xlen], [-xlen,-xlen], psym=-3, thick= 5.0, color= 0, /t3d

plots, [-xlen, xlen], [-xlen,-xlen], [ xlen, xlen], psym=-3, thick= 5.0, color= 0, /t3d
plots, [ xlen, xlen], [-xlen,-xlen], [ xlen,-xlen], psym=-3, thick= 5.0, color= 0, /t3d
plots, [-xlen, xlen], [-xlen,-xlen], [-xlen,-xlen], psym=-3, thick= 5.0, color= 0, /t3d

plots, [ xlen,-xlen], [ xlen, xlen], [ xlen, xlen], psym=-3, thick= 5.0, color= 0, /t3d
plots, [ xlen, xlen], [ xlen,-xlen], [ xlen, xlen], psym=-3, thick= 5.0, color= 0, /t3d


; back sides

plots, [-xlen, xlen], [ xlen, xlen], [-xlen,-xlen], psym=-3, thick= 3.0, color= 0, /t3d, linestyle= 2
plots, [ xlen, xlen], [ xlen, xlen], [ xlen,-xlen], psym=-3, thick= 3.0, color= 0, /t3d, linestyle= 2
plots, [ xlen, xlen], [-xlen, xlen], [-xlen,-xlen], psym=-3, thick= 3.0, color= 0, /t3d, linestyle= 2






;-----------------------------
;
; Mock up the MW disk
;
;-----------------------------
;draw_mw_disk=1
draw_mw_disk=0
if draw_mw_disk eq 1 then begin

	mw_center= [0.307550, 0.179521, 0.00000]
	rd= 0.0025
	hd= 0.0005
   
	;phi = dindgen(101)*2D*!dpi/100
	phi = dindgen(21)*2D*!dpi/20
	; data coord
	re_x = 2.0 * rd * cos(phi)
	re_y = 2.0 * rd * sin(phi)
   
	re_x= re_x + mw_center[0]
	re_y= re_y + mw_center[1]

	z= fltarr(101,101)
	z(*,*)= 0.0

	;shade_surf, wall_z, wall_x, wall_y, color= 200, /noerase, /t3d, $
	surface, mw_z, mw_x, mw_y, color= 200, /noerase, /t3d, $
		;AX = 70, AZ = 25, $
		xrange=[-xlen,xlen], yrange=[-xlen,xlen], zrange=[-xlen,xlen], $
        	;xstyle= 4, ystyle= 4, zstyle= 4, charsize= 1.5, $
        	xstyle= 1, ystyle= 1, zstyle= 1, charsize= 1.5, $
        	;xticklen=0.0, yticklen=0.0, zticklen=0.0, $
        	xtickformat='(a1)', ytickformat='(a1)', ztickformat='(a1)', $
        	position=[0.02, 0.02, 0.98, 0.98, 0.02, 0.98]

	;surface, z, re_x, re_y, color= 50, /noerase, $
        ;xrange=[-xlen,xlen], yrange=[-xlen,xlen], zrange=[-xlen,xlen], $
        ;xstyle= 1, ystyle= 1, zstyle= 1, charsize= 1.5, $
	;xticklen=0, yticklen=0, zticklen=0, $
	;xtickformat='(a1)', ytickformat='(a1)', ztickformat='(a1)', $
        ;position=[0.1, 0.1, 0.95, 0.95, 0.1, 0.95]


endif



;draw_mw_halo=1
draw_mw_halo=0
if draw_mw_halo eq 1 then begin

        mw_center= [0.307550, 0.179521, 0.00000]
        rvir= 0.22 
   
	    nring= 20
            phi = dindgen(nring+1)*2D*!dpi/nring
            ; data coord
            mw_x = rvir * cos(phi)
            mw_y = rvir * sin(phi)

            ;mw_x= mw_x + mw_center[0]
            ;mw_y= mw_y + mw_center[1]

            mw_z= fltarr(nring+1,nring+1)
            mw_z(*,*)= 0.0

        ;shade_surf, wall_z, wall_x, wall_y, color= 200, /noerase, /t3d, $
        surface, mw_z, mw_x, mw_y, color= 200, /noerase, /t3d, $
                xrange=[-xlen,xlen], yrange=[-xlen,xlen], zrange=[-xlen,xlen], $
                xstyle= 4, ystyle= 4, zstyle= 4, charsize= 1.5
                ;xstyle= 1, ystyle= 1, zstyle= 1, charsize= 1.5, $
                ;xticklen=0.0, yticklen=0.0, zticklen=0.0, $
                ;AX = 70, AZ = 25, $
                ;xtickformat='(a1)', ytickformat='(a1)', ztickformat='(a1)', $
                ;position=[0.02, 0.02, 0.98, 0.98, 0.02, 0.98]

endif



draw_mw_halo=1
;draw_mw_halo=0
if draw_mw_halo eq 1 then begin

	N= 10

	; equal to the Euclidean distance from the center:  
	mwZ = SHIFT(DIST(2*N), N, N)

	; Make Sphere with radius of 20
	idx=where(mwZ gt N)
	if idx(0) ne -1 then mwZ(idx)= N
	mwZ = SQRT(N^2 - mwZ^2)

	mwZ= 0.22 * mwZ / max(mwZ)

        mw_center= [0.307550, 0.179521, 0.00000]
        mwX= 2.0 * 0.22 * ((findgen(2*N)+0.5)/float(2*N) - 0.5) + mw_center[0]
        mwY= 2.0 * 0.22 * ((findgen(2*N)+0.5)/float(2*N) - 0.5) + mw_center[1]

	; ditch flat (i.e. zero) values
	help, idx
	idx=where(mwZ eq 0)
	help, idx
	;new_Z= fltarr(n_elements(idx))
	new_X= mwX
	new_Y= mwY
	new_Z= mwZ


        ;surface, mwZ, mwX, mwY, color= 150, /noerase, /t3d, $
        surface, new_Z, new_X, new_Y, color= 150, /noerase, /t3d, $
                xrange=[-xlen,xlen], yrange=[-xlen,xlen], zrange=[-xlen,xlen], $
                xstyle= 4, ystyle= 4, zstyle= 4, charsize= 1.5

	xyouts, 0.7, 0.57, '!6Milky Way', /normal, size=1.5, color= 0, charthick= 4.0

endif


draw_m31_halo=1
;draw_m31_halo=0
if draw_m31_halo eq 1 then begin

	N= 10

        ; equal to the Euclidean distance from the center:  
        m31Z = SHIFT(DIST(2*N), N, N)

        ; Make Sphere with radius of 20
        idx=where(m31Z gt N)
        if idx(0) ne -1 then m31Z(idx)= N
        m31Z = SQRT(N^2 - m31Z^2)

        m31Z= 0.26 * m31Z / max(m31Z)

        ;m31_center= [-0.495632, -0.289306, 0.00000]
        m31_center= [-0.445632, -0.289306, 0.00000]
        m31X= 2.0 * 0.26 * ((findgen(2*N)+0.5)/float(2*N) - 0.5) + m31_center[0] -0.6
        m31Y= 2.0 * 0.26 * ((findgen(2*N)+0.5)/float(2*N) - 0.5) + m31_center[1] -0.3

        surface, m31Z, m31X, m31Y, color= 50, /noerase, /t3d, $
                xrange=[-xlen,xlen], yrange=[-xlen,xlen], zrange=[-xlen,xlen], $
                xstyle= 4, ystyle= 4, zstyle= 4, charsize= 1.5

	xyouts, 0.15, 0.25, '!6Andromeda', /normal, size=1.5, color= 0, charthick= 4.0

	xyouts, 0.05, 0.93, '!6The Local Group', /normal, size=1.5, color= 0, charthick= 4.0

endif




;  Now complete the drawing of the axis box     
     
   ;Axis, !X.Crange[1], !Y.Crange[1], !Z.Crange[1], /YAXIS,   /T3D, COLOR=0, YTICKLEN=0, ytickformat='(a1)'
   ;Axis, !X.Crange[1], !Y.Crange[0], !Z.Crange[0], ZAXIS=-1, /T3D, COLOR=0, ZTICKLEN=0, ytickformat='(a1)'
   ;Axis, !X.Crange[0], !Y.Crange[1], !Z.Crange[0], /XAXIS,   /T3D, COLOR=0, XTICKLEN=0, ytickformat='(a1)'
   ;Axis, !X.Crange[0], !Y.Crange[1], !Z.Crange[1], /XAXIS,   /T3D, COLOR=0, XTICKLEN=0, ytickformat='(a1)'
   ;Axis, !X.Crange[1], !Y.Crange[1], !Z.Crange[0], ZAXIS=-1, /T3D, COLOR=0, ZTICKLEN=0, ytickformat='(a1)'
   ;Axis, !X.Crange[1], !Y.Crange[0], !Z.Crange[1], YAXIS=-1, /T3D, COLOR=0, YTICKLEN=0, ytickformat='(a1)'
   ;Axis, !X.Crange[1], !Y.Crange[0], !Z.Crange[0], /YAXIS,   /T3D, COLOR=0, YTICKLEN=0, ytickformat='(a1)'


; -------------------------------

device, /close


end






;======================================================================




pro test1, junk


if not keyword_set(junk) then begin
	print, "  "
	print, "  showit, junk"
	print, "  "
	print, "  "
	return
endif


filename='lgtest1.eps'


boxsize= 1.4  ; in Mpc
xlen= boxsize/2.0


;  mark today's position
;
t_sep= 0.714 ; in Mpc




; ----------------------------------
;   Try this plot thingy
; ----------------------------------

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=16.0, newysize=16.0


; equal to the Euclidean distance from the center:  
Z = SHIFT(DIST(40), 20, 20)  

; Make Gaussian with a 1/e width of 10:  
;Z = EXP(-(Z/10)^2)   


; Make Sphere with radius of 20
idx=where(z gt 20.0)
if idx(0) ne -1 then Z(idx)= 20.0
Z = SQRT((20.0)^2 - (Z)^2)   

; this ruins the 2D-ness of the array - shit
;idx= where(z gt 0)
;Z= Z(idx)

; Call SURFACE to display plot:  
;SURFACE, Z  


xlen= 50.0


;!p.ticklen=0.0
;
surface, z, /save, $
	xrange=[-xlen,xlen], yrange=[-xlen,xlen], zrange=[-xlen,xlen], $
	;xrange=[-xlen,xlen], yrange=[-xlen,xlen], zrange=[-1,1], $
	xstyle= 1, ystyle= 1, zstyle= 1, charsize= 1.5, color= 0, $
	xtickformat='(a1)', ytickformat='(a1)', ztickformat='(a1)', $
	position=[0.02, 0.02, 0.98, 0.98, 0.02, 0.98]
	;AX = 70, AZ = 25, $


; -------------------------------

device, /close


end






;======================================================================




