;====================================
;
;  Lay slit and measure velocity
;
;====================================
pro slit_velocity_info, frun, snapnum, filename=filename, $
				rotate_theta=rotate_theta, $
				rotate_phi=rotate_phi, $
				pa=pa


if not keyword_set(frun) then begin
        print, " "
        print, " slit_velocity_info, frun, snapnum, filename=filename, $"
	print, "		rotate_theta=rotate_theta, rotate_phi=rotate_phi,"
	print, "		pa=pa"
        print, " "
        print, " "
        return
endif

if not keyword_set(filename) then filename='slit.eps'
if not keyword_set(snapnum) then snapnum= 30


; --------------------------
;  Setup stuff
; --------------------------

; slit dimensions
; -----------------
slit_width = 0.5
slit_bins = 26
slit_len = 5.0

; plot range
slitmax = 1.2*slit_len

vmax= 175.0

; image xlen
xlen= 5.5



; details
; ----------
if not keyword_set(rotate_theta) then rotate_theta= 0.0
if not keyword_set(rotate_phi) then rotate_phi= 0.0

if not keyword_set(pa) then pa= 0.0
pa_orig= pa


; get print stuff ready
; ----------------------
initialize_plotinfo, 1
        
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=13, newysize=12



x0= 0.18
y0= 0.15
x1= 0.98
y1= 0.98




; ---------------
;  Do it Up
; ---------------

;ok=fload_initial_conditions(frun)
ok=fload_snapshot_bh(frun,snapnum)

center= fload_center_alreadycomp(1)
orig_center= center



; STARS - Total          (should always be done!)
;----------------
stellar_slit= 1
if stellar_slit then begin
	;   all stars
	; --------------
	x= fload_allstars_xyz('x', center=[0,0,0])
	y= fload_allstars_xyz('y', center=[0,0,0])
	z= fload_allstars_xyz('z', center=[0,0,0])
	vx= fload_allstars_v('x')
	vy= fload_allstars_v('y')
	vz= fload_allstars_v('z')
	mass= fload_allstars_mass(1)

	center= orig_center
	pa= pa_orig

	process_slit_and_overplot, x, y, z, vx, vy, vz, mass, $
	        x0, y0, x1, y1, $
	        center=center, /labelit, $
	        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	        slit_len=slit_len, slit_width=slit_width, $
	        slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
	        xlen=xlen, pa=pa, /ditchsigma, thiscolor= 150



endif





; STARS - New
;----------------
new_stars_slit= 1
if new_stars_slit eq 1 then begin

        ;   new stars - new, as in new
        ; ------------------------------
        x= fload_newstars_xyz('x', center=[0,0,0])
        y= fload_newstars_xyz('y', center=[0,0,0])
        z= fload_newstars_xyz('z', center=[0,0,0])
        vx= fload_newstars_v('x')
        vy= fload_newstars_v('y')
        vz= fload_newstars_v('z')
	age= fload_newstars_age(1)
        mass= fload_newstars_mass(1)

	idx=where(age gt 0.75)
	if idx(0) ne -1 then begin
		x= x(idx)
		y= y(idx)
		z= z(idx)
		vx= vx(idx)
		vy= vy(idx)
		vz= vz(idx)
		mass= mass(idx)
	endif

        center= orig_center
	pa= pa_orig

        process_slit_and_overplot, x, y, z, vx, vy, vz, mass, $
                x0, y0, x1, y1, /overplot, $
                center=center, /labelit, $
                rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                slit_len=slit_len, slit_width=slit_width, $
                slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
                xlen=xlen, pa=pa, /ditchsigma, thiscolor= 100





        ;   new stars - intermediate age
        ; ---------------------------------
        x= fload_newstars_xyz('x', center=[0,0,0])
        y= fload_newstars_xyz('y', center=[0,0,0])
        z= fload_newstars_xyz('z', center=[0,0,0])
        vx= fload_newstars_v('x')
        vy= fload_newstars_v('y')
        vz= fload_newstars_v('z')
        age= fload_newstars_age(1)
        mass= fload_newstars_mass(1)

        idx=where(age lt 0.75)
        if idx(0) ne -1 then begin
                x= x(idx)
                y= y(idx)
                z= z(idx)
                vx= vx(idx)
                vy= vy(idx)
                vz= vz(idx)
                mass= mass(idx)
        endif

        center= orig_center
        pa= pa_orig

        process_slit_and_overplot, x, y, z, vx, vy, vz, mass, $
                x0, y0, x1, y1, /overplot, $
                center=center, /labelit, $
                rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                slit_len=slit_len, slit_width=slit_width, $
                slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
                xlen=xlen, pa=pa, /ditchsigma, thiscolor= 75




endif




; STARS - Old
;----------------
old_stars_slit= 1
if old_stars_slit eq 1 then begin

        ;   old stars
        ; --------------
        x= fload_oldstars_xyz('x', center=[0,0,0])
        y= fload_oldstars_xyz('y', center=[0,0,0])
        z= fload_oldstars_xyz('z', center=[0,0,0])
        vx= fload_oldstars_v('x')
        vy= fload_oldstars_v('y')
        vz= fload_oldstars_v('z')
        mass= fload_oldstars_mass(1)

        center= orig_center
	pa= pa_orig

        process_slit_and_overplot, x, y, z, vx, vy, vz, mass, $
                x0, y0, x1, y1, /overplot, $
                center=center, /labelit, $
                rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                slit_len=slit_len, slit_width=slit_width, $
                slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
                xlen=xlen, pa=pa, /ditchsigma, thiscolor= 50



endif




; STARS - Add the collisionless stars?
;----------------
;collisionless_version= 0
collisionless_version= 1
if collisionless_version eq 1 then begin

	frun= '/raid4/tcox/collisionless/cvc3vc3f'            ; manually set this for now
	ok=fload_snapshot_bh(frun,snapnum)

        ;   collisionless stars
        ; -----------------------
        x= fload_allstars_xyz('x', center=[0,0,0])
        y= fload_allstars_xyz('y', center=[0,0,0])
        z= fload_allstars_xyz('z', center=[0,0,0])
        vx= fload_allstars_v('x')
        vy= fload_allstars_v('y')
        vz= fload_allstars_v('z')
        mass= fload_allstars_mass(1)

	center= fload_center_alreadycomp(1)
	pa= pa_orig

        process_slit_and_overplot, x, y, z, vx, vy, vz, mass, $
                x0, y0, x1, y1, /overplot, $
                center=center, /labelit, $
                rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                slit_len=slit_len, slit_width=slit_width, $
                slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
                xlen=xlen, pa=pa, /ditchsigma, thiscolor= 5



endif



; Gas
;----------
gaseous_slit= 0
;gaseous_slit= 1
if gaseous_slit eq 1 then begin

        ;   gas
        ; --------------
        x= fload_gas_xyz('x', center=[0,0,0])
        y= fload_gas_xyz('y', center=[0,0,0])
        z= fload_gas_xyz('z', center=[0,0,0])
        vx= fload_gas_v('x')
        vy= fload_gas_v('y')
        vz= fload_gas_v('z')
        mass= fload_gas_mass(1)

	center= orig_center
	pa= pa_orig

	process_slit_and_overplot, x, y, z, vx, vy, vz, mass, $
	        x0, y0, x1, y1, /overplot, $
	        center=center, /labelit, $
	        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	        slit_len=slit_len, slit_width=slit_width, $
	        slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
	        xlen=xlen, pa=pa, /ditchsigma

endif



; -------------
;  Done
; -------------

device, /close


; -------------

end






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Calculate  Velocity & Sigma along Slit
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; returns v_mean, v_disp

pro process_slit, x, y, v, mass, slit, v_mean, v_disp, $
			slit_width=slit_width, $
			slit_len=slit_len, $
			slit_bins=slit_bins, $
			v_mean_unweighted=v_mean_unweighted, $
			v_disp_unweighted=v_disp_unweighted, $
			mags=mags, $
			v_mean_magwt=v_mean_magwt, $
			v_disp_magwt=v_disp_magwt

slit  = fltarr(slit_bins)
slit(*) = (2*slit_len)*findgen(slit_bins)/float(slit_bins) - slit_len
dx = (2*slit_len)/float(slit_bins)

; trim particles
idx= where((y le slit_width/2.0) and (y ge -1.0*slit_width/2.0))
if idx(0) ne -1 then begin
	x_inslit= x(idx)
	mass_inslit= mass(idx)
	v_inslit= v(idx)
	if keyword_set(mags) then mags_inslit= mags(idx)
endif else begin
;	x_inslit= [14000]
endelse

v_mean  = fltarr(slit_bins)
v_mean_unweighted= fltarr(slit_bins)
v_mean_magwt= fltarr(slit_bins)
v_disp = fltarr(slit_bins)
v_disp_unweighted= fltarr(slit_bins)
v_disp_magwt= fltarr(slit_bins)


for i=0,slit_bins-1 do begin
	
	;idx = where((x_inslit gt slit(i)-dx) and (x_inslit le slit(i)+dx))
	idx = where((x_inslit gt slit(i)-0.5*dx) and (x_inslit le slit(i)+0.5*dx))
	if idx(0) ne -1 then begin

		; velocity
		v_mean(i)  = total(mass_inslit(idx)*v_inslit(idx))/total(mass_inslit(idx))
		v_mean_unweighted(i)= mean(v_inslit(idx))
		if keyword_set(mags) then v_mean_magwt(i)  = total(mags_inslit(idx)*v_inslit(idx))/total(mags_inslit(idx))

		; dispersion
		v_disp(i) = sqrt(total(mass_inslit(idx)*(v_inslit(idx)-v_mean(i))^2)/total(mass_inslit(idx)))
		v_disp_unweighted(i)= sqrt(variance(v_inslit(idx)))
		if keyword_set(mags) then v_disp_magwt(i) = sqrt(total(mags_inslit(idx)*(v_inslit(idx)-v_mean_magwt(i))^2)/total(mags_inslit(idx)))

;print, "Slit x= ",slit(i), v_mean(i), mean(v_inslit(idx)), v_disp(i), sqrt(variance(v_inslit(idx)))
	endif else begin
		v_mean(i) = 0.0
		v_disp(i) = 0.0
	endelse

endfor

; average of central three points
; -------------------------------
c_i= slit_bins/2
central_avg= 0.333*(v_disp[c_i-1]+v_disp[c_i]+v_disp[c_i+1])
print, "Central Velocity Disperion: ", central_avg


; rotation velocity
; ------------------
;f = (2.0*!pi/16.0)*findgen(17)
;usersym,cos(f),sin(f),/fill
;oplot,slit,v_rot,psym=-8,color=140


; dispersion
; ------------
;oplot,slit,v_disp,psym=-2,color=40


end









;===================================
;  Given the velocity information,
;  plot it.
;===================================
pro plot_slit_velocity, slitr, slitv, filename=filename, $
				vmax=vmax, slitmax=slitmax, $
				slitsigma=slitsigma


if not keyword_set(slitr) then begin
        print, " "
        print, " plot_slit_velocity, slitr, slitv, filename=filename,"
	print, "                     vmax=vmax, slitmax=slitmax"
        print, " "
        print, " "
        return
endif




; --------------------------
;  Setup stuff
; --------------------------

if not keyword_set(filename) then filename='velocity.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


xaxistitle='R (h!u-1!N kpc)'
yaxistitle='V(R) (km s!U-1!N)'
;yaxistitle='!Ms!N(R) (km s!U-1!N)'
xmax = slitmax
xmin = -slitmax
ymax = vmax
ymin = -vmax


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        ;/xlog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; rotation velocity
; ------------------
f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill
thispsym= 8
thiscolor= 140
oplot,slitr,slitv,psym=-thispsym,color=thiscolor

oplot, [slitr[0],slitr[4]], [ymin*0.62,ymin*0.62], psym=-thispsym,color=thiscolor
xyouts, slitr[5], ymin*0.65, 'velocity', color=thiscolor, charthick=3.0, size=1.33


if keyword_set(slitsigma) then begin
	thispsym= 2
	thiscolor= 50
	oplot, slitr, slitsigma, psym=-thispsym, color=thiscolor

	oplot, [slitr[0],slitr[4]], [ymin*0.82,ymin*0.82], psym=-thispsym,color=thiscolor
	xyouts, slitr[5], ymin*0.85, 'sigma', color=thiscolor, charthick=3.0, size=1.33
endif

; ------------------
;  Print any extras
; ------------------
;xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0
;xyouts, 0.2, 0.90, fload_fid(frun), /normal, charthick=2.0, size=1.33, color=0


; zero
x=[xmin,xmax]
y=[0,0]
oplot, x, y, psym=-3, linestyle=2, color= 0, thick=2.0



; -----------
;   Done
; -----------
device, /close


end








;===================================
;  Given the velocity dipsersion 
;  information, plot it.
;===================================
pro plot_slit_dispersion, slitr, slitsigma, filename=filename, $
				vmax=vmax, slitmax=slitmax


if not keyword_set(slitr) then begin
        print, " "
        print, " plot_slit_dispersion, slitr, slitsigma, filename=filename,"
	print, "                     vmax=vmax, slitmax=slitmax"
        print, " "
        print, " "
        return
endif




; --------------------------
;  Setup stuff
; --------------------------

if not keyword_set(filename) then filename='velocity.eps'


initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable=4


xaxistitle='R (h!u-1!N kpc)'
;yaxistitle='V(R) (km s!U-1!N)'
yaxistitle='!Ms!N(R) (km s!U-1!N)'
xmax = slitmax
xmin = -slitmax
ymax = vmax
ymin = 0.0


;---------------------------

!p.position= [0.18, 0.15, 0.95, 0.95]

plot, [1.0],[1.0], psym=-3, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        color= 0, $
        ;/ylog, $
        ;/xlog, $
        xstyle=1, $
        ystyle=1, $
        ;ystyle=8, $     ; this will suppress the second y axis
        xcharsize=1.50, ycharsize=1.50, $
        xthick=4.0, ythick=4.0, $
        charthick=3.0, $
        ;ytickformat='exp_label', $
        xtitle=xaxistitle, $
        ytitle=yaxistitle, $
        /nodata




; rotation velocity
; ------------------
f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill
thispsym= 8
thiscolor= 140
oplot,slitr,slitsigma,psym=-thispsym,color=thiscolor


; ------------------
;  Print any extras
; ------------------
;xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0
;xyouts, 0.2, 0.90, fload_fid(frun), /normal, charthick=2.0, size=1.33, color=0




; -----------
;   Done
; -----------
device, /close


end




; -----------------------------------------------------------------------------




;=============================================================
;
;  OK, now make an oversweeping procedure to generate a 
; multiplot with images, slits and velocity info.
;
; grid: 4x2
;
; top row, is four images of projected gas and stars, from
; two orthogonal directions.  slits are overlayed, and the 
; bottom row has the corresponding slit velocity info.
;
;=============================================================


pro multi_slit_plot, junk


if not keyword_set(junk) then begin
        print, " "
        print, " multi_slit_plot, junk"
        print, " "
        print, " "
        return
endif


; --------------------------
;  Setup stuff
; --------------------------

filename='slit_yo.eps'

center=[0,0,0]

; view 1
rotate_theta1= 0.0
rotate_phi1= 0.0
; view 2
rotate_theta2= 90.0
rotate_phi2= 0.0


; figure 1
;frun= "/raid4/tcox/vc3vc3h"
;snapnum= 1
;center= [26.7,19.8,0.1]

; figure 2
;frun= "/raid4/tcox/vc3vc3e"
;snapnum= 11
;center= [0.14,0.52,1.25]
;rotate_phi1= 45.0

; figure 3
;frun= "/raid4/tcox/vc3vc3e"
;snapnum= 12
;center= [-0.47,0.88,0.73]

; figure 4
;frun= "/raid4/tcox/vc3vc3e"
;snapnum= 25
;center= [-0.63,-1.91,0.66]

; figure x - whoops this was just testing
frun= "/raid4/tcox/vc3vc3b"
snapnum= 25
center= [0,0,0]


; slit dimensions
; -----------------
slit_width = 0.5
slit_bins = 20
slit_len = 5.0

; plot range
slitmax = 1.2*slit_len

vmax= 350.0

; image xlen
xlen= 1.2*slit_len



; start it up
; -------------

orig_center= center



        x0= 0.004
        xs= 0.248   ; assumes 4 panels
        x1= x0+xs
        x2= x0+xs+xs
        x3= x0+xs+xs+xs
        x4= x0+xs+xs+xs+xs
    
        y0= 0.003
        ys= 0.496   ; assumes 2 panels
        y1= y0+ys
        y2= y0+ys+ys
        y3= y0+ys+ys+ys



; ---------------
;  Do it Up
; ---------------

;ok=fload_initial_conditions(frun)
ok=fload_snapshot_bh(frun,snapnum)


; get print stuff ready
; ----------------------
initialize_plotinfo, 1
        
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=24, newysize=12
;setup_plot_stuff, 'ps', filename=filename, colortable= 1, newxsize=24, newysize=12
;setup_plot_stuff, 'ps', filename=filename, colortable= 0, newxsize=24, newysize=12





; Do Stars first
;=================

	x= fload_allstars_xyz('x', center=[0,0,0])
	y= fload_allstars_xyz('y', center=[0,0,0])
	z= fload_allstars_xyz('z', center=[0,0,0])
	vx= fload_allstars_v('x')
	vy= fload_allstars_v('y')
	vz= fload_allstars_v('z')
	mass= fload_allstars_mass(1)


	; first do orbital (xy, no rotation) plane
	; -----------------------------------------
	x_orig= x
	y_orig= y
	z_orig= z
	center= orig_center
	process_image_and_overplot, x_orig, y_orig, z_orig, mass, xlen, x0, y1, x1, y2, $
			rotate_phi=rotate_phi1, rotate_theta=rotate_theta1, $
			center=center, /showbhs, $
			set_maxden=10.0, set_dynrng=1.0e5, $
			msg1='stars', msg2='top', /drawbox, xslit=slit_len, yslit=0.5*slit_width


	x_orig= x
	y_orig= y
	z_orig= z
	vx_orig= vx
	vy_orig= vy
	vz_orig= vz
	center= orig_center
	process_slit_and_overplot, x_orig, y_orig, z_orig, vx_orig, vy_orig, vz_orig, mass, $
				x0, y0, x1, y1, $
				center=center, /labelit, $
				rotate_phi=rotate_phi1, rotate_theta=rotate_theta1, $
				slit_len=slit_len, slit_width=slit_width, $
				slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
				xlen=xlen
				


	; do 90-degree projection (xz)
	; -----------------------------------------
	x_orig= x
	y_orig= y
	z_orig= z
	center= orig_center
	process_image_and_overplot, x_orig, y_orig, z_orig, mass, xlen, x1, y1, x2, y2, $
			rotate_phi=rotate_phi2, rotate_theta=rotate_theta2, $
			center=center, /showbhs, $
			set_maxden=10.0, set_dynrng=1.0e5, $
			msg1='stars', msg2='side', /drawbox, xslit=slit_len, yslit=0.5*slit_width


	x_orig= x
	y_orig= y
	z_orig= z
	vx_orig= vx
	vy_orig= vy
	vz_orig= vz
	center= orig_center
	process_slit_and_overplot, x_orig, y_orig, z_orig, vx_orig, vy_orig, vz_orig, mass, $
				x1, y0, x2, y1, $
				center=center, $
				rotate_phi=rotate_phi2, rotate_theta=rotate_theta2, $
				slit_len=slit_len, slit_width=slit_width, $
				slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
				xlen=xlen
				



; Now do Gas
;============

        x= fload_gas_xyz('x', center=[0,0,0])
        y= fload_gas_xyz('y', center=[0,0,0])
        z= fload_gas_xyz('z', center=[0,0,0])
        vx= fload_gas_v('x')
        vy= fload_gas_v('y')
        vz= fload_gas_v('z')
        mass= fload_gas_mass(1)
        

	; first do orbital (xy, no rotation) plane
	; -----------------------------------------
	x_orig= x
	y_orig= y
	z_orig= z
	center= orig_center
	process_image_and_overplot, x_orig, y_orig, z_orig, mass, xlen, x2, y1, x3, y2, $
			rotate_phi=rotate_phi1, rotate_theta=rotate_theta1, $
			center=center, /showbhs, $
			set_maxden=10.0, set_dynrng=1.0e5, $
			msg1='gas', msg2='top', /drawbox, xslit=slit_len, yslit=0.5*slit_width


	x_orig= x
	y_orig= y
	z_orig= z
	vx_orig= vx
	vy_orig= vy
	vz_orig= vz
	center= orig_center
	process_slit_and_overplot, x_orig, y_orig, z_orig, vx_orig, vy_orig, vz_orig, mass, $
				x2, y0, x3, y1, $
				center=center, $
				rotate_phi=rotate_phi1, rotate_theta=rotate_theta1, $
				slit_len=slit_len, slit_width=slit_width, $
				slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
				xlen=xlen
				


	; do 90-degree projection (xz)
	; -----------------------------------------
	x_orig= x
	y_orig= y
	z_orig= z
	center= orig_center
	process_image_and_overplot, x_orig, y_orig, z_orig, mass, xlen, x3, y1, x4, y2, $
			rotate_phi=rotate_phi2, rotate_theta=rotate_theta2, $
			center=center, /showbhs, $
			set_maxden=10.0, set_dynrng=1.0e5, $
			msg1='gas', msg2='side', /drawbox, xslit=slit_len, yslit=0.5*slit_width


	x_orig= x
	y_orig= y
	z_orig= z
	vx_orig= vx
	vy_orig= vy
	vz_orig= vz
	center= orig_center
	process_slit_and_overplot, x_orig, y_orig, z_orig, vx_orig, vy_orig, vz_orig, mass, $
				x3, y0, x4, y1, $
				center=center, $
				rotate_phi=rotate_phi2, rotate_theta=rotate_theta2, $
				slit_len=slit_len, slit_width=slit_width, $
				slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
				xlen=xlen
				


; any extras?
; -------------

xyouts, 0.02, 0.45, fload_timelbl(1,2,/noteq), /normal, size= 1.33, charthick=3.0, color= 0


device, /close

; -------------
;  Done
; -------------


end







; process the slit, and overplot it
; -----------------------------------
pro process_slit_and_overplot, x, y, z, vx, vy, vz, mass, $
				x0, y0, x1, y1, $
				center=center, labelit=labelit, $
				rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
				slit_len=slit_len, slit_width=slit_width, $
				slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
				xlen=xlen, pa=pa, mags=mags, overplot=overplot, $
				ditchsigma=ditchsigma, thiscolor= thiscolor



if keyword_set(center) then begin
     print, "process slit, changing center to= ",center
     x=x-center[0]
     y=y-center[1]
     z=z-center[2]

endif


; Now we need to rotate for the slit
; ------------------------------------
;if (rotate_theta gt 0) or (rotate_phi gt 0) or (total(abs(center)) gt 0) then begin
if keyword_set(rotate_theta) or keyword_set(rotate_phi) then begin

	print, "process slit, rotating to theta/phi= ", rotate_theta, rotate_phi

        ; transform to radians
        theta= !PI*rotate_theta/180.0
        phi= !PI*rotate_phi/180.0

        ; 1. first around y axis (theta)
        ; then around z axis (phi)
        ;x_new= x*(cos(theta)*cos(phi)) - y*sin(phi) + z*(cos(phi)*sin(theta))
        ;y_new= x*(cos(theta)*sin(phi)) + y*cos(phi) + z*(sin(phi)*sin(theta))
        ;z_new= -x*sin(theta) + z*cos(theta)

        ; 2. first around z axis (phi)
        ; then around y axis (theta)
        ;x_new= x*(cos(theta)*cos(phi)) + y*(cos(theta)*sin(phi)) + z*sin(theta)
        ;y_new= -x*sin(theta) + y*cos(theta)
        ;z_new= x*(sin(theta)*cos(phi)) - y*(sin(theta)*sin(phi)) + z*cos(theta)

	; 2.5 first around z axis (phi) - clockwise
	; then around y axis (theta)  - counter clockwise
	x_new= x*(cos(theta)*cos(phi)) + y*(cos(theta)*sin(phi)) - z*sin(theta)
	y_new= -x*sin(phi) + y*cos(phi)
	z_new= x*(sin(theta)*cos(phi)) + y*(sin(theta)*sin(phi)) + z*cos(theta)

        ; 3. first around x axis (theta)
        ; then around z axis (phi)
        ;x_new= x*cos(phi) + y*(sin(phi)*cos(theta)) + z*(sin(phi)*sin(theta))
        ;y_new= -x*sin(phi) + y*(cos(phi)*cos(theta)) + z*(cos(phi)*sin(theta))
        ;z_new= -y*sin(theta) + z*cos(theta)

        ; 4. first around y axis (theta)
        ; then around x axis (phi)
        ;x_new= x*cos(theta) + z*sin(theta)
        ;y_new= -x*(sin(phi)*sin(theta)) + y*cos(phi) + z*(cos(theta)*sin(phi))
        ;z_new= -x*(cos(phi)*sin(theta)) - y*sin(phi) + z*(cos(theta)*cos(phi))

        ;x= x_new
        ;y= y_new
        ;z= z_new


        ; rotate the velocities
        ; ----------------------
	; 1.
        ;vx_new= vx*(cos(theta)*cos(phi)) - vy*sin(phi) + vz*(cos(phi)*sin(theta))
        ;vy_new= vx*(cos(theta)*sin(phi)) + vy*cos(phi) + vz*(sin(phi)*sin(theta))
        ;vz_new= -vx*sin(theta) + vz*cos(theta)

	; 2.5 first around z axis (phi) - clockwise
	; then around y axis (theta)  - counter clockwise
	vx_new= vx*(cos(theta)*cos(phi)) + vy*(cos(theta)*sin(phi)) - vz*sin(theta)
	vy_new= -vx*sin(phi) + vy*cos(phi)
	vz_new= vx*(sin(theta)*cos(phi)) + vy*(sin(theta)*sin(phi)) + vz*cos(theta)

	; 3.
        ;vx_new= vx*cos(phi) + vy*(sin(phi)*cos(theta)) + vz*(sin(phi)*sin(theta))
        ;vy_new= -vx*sin(phi) + vy*(cos(phi)*cos(theta)) + vz*(cos(phi)*sin(theta))
        ;vz_new= -vy*sin(theta) + vz*cos(theta)

	; 4. 
        ;vx_new= vx*cos(theta) + vz*sin(theta)
        ;vy_new= -vx*(sin(phi)*sin(theta)) + vy*cos(phi) + vz*(cos(theta)*sin(phi))
        ;vz_new= -vx*(cos(phi)*sin(theta)) - vy*sin(phi) + vz*(cos(theta)*cos(phi))


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


; rotate just x and y so that it's
; along the major axis position angle (pa)
; ----------------------------------------
if keyword_set(pa) then begin
	print, "actual slit pa= ", pa
	par = pa * !PI/180.0
	x_newest= x_new*cos(par) + y_new*sin(par)      ; with this choice of signs, it automatically
	y_newest= -x_new*sin(par) + y_new*cos(par)     ;  rotates clock-wise (hence don't need -par)
	x_new= x_newest
	y_new= y_newest
endif else begin
	print, "pa not set.  fixing to be 0"
	pa= 0.0
endelse



; calculate slit velocities
process_slit, x_new, y_new, vz_new, mass, slit, v_mean, v_disp, $
			slit_width=slit_width, $
			slit_len=slit_len, $
			slit_bins=slit_bins, $
			v_mean_unweighted=v_mean_unweighted, $
			v_disp_unweighted=v_disp_unweighted, $
			mags=mags, $
			v_mean_magwt=v_mean_magwt, $
			v_disp_magwt=v_disp_magwt


; now v_mean and v_disp contain velocity info


slitr= slit
slitv= v_mean
slitsigma= v_disp


; plot it up
; -----------
xaxistitle='!6R (h!u-1!N kpc)'
;yaxistitle='!6V(R) (km s!U-1!N)'
;yaxistitle='!6Velocity (km s!U-1!N)'
yaxistitle='!6(km s!U-1!N)'
;yaxistitle='!Ms!N(R) (km s!U-1!N)'
xmax = slitmax
xmin = -slitmax
ymax = vmax
ymin = -vmax


;---------------------------

if not keyword_set(overplot) then begin
	;!p.position= [0.18, 0.15, 0.95, 0.95]
	!p.position= [x0, y0, x1, y1]
	!p.ticklen=0.03

	plot, [1.0],[1.0], psym=-3, $
	        xrange=[xmin,xmax], yrange=[ymin,ymax], $
	        color= 0, $
	        ;/ylog, $
	        ;/xlog, $
	        xstyle=1, $
	        ystyle=1, $
	        ;ystyle=8, $     ; this will suppress the second y axis
	        xcharsize=1.10, $
		ycharsize=1.10, $
	        ;xcharsize=1.5, $
		;ycharsize=1.5, $
	        xthick=4.0, ythick=4.0, $
	        charthick=3.0, $
	        ;ytickformat='exp_label', $
		;xtickformat='(a1)', $
		;ytickformat='(a1)', $
	        xtitle=xaxistitle, $
	        ytitle=yaxistitle, $
	        /nodata, /noerase

endif




; rotation velocity
; ------------------
f = (2.0*!pi/16.0)*findgen(17)
usersym,cos(f),sin(f),/fill
if not keyword_set(thiscolor) then thiscolor= 150
if thiscolor eq 150 then thispsym= 8
if thiscolor eq 100 then thispsym= 2
if thiscolor eq 75 then thispsym= 7
if thiscolor eq 50 then thispsym= 5
if thiscolor eq 5 then thispsym= 3

;oplot,slitr,slitv,psym=-thispsym,color=thiscolor, thick=5.0, symsize=1.5
oplot,slitr,slitv,psym=-thispsym,color=thiscolor, thick=2.0, symsize=1.00
;oplot,slitr,v_mean_unweighted,psym=-3,linestyle= 1, color=thiscolor

if keyword_set(mags) then begin
	oplot, slitr, v_mean_magwt, psym=-3, linestyle= 2, color=thiscolor, thick=3.0
endif


print, "slit velocity max/min= ",max(slitv),min(slitv)

if keyword_set(labelit) then begin
   yval= 0.65
   ylbl= '!6velocity'
   ;yval= 0.50
   ;ylbl= '!6all stars'
   if thiscolor eq 100 then yval= 0.59
   if thiscolor eq 100 then ylbl= 'young'
   if thiscolor eq 75 then yval= 0.68
   if thiscolor eq 75 then ylbl= 'intermediate'
   if thiscolor eq 50 then yval= 0.77
   if thiscolor eq 50 then ylbl= 'old'
   if thiscolor eq 5 then yval= 0.86
   if thiscolor eq 5 then ylbl= 'dissipationless'
   ;oplot, [slitr[0],slitr[2]], [ymin*yval,ymin*yval], psym=-thispsym,color=thiscolor, symsize= 1.5, thick= 5.0
   oplot, [slitr[0],slitr[2]], [ymin*yval,ymin*yval], psym=-thispsym,color=thiscolor, symsize= 1.00, thick= 2.0
   xyouts, slitr[4], ymin*(yval+0.02), ylbl, color=thiscolor, charthick=3.0, size=1.33
endif


if keyword_set(slitsigma) and not keyword_set(ditchsigma) then begin
	thispsym= 2
	thiscolor= 50
	;oplot, slitr, slitsigma, psym=-thispsym, color=thiscolor, symsize= 1.5, thick= 5.0
	oplot, slitr, slitsigma, psym=-thispsym, color=thiscolor, symsize= 1.33, thick= 2.0
	;oplot, slitr, v_disp_unweighted, psym=-3, linestyle=1, color=thiscolor

	if keyword_set(mags) then begin
	   oplot, slitr, v_disp_magwt, psym=-3, linestyle=2, color=thiscolor, thick=3.0
	endif

	print, "slit dispersion max/min= ", max(slitsigma), min(slitsigma)

	if keyword_set(labelit) then begin
	   oplot, [slitr[0],slitr[2]], [ymin*0.82,ymin*0.82], psym=-thispsym,color=thiscolor
	   xyouts, slitr[4], ymin*0.85, '!7r!6', color=thiscolor, charthick=3.0, size=1.33
	endif
endif


; ------------------
;  Print any extras
; ------------------
;xyouts, 0.2, 0.90, fload_getid(frun), /normal, charthick=1, size=1.33, color=0
;xyouts, 0.2, 0.90, fload_fid(frun), /normal, charthick=2.0, size=1.33, color=0


if not keyword_set(overplot) then begin
	; zero
	x=[xmin,xmax]
	y=[0,0]
	oplot, x, y, psym=-3, linestyle=2, color= 0, thick=2.0
endif



end








; overplot the image
; --------------------
pro process_image_and_overplot, x, y, z, mass, xlen, $
			x0, y0, x1, y1, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			center=center, showbhs=showbhs, $
			show_re=show_re, $
			msg1=msg1, msg2=msg2, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			ReturnedMap=ReturnedMap, $
			msg=msg, drawbox=drawbox, xslit=xslit, yslit=yslit, $
			pa=pa, HalfMassSB=HalfMassSB

	m= mass
	if keyword_set(center) then orig_center= center else orig_center=[0,0,0]

	print, "=================================="

	xin= x
	yin= y
	zin= z
        contour_makepic, x, y, z, m, xlen, xz=xz, yz=yz, $
                        filename=filename, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, $
                        pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
                        NxNImage=NxNImage, HalfMassSB=HalfMassSB
        
        Pic= NxNImage
	ReturnedMap= NxNImage

    ;print, "center= ", center
    ;print, "center (original)= ", orig_center
    center= [0,0,0]

    xmin= -xlen+center[0]
    xmax=  xlen+center[0]
    ymin= -xlen+center[1]
    ymax=  xlen+center[1]


    tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


    !p.position=[x0,y0,x1,y1]
    !p.ticklen=0.03

                ; creates axes and plot style
                ; ---------------------------
                plot,[0],[0], psym=3,  $
                      xrange=[xmin,xmax], $
                      yrange=[ymin,ymax], $
                      /NOERASE, $
                      color= 0, $
                      xcharsize=0.01, $
		      ycharsize=0.01, $
                      xthick=4.0, $
		      ythick=4.0, $
                      charthick=3.0, $
		      ;xticks=10, $
		      ;yticks=10, $
		      xtickformat='(a1)', ytickformat='(a1)', $
                      xstyle=1, ystyle=1, $
                      xtitle=xtit, $
                      ytitle=ytit, $ ;/normal, $
                      /nodata


     ;  key
     ; ------

     ; white box
     bx0=xmin+1.0
     bx1=xmin+6.0
     by0=ymax-1.0
     by1=ymax-3.6
     polyfill, [bx0,bx0,bx1,bx1], [by0,by1,by1,by0], color= 1
     oplot, [bx0,bx0,bx1,bx1,bx0], [by0,by1,by1,by0,by0] , color= 0, thick=3.0

     ;xyouts, x0+0.04, y1-0.08, fload_timelbl(1,2,/noteq), /normal, size= 1.7, charthick=3.0, color= 0
     xyouts, x0+0.02, y1-0.05, msg1, /normal, size= 1.3, charthick=4.0, color= 0
     xyouts, x0+0.02, y1-0.10, msg2, /normal, size= 1.3, charthick=4.0, color= 0




                ; -----------------------------
                ;  X marks the black hole spot
                ; -----------------------------
                if keyword_set(showbhs) and fload_npart(5) ge 1 then begin

                        loadct, 4          ; blue/green/red/yellow (std)
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0

                        bhpsym= 7
                        bhsymsize= 2.0
                        bhcolor= 220

                        n_bh= fload_npart(5)
                        bhid= fload_blackhole_id(1)
                        bhid1= bhid[0]
                        bh_center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
                        
                        bh_new= bh_center_1 - orig_center
                        if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin                     
                           bh_c= bh_new
                           theta= !PI*rotate_theta/180.0
                           phi= !PI*rotate_phi/180.0
                           ; 4. rotate y (theta) and then x (phi)
                           ;bh_new[0]= bh_c[0]*cos(theta) + bh_c[2]*sin(theta)
                           ;bh_new[1]= -bh_c[0]*(sin(phi)*sin(theta)) + bh_c[1]*cos(phi) + bh_c[2]*(cos(theta)*sin(phi))
                           ;bh_new[2]= -bh_c[0]*(cos(phi)*sin(theta)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(theta)*cos(phi))
                           ; 3. rotate x (theta) and then z (phi)
                           ;bh_new[0]= bh_c[0]*cos(phi) + bh_c[1]*(sin(phi)*cos(theta)) + bh_c[2]*(sin(phi)*sin(theta))
                           ;bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*(cos(phi)*cos(theta)) + bh_c[2]*(cos(phi)*sin(theta))
                           ;bh_new[2]= -bh_c[1]*sin(theta) + bh_c[2]*cos(theta)
			   ; 1. rotate y (theta) and then z (phi)
			   ;bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(phi)*sin(theta))
			   ;bh_new[1]= bh_c[0]*(cos(theta)*sin(phi)) + bh_c[1]*cos(phi) + bh_c[2]*(sin(phi)*sin(theta))
			   ;bh_new[2]= -bh_c[0]*sin(theta) + bh_c[2]*cos(theta)
	        	   ; 2.5 first around z axis (phi) - clockwise
	        	   ; then around y axis (theta)  - counter clockwise
	        	   bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) + bh_c[1]*(cos(theta)*sin(phi)) - bh_c[2]*sin(theta)
	        	   bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*cos(phi)
	        	   bh_new[2]= bh_c[0]*(sin(theta)*cos(phi)) + bh_c[1]*(sin(theta)*sin(phi)) + bh_c[2]*cos(theta)
                        endif
                        
                        ; for the moment, just plot xy black hole position
                        oplot, [bh_new[0]], [bh_new[1]], psym=bhpsym, thick=6.0, $
                                        color=bhcolor, symsize=bhsymsize
                        
                        if n_bh gt 1 then begin
                                bhid2= bhid[1]
                                bh_center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
                                bh_new= bh_center_2 - orig_center
                                if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin
                                   bh_c= bh_new
                                   theta= !PI*rotate_theta/180.0
                                   phi= !PI*rotate_phi/180.0
                                   ; 4. rotate y (theta) and then x (phi)
                                   ;bh_new[0]= bh_c[0]*cos(theta) + bh_c[2]*sin(theta)                                 
                                   ;bh_new[1]= -bh_c[0]*(sin(phi)*sin(theta)) + bh_c[1]*cos(phi) + bh_c[2]*(cos(theta)*sin(phi))
                                   ;bh_new[2]= -bh_c[0]*(cos(phi)*sin(theta)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(theta)*cos(phi))
                                   ; 3. rotate x (theta) and then z (phi)
                                   ;bh_new[0]= bh_c[0]*cos(phi) + bh_c[1]*(sin(phi)*cos(theta)) + bh_c[2]*(sin(phi)*sin(theta))
                                   ;bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*(cos(phi)*cos(theta)) + bh_c[2]*(cos(phi)*sin(theta))
                                   ;bh_new[2]= -bh_c[1]*sin(theta) + bh_c[2]*cos(theta)                              
				   ; 1. rotate y (theta) and then z (phi)
				   ;bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(phi)*sin(theta))
				   ;bh_new[1]= bh_c[0]*(cos(theta)*sin(phi)) + bh_c[1]*cos(phi) + bh_c[2]*(sin(phi)*sin(theta))
				   ;bh_new[2]= -bh_c[0]*sin(theta) + bh_c[2]*cos(theta)
				   ; 2.5 first around z axis (phi) - clockwise
				   ; then around y axis (theta)  - counter clockwise
				   bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) + bh_c[1]*(cos(theta)*sin(phi)) - bh_c[2]*sin(theta)
				   bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*cos(phi)
				   bh_new[2]= bh_c[0]*(sin(theta)*cos(phi)) + bh_c[1]*(sin(theta)*sin(phi)) + bh_c[2]*cos(theta)
                                endif
                                
                                print, "original center= ",orig_center
                                print, "Blackhole #1: ",bhid1, bh_center_1
                                print, "Blackhole #2: ",bhid2, bh_center_2
                                
                                oplot, [bh_new[0]], [bh_new[1]], psym=bhpsym, thick=6.0, $                                  
                                        color=bhcolor, symsize=bhsymsize
                                
                        endif else begin
                                print, "original center= ",orig_center
                                print, "Blackhole: ", bhid1, bh_center_1
                        endelse
                
                endif

                loadct, 4          ; blue/green/red/yellow (std)
                tvlct,r,g,b,/get
                v1=[0,255]
                v2=[0,255]
                v3=[0,255]
                tvlct,v1,v2,v3,0


		if keyword_set(pa) then pa_orig= pa else pa_orig= -1

		; overplot contours, and ellipse
		overplot_contours, Pic, x0, y0, x1, y1, hmsb=HalfMassSB, /fitellipse, xlen=xlen, pa=pa

		if pa_orig ne -1 then pa= -1.0*pa_orig   ; minus sign so manual setting
							 ; goes clockwise


                ; -------------------------------
                ;  Draw Box (slit, for instance)
                ; -------------------------------
                ;drawbox= 1
                if keyword_set(drawbox) and not keyword_set(pa) then begin
                   slitcolor= 1
                   slitthick= 4.0
                   oplot, [-xslit, -xslit], [-yslit,yslit], psym=-3, color= slitcolor, thick= slitthick
                   oplot, [-xslit, xslit], [yslit,yslit], psym=-3, color= slitcolor, thick= slitthick  
                   oplot, [xslit, xslit], [-yslit,yslit], psym=-3, color= slitcolor, thick= slitthick
                   oplot, [-xslit, xslit], [-yslit,-yslit], psym=-3, color= slitcolor, thick= slitthick
                endif

		print, "xyslit= ",xslit, yslit


		; rotate just x and y so that it's
		; along the major axis position angle (pa)
		; ----------------------------------------
		if keyword_set(pa) then begin
		   print, "slit pa= ",pa
		   ; actually the pa passed is clockwise rotation, so we 
		   ; need to flips this around
		   if pa gt 0.0 then pa= 180.0 - pa else pa= -pa
		   par = pa * !PI/180.0
		   print, "new pa= ", pa

		   xtr= xslit*cos(par) - yslit*sin(par)       ; this choice of signs rotates
		   ytr= xslit*sin(par) + yslit*cos(par)       ; counter-clockwise

		   xbr= xslit*cos(par) - (-yslit)*sin(par)
		   ybr= xslit*sin(par) +(-yslit)*cos(par)

		   xtl= (-xslit)*cos(par) - yslit*sin(par)
		   ytl= (-xslit)*sin(par) + yslit*cos(par)

		   xbl= (-xslit)*cos(par) - (-yslit)*sin(par)
		   ybl= (-xslit)*sin(par) + (-yslit)*cos(par)

		   ; need axes to be reset

    		   !p.position=[x0,y0,x1,y1]
    		   !p.ticklen=0.03

                   ; creates axes and plot style
                   ; ---------------------------
                   plot,[0],[0], psym=3,  $
                      xrange=[xmin,xmax], $
                      yrange=[ymin,ymax], $
                      /NOERASE, $
                      color= 0, $
                      xcharsize=0.01, $
		      ycharsize=0.01, $
                      xthick=4.0, $
		      ythick=4.0, $
                      charthick=3.0, $
		      ;xticks=10, $
		      ;yticks=10, $
		      xtickformat='(a1)', ytickformat='(a1)', $
                      xstyle=1, ystyle=1, $
                      xtitle=xtit, $
                      ytitle=ytit, $ ;/normal, $
                      /nodata


                   ;slitcolor= 1    ; white
                   slitcolor= 240     ; yellow
                   slitthick= 4.0
                   oplot, [xtr,xbr,xbl,xtl,xtr], [ytr,ybr,ybl,ytl,ytr], psym=-3, color= slitcolor, thick= slitthick


		endif


	if keyword_set(show_re) then begin
	       ; 2.5 first around z axis (phi) - clockwise
	       ; then around y axis (theta)  - counter clockwise
		if keyword_set(rotate_theta) then theta= !PI*rotate_theta/180.0 else theta= 0.0
		if keyword_set(rotate_phi) then phi= !PI*rotate_phi/180.0 else phi= 0.0
		xin= xin-orig_center[0]
		yin= yin-orig_center[1]
		zin= zin-orig_center[2]
		x_new= xin*(cos(theta)*cos(phi)) + yin*(cos(theta)*sin(phi)) - zin*sin(theta)
		y_new= -xin*sin(phi) + yin*cos(phi)
		;R_e= 2.4384
		R_e= determine_reff(x_new,y_new,m)

		phi = dindgen(101)*2D*!dpi/100
		; data coord
		re_x = R_e * cos(phi)
		re_y = R_e * sin(phi)

		re_x= re_x + center[0]
		re_y= re_y + center[1]

		oplot, re_x, re_y, color= 50, thick=4.0, linestyle=0, psym=-3
		;xyouts, 0.22, 0.68, 'R!De!N=2.4', /normal, size=1.3, color= 50, charthick= 4.0
	endif


end





; -----------------------------------------------------------------------------




;=============================================================
;
; grid: 2x2
;
; top row: 1) project gas or star density
;          2) velocity field
; bottom row: 1) slit info
;             2) dispersion field
;
;=============================================================


;pro multi_slit_plot_2, junk
pro multi_slit_plot_2, frun, rt, rp, pa=pa, snapnum=snapnum, $
					filename=filename


;if not keyword_set(junk) then begin
if not keyword_set(frun) then begin
        print, " "
        ;print, " multi_slit_plot_2, junk"
        print, " multi_slit_plot_2, frun, rotate_theta, rotate_phi,"
        print, "                          pa=pa, snapnum=snapnum,"
        print, "                          filename=filename"
        print, " "
        return
endif


; --------------------------
;  Setup stuff
; --------------------------

if not keyword_set(filename) then filename='slit_ho.eps'


; view 0
if keyword_set(rt) then rotate_theta= rt else rotate_theta= 0.0
if keyword_set(rp) then rotate_phi= rp else rotate_phi= 0.0


if keyword_set(pa) then pa_fixed= pa else pa_fixed= -1


;frun= "/raid4/tcox/vc3vc3b"
;frun= "/raid4/tcox/vc3vc3d"
;frun= "/raid4/tcox/vc3vc3e"
;frun= "/raid4/tcox/vc3vc3f"
;frun= "/raid4/tcox/cvc3vc3f"
;frun= "/raid4/tcox/cvc3vc3h"
;frun= "/raid4/tcox/vc3vc3k"
;frun= "/raid4/tcox/vc3vc3l"
;frun= "/raid4/tcox/vc3vc3o"
if not keyword_set(frun) then frun= "/raid4/tcox/vc3vc3b"

if not keyword_set(snapnum) then snapnum= 30



; slit dimensions
; -----------------
slit_width = 0.5
slit_bins = 26
slit_len = 5.0

; plot range
slitmax = 1.2*slit_len

vmax= 250.0
vscale= 200.0

; image size
;xlen= 1.2*slit_len
xlen= 8.0



; start it up
; -------------



; these are for 4x4 image that takes up ALL the room
;x0= 0.003
;xs= 0.496   ; assumes 2 panels
;y0= 0.003
;ys= 0.496   ; assumes 2 panels

; these are for 4x4 image that leaves room for labels
x0= 0.1111111
xs= 0.3888889   ; assumes 2 panels
y0= 0.0666667
ys= 0.4666667   ; assumes 2 panels

x1= x0+xs
x2= x0+xs+xs
x3= x0+xs+xs+xs
    
y1= y0+ys
y2= y0+ys+ys
y3= y0+ys+ys+ys



; ---------------
;  Do it Up
; ---------------

;ok=fload_initial_conditions(frun)
ok=fload_snapshot_bh(frun,snapnum)


; get print stuff ready
; ----------------------
initialize_plotinfo, 1
        
;setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=12, newysize=12
;setup_plot_stuff, 'ps', filename=filename, colortable= 1, newxsize=12, newysize=12
setup_plot_stuff, 'ps', filename=filename, colortable= 0, newxsize=18, newysize=15


center= fload_center_alreadycomp(1)
orig_center= center



; Do Stars
;===========
;do_stars_slit= 0
do_stars_slit= 1
if do_stars_slit eq 1 then begin
	x= fload_allstars_xyz('x', center=[0,0,0])
	y= fload_allstars_xyz('y', center=[0,0,0])
	z= fload_allstars_xyz('z', center=[0,0,0])
	vx= fload_allstars_v('x')
	vy= fload_allstars_v('y')
	vz= fload_allstars_v('z')
	mass= fload_allstars_mass(1)
endif


; Do Gas
;===========
do_gas_slit= 0
;do_gas_slit= 1
if do_gas_slit eq 1 then begin
        x= fload_gas_xyz('x', center=[0,0,0])
        y= fload_gas_xyz('y', center=[0,0,0])
        z= fload_gas_xyz('z', center=[0,0,0])
        vx= fload_gas_v('x')
        vy= fload_gas_v('y')
        vz= fload_gas_v('z')
        mass= fload_gas_mass(1)
endif

        

;  top-left, projected surface density + slit
; ---------------------------------------------
x_orig= x
y_orig= y
z_orig= z
center= orig_center
process_image_and_overplot, x_orig, y_orig, z_orig, mass, xlen, x0, y1, x1, y2, $
	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	center=center, /show_re, $ ;/showbhs, $
	set_maxden=10.0, set_dynrng=1.0e5, $
	;msg1='stars', msg2=' ', /drawbox, xslit=slit_len, yslit=0.5*slit_width, $
	msg1=' ', msg2=' ', /drawbox, xslit=slit_len, yslit=0.5*slit_width, $
	ReturnedMap=ReturnedMap, pa=pa, HalfMassSB=HalfMassSB



if pa_fixed ne -1 then pa= pa_fixed


;  bottom-left, velocity field within slit
; ------------------------------------------
x_orig= x
y_orig= y
z_orig= z
vx_orig= vx
vy_orig= vy
vz_orig= vz
center= orig_center
process_slit_and_overplot, x_orig, y_orig, z_orig, vx_orig, vy_orig, vz_orig, mass, $
        x0, y0, x1, y1, $
        center=center, /labelit, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        slit_len=slit_len, slit_width=slit_width, $
        slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
        xlen=xlen, pa=pa




;  top-right, velocity field
; ----------------------------
x_orig= x
y_orig= y
z_orig= z
vx_orig= vx
vy_orig= vy
vz_orig= vz
center= orig_center
process_velimage_and_overplot, x_orig, y_orig, z_orig, vx_orig, vy_orig, vz_orig, mass, $
	xlen, x1, y1, x2, y2, /velocitymap, $
	rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
	center=center, $ ;/showbhs, $
	set_maxden=10.0, set_dynrng=1.0e5, $
	msg1=' ', msg2=' ', $ ; /drawbox, xslit=slit_len, yslit=0.5*slit_width, $
	ContourMap=ReturnedMap, HalfMassSB=HalfMassSB






;  bottom-right, dispersion field
; ---------------------------------
x_orig= x
y_orig= y
z_orig= z
vx_orig= vx
vy_orig= vy
vz_orig= vz
center= orig_center
process_velimage_and_overplot, x_orig, y_orig, z_orig, vx_orig, vy_orig, vz_orig, mass, $
        xlen, x1, y0, x2, y1, /dispersionmap, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        center=center, $ ;/showbhs, $
        set_maxden=10.0, set_dynrng=1.0e5, $
        msg1=' ', msg2=' ', $ ; /drawbox, xslit=slit_len, yslit=0.5*slit_width, $
        ContourMap=ReturnedMap, HalfMassSB=HalfMassSB





; any extras?
; -------------

;xyouts, 0.02, 0.45, fload_timelbl(1,2,/noteq), /normal, size= 1.33, charthick=3.0, color= 0


device, /close

; -------------
;  Done
; -------------


end











; overplot the velocity image
; -----------------------------
pro process_velimage_and_overplot, x, y, z, vx, vy, vz, mass, xlen, $
			x0, y0, x1, y1, $
			rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
			center=center, showbhs=showbhs, $
			msg1=msg1, msg2=msg2, $
			set_maxden=set_maxden, set_dynrng=set_dynrng, $
			msg=msg, drawbox=drawbox, xslit=xslit, yslit=yslit, $
			velocitymap=velocitymap, dispersionmap=dispersionmap, $
			ContourMap=ContourMap, HalfMassSB=HalfMassSB

	m= mass
	if keyword_set(center) then orig_center= center else orig_center=[0,0,0]

	print, "=================================="


	if keyword_set(velocitymap) then begin
		set_maxden= 150.0
		set_dynrng= 2.0*set_maxden
                        loadct, 4          ; blue/green/red/yellow (std)
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0

	endif
	if keyword_set(dispersionmap) then begin
		set_maxden= 200.0
		set_dynrng= 100.0
                        loadct, 5
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0

	endif

       	contour_makevelpic, x, y, z, vx, vy, vz, m, xlen, xz=xz, yz=yz, $
                        filename=filename, fitstoo=fitstoo, $
                        xthickness=xthickness, ythickness=ythickness, $
                        pixels=pixels, zthickness=zthickness, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        crude=crude, center=center, $
                        set_maxden=set_maxden, set_dynrng=set_dynrng, $
			velocitymap=velocitymap, dispersionmap=dispersionmap, $
                        NxNImage=NxNImage
        
       	Pic= NxNImage



    if n_elements(Pic) le 0 then begin
	print, "  "
	print, "  WARNING: no map set"
	print, "  "
	return
    endif

    ;loadct, 33          ; temp
    ;tvlct,r,g,b,/get
    ;vv1=[0,255]
    ;vv2=[0,255]
    ;vv3=[0,255]
    ;tvlct,vv1,vv2,vv3,0


    ;print, "center= ", center
    ;print, "center (original)= ", orig_center
    center= [0,0,0]

    xmin= -xlen+center[0]
    xmax=  xlen+center[0]
    ymin= -xlen+center[1]
    ymax=  xlen+center[1]


    tv, Pic, x0, y0, xsize=(x1-x0), ysize=(y1-y0), /normal


    !p.position=[x0,y0,x1,y1]
    !p.ticklen=0.03

    ; creates axes and plot style
    ; ---------------------------
    if x0 lt 0.2 then begin    ; in this case, print x tick labels
	plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, color= 0,  $
		      xstyle=1, ystyle=1, xtitle=xtit, ytitle=ytit, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
		      ;xticks=10, $
		      ;yticks=10, $
		      ytickformat='(a1)'
    endif else begin
	plot,[0],[0], psym=3, xrange=[xmin,xmax], yrange=[ymin,ymax], /noerase, /nodata, color= 0,  $
		      xstyle=1, ystyle=1, xtitle=xtit, ytitle=ytit, $
                      xcharsize=0.01, ycharsize=0.01, xthick=4.0, ythick=4.0, charthick=3.0, $
		      ;xticks=10, $
		      ;yticks=10, $
		      xtickformat='(a1)', ytickformat='(a1)'
    endelse


     ;xyouts, x0+0.04, y1-0.08, fload_timelbl(1,2,/noteq), /normal, size= 1.7, charthick=3.0, color= 0
     xyouts, x0+0.02, y1-0.05, msg1, /normal, size= 1.3, charthick=4.0, color= 0
     xyouts, x0+0.02, y1-0.10, msg2, /normal, size= 1.3, charthick=4.0, color= 0




                ; -----------------------------
                ;  X marks the black hole spot
                ; -----------------------------
                if keyword_set(showbhs) and fload_npart(5) ge 1 then begin

                        loadct, 4          ; blue/green/red/yellow (std)
                        tvlct,r,g,b,/get
                        v1=[0,255]
                        v2=[0,255]
                        v3=[0,255]
                        tvlct,v1,v2,v3,0

                        bhpsym= 7
                        bhsymsize= 2.0
                        bhcolor= 220

                        n_bh= fload_npart(5)
                        bhid= fload_blackhole_id(1)
                        bhid1= bhid[0]
                        bh_center_1= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid1)
                        
                        bh_new= bh_center_1 - orig_center
                        if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin                     
                           bh_c= bh_new
                           theta= !PI*rotate_theta/180.0
                           phi= !PI*rotate_phi/180.0
                           ; 4. rotate y (theta) and then x (phi)
                           ;bh_new[0]= bh_c[0]*cos(theta) + bh_c[2]*sin(theta)
                           ;bh_new[1]= -bh_c[0]*(sin(phi)*sin(theta)) + bh_c[1]*cos(phi) + bh_c[2]*(cos(theta)*sin(phi))
                           ;bh_new[2]= -bh_c[0]*(cos(phi)*sin(theta)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(theta)*cos(phi))
                           ; 3. rotate x (theta) and then z (phi)
                           ;bh_new[0]= bh_c[0]*cos(phi) + bh_c[1]*(sin(phi)*cos(theta)) + bh_c[2]*(sin(phi)*sin(theta))
                           ;bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*(cos(phi)*cos(theta)) + bh_c[2]*(cos(phi)*sin(theta))
                           ;bh_new[2]= -bh_c[1]*sin(theta) + bh_c[2]*cos(theta)
                           ; 1. rotate y (theta) and then z (phi)
                           ;bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(phi)*sin(theta))
                           ;bh_new[1]= bh_c[0]*(cos(theta)*sin(phi)) + bh_c[1]*cos(phi) + bh_c[2]*(sin(phi)*sin(theta))
                           ;bh_new[2]= -bh_c[0]*sin(theta) + bh_c[2]*cos(theta)
	        	   ; 2.5 first around z axis (phi) - clockwise
	        	   ; then around y axis (theta)  - counter clockwise
	        	   bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) + bh_c[1]*(cos(theta)*sin(phi)) - bh_c[2]*sin(theta)
	        	   bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*cos(phi)
	        	   bh_new[2]= bh_c[0]*(sin(theta)*cos(phi)) + bh_c[1]*(sin(theta)*sin(phi)) + bh_c[2]*cos(theta)
                        endif
                        
                        ; for the moment, just plot xy black hole position
                        oplot, [bh_new[0]], [bh_new[1]], psym=bhpsym, thick=6.0, $
                                        color=bhcolor, symsize=bhsymsize
                        
                        if n_bh gt 1 then begin
                                bhid2= bhid[1]
                                bh_center_2= fload_blackhole_xyz('xyz',center=[0,0,0],idtofollow=bhid2)
                                bh_new= bh_center_2 - orig_center
                                if keyword_set(rotate_phi) or keyword_set(rotate_theta) then begin
                                   bh_c= bh_new
                                   theta= !PI*rotate_theta/180.0
                                   phi= !PI*rotate_phi/180.0
                                   ; 4. rotate y (theta) and then x (phi)
                                   ;bh_new[0]= bh_c[0]*cos(theta) + bh_c[2]*sin(theta)                                 
                                   ;bh_new[1]= -bh_c[0]*(sin(phi)*sin(theta)) + bh_c[1]*cos(phi) + bh_c[2]*(cos(theta)*sin(phi))
                                   ;bh_new[2]= -bh_c[0]*(cos(phi)*sin(theta)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(theta)*cos(phi))
                                   ; 3. rotate x (theta) and then z (phi)
                                   ;bh_new[0]= bh_c[0]*cos(phi) + bh_c[1]*(sin(phi)*cos(theta)) + bh_c[2]*(sin(phi)*sin(theta))
                                   ;bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*(cos(phi)*cos(theta)) + bh_c[2]*(cos(phi)*sin(theta))
                                   ;bh_new[2]= -bh_c[1]*sin(theta) + bh_c[2]*cos(theta)                              
				   ; 1. rotate y (theta) and then z (phi)
				   ;bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) - bh_c[1]*sin(phi) + bh_c[2]*(cos(phi)*sin(theta))
				   ;bh_new[1]= bh_c[0]*(cos(theta)*sin(phi)) + bh_c[1]*cos(phi) + bh_c[2]*(sin(phi)*sin(theta))
				   ;bh_new[2]= -bh_c[0]*sin(theta) + bh_c[2]*cos(theta)
				   ; 2.5 first around z axis (phi) - clockwise
				   ; then around y axis (theta)  - counter clockwise
				   bh_new[0]= bh_c[0]*(cos(theta)*cos(phi)) + bh_c[1]*(cos(theta)*sin(phi)) - bh_c[2]*sin(theta)
				   bh_new[1]= -bh_c[0]*sin(phi) + bh_c[1]*cos(phi)
				   bh_new[2]= bh_c[0]*(sin(theta)*cos(phi)) + bh_c[1]*(sin(theta)*sin(phi)) + bh_c[2]*cos(theta)
                                endif
                                
                                print, "original center= ",orig_center
                                print, "Blackhole #1: ",bhid1, bh_center_1
                                print, "Blackhole #2: ",bhid2, bh_center_2
                                
                                oplot, [bh_new[0]], [bh_new[1]], psym=bhpsym, thick=6.0, $                                  
                                        color=bhcolor, symsize=bhsymsize
                                
                        endif else begin
                                print, "original center= ",orig_center
                                print, "Blackhole: ", bhid1, bh_center_1
                        endelse
                
                endif



                ; -------------------------------
                ;  Draw Box (slit, for instance)
                ; -------------------------------
                ;drawbox= 1
                if keyword_set(drawbox) then begin
                   slitcolor= 1
                   slitthick= 4.0
                   oplot, [-xslit, -xslit], [-yslit,yslit], psym=-3, color= slitcolor, thick= slitthick
                   oplot, [-xslit, xslit], [yslit,yslit], psym=-3, color= slitcolor, thick= slitthick  
                   oplot, [xslit, xslit], [-yslit,yslit], psym=-3, color= slitcolor, thick= slitthick
                   oplot, [-xslit, xslit], [-yslit,-yslit], psym=-3, color= slitcolor, thick= slitthick
                endif

		;print, "xyslit= ",xslit, yslit



		; ---------------------------
		;  Stellar Density Contour
		; ---------------------------
		overplot_contours, ContourMap, x0, y0, x1, y1, hmsb=HalfMassSB

		; now turn this off!
		;ContourMap= 0   - bad, we loose our map information
		;if keyword_set(ContourMap) then begin
                ;    levels = 12
                ;    ;step = (Max(XrayPic) - Min(XrayPic)) / levels
                ;    step = (256)/levels
                ;    ;userLevels = IndGen(levels) * step + Min(XrayPic)
                ;    userLevels = IndGen(levels) * step
                ;    ;userLevels = [0,userLevels[2:5]]
                ;    ;userLevels = [userLevels[3],userLevels[5:11]]
                ;    userLevels = userLevels[2:11]
                ;    print, "userLevels= ",userLevels
;
		;    ; parse map
		;    FixedMapIdx= where(ContourMap eq 1)
		;print, "Max ContourMap= ", max(ContourMap)
		;    FixedMap= ContourMap * (256.0/max(ContourMap))
		;    FixedMap= 256-FixedMap
		;    if FixedMapIdx(0) ne -1 then FixedMap(FixedMapIdx)= 1

		;    ; smooth the map
		;    ;NewMap= min_curve_surf(ContourMap, /regular)
		;    NewMap= smooth(FixedMap, 10, /edge_truncate)

                ;    ; load contours into variables (this doesn't print)
                ;    ;contour, ContourMap, path_xy=xy_contours, path_info=info, levels=userLevels, $
                ;    contour, NewMap, path_xy=xy_contours, path_info=info, levels=userLevels, $
                ;                        min_value=2, /PATH_DATA_COORDS

                    ;loadct, 4
                    ;tvlct,r1,g1,b1,/get
                    ;v12=[0,255]
                    ;v22=[0,255]
                    ;v32=[0,255]
                    ;tvlct,v12,v22,v32,0        ; define colors 0 and 1 as black and white

                ;    print, "N_contours= ", n_elements(info)
                
                ;    for i=0, (n_elements(info)-1) do begin
                ;        cnt_index= [indgen(info(i).N),0]
                ;        xycnt= xy_contours(*,info(i).offset+cnt_index)
                ;        n_level= info(i).level
;print, "i= ",i,"  n_level= ",n_level
                ;        xycnt[0,*]= xycnt[0,*]*(x1-x0)/480.0 + x0
                ;        xycnt[1,*]= xycnt[1,*]*(y1-y0)/480.0 + y0
                ;        idx= where(xycnt[0,*] eq x0)
                ;        if idx(0) eq -1 then begin
                ;           ;if n_level eq 0 then plots, xycnt, /normal, color= 150, thick=4.0
                ;           ;if n_level eq 1 then plots, xycnt, /normal, color= 150, thick=3.0
                ;           ;if n_level eq 2 then plots, xycnt, /normal, color= 150, thick=2.0
                ;           ;if n_level eq 3 then plots, xycnt, /normal, color= 150, thick=1.0, linestyle= 1
                ;           ;if n_level eq 4 then plots, xycnt, /normal, color= 130
		;	   ;if (n_level ge 0 and n_level le 4) then plots, xycnt, /normal, color= 0, thick=2.0
		;	   if n_level gt 0 then plots, xycnt, /normal, color= 0, thick=2.0
                ;           ;plots, xy_contours(*,info(i).offset+cnt_index), /normal, color= 150
                ;        endif
                ;    endfor

		;endif


                ; -----------------------
                ;  ColorBar #1
                ; 
                ; The original colorbar.
                ; Oriented on the right
                ; side of the plot, keyword
                ; also affects the plot
                ; positioning above.
                ; ----------------
		colorbar= 1
                if keyword_set(colorbar) then begin

                invertcolors= 1    ; warning: this should match the value in contour_makepic
		if keyword_set(velocitymap) then begin
		   set_maxden= 150.0
		   set_dynrng= 2.0*set_maxden
		endif
		if keyword_set(dispersionmap) then begin
		   set_maxden= 200.0
		   set_dynrng= 100.0
		endif
                MaxDensXY= set_maxden
                DynRange= set_dynrng
                ma=MaxDensXY
                mi=MaxDensXY-DynRange

                  bar= REPLICATE(1B, 10)#BINDGEN(256)
                  if invertcolors eq 1 then bar= 255-bar
                  idx= where(bar eq 0 or bar eq 1)
                  if idx(0) ne -1 then begin
                        ;print, "idx has ", n_elements(idx), idx(0)
                        ;colornumber= idx(0)-1
                        ;if colornumber lt 0 then colornumber= 11     ; if it's bar[*,0] then set to bar[*,1]
                        bar(idx)= 2   ; bar(colornumber)
                  endif
                  barwidth= 0.025

                  tv, bar, x1+0.01, y0+0.04, xsize=barwidth, ysize=(y1-y0-0.08), /normal

                  !p.position= [x1+0.01,y0+0.04,x1+0.01+barwidth,y1-0.04]
		  !p.ticklen= 0.20

                  plot, [0],[0], psym=3, xrange=[0,1], yrange=[mi,ma], $
                        /normal, $
                        /noerase, color= 0, xstyle=1, ystyle=1, $
                        xthick=4.0, ythick=4.0, $
                        xticks=1, xtickformat='(a1)', ytickformat='(a1)', $
                        /nodata

		  if keyword_set(velocitymap) then ytit= '!6Velocity (km s!E-1!N)' 
		  if keyword_set(dispersionmap) then ytit= '!7r!6 (km s!E-1!N)'
                  axis, yaxis=1, yrange=[mi,ma], ystyle=1, /normal, $
                        ycharsize=0.80, charthick=3.0, ytitle= ytit   ;, $ 
                        ;ytickformat='(a1)', ymargin= [1,6]
                  ;axis, yaxis=1, yrange=[mi,ma], ystyle=1, /normal, $
                        ;ycharsize=0.80, charthick=2.0
                endif




end



; =======================================================================================
;
;  Add procedure to control contours and ellipse fitting.
;
;
pro overplot_contours, ContourMap, x0, y0, x1, y1, hmsb=hmsb, $
				xlen=xlen, $
				fitellipse=fitellipse, $
				pa=pa


	Pic= ContourMap
	HalfMassSB= hmsb

                idx=where(Pic eq 1)
                Pic= 256-Pic
                if idx(0) ne -1 then Pic(idx)= 1
                print, "before: min/max= ",min(Pic), max(Pic)
                print, "n > 240= ", n_elements(where(Pic gt 240))
                print, "n < 20= ", n_elements(where(Pic lt 20))
                ; smooth?
                width= 10      ; try other things?
                newMap= smooth(Pic, width, /edge_truncate)
                Pic= newMap
                print, "after: min/max= ",min(Pic), max(Pic)

                levels = 8
                ;levels = 16
                ;step = (Max(Pic) - Min(Pic)) / levels
                ;userLevels = IndGen(levels) * step + Min(Pic)
                step = (256)/levels
                userLevels = IndGen(levels) * step
                ;userLevels = userLevels[2:7]
                ;userLevels = userLevels[10:15]
                print, "userLevels= ",userLevels
                print, '256-HalfMassSB= ', 256-HalfMassSB
                ; this is 256 minus because we flip Pic around above
                for i=2,n_elements(userLevels)-1 do begin
                   if userLevels[i] gt (256-HalfMassSB) then begin
                        nuL=[nuL,256-HalfMassSB,userLevels[i:7]]
                        level_to_fit= i-2
                        break
                   endif else begin
                        if i eq 2 then nuL= [userLevels[i]] else nuL=[nuL,userLevels[i]]
                   endelse
                endfor
                print, "new userLevels= ",nuL
                userLevels= nuL
                contour_to_fit= -1

                ; load contours into variables (this doesn't print)
                contour, Pic, path_xy=xy_contours, path_info=info, levels=userLevels, $
                                min_value=2, /PATH_DATA_COORDS


                for i=0, (n_elements(info)-1) do begin
                        cnt_index= [indgen(info(i).N),0]
                        xycnt= xy_contours(*,info(i).offset+cnt_index)
                        n_level= info(i).level
print, "i= ", i, "  level= ",n_level, "       N= ", info(i).N
                        xycnt[0,*]= xycnt[0,*]*(x1-x0)/480.0 + x0
                        xycnt[1,*]= xycnt[1,*]*(y1-y0)/480.0 + y0
                        idx= where(xycnt[0,*] eq x0)
                        if idx(0) eq -1 then begin
                           ;if n_level eq 0 then plots, xycnt, /normal, color= 0, thick=2.0
                           ;if n_level eq 1 then plots, xycnt, /normal, color= 0, thick=2.0
                           ;if n_level eq 2 then plots, xycnt, /normal, color= 0, thick=2.0
                           ;if n_level eq 3 then plots, xycnt, /normal, color= 0, thick=2.0
                           ;if n_level eq 4 then plots, xycnt, /normal, color= 130
                           ;plots, xy_contours(*,info(i).offset+cnt_index), /normal, color= 150
                           ;plots, xycnt, /normal, color= 0, thick=2.0
                           plots, xycnt, /normal, color= 1, thick=2.0
                        endif

                        ;if info(i).N gt 200 then contour_to_fit= i
                        if info(i).level eq level_to_fit and contour_to_fit eq -1 then contour_to_fit= i
                endfor




        ; -------------------------
        ;
        ;  Fit Ellipse to Contour
        ;
        ;  Use the MPFit routine
        ; to fit and ellipse to one
        ; of the above contours.
        ; -------------------------
        ;fitellipse= 1
        ;fitellipse= 0
        if keyword_set(fitellipse) then begin

           ; --------------------------
           if contour_to_fit eq -1 then i= 3 else i= contour_to_fit
           print, "fit i,N,levle= ",i,info(i).N,info(i).level
           cnt_index= [indgen(info(i).N),0]
           xy= xy_contours(*,info(i).offset+cnt_index)

           x=fltarr((size(xy))[2])
           y=fltarr((size(xy))[2])

           x(*)= xy[0,*]
           y(*)= xy[1,*]

           ; x,y are in pixel units at this point

           ; here's where we actually fit to ellipse
           params=mpfitellipse(x,y, /tilt)

           ; a and b are in normalized coordinates
           if params(0) GT params(1) then begin
                a= params(0)
                b= params(1)
		pa_rad= params(4)
           endif else begin
                a= params(1)
                b= params(0)
		if params(4) gt 0.0 then pa_rad= !PI*0.5 + params(4)
		if params(4) lt 0.0 then pa_rad= -0.5*!PI + params(4)
           endelse
           ellipticity= 1-b/a
           print, " xxxxxxxxxxxxxxxx "
           print, "  Fit Results "
           print, "Semiaxes:   a=", a, "   b=",b
           print, "ellipticity(1-b/a)= ", ellipticity

           ; convert to physics coordinates: normal= (2.0*xlen)/(x1-x0)
           ; CAREFUL: we need the image to be square to make this easy
           ; conversion
           ;a= a*2.0*xlen/(x1-x0)
           ;b= b*2.0*xlen/(x1-x0)   -> old methods

           ; new methods has x,y in pixel coordinates
           a= a*2.*xlen/480.0
           b= b*2.*xlen/480.0

           print, "(phy units) a=", a,"   b=",b
	   print, "  params(4, in deg)=", params(4)*180.0/!PI
           print, "       phi (in deg)=", pa_rad*180.0/!PI
           print, "             params=", params

;pa_rad= 110.0
;print, "MANUALLY fixing pa_rad to ", pa_rad
;pa_rad= 110.0*!PI/180.0

           ; replot the contour we're fitting to.
;          plots, x, y, psym=-3, color=getcolor('red'), /normal

           ; plot the ellipse 
           ; -----------------
           phi = dindgen(101)*2D*!dpi/100
           ; data coord
           ;e_x = 2.*xlen/480.0 * params(0)*cos(phi)
           ;e_y = 2.*xlen/480.0 * params(1)*sin(phi)
           ; normal coord
           e_x = (x1-x0)*1./480.0 * params(0)*cos(phi)
           e_y = (y1-y0)*1./480.0 * params(1)*sin(phi)

           ; rotate to tilted frame
           ; ------------------------
           if params(4) NE 0 then begin
                e_x_prime = e_x*cos(params(4)) + e_y*sin(params(4))
                e_y_prime = -e_x*sin(params(4)) + e_y*cos(params(4))
                e_x = e_x_prime
                e_y = e_y_prime
		pa= pa_rad*180.0/!PI
           endif

           ; relocate to center
           ; ------------------
           ; data coord norm.
           ;e_x = 2.*xlen/480.0 * params(2) + e_x
           ;e_y = 2.*xlen/480.0 * params(3) + e_y
           ; normal coordinates
           e_x = x0 + (x1-x0)*1./480.0 * params(2) + e_x
           e_y = y0 + (y1-y0)*1./480.0 * params(3) + e_y

           loadct, 4
           tvlct,r,g,b,/get
           v1=[0,255]
           v2=[0,255]
           v3=[0,255]
           tvlct,v1,v2,v3,0        ; define colors 0 and 1 as black and white

           plots, e_x, e_y, color= 150, thick= 4.0, /normal

           ;oplot, e_x, e_y, color= 150, thick= 4.0

           ; labels, if we want 'em
           ecclbl= strcompress(string(ellipticity),/remove_all)
           ecclbl= '!7e!6='+strmid(ecclbl,0,4)
           xyouts, x0+0.03, y1-0.10, ecclbl, /normal, size=1.4, color= 0, charthick= 3.0
           albl= strcompress(string(a),/remove_all)
           ;albl= '!6R!D!8a!6!N='+strmid(albl,0,4)
           albl= '!8a!6='+strmid(albl,0,4)
           xyouts, x0+0.03, y1-0.06, albl, /normal, size=1.4, color= 0, charthick= 3.0

        endif




end











;#########################################################################
;
;
;   Slit Time Evolution
;
;
;
;
;#########################################################################


pro slit_time_vel, junk


if not keyword_set(junk) then begin
        print, " "
        print, " slit_time_vel, junk"
        print, " "
        print, " "
        return
endif


; --------------------------
;  Setup stuff
; --------------------------

center=[0,0,0]

; view 1
rotate_theta= 90.0
;rotate_theta= 0.0
rotate_phi= 0.0
;rotate_phi= 90.0
;rotate_phi= 105.0


; figure 1
;frun= "/raid4/tcox/vc3vc3b"
;frun= "/raid4/tcox/vc3vc3e"
frun= "/raid4/tcox/vc3vc3f"
;frun= "/raid4/tcox/cvc3vc3f"
;frun= "/raid4/tcox/vc3vc3l"
;frun= "/raid4/tcox/vc3vc3o"

startsnap= 13
endsnap= 30
nsnaps= endsnap-startsnap

;center= [-0.6,-1.9,0.6]
center= [0.0,0.0,0.0]


manual_entry= 1
;manual_entry= 0

Ttime= [1.3, 1.4, 1.5, 1.6, 1.7, 1.8, $
	1.9, 2.0, 2.1, 2.2, 2.3, 2.4, $
	2.5, 2.6, 2.7, 2.8, 2.9, 3.0]

a_t=   [3.36, 3.44, 3.33, 3.23, 3.22, 3.18, $
	3.16, 3.07, 3.07, 3.18, 3.16, 3.19, $
	3.17, 3.16, 3.20, 3.14, 3.16, 3.21]

Re_t=  [3.02, 2.89, 2.82, 2.83, 2.76, 2.77, $
	2.74, 2.73, 2.73, 2.75, 2.76, 2.77, $
	2.74, 2.78, 2.77, 2.77, 2.79, 2.79]

ellip_t= [0.20, 0.27, 0.29, 0.28, 0.28, 0.25, $
	  0.27, 0.24, 0.24, 0.27, 0.26, 0.27, $
	  0.31, 0.29, 0.30, 0.29, 0.29, 0.31]

V_max_t= [111.3, 111.2, 113.8, 115.5, 107.2, 100.5, $
	  101.2, 107.2, 116.2, 109.2, 107.2, 112.7, $
	  121.2, 119.1, 109.5, 114.1, 105.2, 113.1]

Sigma_e_t= [164.3, 149.3, 152.2, 148.1, 144.8, 143.5, $
	    143.7, 146.0, 148.1, 147.5, 147.4, 144.1, $
	    144.8, 143.8, 144.4, 143.5, 144.7, 145.5]



; slit dimensions
; -----------------
slit_width = 0.5
slit_bins = 20
slit_len = 5.0

; plot range
slitmax = 1.2*slit_len

vmax= 150.0       ; also set below
vscale= 200.0

; image size
;xlen= 1.2*slit_len
xlen= 20.0



slitr_t= fltarr(nsnaps+1,slit_bins)
slitv_t= fltarr(nsnaps+1,slit_bins)
slitsigma_t= fltarr(nsnaps+1,slit_bins)
;V_max_t= fltarr(nsnaps+1)
;Sigma_e_t= fltarr(nsnaps+1)
;Ttime= fltarr(nsnaps+1)



; ---------------
;  Do it Up
; ---------------

if manual_entry eq 0 then begin
    for ii=0,nsnaps do begin

	snapnum= startsnap+ii

	ok=fload_snapshot_bh(frun,snapnum)

	center= fload_center_alreadycomp(1)
	orig_center= center



	; Do Stars
	;===========
	;do_stars_slit= 0
	do_stars_slit= 1
	if do_stars_slit eq 1 then begin
		x= fload_allstars_xyz('x', center=[0,0,0])
		y= fload_allstars_xyz('y', center=[0,0,0])
		z= fload_allstars_xyz('z', center=[0,0,0])
		vx= fload_allstars_v('x')
		vy= fload_allstars_v('y')
		vz= fload_allstars_v('z')
		mass= fload_allstars_mass(1)
	endif


	; Do Gas
	;===========
	do_gas_slit= 0
	;do_gas_slit= 1
	if do_gas_slit eq 1 then begin
	        x= fload_gas_xyz('x', center=[0,0,0])
	        y= fload_gas_xyz('y', center=[0,0,0])
	        z= fload_gas_xyz('z', center=[0,0,0])
	        vx= fload_gas_v('x')
	        vy= fload_gas_v('y')
	        vz= fload_gas_v('z')
	        mass= fload_gas_mass(1)
	endif

        


	if keyword_set(center) then begin
	     print, "process slit, changing center to= ",center
	     x=x-center[0]
	     y=y-center[1]
	     z=z-center[2]
	endif


	; Now we need to rotate for the slit
	; ------------------------------------
	;if (rotate_theta gt 0) or (rotate_phi gt 0) or (total(abs(center)) gt 0) then begin
	if (rotate_theta gt 0) or (rotate_phi gt 0) then begin
	;if keyword_set(rotate_theta) or keyword_set(rotate_phi) then begin

	        ; transform to radians
	        theta= !PI*rotate_theta/180.0
	        phi= !PI*rotate_phi/180.0

	        ; 1. first around y axis (theta)
	        ; then around z axis (phi)
	        ;x_new= x*(cos(theta)*cos(phi)) - y*sin(phi) + z*(cos(phi)*sin(theta))
	        ;y_new= x*(cos(theta)*sin(phi)) + y*cos(phi) + z*(sin(phi)*sin(theta))
	        ;z_new= -x*sin(theta) + z*cos(theta)

	        ; 2. first around z axis (phi)
	        ; then around y axis (theta)  - both counter clockwise
	        ;x_new= x*(cos(theta)*cos(phi)) + y*(cos(theta)*sin(phi)) + z*sin(theta)
	        ;y_new= -x*sin(phi) + y*cos(theta)
	        ;z_new= x*(sin(theta)*cos(phi)) - y*(sin(theta)*sin(phi)) + z*cos(theta)

	        ; 2.5 first around z axis (phi) - clockwise
	        ; then around y axis (theta)  - counter clockwise
	        x_new= x*(cos(theta)*cos(phi)) + y*(cos(theta)*sin(phi)) - z*sin(theta)
	        y_new= -x*sin(phi) + y*cos(phi)
	        z_new= x*(sin(theta)*cos(phi)) + y*(sin(theta)*sin(phi)) + z*cos(theta)

	        ; 3. first around x axis (theta)
	        ; then around z axis (phi)
	        ;x_new= x*cos(phi) + y*(sin(phi)*cos(theta)) + z*(sin(phi)*sin(theta))
	        ;y_new= -x*sin(phi) + y*(cos(phi)*cos(theta)) + z*(cos(phi)*sin(theta))
	        ;z_new= -y*sin(theta) + z*cos(theta)

	        ; 4. first around y axis (theta)
	        ; then around x axis (phi)
	        ;x_new= x*cos(theta) + z*sin(theta)
	        ;y_new= -x*(sin(phi)*sin(theta)) + y*cos(phi) + z*(cos(theta)*sin(phi))
	        ;z_new= -x*(cos(phi)*sin(theta)) - y*sin(phi) + z*(cos(theta)*cos(phi))

	        x= x_new
	        y= y_new
	        z= z_new


	        ; rotate the velocities
	        ; ----------------------
		; 1.
	        ;vx_new= vx*(cos(theta)*cos(phi)) - vy*sin(phi) + vz*(cos(phi)*sin(theta))
	        ;vy_new= vx*(cos(theta)*sin(phi)) + vy*cos(phi) + vz*(sin(phi)*sin(theta))
	        ;vz_new= -vx*sin(theta) + vz*cos(theta)

	        ; 2.5
	        vx_new= vx*(cos(theta)*cos(phi)) + vy*(cos(theta)*sin(phi)) - vz*sin(theta)
	        vy_new= -vx*sin(phi) + vy*cos(phi)
	        vz_new= vx*(sin(theta)*cos(phi)) + vy*(sin(theta)*sin(phi)) + vz*cos(theta)

		; 3.
	        ;vx_new= vx*cos(phi) + vy*(sin(phi)*cos(theta)) + vz*(sin(phi)*sin(theta))
	        ;vy_new= -vx*sin(phi) + vy*(cos(phi)*cos(theta)) + vz*(cos(phi)*sin(theta))
	        ;vz_new= -vy*sin(theta) + vz*cos(theta)

		; 4. 
	        ;vx_new= vx*cos(theta) + vz*sin(theta)
	        ;vy_new= -vx*(sin(phi)*sin(theta)) + vy*cos(phi) + vz*(cos(theta)*sin(phi))
	        ;vz_new= -vx*(cos(phi)*sin(theta)) - vy*sin(phi) + vz*(cos(theta)*cos(phi))


	        vx= vx_new
	        vy= vy_new
	        vz= vz_new

	endif



	; calculate slit velocities
	process_slit, x, y, vz, mass, slit, v_mean, v_disp, $
			slit_width=slit_width, $
			slit_len=slit_len, $
			slit_bins=slit_bins, $
			v_mean_unweighted=v_mean_unweighted, $
			v_disp_unweighted=v_disp_unweighted, $
			mags=mags, $
			v_mean_magwt=v_mean_magwt, $
			v_disp_magwt=v_disp_magwt



	; now v_mean and v_disp contain velocity info
	slitr_t[ii,*]= slit
	slitv_t[ii,*]= v_mean
	slitsigma_t[ii,*]= v_disp


	; Suvendra's procedure to determine V_max and Sigma_e
	; ----------------------------------------------------
	; determine V_max
	Vmax= max(v_mean) & imax= where(v_mean eq Vmax)
	Vmin= min(v_mean) & imin= where(v_mean eq Vmin)

	Vmin= 0.0
	Vmax= 0.0

	imincnt= 0
	imaxcnt= 0
	
	for i=-1,1 do begin
	    if (((imin+i) ge 0) and ((imin+i) lt n_elements(v_mean))) then begin
		Vmin= Vmin + v_mean[imin+i]
		imincnt++
	    endif
	    if (((imax+i) ge 0) and ((imax+i) lt n_elements(v_mean))) then begin
		Vmax= Vmax + v_mean[imax+i]
		imaxcnt++
	    endif
	endfor

	if imaxcnt ne 0 then Vmax=Vmax/double(imaxcnt)
	if imincnt ne 0 then Vmin=Vmin/double(imincnt)
	V_max= 0.5*(abs(Vmax)+abs(Vmin))

	Sigma_e= 0
	ire= where(abs(slit) le 1.0)
	if ire(0) ne -1 then Sigma_e= total(v_disp(ire))/double(n_elements(ire))

	V_max_t[ii]= V_max
	Sigma_e_t[ii]= Sigma_e
	Ttime[ii]= fload_time(1)

    endfor
endif



; =====================
; now plot things
; =====================


plot_lines_versus_time= 0
;plot_lines_versus_time= 1

if plot_lines_versus_time eq 1 then begin
	; -------------
	x0= 0.21
	x1= 0.98
    
	y0= 0.15
	y1= 0.98


	vmax= 150.0       ; also set above
	vscale= 250.0


	; Velocity vs. Time
	; -------------------

	; plot it up
	; -----------
	filename='slit_vel_time.eps'
	xaxistitle='R (h!u-1!N kpc)'
	yaxistitle='V(R) (km s!E-1!N)'
	xmax = slitmax
	xmin = -slitmax
	ymax = vmax
	ymin = -vmax


	; get print stuff ready
	; ----------------------
	initialize_plotinfo, 1
	setup_plot_stuff, 'ps', filename=filename, colortable= 1, newxsize=13, newysize=12

	thiscolor= 30

	!p.position= [x0, y0, x1, y1]

	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
	        xrange=[xmin,xmax], yrange=[ymin,ymax], $
	        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	        xtitle=xaxistitle, ytitle=yaxistitle, /nodata



	for ii=0,nsnaps do begin
		oplot,slitr_t[ii,*],slitv_t[ii,*],psym=-3,color=thiscolor, thick=3.0
		thiscolor= thiscolor+10
	endfor


	; zero
	; ------ 
	x=[xmin,xmax]
	y=[0,0]
	oplot, x, y, psym=-3, linestyle=2, color= 0, thick=2.0


	; general info
	xyouts, 0.65, 0.85, fload_fid(frun), /normal, charthick=3.0, size=1.33, color=0
	xyouts, 0.65, 0.80, 'xy projection', /normal, charthick=3.0, size=1.33, color=0

	device, /close






	; Dispersion vs. Time
	; ----------------------

	filename='slit_sigma_time.eps'
	xaxistitle='R (h!u-1!N kpc)'
	yaxistitle='!7r!3(R) (km s!E-1!N)'
	xmax = slitmax
	xmin = -slitmax
	ymax = vscale
	ymin = 0.0


	; get print stuff ready
	; ----------------------
	initialize_plotinfo, 1
	setup_plot_stuff, 'ps', filename=filename, colortable= 1, newxsize=13, newysize=12

	thiscolor= 30

	!p.position= [x0, y0, x1, y1]

	plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
	        xrange=[xmin,xmax], yrange=[ymin,ymax], $
	        xcharsize=1.50, ycharsize=1.50, xthick=4.0, ythick=4.0, charthick=4.0, $
	        xtitle=xaxistitle, ytitle=yaxistitle, /nodata



	for ii=0,nsnaps do begin
		oplot, slitr_t[ii,*], slitsigma_t[ii,*], psym=-3, color=thiscolor
		thiscolor= thiscolor+10
	endfor



	; zero
	; ------ 
	x=[xmin,xmax]
	y=[0,0]
	oplot, x, y, psym=-3, linestyle=2, color= 0, thick=2.0


	xyouts, 0.65, 0.45, fload_fid(frun), /normal, charthick=3.0, size=1.33, color=0
	xyouts, 0.65, 0.40, 'xy projection', /normal, charthick=3.0, size=1.33, color=0

	device, /close

endif



; -------------------------------------------
; -------------------------------------------
; now do the V_max and Sigma_e versus time
; -------------------------------------------
; -------------------------------------------


filename='slit_time.eps'

initialize_plotinfo, 1
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=13, newysize=20

x0= 0.15 & x1= 0.99
y0= 0.08 & ysize=0.182
y1= y0+ysize 
y2= y0+ysize+ysize
y3= y0+ysize+ysize+ysize
y4= y0+ysize+ysize+ysize+ysize
y5= y0+ysize+ysize+ysize+ysize+ysize

xmax= 1.9 & xmin= 0.0
xaxistitle='!6Time After Merger (h!E-1!N Gyr)'
TaM= Ttime-1.2



; zero: ellipticity
; --------------------
yaxistitle='!6(h!E-1!N kpc)'
ymax = 3.6
ymin= 2.5
!p.position= [x0, y4, x1, y5]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, TaM, a_t, psym=-2, color=150, thick= 3.0
oplot, TaM, Re_t, psym=-5, color=50, thick= 3.0

xyouts, 0.6, 0.95, '!16a!6', /normal, color=150, size=1.2
xyouts, 0.6, 0.87, 'R!De!N', /normal, color=50, size=1.2



; zero: ellipticity
; --------------------
yaxistitle='!7e!6'
ymax = 0.8
ymin= 0.0
!p.position= [x0, y3, x1, y4]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, TaM, ellip_t, psym=-2, color=150, thick= 3.0




; first: V_max
; -------------
yaxistitle='!6V!Dmaj!N (km s!E-1!N)'
;ymax = 60.0
;ymin= -10.0
ymax = 145.0
ymin= 85.0
!p.position= [x0, y2, x1, y3]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, TaM, V_max_t, psym=-2, color=150, thick= 3.0



; second: Sigma_e
; ----------------
yaxistitle='!7r!6 (km s!E-1!N)'
ymax = 185.0
ymin= 125.0
!p.position= [x0, y1, x1, y2]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtickformat='(a1)', ytitle=yaxistitle, /nodata, /noerase

oplot, TaM, Sigma_e_t, psym=-2, color=150, thick= 3.0



; second: Sigma_e
; ----------------
yaxistitle='!6V!Dmaj!N / !7r!6'
ymax = 1.1
ymin= -0.1
!p.position= [x0, y0, x1, y1]
plot, [0], [0], psym=-3, linestyle=0, color= 0, xstyle=1, ystyle=1, $
        xrange=[xmin,xmax], yrange=[ymin,ymax], $
        xcharsize=1.20, ycharsize=1.20, xthick=4.0, ythick=4.0, charthick=4.0, $
        xtitle=xaxistitle, ytitle=yaxistitle, /nodata, /noerase

oplot, TaM, V_max_t/Sigma_e_t, psym=-2, color=150, thick= 3.0





;xyouts, 0.65, 0.95, fload_fid(frun), /normal, charthick=3.0, size=1.33, color=0
;xyouts, 0.65, 0.91, 'xy projection', /normal, charthick=3.0, size=1.33, color=0

device, /close


; -------------
;  Done
; -------------


end





;------------------------------------------------------------------------------------






;=============================================================
;
;  Velocity along a slit
;
;=============================================================


pro slit, frun, rt, rp, pa=pa, snapnum=snapnum, $
				filename=filename


;if not keyword_set(junk) then begin
if not keyword_set(frun) then begin
        print, " "
        print, " slit, frun, rt, rp, pa=pa, snapnum=snapnum,"
        print, "                      filename=filename"
        print, " "
        return
endif


; --------------------------
;  Setup stuff
; --------------------------

if not keyword_set(filename) then filename='slit_yo.eps'


; view 0
if keyword_set(rt) then rotate_theta= rt else rotate_theta= 0.0
if keyword_set(rp) then rotate_phi= rp else rotate_phi= 0.0


if keyword_set(pa) then pa_fixed= pa else pa_fixed= -1


;frun= "/raid4/tcox/vc3vc3b"
;frun= "/raid4/tcox/vc3vc3d"
;frun= "/raid4/tcox/vc3vc3e"
;frun= "/raid4/tcox/vc3vc3f"
;frun= "/raid4/tcox/cvc3vc3f"
;frun= "/raid4/tcox/cvc3vc3h"
;frun= "/raid4/tcox/vc3vc3k"
;frun= "/raid4/tcox/vc3vc3l"
;frun= "/raid4/tcox/vc3vc3o"
if not keyword_set(frun) then frun= "/raid4/tcox/vc3vc3b"

if not keyword_set(snapnum) then snapnum= 30



; slit dimensions
; -----------------
slit_width = 0.5
slit_bins = 26
slit_len = 5.0

; plot range
slitmax = 1.2*slit_len

vmax= 250.0
vscale= 200.0

; image size
;xlen= 1.2*slit_len
xlen= 8.0



; start it up
; -------------

; these are for 4x4 image that leaves room for labels
x0= 0.15
xs= 0.83
y0= 0.15
ys= 0.83

x1= x0+xs
y1= y0+ys



; ---------------
;  Do it Up
; ---------------

;ok=fload_initial_conditions(frun)
ok=fload_snapshot_bh(frun,snapnum)


; get print stuff ready
; ----------------------
initialize_plotinfo, 1
        
setup_plot_stuff, 'ps', filename=filename, colortable= 4, newxsize=13, newysize=12


center= fload_center_alreadycomp(1)
orig_center= center



; Stars
;========
x= fload_allstars_xyz('x', center=[0,0,0])
y= fload_allstars_xyz('y', center=[0,0,0])
z= fload_allstars_xyz('z', center=[0,0,0])
vx= fload_allstars_v('x')
vy= fload_allstars_v('y')
vz= fload_allstars_v('z')
mass= fload_allstars_mass(1)


; load magnitudes
TTime= float(fload_time(1))

ndisk= fload_npart(2)                               ; disk particles
nbulge= fload_npart(3)                              ; bulge particles
nstars= fload_npart(4)
npart= long(ndisk) + long(nbulge) + long(nstars)
N= npart

m= 1.0e+10*mass
age=fload_allstars_age(1)
age=float(TTime-age) 
zmets=fload_allstars_z(1)


; get the luminosities
print, "load luminosities"
load_all_stellar_luminosities, N, TTime, m, age, zmets, $
                                Lum_Bol, Lum_U, Lum_B, Lum_V, Lum_R, Lum_I, Lum_J, Lum_H, Lum_K, $
                                Lum_sdss_u, Lum_sdss_g, Lum_sdss_r, Lum_sdss_i, Lum_sdss_z   ;, $
                                ;/notsolar

; trap for any NaN's
idx=where(finite(Lum_B) eq 0)
if idx(0) ne -1 then begin
        Lum_Bol(idx)= 100.0 & Lum_U(idx)= 100.0 & Lum_B(idx)= 100.0 & Lum_V(idx)= 100.0 & Lum_R(idx)= 100.0
        Lum_I(idx)= 100.0 & Lum_J(idx)= 100.0 & Lum_H(idx)= 100.0 & Lum_K(idx)= 100.0
        Lum_sdss_u(idx)= 100.0 & Lum_sdss_g(idx)= 100.0 & Lum_sdss_r(idx)= 100.0
        Lum_sdss_i(idx)= 100.0 & Lum_sdss_z(idx)= 100.0
	print, "Fixed ",n_elements(idx), " luminosities when trapping for NaN values."
endif


mags= Lum_B


        

;  bottom-left, velocity field within slit
; ------------------------------------------
x_orig= x
y_orig= y
z_orig= z
vx_orig= vx
vy_orig= vy
vz_orig= vz
center= orig_center
process_slit_and_overplot, x_orig, y_orig, z_orig, vx_orig, vy_orig, vz_orig, mass, $
        x0, y0, x1, y1, $
        center=center, /labelit, $
        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
        slit_len=slit_len, slit_width=slit_width, $
        slit_bins=slit_bins, slitmax=slitmax, vmax=vmax, $
        xlen=xlen, pa=pa, mags=mags



device, /close

; -------------
;  Done
; -------------


end









;====================================================================================







